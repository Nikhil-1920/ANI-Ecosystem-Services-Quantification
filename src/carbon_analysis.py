"""
ANI Ecosystem Services — Week 6: Carbon Storage & Loss Estimation
=================================================================
Inputs  : data/processed/  (clipped GeoTIFFs)
Outputs : results/carbon_annual_loss_by_year.csv
          figures/carbon/agb_gedi_baseline_map.png
          figures/carbon/agb_cross_validation_scatter.png
          figures/carbon/carbon_annual_loss_timeseries.png
          figures/carbon/carbon_loss_hotspots_map.png

Run with: venv/bin/python src/carbon_analysis.py
"""

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import rasterio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy import stats
from pathlib import Path

# ── Directory Paths ────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).parent
PROC_DIR   = SCRIPT_DIR.parent / 'data' / 'processed'
RES_DIR    = SCRIPT_DIR.parent / 'results'
FIG_DIR    = SCRIPT_DIR.parent / 'figures'
RES_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from scipy.ndimage import binary_dilation

# ── Print Separator ────────────────────────────────────────────────
SEP = '=' * 60

# ── Split-map row boundaries (Ten Degree Channel) ────────────────────
ANDAMAN_ROW_END   = 11750
NICOBAR_ROW_START = 16200
_PAD = 200


def _split_extent(arr):
    h = arr.shape[0]
    return arr[:min(ANDAMAN_ROW_END + _PAD, h), :], \
           arr[max(0, NICOBAR_ROW_START - _PAD):, :]


def _tight_crop(arr, ref=None):
    """Crop arr to the bounding box of valid (non-NaN) pixels.
    If ref is provided, use its valid mask to determine the bbox
    (useful for co-registered arrays like carbon_map vs agb_grid)."""
    src = ref if ref is not None else arr
    if np.ma.is_masked(src):
        valid = ~src.mask
    else:
        valid = ~np.isnan(src)
    if not valid.any():
        return arr
    rows = np.any(valid, axis=1); cols = np.any(valid, axis=0)
    r0, r1 = np.where(rows)[0][[0, -1]]
    c0, c1 = np.where(cols)[0][[0, -1]]
    pad = 60
    r0 = max(0, r0-pad); r1 = min(arr.shape[0], r1+pad)
    c0 = max(0, c0-pad); c1 = min(arr.shape[1], c1+pad)
    return arr[r0:r1, c0:c1]


def split_imshow(fig, gs_l, gs_r, data, cmap, vmin, vmax,
                 label_left='Andaman Islands', label_right='Nicobar Islands'):
    """Render data as Andaman (left) + Nicobar (right) sub-axes."""
    and_r, nic_r = _split_extent(data)
    ax_a = fig.add_subplot(gs_l); ax_n = fig.add_subplot(gs_r)
    for ax, crop, title in [(ax_a, _tight_crop(and_r), label_left),
                             (ax_n, _tight_crop(nic_r), label_right)]:
        ax.set_facecolor('white')
        ax.imshow(crop, cmap=cmap, vmin=vmin, vmax=vmax, interpolation='nearest')
        ax.set_title(title, color='#1a1a2e', fontsize=11, fontweight='bold', pad=8)
        ax.axis('off')
    return ax_a, ax_n

# ── Biophysical Constants ──────────────────────────────────────────────
PIXEL_AREA_HA  = 0.09    # 30 m × 30 m = 900 m² = 0.09 ha
IPCC_CARBON_F  = 0.47    # Biomass-to-carbon fraction   (IPCC Tier 1)
CO2E_MW_RATIO  = 44/12   # Carbon-to-CO₂ molecular weight ratio
ROOT_SHOOT_R   = 0.24    # Root-to-shoot ratio, tropical (IPCC Tier 1)

# ── ESA WorldCover Class Codes ─────────────────────────────────────────
ESA_FOREST_CLS  = 10     # Tropical evergreen tree cover
ESA_MANGROVE_CLS = 95    # Mangrove forest


# ══════════════════════════════════════════════════════════════════════
# UTILITY
# ══════════════════════════════════════════════════════════════════════
def load_clipped_raster(filename: str):
    """Load a clipped GeoTIFF from PROC_DIR as a float array.

    Returns (data_array, nodata_value, raster_profile).
    Returns (None, None, None) if the file does not exist.
    """
    file_path = PROC_DIR / filename
    if not file_path.exists():
        print(f"  ⚠️  File not found: {filename}")
        return None, None, None
    with rasterio.open(file_path) as src:
        data_array  = src.read(1).astype(float)
        nodata_val  = src.nodata
        rast_profile = src.profile
    if nodata_val is not None:
        data_array = np.where(data_array == nodata_val, np.nan, data_array)
    data_array = np.where(data_array < -1e9, np.nan, data_array)
    return data_array, nodata_val, rast_profile


# ══════════════════════════════════════════════════════════════════════
# PART 1 — Baseline Aboveground Biomass from GEDI L4B
# ══════════════════════════════════════════════════════════════════════
def estimate_gedi_agb_baseline():
    """Map GEDI L4B aboveground biomass and stratify by forest type."""
    print(f"\n{SEP}")
    print("PART 1 — GEDI L4B Baseline Aboveground Biomass")
    print(SEP)

    agb_grid, _, _  = load_clipped_raster('ANI_GEDI_Biomass_Density_clipped.tif')
    landcover_grid, _, _ = load_clipped_raster('ANI_ESA_WorldCover_mosaic_clipped.tif')
    if agb_grid is None:
        print("  ❌  GEDI layer missing — aborting Part 1.")
        return None, None

    # Remove physically implausible GEDI outliers
    agb_grid = np.where(agb_grid < 0,   np.nan, agb_grid)
    agb_grid = np.where(agb_grid > 800, np.nan, agb_grid)

    valid_pixels = ~np.isnan(agb_grid)

    print(f"  Valid pixels     : {valid_pixels.sum():,}")
    print(f"  Mean AGB         : {np.nanmean(agb_grid):.2f} Mg/ha")
    print(f"  Max AGB          : {np.nanmax(agb_grid):.2f} Mg/ha")
    total_agb_tg = np.nansum(agb_grid) * PIXEL_AREA_HA / 1e6
    print(f"  Total AGB stock  : {total_agb_tg:.4f} Tg")

    # ── Stratify by forest ecosystem type ─────────────────────────────
    if landcover_grid is not None:
        strata = [
            ("Tropical Evergreen Forest", ESA_FOREST_CLS),
            ("Mangrove Forest",           ESA_MANGROVE_CLS),
        ]
        for stratum_label, class_code in strata:
            stratum_agb = np.where(landcover_grid == class_code, agb_grid, np.nan)
            stratum_ha  = np.sum(~np.isnan(stratum_agb)) * PIXEL_AREA_HA
            print(f"\n  ▸ {stratum_label}:")
            print(f"    Area           : {stratum_ha:,.0f} ha")
            if stratum_ha > 0:
                print(f"    Mean AGB       : {np.nanmean(stratum_agb):.2f} Mg/ha")
                print(f"    Total AGB      : {np.nansum(stratum_agb) * PIXEL_AREA_HA / 1e6:.4f} Tg")

    # ── Figure: GEDI AGB — split Andaman / Nicobar ───────────────────
    forest_gain_grid, _, _ = load_clipped_raster('ANI_GFW_Forest_Gain_clipped.tif')

    # Build a proper binary land mask so we can NaN-out ocean pixels.
    # Without this, the cream low-end of YlGn paints the ocean and
    # makes it indistinguishable from low-biomass land.
    land_binary = None
    for _mp in [PROC_DIR / 'ANI_GFW_DataMask_Land_Water_clipped.tif',
                PROC_DIR / 'ANI_land_mask.tif']:
        if _mp.exists():
            with rasterio.open(_mp) as _src:
                _raw = _src.read(1).astype(float)
                if _src.nodata is not None:
                    _raw[_raw == _src.nodata] = 0
                if _raw.shape != agb_grid.shape:
                    from skimage.transform import resize as sk_resize
                    _raw = sk_resize(_raw, agb_grid.shape,
                                     order=0, preserve_range=True,
                                     anti_aliasing=False)
                land_binary = (_raw > 0)
            break
    if land_binary is None:
        land_binary = ~np.isnan(agb_grid)

    agb_float = agb_grid.astype(float).copy()
    agb_float[~land_binary] = np.nan        # ocean → NaN → set_bad colour

    cmap_agb = plt.cm.YlGn.copy()
    cmap_agb.set_bad(color='#e6eef5')        # ocean = pale blue-grey

    fig = plt.figure(figsize=(18, 13), facecolor='white')
    gs_agb = gridspec.GridSpec(1, 3, figure=fig,
                               width_ratios=[2.0, 1.4, 0.06], wspace=0.05,
                               left=0.04, right=0.94, top=0.92, bottom=0.07)

    # Dilate single-pixel forest-gain dots so they're legible at this zoom.
    if forest_gain_grid is not None:
        h, w = agb_grid.shape
        g_full = (forest_gain_grid[:h, :w] == 1) & land_binary
        g_full = binary_dilation(g_full, iterations=2)
    else:
        g_full = None

    panel_axes = []
    for col, (region, row_sl) in enumerate([
            ('Andaman Islands',
             slice(0, min(ANDAMAN_ROW_END + _PAD, agb_float.shape[0]))),
            ('Nicobar Islands',
             slice(max(0, NICOBAR_ROW_START - _PAD), agb_float.shape[0]))
    ]):
        ax = fig.add_subplot(gs_agb[0, col])
        ax.set_facecolor('#e6eef5')

        # Crop AGB and the LAND mask to the same bbox so the coastline
        # is drawn from the true land/water boundary rather than the
        # (sparse) GEDI sampling pattern. Use agb_float (which now
        # already has ocean NaN'd via land_binary) as the ref so both
        # crops have identical extents.
        agb_region = agb_float[row_sl]
        crop       = _tight_crop(agb_region)
        land_crop  = _tight_crop(land_binary[row_sl].astype(float),
                                 ref=agb_region)
        ax.imshow(crop, cmap=cmap_agb, vmin=0, vmax=400,
                  interpolation='bilinear')

        # Bold coastline drawn from the actual land mask (not the AGB
        # validity mask). Two-pass stroke: a wider light halo behind a
        # darker line so the boundary stays visible on both ocean-side
        # and forest-side of the contour.
        ax.contour(land_crop, levels=[0.5], colors=['#ffffff'],
                   linewidths=2.4, alpha=0.95)
        ax.contour(land_crop, levels=[0.5], colors=['#0d1b2a'],
                   linewidths=1.4, alpha=1.0)

        # Forest gain overlay (vivid magenta — contrasts both green AGB
        # and pale ocean; the old gold blended into the cream background)
        if g_full is not None:
            g_rgn  = g_full[row_sl]
            g_crop = _tight_crop(g_rgn, ref=agb_float[row_sl])
            rgba   = np.zeros((*g_crop.shape, 4), dtype=float)
            rgba[g_crop, 0] = 0.85
            rgba[g_crop, 1] = 0.10
            rgba[g_crop, 2] = 0.55
            rgba[g_crop, 3] = 0.95
            ax.imshow(rgba, interpolation='nearest')

        # Per-panel headline: pixel-count + mean AGB
        region_agb_vals = agb_grid[row_sl][land_binary[row_sl]]
        region_agb_vals = region_agb_vals[~np.isnan(region_agb_vals)]
        if region_agb_vals.size:
            ax.text(
                0.02, 0.02,
                f'mean AGB: {region_agb_vals.mean():.0f} Mg/ha\n'
                f'P95: {np.percentile(region_agb_vals, 95):.0f} Mg/ha',
                transform=ax.transAxes, ha='left', va='bottom',
                fontsize=9, color='#1a1a2e',
                bbox=dict(boxstyle='round,pad=0.35', facecolor='white',
                          edgecolor='#cccccc', linewidth=0.7, alpha=0.92),
            )

        ax.set_title(region, color='#1a1a2e',
                     fontsize=13, fontweight='bold', pad=8)
        ax.axis('off')
        panel_axes.append(ax)

    # Shared colorbar — sits to the right of the Nicobar panel, not
    # wedged between the two maps.
    cax_agb = fig.add_subplot(gs_agb[0, 2])
    sm_agb  = plt.cm.ScalarMappable(cmap=cmap_agb, norm=plt.Normalize(0, 400))
    sm_agb.set_array([])
    cbar_agb = fig.colorbar(sm_agb, cax=cax_agb)
    cbar_agb.set_label('Aboveground Biomass (Mg/ha)', color='#1a1a2e',
                       fontsize=10)
    cbar_agb.ax.yaxis.set_tick_params(color='#1a1a2e', labelsize=8)
    plt.setp(cbar_agb.ax.yaxis.get_ticklabels(), color='#1a1a2e')

    legend_handles = [
        mpatches.Patch(facecolor='#e6eef5', edgecolor='#37474f',
                       linewidth=0.5, label='Ocean / Non-Land'),
    ]
    if g_full is not None:
        legend_handles.append(
            mpatches.Patch(facecolor='#d91a8c', alpha=0.95,
                           label='Forest Gain (GFW 2000–2012, dilated ×2)')
        )
    fig.legend(handles=legend_handles, loc='lower center',
               bbox_to_anchor=(0.5, 0.005), ncol=len(legend_handles),
               facecolor='white', edgecolor='#cccccc',
               labelcolor='#1a1a2e', fontsize=10)

    fig.suptitle('GEDI L4B Aboveground Biomass Baseline — '
                 'Andaman & Nicobar Islands',
                 color='#1a1a2e', fontsize=14, fontweight='bold', y=0.97)
    (FIG_DIR / 'carbon').mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_DIR / 'carbon' / 'agb_gedi_baseline_map.png',
                dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  ✅  Figure saved → figures/carbon/agb_gedi_baseline_map.png")
    return agb_grid, landcover_grid



# ══════════════════════════════════════════════════════════════════════
# PART 2 — Cross-Validation: GEDI vs. Saatchi Independent AGB
# ══════════════════════════════════════════════════════════════════════
def cross_validate_agb_datasets(agb_grid):
    """Cross-validate GEDI AGB against Saatchi independent AGB product."""
    print(f"\n{SEP}")
    print("PART 2 — Cross-Validation: GEDI vs. Saatchi AGB")
    print(SEP)

    saatchi_grid, _, _ = load_clipped_raster('ANI_Saatchi_AGB_CrossValidation_clipped.tif')
    if saatchi_grid is None or agb_grid is None:
        print("  ⚠️  Skipping cross-validation (missing layer).")
        return None, None

    # Remove background/ocean artefacts (Saatchi background < 1 Mg/ha)
    saatchi_grid = np.where(saatchi_grid <   1.0, np.nan, saatchi_grid)
    saatchi_grid = np.where(saatchi_grid > 800.0, np.nan, saatchi_grid)

    # Both grids are now at 30m — align shapes to smallest common extent
    h = min(agb_grid.shape[0], saatchi_grid.shape[0])
    w = min(agb_grid.shape[1], saatchi_grid.shape[1])
    gedi_crop    = agb_grid[:h, :w]
    saatchi_crop = saatchi_grid[:h, :w]

    overlap_mask = ~np.isnan(gedi_crop) & ~np.isnan(saatchi_crop)

    if overlap_mask.sum() < 5:
        print("  ⚠️  Insufficient overlapping pixels for cross-validation.")
        return None, None

    gedi_vals    = gedi_crop[overlap_mask]
    saatchi_vals = saatchi_crop[overlap_mask]

    pearson_r, p_value = stats.pearsonr(gedi_vals, saatchi_vals)
    rmse_val = np.sqrt(np.mean((gedi_vals - saatchi_vals) ** 2))
    bias_val = np.mean(gedi_vals - saatchi_vals)
    mae_val  = np.mean(np.abs(gedi_vals - saatchi_vals))

    print(f"  Compared pixels    : {overlap_mask.sum():,}")
    print(f"  Pearson r          : {pearson_r:.4f}  (p = {p_value:.2e})")
    print(f"  RMSE               : {rmse_val:.2f} Mg/ha")
    print(f"  Bias (GEDI−Saatchi): {bias_val:+.2f} Mg/ha  "
          f"({'GEDI underestimates' if bias_val < 0 else 'GEDI overestimates'})")
    print(f"  MAE                : {mae_val:.2f} Mg/ha")

    # Figure rendering is owned by src/validation_stats.py (canonical
    # forest-only cross-validation panel: scatter + AGB histograms +
    # CCC/RMA/CI table). This function now only returns the metrics that
    # downstream carbon-loss accounting needs.
    return rmse_val, bias_val


# ══════════════════════════════════════════════════════════════════════
# PART 3 — Annual Carbon & CO₂e Loss (GFW lossyear × GEDI AGB)
# ══════════════════════════════════════════════════════════════════════
def compute_annual_carbon_loss(agb_grid, rmse_val, bias_val):
    """Calculate per-year carbon and CO₂e loss using GFW loss-year layer."""
    print(f"\n{SEP}")
    print("PART 3 — Annual Carbon & CO₂e Loss (2001–2024)")
    print(SEP)

    import pandas as pd

    deforestation_yr, _, _ = load_clipped_raster('ANI_GFW_Forest_Loss_2001_2023_clipped.tif')
    forest_gain_grid, _, _ = load_clipped_raster('ANI_GFW_Forest_Gain_clipped.tif')

    if deforestation_yr is None or agb_grid is None:
        print("  ❌  Missing GFW or GEDI layer — aborting.")
        return None

    # GFW nodata → 0 (no loss year recorded)
    defor_yr_int = np.nan_to_num(deforestation_yr, nan=0).astype(int)

    year_records = []

    for yr_code in range(1, 25):          # code 1 = 2001 … 24 = 2024
        calendar_yr  = 2000 + yr_code
        loss_px_mask = defor_yr_int == yr_code

        if loss_px_mask.sum() == 0:
            year_records.append({
                'year': calendar_yr, 'area_ha': 0.0,
                'total_biomass_mg': 0.0, 'carbon_mgc': 0.0,
                'co2e_mgco2e': 0.0,
            })
            continue

        agb_in_loss  = agb_grid[loss_px_mask]
        agb_in_loss  = agb_in_loss[~np.isnan(agb_in_loss)]
        area_ha      = loss_px_mask.sum() * PIXEL_AREA_HA
        agb_mg       = agb_in_loss.sum() * PIXEL_AREA_HA
        total_bio_mg = agb_mg * (1 + ROOT_SHOOT_R)   # AGB + BGB
        carbon_mgc   = total_bio_mg * IPCC_CARBON_F
        co2e_mgco2e  = carbon_mgc   * CO2E_MW_RATIO

        year_records.append({
            'year':           calendar_yr,
            'area_ha':        round(area_ha,      2),
            'total_biomass_mg': round(total_bio_mg, 2),
            'carbon_mgc':     round(carbon_mgc,   2),
            'co2e_mgco2e':    round(co2e_mgco2e,  2),
        })

    annual_df = pd.DataFrame(year_records)
    annual_df['carbon_ggc']   = annual_df['carbon_mgc']  / 1e3
    annual_df['co2e_ggco2e']  = annual_df['co2e_mgco2e'] / 1e3

    print(f"\n{'Year':>5} {'Area_ha':>10} {'Biomass_Mg':>12} {'C_GgC':>10} {'CO2e_Gg':>10}")
    print('-' * 55)
    for _, row in annual_df.iterrows():
        if row['area_ha'] > 0:
            print(f"{int(row['year']):>5} {row['area_ha']:>10,.1f} "
                  f"{row['total_biomass_mg']:>12,.1f} {row['carbon_ggc']:>10.4f} "
                  f"{row['co2e_ggco2e']:>10.4f}")

    gross_area_ha  = annual_df['area_ha'].sum()
    gross_carbon   = annual_df['carbon_ggc'].sum()
    gross_co2e     = annual_df['co2e_ggco2e'].sum()
    print(f"\n  GROSS DEFORESTATION (2001–2024):")
    print(f"    Total area deforested : {gross_area_ha:,.0f} ha")
    print(f"    Total carbon lost     : {gross_carbon:.4f} GgC")
    print(f"    Total CO₂e emitted    : {gross_co2e:.4f} GgCO₂e")

    # ── Reforestation / Afforestation net balance ──────────────────────
    gain_stats = None
    gain_px_mask = None
    if forest_gain_grid is not None:
        gain_px_mask = (forest_gain_grid == 1)
        gain_agb     = agb_grid[gain_px_mask]
        gain_agb     = gain_agb[~np.isnan(gain_agb)]

        gain_area_ha = gain_px_mask.sum() * PIXEL_AREA_HA
        gain_agb_mg  = gain_agb.sum() * PIXEL_AREA_HA
        gain_bio_mg  = gain_agb_mg * (1 + ROOT_SHOOT_R)
        gain_carbon  = (gain_bio_mg * IPCC_CARBON_F) / 1e3   # GgC
        gain_co2e    = (gain_bio_mg * IPCC_CARBON_F * CO2E_MW_RATIO) / 1e3

        gain_stats = dict(
            area_ha     = gain_area_ha,
            carbon_ggc  = gain_carbon,
            co2e_ggco2e = gain_co2e,
        )

        print(f"\n  REFORESTATION & AFFORESTATION (GFW Gain Layer):")
        print(f"    Total area gained     : {gain_area_ha:,.0f} ha")
        print(f"    Carbon sequestered    : {gain_carbon:.4f} GgC")
        print(f"    CO₂e sequestered      : {gain_co2e:.4f} GgCO₂e")
        print(f"\n  NET BALANCE:")
        print(f"    Net Area Loss         : {gross_area_ha - gain_area_ha:,.0f} ha")
        print(f"    Net Carbon Flux       : {gross_carbon - gain_carbon:.4f} GgC")
        print(f"    Net Climate Impact    : {gross_co2e - gain_co2e:.4f} GgCO₂e")

    # ── Uncertainty estimate from cross-validation RMSE ───────────────
    if rmse_val is not None and np.nanmean(agb_grid) > 0:
        uncertainty_frac = rmse_val / np.nanmean(agb_grid)
        print(f"\n  Uncertainty (RMSE/mean AGB) = ±{uncertainty_frac*100:.1f}%")
        print(f"  Gross Carbon range: "
              f"{gross_carbon*(1-uncertainty_frac):.4f} – "
              f"{gross_carbon*(1+uncertainty_frac):.4f} GgC")

    # ── Append gain summary to DataFrame — cumulative total only ──────
    if gain_stats is not None:
        annual_df['cumulative_gain_co2e_ggco2e'] = gain_stats['co2e_ggco2e']
        annual_df['cumulative_gain_area_ha']     = gain_stats['area_ha']
        gross_co2e_total = annual_df['co2e_ggco2e'].sum()
        net_total        = gross_co2e_total - gain_stats['co2e_ggco2e']
        annual_df['net_co2e_total_ggco2e']       = round(net_total, 6)

    # ── Save result CSV ────────────────────────────────────────────────
    csv_out = RES_DIR / 'carbon_annual_loss_by_year.csv'
    annual_df.to_csv(csv_out, index=False)
    print(f"\n  ✅  CSV saved → {csv_out}")
    return annual_df, defor_yr_int, gain_stats, gain_px_mask


# ══════════════════════════════════════════════════════════════════════
# PART 4 — Publication-Quality Figures
# ══════════════════════════════════════════════════════════════════════
def render_carbon_figures(annual_df, agb_grid, defor_yr_int,
                          gain_stats=None, gain_px_mask=None):
    """Render annual time-series and spatial hotspot figures."""
    print(f"\n{SEP}")
    print("PART 4 — Generating Figures")
    print(SEP)

    if annual_df is None or agb_grid is None:
        print("  ⚠️  Skipping figures (missing data).")
        return

    # ── 4a: Annual + cumulative CO₂e + forest-gain reference band ───────
    # GFW Forest-Gain is a single period layer (2000–2012) with no annual
    # resolution. The previous render faked an annual gain bar by spreading
    # the total evenly across 24 years, which visually implied an ongoing
    # annual offset that doesn't exist. We now show gain as a period-band
    # at the top of the panel with a single total-budget annotation.
    df_active   = annual_df[annual_df['area_ha'] > 0].copy().sort_values('year')
    years       = df_active['year'].to_numpy()
    losses      = df_active['co2e_ggco2e'].to_numpy()
    cumulative  = losses.cumsum()

    fig, ax_bar = plt.subplots(figsize=(14, 5.8), facecolor='white')
    ax_bar.set_facecolor('#f5f7fa')

    bar_w = 0.72
    bars  = ax_bar.bar(years, losses, width=bar_w,
                       color='#d32f2f', alpha=0.88, edgecolor='white',
                       linewidth=0.4,
                       label='Annual CO₂e Loss (deforestation)', zorder=3)

    # Mean-annual reference (dashed) so individual years can be read
    # against the long-run rate.
    mean_loss = losses.mean()
    ax_bar.axhline(mean_loss, color='#555555', linestyle='--',
                   linewidth=1.0, alpha=0.7, zorder=2,
                   label=f'Mean annual loss ({mean_loss:.0f} Gg CO₂e)')

    # Period band: GFW Forest Gain 2000–2012, total budget annotation.
    if gain_stats is not None:
        ax_bar.axvspan(2000.5, 2012.5, color='#388e3c', alpha=0.10,
                       zorder=1)
        gain_total = gain_stats['co2e_ggco2e']
        gain_mean  = gain_total / 12.0   # period spans 12 years
        ax_bar.annotate(
            f'GFW Forest Gain period\n'
            f'2000–2012 total: {gain_total:.0f} Gg CO₂e\n'
            f'(period-mean ≈ {gain_mean:.0f} Gg/yr)',
            xy=(2006.5, max(losses) * 0.94),
            ha='center', va='top', fontsize=9, color='#1b5e20',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='#e8f5e9',
                      edgecolor='#388e3c', linewidth=0.8, alpha=0.95),
            zorder=6,
        )

    # Peak-year annotations (top 3 loss years) — these are the story.
    top_idx = np.argsort(losses)[-3:]
    for i in top_idx:
        ax_bar.annotate(
            f'{years[i]}\n{losses[i]:.0f} Gg',
            xy=(years[i], losses[i]),
            xytext=(0, 8), textcoords='offset points',
            ha='center', va='bottom', fontsize=8.5,
            color='#7f1d1d', fontweight='bold', zorder=7,
        )

    ax_bar.set_xlabel('Year', color='#1a1a2e', fontsize=11)
    ax_bar.set_ylabel('Annual CO₂e Loss (Gg CO₂e)',
                      color='#1a1a2e', fontsize=11)
    ax_bar.tick_params(colors='#1a1a2e', axis='both')
    ax_bar.spines[:].set_color('#cccccc')
    ax_bar.set_xlim(years.min() - 0.8, years.max() + 0.8)
    ax_bar.set_ylim(0, max(losses) * 1.18)
    ax_bar.grid(True, axis='y', color='#dddddd', linewidth=0.6, alpha=0.7)
    ax_bar.set_axisbelow(True)

    # Cumulative line on a secondary axis — muted so the bars remain
    # the visual focus.
    ax_cum = ax_bar.twinx()
    ax_cum.plot(years, cumulative,
                color='#1565c0', lw=1.8, marker='o', ms=3.5,
                alpha=0.85, label='Cumulative loss', zorder=5)
    ax_cum.fill_between(years, 0, cumulative,
                        color='#1565c0', alpha=0.06, zorder=2)
    ax_cum.set_ylabel('Cumulative CO₂e Loss (Gg CO₂e)',
                      color='#1565c0', fontsize=11)
    ax_cum.tick_params(colors='#1565c0', axis='y')
    ax_cum.spines[:].set_color('#cccccc')
    ax_cum.set_ylim(0, cumulative.max() * 1.08)

    # Endpoint label on the cumulative line — placed above the last
    # marker (offset in points) so it can't collide with the right-axis
    # tick labels.
    ax_cum.annotate(
        f'{cumulative[-1]:,.0f} Gg',
        xy=(years[-1], cumulative[-1]),
        xytext=(-2, 10), textcoords='offset points',
        ha='right', va='bottom',
        color='#1565c0', fontsize=9, fontweight='bold',
    )

    h1, l1 = ax_bar.get_legend_handles_labels()
    h2, l2 = ax_cum.get_legend_handles_labels()
    ax_bar.legend(h1 + h2, l1 + l2,
                  facecolor='white', edgecolor='#cccccc',
                  labelcolor='#1a1a2e', fontsize=9,
                  loc='upper left', bbox_to_anchor=(0.30, 0.99),
                  ncol=3, framealpha=0.95)

    fig.suptitle('Annual CO₂e Loss (2001–2024) & GFW Forest-Gain Period (2000–2012)\n'
                 'Andaman & Nicobar Islands',
                 color='#1a1a2e', fontsize=13, fontweight='bold', y=1.02)
    fig.tight_layout()
    fig.savefig(FIG_DIR / 'carbon' / 'carbon_annual_loss_timeseries.png',
                dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print("  ✅  Figure saved → figures/carbon/carbon_annual_loss_timeseries.png")


    # ── 4b: Net Carbon Balance figure (cumulative loss / gain / net) ──
    # Single-panel design: Net Flux is the headline. The old left-panel
    # annual-bar view duplicated carbon_annual_loss_timeseries.png, so
    # we dropped it and let the cumulative balance be the focus.
    if gain_stats is not None:
        df_active  = annual_df[annual_df['area_ha'] > 0].copy().sort_values('year')
        total_loss = df_active['co2e_ggco2e'].sum()
        total_gain = gain_stats['co2e_ggco2e']
        net_flux   = total_loss - total_gain
        x2         = df_active['year'].to_numpy()

        # GFW gain data: 2000-2012 cumulative only — no annual breakdown.
        GAIN_START, GAIN_END = 2000, 2012
        gain_n_yrs  = GAIN_END - GAIN_START + 1   # 13 years

        # Cumulative series
        cum_loss = df_active['co2e_ggco2e'].cumsum().values
        cum_gain = np.array([
            total_gain * max(0, min(1.0, (yr - GAIN_START + 1) / gain_n_yrs))
            for yr in x2
        ])
        cum_net  = cum_loss - cum_gain

        fig2, ax_nb = plt.subplots(figsize=(13, 6.2), facecolor='white')
        ax_nb.set_facecolor('#f5f8fc')
        fig2.subplots_adjust(left=0.07, right=0.95, top=0.86, bottom=0.20)

        # Cumulative Loss — pale red ramp + thin outline (context only).
        ax_nb.fill_between(x2, cum_loss, alpha=0.10, color='#ef5350',
                           zorder=1)
        ax_nb.plot(x2, cum_loss, color='#c62828', lw=1.8, alpha=0.75,
                   label='Cumulative Loss (GFW 2001–2024)', zorder=3)

        # Cumulative Gain — green dashed line + faint fill.
        ax_nb.fill_between(x2, cum_gain, alpha=0.10, color='#2e7d32',
                           zorder=1)
        ax_nb.plot(x2, cum_gain, color='#2e7d32', lw=1.8, alpha=0.85,
                   linestyle='--',
                   label='Cumulative Gain (GFW 2000–2012, plateau after)',
                   zorder=3)

        # Net Flux — the headline. Thick, dark blue, on top of everything.
        ax_nb.plot(x2, cum_net, color='#0d47a1', lw=3.0, marker='o',
                   ms=4.5, label='Net Flux (Loss − Gain)', zorder=6)
        ax_nb.axhline(0, color='#888', lw=0.8, linestyle='--', alpha=0.5,
                      zorder=2)

        # End-of-gain-period marker.
        ax_nb.axvline(GAIN_END, color='#2e7d32', lw=1.0, linestyle=':',
                      alpha=0.7, zorder=2)
        ax_nb.annotate(
            'Gain period ends',
            xy=(GAIN_END, cum_gain.max()),
            xytext=(GAIN_END - 0.3, cum_gain.max() * 1.05),
            ha='right', va='bottom', fontsize=9, color='#2e7d32',
            fontweight='bold',
            arrowprops=dict(arrowstyle='->', color='#2e7d32', lw=0.8,
                            alpha=0.7),
            zorder=7,
        )

        # Endpoint labels for all three series.
        last_x = x2[-1]
        ax_nb.annotate(f'{cum_loss[-1]:,.0f} Gg',
                       xy=(last_x, cum_loss[-1]),
                       xytext=(6, 0), textcoords='offset points',
                       ha='left', va='center', fontsize=9,
                       color='#c62828', fontweight='bold')
        ax_nb.annotate(f'{cum_gain[-1]:,.0f} Gg',
                       xy=(last_x, cum_gain[-1]),
                       xytext=(6, 0), textcoords='offset points',
                       ha='left', va='center', fontsize=9,
                       color='#2e7d32', fontweight='bold')

        # Net Flux headline callout pointing at the endpoint.
        ax_nb.annotate(
            f'Net release: {net_flux:,.0f} Gg CO₂e',
            xy=(last_x, cum_net[-1]),
            xytext=(-110, -55), textcoords='offset points',
            ha='left', va='center', fontsize=11, fontweight='bold',
            color='#0d47a1',
            bbox=dict(boxstyle='round,pad=0.45', facecolor='#e3f2fd',
                      edgecolor='#0d47a1', linewidth=1.0),
            arrowprops=dict(arrowstyle='->', color='#0d47a1', lw=1.2),
            zorder=10,
        )

        ax_nb.set_xlabel('Year', color='#1a1a2e', fontsize=11)
        ax_nb.set_ylabel('Cumulative CO₂e Flux (Gg CO₂e)',
                         color='#1a1a2e', fontsize=11)
        ax_nb.tick_params(colors='#1a1a2e')
        ax_nb.spines[:].set_color('#cccccc')
        ax_nb.grid(True, color='#dddddd', linewidth=0.6, alpha=0.7)
        ax_nb.set_axisbelow(True)
        ax_nb.set_xlim(x2.min() - 0.5, x2.max() + 2.0)

        leg = ax_nb.legend(facecolor='white', edgecolor='#cccccc',
                           labelcolor='#1a1a2e', fontsize=9.5,
                           loc='upper left')
        leg.set_zorder(11)

        # Single figure-level caveat footnote.
        fig2.text(
            0.5, 0.03,
            '★ GFW Forest-Gain data covers 2000–2012 only (no annual '
            'breakdown). The cumulative-gain curve rises linearly across '
            'the period and plateaus thereafter.',
            ha='center', va='bottom', fontsize=8.5, color='#4a5568',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='#fffde7',
                      edgecolor='#f0c000', alpha=0.9),
        )

        fig2.suptitle(
            'Cumulative Carbon Balance (2001–2024) — '
            'Andaman & Nicobar Islands',
            color='#1a1a2e', fontsize=14, fontweight='bold', y=0.96,
        )

        fig2.savefig(FIG_DIR / 'carbon' / 'carbon_net_balance_figure.png',
                     dpi=200, bbox_inches='tight', facecolor='white')
        plt.close()
        print("  ✅  Figure saved → figures/carbon/carbon_net_balance_figure.png")


    # ── 4c: Spatial carbon loss hotspot map — split Andaman / Nicobar ──
    if defor_yr_int is not None:
        loss_pixels = (defor_yr_int >= 1) & (defor_yr_int <= 24)

        # ── Load proper binary land mask (GFW DataMask) ────────────────
        # Using a true land/water raster prevents loss dots appearing in ocean.
        # The GFW DataMask is a full-coverage binary: 1=land, 0=water.
        _datamask_path = PROC_DIR / 'ANI_GFW_DataMask_Land_Water_clipped.tif'
        _landmask_path = PROC_DIR / 'ANI_land_mask.tif'
        land_binary = None
        for _mp in [_datamask_path, _landmask_path]:
            if _mp.exists():
                with rasterio.open(_mp) as _src:
                    _raw = _src.read(1).astype(float)
                    _raw[_raw == _src.nodata] = 0 if _src.nodata is not None else 0
                    # Resize to match loss_pixels shape if needed
                    if _raw.shape != loss_pixels.shape:
                        from skimage.transform import resize as sk_resize
                        _raw = sk_resize(_raw, loss_pixels.shape,
                                         order=0, preserve_range=True, anti_aliasing=False)
                    land_binary = (_raw > 0)
                print(f"  ℹ️   Land mask loaded: {_mp.name}  shape={land_binary.shape}")
                break
        if land_binary is None:
            # Fallback: derive from agb_grid (sparse but better than nothing)
            land_binary = ~np.isnan(agb_grid)
            print("  ⚠️   Using GEDI AGB NaN pattern as land mask (sparse).")

        # ── Mask loss pixels to LAND ONLY before block-averaging ───────
        # This ensures coastal 500m blocks only count land-area loss,
        # eliminating the bleed of hotspot colour into the ocean.
        loss_land = loss_pixels & land_binary

        # ── Hotspot intensity: forest loss density per ~500 m cell ─────
        BLKSZ = 17          # 17 × 30 m ≈ 510 m
        h0, w0 = loss_land.shape
        ht = h0 - h0 % BLKSZ;  wt = w0 - w0 % BLKSZ
        lp_trim = loss_land[:ht, :wt].astype(float)
        lb_trim = land_binary[:ht, :wt].astype(float)
        nb_h = ht // BLKSZ;  nb_w = wt // BLKSZ
        # Sum loss and land pixels per block separately, then divide
        # → loss % relative to land area within each block (not whole block area)
        loss_sum = lp_trim.reshape(nb_h, BLKSZ, nb_w, BLKSZ).sum(axis=(1, 3))
        land_sum = lb_trim.reshape(nb_h, BLKSZ, nb_w, BLKSZ).sum(axis=(1, 3))
        land_sum = np.where(land_sum == 0, np.nan, land_sum)   # avoid /0
        loss_density = (loss_sum / land_sum) * 100             # % of land lost

        # Upsample back to pixel grid
        loss_density_full = np.repeat(np.repeat(loss_density, BLKSZ, axis=0),
                                       BLKSZ, axis=1)
        pad_h = h0 - loss_density_full.shape[0]
        pad_w = w0 - loss_density_full.shape[1]
        loss_density_full = np.pad(loss_density_full, ((0, pad_h), (0, pad_w)),
                                    constant_values=np.nan)

        # Mask to land pixels using the SAME binary land mask (not GEDI AGB).
        # Both the loss density and the contour now share identical land coverage.
        ld_masked = np.where(land_binary, loss_density_full, np.nan)

        # ── Colormaps ──────────────────────────────────────────────────
        from matplotlib.colors import LinearSegmentedColormap
        cmap_hot = LinearSegmentedColormap.from_list(
            'carbon_hot', ['#ffffcc', '#fd8d3c', '#e31a1c', '#67000d'], N=256)
        cmap_hot.set_bad(alpha=0)

        # Stronger green ramp so the land base actually reads against a
        # pale ocean. Previous c8e6c9→388e3c was too washed out at low AGB.
        cmap_land_base = LinearSegmentedColormap.from_list(
            'land_base', ['#a5c8a9', '#1b5e20'], N=256)
        cmap_land_base.set_bad(alpha=0)

        # For land base shading use AGB (NaN where GEDI missing → fade there)
        agb_norm = np.where(land_binary, agb_grid / np.nanmax(agb_grid), np.nan)

        vmax_ld = float(np.nanmax(ld_masked))

        fig = plt.figure(figsize=(18, 13), facecolor='white')
        gs_cl = gridspec.GridSpec(1, 3, figure=fig,
                                   width_ratios=[2.0, 1.4, 0.06], wspace=0.05,
                                   left=0.04, right=0.94, top=0.92, bottom=0.07)

        for col, (region, row_sl) in enumerate([
                ('Andaman Islands',
                 slice(0, min(ANDAMAN_ROW_END + _PAD, h0))),
                ('Nicobar Islands',
                 slice(max(0, NICOBAR_ROW_START - _PAD), h0))
        ]):
            ax = fig.add_subplot(gs_cl[0, col])
            ax.set_facecolor('#e6eef5')         # pale blue-grey ocean

            lb_rgn   = land_binary[row_sl].astype(float)
            agb_rgn  = agb_norm[row_sl]
            ld_rgn   = ld_masked[row_sl]

            # Use land_binary as the crop reference (full coverage, not sparse GEDI)
            lb_crop  = _tight_crop(lb_rgn,  ref=lb_rgn)
            agb_crop = _tight_crop(agb_rgn, ref=lb_rgn)
            ld_crop  = _tight_crop(ld_rgn,  ref=lb_rgn)

            # Turn 0 (ocean) in lb_crop to NaN for colormaps
            agb_for_base = np.where(lb_crop > 0, agb_crop, np.nan)
            ld_for_hot   = np.where(lb_crop > 0, ld_crop,  np.nan)

            # 1) Flat pale-green land underlay (so islands without GEDI
            #    sampling still appear as "land", not as ocean).
            land_only = np.where(lb_crop > 0, 1.0, np.nan)
            ax.imshow(land_only,
                      cmap=LinearSegmentedColormap.from_list('flatland',
                                                              ['#cfe3d0',
                                                               '#cfe3d0']),
                      vmin=0, vmax=1, interpolation='nearest')

            # 2) Green AGB-shaded base on top of the flat underlay
            ax.imshow(agb_for_base, cmap=cmap_land_base, vmin=0, vmax=1,
                      interpolation='bilinear')

            # 3) Loss-density overlay — strictly inside land boundary
            ax.imshow(ld_for_hot, cmap=cmap_hot, vmin=0, vmax=vmax_ld,
                      alpha=0.88, interpolation='nearest')

            # 4) Bold coastline (two-pass: white halo + dark line) so
            #    every island reads clearly against ocean and forest.
            ax.contour(lb_crop, levels=[0.5], colors=['#ffffff'],
                       linewidths=2.4, alpha=0.95)
            ax.contour(lb_crop, levels=[0.5], colors=['#0d1b2a'],
                       linewidths=1.4, alpha=1.0)

            # Region label + loss-area summary, rolled into a single
            # title so it isn't floating in aspect-ratio padding.
            r_area = loss_land[row_sl].sum() * PIXEL_AREA_HA
            ax.set_title(
                f'{region}\nForest loss area: {r_area:,.0f} ha',
                color='#1a1a2e', fontsize=13, fontweight='bold', pad=10,
            )
            ax.axis('off')

        cax_cl  = fig.add_subplot(gs_cl[0, 2])
        sm_cl   = plt.cm.ScalarMappable(cmap=cmap_hot,
                                         norm=plt.Normalize(0, vmax_ld))
        sm_cl.set_array([])
        cbar_cl = fig.colorbar(sm_cl, cax=cax_cl)
        cbar_cl.set_label('Forest Loss Density (% land area lost per ~500 m cell)',
                          color='#1a1a2e', fontsize=10)
        cbar_cl.ax.yaxis.set_tick_params(color='#1a1a2e', labelsize=8)
        plt.setp(cbar_cl.ax.yaxis.get_ticklabels(), color='#1a1a2e')

        # Legend explaining the non-hotspot layers
        legend_handles = [
            mpatches.Patch(facecolor='#e6eef5', edgecolor='#0d1b2a',
                           linewidth=0.8, label='Ocean / Non-Land'),
            mpatches.Patch(facecolor='#1b5e20',
                           label='Intact forest (darker = higher AGB)'),
            mpatches.Patch(facecolor='#cfe3d0',
                           label='Land without GEDI sampling'),
        ]
        fig.legend(handles=legend_handles, loc='lower center',
                   bbox_to_anchor=(0.5, 0.005), ncol=3,
                   facecolor='white', edgecolor='#cccccc',
                   labelcolor='#1a1a2e', fontsize=9.5)

        fig.suptitle('Carbon Loss Hotspots (2001–2024) — '
                     'Andaman & Nicobar Islands',
                     color='#1a1a2e', fontsize=14, fontweight='bold', y=0.97)
        (FIG_DIR / 'carbon').mkdir(parents=True, exist_ok=True)
        fig.savefig(FIG_DIR / 'carbon' / 'carbon_loss_hotspots_map.png',
                    dpi=200, bbox_inches='tight', facecolor='white')
        plt.close()
        print("  ✅  Figure saved → figures/carbon/carbon_loss_hotspots_map.png")



        if gain_px_mask is not None:
            h2, w2 = loss_pixels.shape
            gain_aligned = gain_px_mask[:h2, :w2]

            # ── Reuse land_binary from the hotspot section (GFW DataMask) ──
            # Ensures loss/gain pixels and the boundary contour all use the
            # same true coastline — no coloured dots outside the boundary.
            lb2 = land_binary[:h2, :w2]           # binary: True = land

            # True per-pixel loss / gain masks (no dilation) — used for
            # the area-in-hectares tallies in the titles / legend.
            loss_on_land = loss_pixels[:h2, :w2] & lb2
            gain_on_land = gain_aligned & lb2

            # Dilated versions for *rendering* only. Single 30 m pixels
            # disappear at full-island zoom; dilating ×2 makes them
            # legible without changing the underlying statistics.
            loss_vis = binary_dilation(loss_on_land, iterations=2)
            gain_vis = binary_dilation(gain_on_land, iterations=2)

            # ── Build RGBA composite ────────────────────────────────────────
            # Layer order: ocean (transparent) → land base → loss → gain
            composite = np.zeros((h2, w2, 4), dtype=float)

            # 1) Land base — pale sage so islands read as land even where
            #    there's no loss/gain pixel.
            LAND_R, LAND_G, LAND_B = 0.87, 0.93, 0.80
            composite[lb2, 0] = LAND_R
            composite[lb2, 1] = LAND_G
            composite[lb2, 2] = LAND_B
            composite[lb2, 3] = 1.00

            # 2) Loss pixels → vivid crimson  (strictly inside land_binary)
            composite[loss_vis, 0] = 0.92
            composite[loss_vis, 1] = 0.10
            composite[loss_vis, 2] = 0.10
            composite[loss_vis, 3] = 1.00

            # 3) Gain pixels → bright gold  (wins over loss where overlap)
            composite[gain_vis, 0] = 1.00
            composite[gain_vis, 1] = 0.84
            composite[gain_vis, 2] = 0.00
            composite[gain_vis, 3] = 1.00

            fig = plt.figure(figsize=(18, 13), facecolor='white')
            gs_lg = gridspec.GridSpec(1, 2, figure=fig,
                                       width_ratios=[2.0, 1.4], wspace=0.06,
                                       left=0.04, right=0.97, top=0.92,
                                       bottom=0.08)

            for col, (region_label, grp_rows) in enumerate([
                    ('Andaman Islands', slice(0, min(ANDAMAN_ROW_END + _PAD, h2))),
                    ('Nicobar Islands', slice(max(0, NICOBAR_ROW_START - _PAD), h2))
            ]):
                ax_p = fig.add_subplot(gs_lg[0, col])
                ax_p.set_facecolor('#e6eef5')    # pale blue-grey ocean

                lb_rgn   = lb2[grp_rows].astype(float)
                comp_rgn = composite[grp_rows]

                lb_crop_lg   = _tight_crop(lb_rgn,  ref=lb_rgn)
                comp_crop_lg = _tight_crop(comp_rgn, ref=lb_rgn)

                ax_p.imshow(comp_crop_lg, interpolation='nearest')

                # Bold coastline — white halo + dark stroke. Drawn from
                # the binary land mask so every island shows up.
                ax_p.contour(lb_crop_lg, levels=[0.5], colors=['#ffffff'],
                             linewidths=2.4, alpha=0.95)
                ax_p.contour(lb_crop_lg, levels=[0.5], colors=['#0d1b2a'],
                             linewidths=1.4, alpha=1.0)

                # Per-panel stats rolled into the title so they don't
                # float in aspect-ratio padding.
                r_loss = loss_on_land[grp_rows].sum() * PIXEL_AREA_HA
                r_gain = gain_on_land[grp_rows].sum() * PIXEL_AREA_HA
                ax_p.set_title(
                    f'{region_label}\n'
                    f'Loss: {r_loss:,.0f} ha   •   Gain: {r_gain:,.0f} ha',
                    color='#1a1a2e', fontsize=13, fontweight='bold', pad=10,
                )
                ax_p.axis('off')

            from matplotlib.patches import Patch
            fig.legend(handles=[
                Patch(facecolor='#EB1A1A', edgecolor='#666',
                      label=f'Forest Loss 2001–2024 — '
                            f'{loss_on_land.sum()*PIXEL_AREA_HA:,.0f} ha '
                            f'(dilated ×2 for visibility)'),
                Patch(facecolor='#FFD600', edgecolor='#666',
                      label=f'Forest Gain 2000–2012 — '
                            f'{gain_on_land.sum()*PIXEL_AREA_HA:,.0f} ha '
                            f'(dilated ×2 for visibility)'),
                Patch(facecolor='#DDECC5', edgecolor='#33691e',
                      label='Intact Land / No Change'),
                Patch(facecolor='#e6eef5', edgecolor='#0d1b2a',
                      linewidth=0.8, label='Ocean / Non-Land'),
            ], loc='lower center', bbox_to_anchor=(0.5, 0.005), ncol=4,
               facecolor='white', edgecolor='#cccccc',
               labelcolor='#1a1a2e', fontsize=9.5)
            fig.suptitle('Forest Loss vs Forest Gain — '
                         'Andaman & Nicobar Islands',
                         color='#1a1a2e', fontsize=14, fontweight='bold',
                         y=0.97)
            fig.savefig(FIG_DIR / 'carbon' / 'carbon_loss_vs_gain_spatial.png',
                        dpi=200, bbox_inches='tight', facecolor='white')
            plt.close()
            print('  ✅  Figure saved → figures/carbon/carbon_loss_vs_gain_spatial.png')



# ══════════════════════════════════════════════════════════════════════
# ENTRY POINT
# ══════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    print(SEP)
    print("  ANI Ecosystem Services — Week 6: Carbon Analysis")
    print(f"  Input : {PROC_DIR}")
    print(f"  Output: {RES_DIR}  |  {FIG_DIR}")
    print(SEP)

    agb_grid, landcover_grid  = estimate_gedi_agb_baseline()
    rmse_val, bias_val        = cross_validate_agb_datasets(agb_grid)
    result_tuple              = compute_annual_carbon_loss(agb_grid, rmse_val, bias_val)

    if result_tuple is not None and len(result_tuple) == 4:
        annual_df, defor_yr_int, gain_stats, gain_px_mask = result_tuple
    else:
        annual_df, defor_yr_int, gain_stats, gain_px_mask = None, None, None, None

    render_carbon_figures(annual_df, agb_grid, defor_yr_int,
                          gain_stats=gain_stats, gain_px_mask=gain_px_mask)

    print(f"\n{SEP}")
    print("  ✅  Week 6 Carbon Analysis Complete!")
    print(f"      Results : {RES_DIR}/carbon_annual_loss_by_year.csv")
    print(f"      Figures : {FIG_DIR}/carbon/")
    print(f"  Next → run: venv/bin/python src/habitat_quality.py")
    print(SEP)
