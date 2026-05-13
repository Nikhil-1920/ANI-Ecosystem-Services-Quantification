"""
ANI Ecosystem Services — Multifunctional Hotspot Synthesis (ECI)
================================================================
Fuses Carbon, Habitat Quality, and Soil Retention loss layers into a
Unified Ecosystem Collapse Index (ECI) to identify pixels where all
three ecosystem services are simultaneously degraded.

    ECI = (habitat_harm_norm + soil_harm_norm + deforestation_flag) / 3

Inputs  : data/processed/ANI_GFW_Forest_Loss_2001_2023_clipped.tif
          results/habitat_quality_delta.tif
          results/rusle_soil_loss_delta.tif
Outputs : results/eci_collapse_hotspots.tif
          results/eci_collapse_hotspots_summary.csv
          figures/synthesis/eci_collapse_hotspots_map.png

Run with: venv/bin/python src/synthesis_hotspots.py
"""

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import rasterio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap as LSC
import pandas as pd
from pathlib import Path

# ── Directory Paths ────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).parent
PROC_DIR   = SCRIPT_DIR.parent / 'data' / 'processed'
RES_DIR    = SCRIPT_DIR.parent / 'results'
FIG_DIR    = SCRIPT_DIR.parent / 'figures'
RES_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ── Print Separator ────────────────────────────────────────────────────
SEP = '=' * 60

# ── Split-map row boundaries (Ten Degree Channel) ─────────────────────
ANDAMAN_ROW_END   = 11750
NICOBAR_ROW_START = 16200
_PAD = 200


def _split_extent(arr):
    h = arr.shape[0]
    return arr[:min(ANDAMAN_ROW_END + _PAD, h), :], \
           arr[max(0, NICOBAR_ROW_START - _PAD):, :]


def _tight_crop(arr, valid=None):
    """Tight crop to non-NaN bounding box."""
    if valid is None:
        if np.ma.is_masked(arr):
            valid = ~arr.mask
        else:
            valid = ~np.isnan(arr)
    if not valid.any():
        return arr, valid
    rows = np.any(valid, axis=1); cols = np.any(valid, axis=0)
    r0, r1 = np.where(rows)[0][[0, -1]]
    c0, c1 = np.where(cols)[0][[0, -1]]
    pad = 60
    r0 = max(0, r0-pad); r1 = min(arr.shape[0], r1+pad)
    c0 = max(0, c0-pad); c1 = min(arr.shape[1], c1+pad)
    return arr[r0:r1, c0:c1], valid[r0:r1, c0:c1]

# ── ECI Percentile Thresholds ──────────────────────────────────────────
HOTSPOT_PERCENTILE = 95   # Top 5 % of impacted pixels = critical hotspot


# ══════════════════════════════════════════════════════════════════════
# UTILITY
# ══════════════════════════════════════════════════════════════════════
def load_result_raster(file_path: Path):
    """Load a results GeoTIFF as a float array.

    Returns (data_array, raster_profile).
    Returns (None, None) if the file does not exist.
    """
    if not file_path.exists():
        return None, None
    with rasterio.open(file_path) as src:
        data_array   = src.read(1)
        rast_profile = src.profile
    return data_array, rast_profile


# ══════════════════════════════════════════════════════════════════════
# ENTRY POINT
# ══════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    print(SEP)
    print("  ANI Ecosystem Services — ECI Hotspot Synthesis")
    print(SEP)

    # 1. Load all three independent loss layers ────────────────────────
    print("  Loading loss layers …")

    gfw_path  = PROC_DIR / 'ANI_GFW_Forest_Loss_2001_2023_clipped.tif'
    hab_path  = RES_DIR  / 'habitat_quality_delta.tif'
    soil_path = RES_DIR  / 'rusle_soil_loss_delta.tif'

    gfw_grid,  gfw_profile  = load_result_raster(gfw_path)
    hab_delta, hab_profile  = load_result_raster(hab_path)
    soil_delta, _           = load_result_raster(soil_path)

    if gfw_grid is None or hab_delta is None or soil_delta is None:
        print("  ❌  Missing required raster layers — run full pipeline first.")
        exit(1)

    print("  All layers loaded. Building valid-land mask …")

    # 2. Build valid-land mask ─────────────────────────────────────────
    valid_land = (
        (gfw_grid  != gfw_profile['nodata'])
        & (~np.isnan(hab_delta))
        & (~np.isnan(soil_delta))
        & (gfw_grid >= 0)
    )
    print(f"  Valid land pixels : {valid_land.sum():,}")

    # 3. Normalise each harm layer to [0, 1] ───────────────────────────

    # HABITAT HARM  — delta Q is negative where quality fell; invert → positive harm
    hab_harm     = np.zeros_like(hab_delta, dtype=float)
    hab_harm[valid_land] = -1.0 * hab_delta[valid_land]
    hab_harm     = np.clip(hab_harm, 0.0, None)                 # ignore quality gains
    hab_max      = np.percentile(hab_harm[valid_land], 99.9) if hab_harm.max() > 0 else 1.0
    hab_norm     = np.clip(hab_harm / hab_max, 0, 1)

    # SOIL HARM     — delta A is positive where erosion increased
    soil_harm    = np.zeros_like(soil_delta, dtype=float)
    soil_harm[valid_land] = soil_delta[valid_land]
    soil_harm    = np.clip(soil_harm, 0.0, None)                # ignore erosion reductions
    soil_log     = np.log1p(soil_harm)                          # log1p for right-skewed data
    soil_max     = np.percentile(soil_log[valid_land], 99.9) if soil_log.max() > 0 else 1.0
    soil_norm    = np.clip(soil_log / soil_max, 0, 1)

    # CARBON DRIVER — binary: 1 where GFW recorded any deforestation
    defor_flag   = np.zeros_like(gfw_grid, dtype=float)
    defor_flag[valid_land & (gfw_grid > 0)] = 1.0

    print("  Layers normalised. Computing ECI …")

    # 4. Ecosystem Collapse Index ──────────────────────────────────────
    eci_grid     = np.zeros_like(hab_delta, dtype=np.float32)
    eci_grid[valid_land] = (
        hab_norm[valid_land] + soil_norm[valid_land] + defor_flag[valid_land]
    ) / 3.0
    eci_grid[~valid_land] = np.nan

    # Analyse impacted pixels only (ECI > 0)
    valid_eci    = eci_grid[valid_land]
    impacted_eci = valid_eci[valid_eci > 0]

    print(f"\n{SEP}")
    print("  ECI ANALYTICAL METRICS")
    print(f"{SEP}")

    if len(impacted_eci) == 0:
        print("  No ecosystem harm detected.")
        exit(0)

    hotspot_threshold = np.percentile(impacted_eci, HOTSPOT_PERCENTILE)
    hotspot_px_mask   = (eci_grid >= hotspot_threshold) & valid_land
    hotspot_area_ha   = hotspot_px_mask.sum() * 0.09

    print(f"  Total degraded footprint : {len(impacted_eci) * 0.09:,.1f} ha")
    print(f"  Mean ECI                 : {impacted_eci.mean():.4f}")
    print(f"  {HOTSPOT_PERCENTILE}th Percentile Threshold : {hotspot_threshold:.4f}")
    print(f"  Critical Hotspot Area    : {hotspot_area_ha:,.1f} ha")

    # 5. Save GeoTIFF ──────────────────────────────────────────────────
    print(f"\n{SEP}")
    print("  SAVING OUTPUTS")
    print(f"{SEP}")

    export_profile = hab_profile.copy()
    export_profile.update(dtype=rasterio.float32, nodata=-9999.0, compress='lzw')
    eci_export = eci_grid.copy()
    eci_export[~valid_land] = -9999.0

    tif_out = RES_DIR / 'eci_collapse_hotspots.tif'
    with rasterio.open(tif_out, 'w', **export_profile) as dst:
        dst.write(eci_export.astype(np.float32), 1)
    print(f"  ✅  GeoTIFF saved → {tif_out.name}")

    # Save summary CSV ─────────────────────────────────────────────────
    summary_df = pd.DataFrame([{
        'total_degraded_ha':    round(len(impacted_eci) * 0.09, 1),
        'mean_eci':             round(float(impacted_eci.mean()), 4),
        'hotspot_threshold':    round(float(hotspot_threshold),   4),
        'critical_hotspot_ha':  round(float(hotspot_area_ha),     1),
    }])
    csv_out = RES_DIR / 'eci_collapse_hotspots_summary.csv'
    summary_df.to_csv(csv_out, index=False)
    print(f"  ✅  CSV saved → {csv_out.name}")

    # 6. Generate Figure — redesigned: Andaman | Nicobar | colorbar | histogram
    # ─────────────────────────────────────────────────────────────────────────
    # Root-cause fixes vs old version:
    #   (a) Land pixels without ECI data showed as ocean-blue (masked=bad colour)
    #       → Now we load the LC raster and build an RGBA composite so all land
    #         is always visible in sage-green regardless of ECI coverage.
    #   (b) Colorbar spanned 0–1 but all data ≥ 0.667, washing out contrast.
    #       → Now vmin is set to the 5th-percentile of actual ECI values.
    #   (c) Hotspot dots were s=3 (invisible at report scale).
    #       → Large bright-yellow filled circles with a dark-red edge.
    #   (d) _tight_crop used the ECI-valid mask (sparse), not the full land mask,
    #       so panels were zoomed to tiny clusters instead of the whole chain.
    #       → Tight-crop is now driven by the full LC land mask.
    # ─────────────────────────────────────────────────────────────────────────

    import matplotlib.patches as mpatches

    BG_COLOR    = 'white'
    TEXT_COLOR  = '#1a1a2e'
    OCEAN_COLOR = '#b3d9f2'   # light sky-blue ocean
    LAND_BASE   = (220, 237, 200)   # sage-green RGB (0–255) for unaffected land

    # ── Load full LC land mask for geographic context ──────────────────────
    lc_path = PROC_DIR / 'ANI_ESA_WorldCover_mosaic_clipped.tif'
    try:
        with rasterio.open(lc_path) as _lc:
            _lc_data = _lc.read(1).astype(float)
            _lc_nd   = _lc.nodata
        if _lc_nd is not None:
            _lc_data[_lc_data == _lc_nd] = np.nan
        # ESA WorldCover: value=0 means ocean/no-data (nodata attr may be None)
        land_full = (~np.isnan(_lc_data)) & (_lc_data > 0)
    except Exception:
        land_full = valid_land.copy()   # fallback

    # Align land_full to ECI grid shape
    def _align_grid(arr, ref_shape):
        if arr.shape == ref_shape:
            return arr
        out = np.zeros(ref_shape, dtype=arr.dtype)
        r = min(arr.shape[0], ref_shape[0])
        c = min(arr.shape[1], ref_shape[1])
        out[:r, :c] = arr[:r, :c]
        return out

    land_full = _align_grid(land_full, eci_grid.shape)

    # ── ECI normalisation ──────────────────────────────────────────────────
    # ECI = (habitat + soil + deforest) / 3  →  values are quantized to
    # {0, 0.33, 0.67, 1.0}. Anchor the colorbar to the meaningful category
    # range [0.33, 1.0] (one-collapse → triple-collapse) so the three risk
    # tiers in the legend map directly to three colors in the colorbar.
    eci_data_vals = eci_grid[valid_land]
    eci_vmin = 1.0 / 3.0
    eci_vmax = 1.0

    fire_colors = ['#ffffcc', '#fed976', '#fd8d3c', '#f03b20', '#bd0026', '#800026']
    cmap_fire   = LSC.from_list('eci_fire', fire_colors, N=256)

    # ── Per-pixel RGBA composite builder ──────────────────────────────────
    def build_rgba(eci_arr, land_arr, vmin, vmax):
        """
        Returns an RGBA float32 image (fully opaque):
          - Non-land pixels      → ocean blue (opaque)
          - Land with no ECI     → sage-green (opaque)
          - Land with ECI data   → fire colormap (opaque)
        """
        H, W = land_arr.shape
        # Start with ocean blue everywhere (fully opaque)
        oc_r, oc_g, oc_b = [int(OCEAN_COLOR.lstrip('#')[i:i+2], 16) / 255.0
                             for i in (0, 2, 4)]
        rgba = np.ones((H, W, 4), dtype=np.float32)
        rgba[:, :, 0] = oc_r
        rgba[:, :, 1] = oc_g
        rgba[:, :, 2] = oc_b

        # Sage-green land base
        r_f, g_f, b_f = [x / 255.0 for x in LAND_BASE]
        rgba[land_arr, 0] = r_f
        rgba[land_arr, 1] = g_f
        rgba[land_arr, 2] = b_f

        # ECI overlay (fire colourmap, dynamic range)
        eci_valid = ~np.ma.getmaskarray(eci_arr) if hasattr(eci_arr, 'mask') \
                    else ~np.isnan(np.ma.filled(eci_arr, np.nan))
        if eci_valid.any():
            filled = np.ma.filled(eci_arr, 0.0).astype(float)
            norm   = np.clip((filled - vmin) / max(vmax - vmin, 1e-8), 0.0, 1.0)
            colours = cmap_fire(norm)   # H×W×4
            rgba[eci_valid, :3] = colours[eci_valid, :3]

        return rgba

    # ── Tight-crop helper driven by FULL LAND extent ───────────────────────
    def land_tight_crop(*arrays, ref_mask, pad=80):
        """Crop all arrays to the bounding box of ref_mask with padding."""
        if not ref_mask.any():
            return arrays
        rows = np.any(ref_mask, axis=1)
        cols = np.any(ref_mask, axis=0)
        r0, r1 = np.where(rows)[0][[0, -1]]
        c0, c1 = np.where(cols)[0][[0, -1]]
        r0 = max(0, r0 - pad);  r1 = min(ref_mask.shape[0], r1 + pad)
        c0 = max(0, c0 - pad);  c1 = min(ref_mask.shape[1], c1 + pad)
        return tuple(a[r0:r1, c0:c1] for a in arrays)

    # ── Split arrays into Andaman / Nicobar slices ─────────────────────────
    eci_masked = np.ma.masked_where(~valid_land, eci_grid)

    (and_eci_r, nic_eci_r)   = _split_extent(eci_masked)
    (and_land_r, nic_land_r) = _split_extent(land_full)
    (and_hot_r,  nic_hot_r)  = _split_extent(hotspot_px_mask)
    (and_vl_r,   nic_vl_r)   = _split_extent(valid_land)

    # Crop to full land extent (not just ECI-valid) so the whole island chain
    # is visible. Using and_land_r/nic_land_r as ref_mask draws the real
    # geography; ECI-valid mask alone would zoom to sparse hotspot clusters.
    and_eci, and_land, and_hot, and_vl = land_tight_crop(
        and_eci_r, and_land_r, and_hot_r, and_vl_r,
        ref_mask=and_land_r, pad=80)
    nic_eci, nic_land, nic_hot, nic_vl = land_tight_crop(
        nic_eci_r, nic_land_r, nic_hot_r, nic_vl_r,
        ref_mask=nic_land_r, pad=80)

    # ── Clean gridspec layout: [Andaman | Nicobar | colorbar] ─────────────
    # Direct two-layer imshow: ocean facecolor → sage land base → ECI overlay.
    # (Histogram dropped because ECI is near-binary — almost all impacted
    #  pixels sit at 1.0, so the distribution conveys no information.)
    fig = plt.figure(figsize=(16, 12), facecolor='white')
    gs  = gridspec.GridSpec(1, 3, figure=fig,
                            width_ratios=[2, 2, 0.08],
                            wspace=0.04)

    land_rgb = tuple(x / 255.0 for x in LAND_BASE)
    cmap_land = LSC.from_list('land_only', [land_rgb, land_rgb], N=2)

    def _draw_panel(ax, eci_arr, land_arr, hot_mask, title):
        ax.set_facecolor(OCEAN_COLOR)
        # Sage land base — masked where there's no land (ocean shows through)
        land_only = np.ma.masked_where(~land_arr, np.ones_like(land_arr, dtype=float))
        ax.imshow(land_only, cmap=cmap_land, vmin=0, vmax=1,
                  interpolation='nearest', origin='upper', aspect='equal')
        # ECI overlay — non-ECI land stays sage; ocean stays blue
        cmap_fire.set_bad(alpha=0.0)
        ax.imshow(eci_arr, cmap=cmap_fire, vmin=eci_vmin, vmax=eci_vmax,
                  interpolation='nearest', origin='upper', aspect='equal')
        hr, hc = np.where(hot_mask)
        if len(hr) > 0:
            # Small open black circles — visible against red ECI without
            # blotting it out. With ~2k hotspot pixels per island, filled
            # markers would obscure the underlying overlay.
            ax.scatter(hc, hr, s=6, facecolors='none', marker='o', alpha=0.9,
                       linewidths=0.5, edgecolors='#0d0d0d', zorder=6)
        ax.set_title(title, color=TEXT_COLOR, fontsize=11,
                     fontweight='bold', pad=6)
        ax.axis('off')
        n_land_px = int(land_arr.sum())
        ax.text(0.03, 0.03,
                f'Land: {n_land_px * 0.09:,.0f} ha\n'
                f'Critical hotspots: {hot_mask.sum() * 0.09:,.0f} ha',
                transform=ax.transAxes, color='#333333', fontsize=9,
                bbox=dict(boxstyle='round,pad=0.4',
                          facecolor='white', alpha=0.9, edgecolor='#cccccc'))

    ax_a = fig.add_subplot(gs[0, 0])
    ax_n = fig.add_subplot(gs[0, 1])
    cax  = fig.add_subplot(gs[0, 2])

    _draw_panel(ax_a, and_eci, and_land, and_hot,
                'Andaman Islands — Ecosystem Collapse Index (ECI)')
    _draw_panel(ax_n, nic_eci, nic_land, nic_hot,
                'Nicobar Islands — Ecosystem Collapse Index (ECI)')

    # ── Legend on Andaman panel ────────────────────────────────────────────
    legend_patches = [
        mpatches.Patch(facecolor=tuple(x/255 for x in LAND_BASE),
                       edgecolor='#888', label='Unaffected land'),
        mpatches.Patch(facecolor='#fed976', edgecolor='none',
                       label='Moderate risk (ECI ~0.67)'),
        mpatches.Patch(facecolor='#f03b20', edgecolor='none',
                       label='High risk (ECI ~0.85)'),
        mpatches.Patch(facecolor='#800026', edgecolor='none',
                       label='Critical risk (ECI = 1.0)'),
        mpatches.Patch(facecolor='none', edgecolor='#0d0d0d',
                       label='Top 5% hotspot (triple-collapse)'),
        mpatches.Patch(facecolor=OCEAN_COLOR, edgecolor='#aaa',
                       label='Ocean / no data'),
    ]
    ax_a.legend(handles=legend_patches, loc='upper right', fontsize=8,
                framealpha=0.95, edgecolor='#cccccc', title='Legend',
                title_fontsize=8.5)

    # ── Colorbar — explicit ticks at the three ECI categories ────────────
    sm = plt.cm.ScalarMappable(cmap=cmap_fire,
                                norm=plt.Normalize(vmin=eci_vmin, vmax=eci_vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label('ECI  (Soil + Habitat + Carbon harm)',
                   color=TEXT_COLOR, fontsize=9, labelpad=6)
    cbar.set_ticks([1/3, 2/3, 1.0])
    cbar.set_ticklabels(['0.33\nModerate', '0.67\nHigh', '1.00\nCritical'])
    cbar.ax.tick_params(colors=TEXT_COLOR, labelsize=8)
    # Force the colorbar's y-axis to exactly span the data range — prevents
    # matplotlib's auto-ticker from stepping past 1.0 (ECI is bounded [0, 1]).
    cbar.ax.set_ylim(eci_vmin, eci_vmax)

    fig.suptitle('Triple-Collapse Ecosystem Risk — Andaman & Nicobar Islands (2000–2024)',
                 color=TEXT_COLOR, fontsize=14, fontweight='bold', y=0.995)

    (FIG_DIR / 'synthesis').mkdir(parents=True, exist_ok=True)
    fig_out = FIG_DIR / 'synthesis' / 'eci_collapse_hotspots_map.png'
    fig.savefig(fig_out, dpi=180, bbox_inches='tight',
                facecolor='white', pad_inches=0.15)
    plt.close()
    print(f"  ✅  Figure saved → {fig_out.name}")






    print(f"\n{SEP}")
    print("  ✅  ECI Hotspot Synthesis Complete!")
    print(f"      GeoTIFF : results/eci_collapse_hotspots.tif")
    print(f"      CSV     : results/eci_collapse_hotspots_summary.csv")
    print(f"      Figure  : figures/synthesis/eci_collapse_hotspots_map.png")
    print(SEP)
