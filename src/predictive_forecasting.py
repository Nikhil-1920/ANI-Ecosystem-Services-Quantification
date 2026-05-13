"""
ANI Ecosystem Services — Predictive BAU Scenario Forecasting (2040)
====================================================================
Projects the historical 2000–2024 deforestation rate to Year 2040
using morphological edge-dilation: deforestation expands radially
from existing infrastructure / clearings at the observed annual rate.
Computes the Economic Value of Foregone Ecosystem Services (ESV).

Inputs  : data/processed/ANI_ESA_WorldCover_mosaic_clipped.tif
          data/processed/ANI_GFW_Forest_Loss_2001_2023_clipped.tif
          data/processed/ANI_GEDI_Biomass_Density_clipped.tif
Outputs : results/economic_damage_2040_forecast.csv
          figures/predictive/forecast_forest_cover_2040.png

Run with: venv/bin/python src/predictive_forecasting.py
"""

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import rasterio
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
from scipy.ndimage import binary_dilation
from pathlib import Path

# ── Directory Paths ────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).parent
PROC_DIR   = SCRIPT_DIR.parent / 'data' / 'processed'
RES_DIR    = SCRIPT_DIR.parent / 'results'
FIG_DIR    = SCRIPT_DIR.parent / 'figures'

# ── Print Separator ────────────────────────────────────────────────────
SEP = '=' * 60

# ── Temporal Constants ─────────────────────────────────────────────────
HISTORIC_YEARS  = 24     # 2000–2024
FORECAST_YEARS  = 16     # 2024–2040

# ── Physical Constants ─────────────────────────────────────────────────
PIXEL_AREA_HA   = 0.09   # 30 m × 30 m = 0.09 ha
ROOT_SHOOT_R    = 1.24   # Root-to-shoot ratio (IPCC tropical)
IPCC_CARBON_F   = 0.47   # Biomass to carbon fraction
CO2E_MW_RATIO   = 44/12  # Carbon to CO₂-equivalent molecular weight ratio

# ── Economic Valuation Constants (USD) ────────────────────────────────
SCC_PER_TONNE_CO2E  = 51.00        # EPA Social Cost of Carbon
DREDGE_PER_TONNE    = 5.00         # Coastal sediment dredging unit cost
HABITAT_RESTORE_HA  = 12_000.00    # Tropical habitat replacement cost per ha
SOIL_EROSION_SPIKE  = 50.0         # Induced erosion uplift vs. forest (t/ha/yr)

# ── Display Constants ──────────────────────────────────────────────────
DISPLAY_DOWNSAMPLE  = 10    # Spatial downsample factor for matplotlib imshow
DILATION_MAX_ITER   = 500   # Safety cap on while-loop iterations

# ── Map Colours ────────────────────────────────────────────────────────
MAP_PALETTE = ['#dce3ed', '#2e7d32', '#ffb300', '#e53935']
# Index:         0=ocean    1=forest   2=hist_loss  3=future_loss


# ══════════════════════════════════════════════════════════════════════
# UTILITY
# ══════════════════════════════════════════════════════════════════════
def load_raster(file_path: Path):
    """Load band 1 of a GeoTIFF with nodata replaced by nan.

    Returns (data_array, raster_profile).
    """
    with rasterio.open(file_path) as src:
        data_array   = src.read(1).astype(float)
        rast_profile = src.profile
        nodata_val   = src.nodata
    if nodata_val is not None:
        data_array = np.where(data_array == nodata_val, np.nan, data_array)
    return data_array, rast_profile


# ── Split-map row boundaries (Ten Degree Channel) ─────────────────────
ANDAMAN_ROW_END   = 11750
NICOBAR_ROW_START = 16200
_PAD = 200


def _split_extent(arr, downsample=1):
    """Return (andaman_slice, nicobar_slice) of a 2-D array.

    Scales the row boundaries by `downsample` so the split is correct on
    arrays that have already been spatially downsampled for display.
    """
    h = arr.shape[0]
    and_end = (ANDAMAN_ROW_END   + _PAD) // downsample
    nic_beg = (NICOBAR_ROW_START - _PAD) // downsample
    return arr[:min(and_end, h), :], \
           arr[max(0, nic_beg):, :]


def _crop(arr, dom):
    """Tight-crop arr to the bounding box of dom (bool mask) + padding."""
    if not dom.any():
        return arr
    rows = np.any(dom, axis=1)
    cols = np.any(dom, axis=0)
    r0, r1 = np.where(rows)[0][[0, -1]]
    c0, c1 = np.where(cols)[0][[0, -1]]
    pad = 60
    r0 = max(0, r0 - pad); r1 = min(arr.shape[0], r1 + pad)
    c0 = max(0, c0 - pad); c1 = min(arr.shape[1], c1 + pad)
    return arr[r0:r1, c0:c1]




# ══════════════════════════════════════════════════════════════════════
# SPATIAL EDGE-DILATION ENGINE
# ══════════════════════════════════════════════════════════════════════
def expand_deforestation_frontier(target_pixel_count: int,
                                  intact_forest_mask: np.ndarray,
                                  historic_loss_mask: np.ndarray) -> np.ndarray:
    """Simulate future deforestation by dilating the historic loss frontier.

    Each iteration expands the clearing boundary by one 30 m pixel,
    consuming only currently intact forest, until the target pixel quota
    (derived from the historic annual rate × forecast years) is reached.

    Returns a boolean array of predicted future loss pixels.
    """
    future_loss_mask = np.zeros_like(historic_loss_mask, dtype=bool)
    active_frontier  = historic_loss_mask.copy()
    n_consumed       = 0

    for iteration in range(DILATION_MAX_ITER):
        dilated_mask   = binary_dilation(active_frontier)
        new_frontier   = dilated_mask & intact_forest_mask & (~future_loss_mask)

        if new_frontier.sum() == 0:
            print(f"  ⚠️  Spatial saturation at iteration {iteration}. "
                  f"Target {target_pixel_count:,} px unattainable.")
            break

        future_loss_mask = future_loss_mask | new_frontier
        active_frontier  = new_frontier
        n_consumed       = future_loss_mask.sum()

        if n_consumed >= target_pixel_count:
            break

    return future_loss_mask


# ══════════════════════════════════════════════════════════════════════
# ENTRY POINT
# ══════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    print(SEP)
    print("  ANI Ecosystem Services — Predictive 2040 BAU Forecasting")
    print(SEP)

    # 1. Load baseline spatial layers ──────────────────────────────────
    print("  Loading baseline spatial layers …")
    landcover_grid, _ = load_raster(PROC_DIR / 'ANI_ESA_WorldCover_mosaic_clipped.tif')
    defor_yr_grid,  _ = load_raster(PROC_DIR / 'ANI_GFW_Forest_Loss_2001_2023_clipped.tif')
    agb_density,    _ = load_raster(PROC_DIR / 'ANI_GEDI_Biomass_Density_clipped.tif')

    # 2. Derive historic velocity ───────────────────────────────────────
    valid_domain     = (~np.isnan(landcover_grid)) & (~np.isnan(defor_yr_grid))
    intact_forest    = (landcover_grid == 10) & (defor_yr_grid == 0) & valid_domain
    historic_loss    = (defor_yr_grid > 0) & valid_domain

    px_lost_historic  = int(np.sum(historic_loss))
    ha_lost_historic  = px_lost_historic * PIXEL_AREA_HA
    annual_px_rate    = px_lost_historic / HISTORIC_YEARS
    annual_ha_rate    = ha_lost_historic / HISTORIC_YEARS

    print(f"\n{SEP}")
    print("  HISTORICAL VELOCITY ANALYSIS (2000–2024)")
    print(f"{SEP}")
    print(f"  Total Forest Cleared : {ha_lost_historic:,.1f} ha")
    print(f"  Annual Velocity      : {annual_ha_rate:,.1f} ha / year")

    # 3. BAU forecast to 2040 ──────────────────────────────────────────
    target_px_2040    = int(annual_px_rate * FORECAST_YEARS)
    target_ha_2040    = target_px_2040 * PIXEL_AREA_HA

    print(f"\n{SEP}")
    print("  BUSINESS-AS-USUAL FORECAST TARGET (2040)")
    print(f"{SEP}")
    print(f"  Projected Additional Deforestation : +{target_ha_2040:,.1f} ha")
    print(f"  Target Expansion Pixel Quota       : {target_px_2040:,} pixels")

    print("  Running morphological edge-dilation simulation …")
    future_loss_mask  = expand_deforestation_frontier(
        target_px_2040, intact_forest, historic_loss
    )

    px_achieved       = int(future_loss_mask.sum())
    ha_achieved       = px_achieved * PIXEL_AREA_HA
    print(f"  Simulation achieved  : {ha_achieved:,.1f} ha projected loss by 2040")

    # 4. Economic Services Valuation (ESV) ────────────────────────────
    print(f"\n{SEP}")
    print("  ECOSYSTEM SERVICES VALUATION (ESV) — 2040 ECONOMIC DAMAGES")
    print(f"{SEP}")

    # Carbon damage — GEDI biomass at predicted loss pixels → CO₂e → USD
    bio_at_loss      = agb_density[future_loss_mask]
    bio_at_loss      = np.where(np.isnan(bio_at_loss), 0, bio_at_loss)
    total_biomass_mg  = np.sum(bio_at_loss) * PIXEL_AREA_HA * ROOT_SHOOT_R
    total_carbon_mg   = total_biomass_mg * IPCC_CARBON_F
    total_co2e_mg     = total_carbon_mg  * CO2E_MW_RATIO
    carbon_damage_usd = total_co2e_mg * SCC_PER_TONNE_CO2E

    # Sediment / dredging damage — 16 years of induced erosion uplift
    induced_erosion_t  = ha_achieved * SOIL_EROSION_SPIKE * FORECAST_YEARS
    dredging_cost_usd  = induced_erosion_t * DREDGE_PER_TONNE

    # Habitat replacement cost
    habitat_cost_usd   = ha_achieved * HABITAT_RESTORE_HA

    total_damage_usd   = carbon_damage_usd + dredging_cost_usd + habitat_cost_usd

    print(f"  Projected CO₂e Emissions        : {total_co2e_mg:,.1f} t CO₂e")
    print(f"    → Atmospheric Damage Cost     : ${carbon_damage_usd:,.2f}")
    print(f"  Projected Sediment Loss         : {induced_erosion_t:,.1f} t")
    print(f"    → Dredging Cost               : ${dredging_cost_usd:,.2f}")
    print(f"  Habitat Replacement Cost        : ${habitat_cost_usd:,.2f}")
    print(f"\n  💵  TOTAL 2040 DAMAGES          : ${total_damage_usd:,.2f}")

    # 5. Save results CSV ──────────────────────────────────────────────
    forecast_df = pd.DataFrame([{
        'historic_deforestation_ha':    round(ha_lost_historic,       1),
        'annual_velocity_ha_per_yr':    round(annual_ha_rate,         1),
        'projected_loss_ha_2040':       round(ha_achieved,            1),
        'co2e_emissions_tonnes':        round(total_co2e_mg,          1),
        'carbon_damage_usd':            round(carbon_damage_usd,      2),
        'sediment_loss_tonnes':         round(induced_erosion_t,      1),
        'dredging_cost_usd':            round(dredging_cost_usd,      2),
        'habitat_replacement_usd':      round(habitat_cost_usd,       2),
        'total_economic_damage_usd':    round(total_damage_usd,       2),
    }])
    csv_out = RES_DIR / 'economic_damage_2040_forecast.csv'
    forecast_df.to_csv(csv_out, index=False)
    print(f"\n  ✅  CSV saved → {csv_out.name}")

    # 6. Render 2024 vs 2040 comparison map — split Andaman / Nicobar ──
    print(f"\n{SEP}")
    print("  RENDERING 2040 PREDICTIVE MAP")
    print(f"{SEP}")

    # Downsample for display stability (full 25 k × 7 k grid is too large)
    ds = DISPLAY_DOWNSAMPLE
    domain_ds      = valid_domain[::ds, ::ds]
    forest_ds      = intact_forest[::ds, ::ds]
    hist_loss_ds   = historic_loss[::ds, ::ds]
    fut_loss_ds    = future_loss_mask[::ds, ::ds]

    # ── Top row: 2024 baseline class map (ocean / forest / historic loss)
    map_2024 = np.zeros_like(domain_ds, dtype=int)
    map_2024[domain_ds]    = 0
    map_2024[forest_ds]    = 1
    map_2024[hist_loss_ds] = 2

    # ── Bottom row: delta map. Only the *new* 2024→2040 loss is shown,
    # dilated a few pixels so it's legible at this zoom, on a faded
    # backdrop of remaining forest so the red signal pops.
    delta_only = fut_loss_ds & ~hist_loss_ds
    delta_vis  = binary_dilation(delta_only, iterations=2)

    map_delta = np.zeros_like(domain_ds, dtype=int)
    map_delta[domain_ds] = 0     # ocean / non-forest
    map_delta[forest_ds] = 1     # surviving forest (faded backdrop)
    map_delta[delta_vis] = 2     # newly predicted collapse 2024→2040

    cmap_baseline = ListedColormap(['#dce3ed', '#2e7d32', '#ffb300'])
    cmap_delta    = ListedColormap(['#dce3ed', '#b8d4ba', '#e53935'])

    # Split each map by island group (Ten Degree Channel boundary)
    and_2024, nic_2024  = _split_extent(map_2024,  downsample=ds)
    and_delta, nic_delta = _split_extent(map_delta, downsample=ds)
    dom_and,  dom_nic   = _split_extent(domain_ds, downsample=ds)

    and_2024_c  = _crop(and_2024,  dom_and)
    nic_2024_c  = _crop(nic_2024,  dom_nic)
    and_delta_c = _crop(and_delta, dom_and)
    nic_delta_c = _crop(nic_delta, dom_nic)

    # 2-row × 2-col figure: top row = 2024 baseline, bottom row = 2024→2040 delta
    fig = plt.figure(figsize=(18, 22), facecolor='#1a1a2e')
    gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.10, wspace=0.04)
    panels = [
        (gs[0, 0], and_2024_c,  cmap_baseline, 2, 'Andaman Islands — Baseline (2024)'),
        (gs[0, 1], nic_2024_c,  cmap_baseline, 2, 'Nicobar Islands — Baseline (2024)'),
        (gs[1, 0], and_delta_c, cmap_delta,    2, 'Andaman Islands — New Predicted Loss (2024→2040)'),
        (gs[1, 1], nic_delta_c, cmap_delta,    2, 'Nicobar Islands — New Predicted Loss (2024→2040)'),
    ]
    for gs_pos, arr, cmap_arg, vmax, title in panels:
        ax = fig.add_subplot(gs_pos)
        ax.set_facecolor('white')
        ax.imshow(arr, cmap=cmap_arg, vmin=0, vmax=vmax, interpolation='nearest')
        ax.set_title(title, color='#1a1a2e', fontsize=12, fontweight='bold', pad=10)
        ax.axis('off')

    legend_patches = [
        mpatches.Patch(color='#d0d8e4', label='Ocean / Non-Forest'),
        mpatches.Patch(color='#2e7d32', label='Intact Forest Sanctuary (2024)'),
        mpatches.Patch(color='#ffb300', label='Anthropogenic Loss (2000–24)'),
        mpatches.Patch(color='#b8d4ba', label='Surviving Forest (delta backdrop)'),
        mpatches.Patch(color='#e53935', label='New Predicted Collapse (2024–40)'),
    ]
    fig.legend(handles=legend_patches, loc='lower center',
               bbox_to_anchor=(0.5, 0.01), ncol=5,
               facecolor='#f0f2f5', edgecolor='#cccccc',
               labelcolor='#1a1a2e', fontsize=11)

    plt.suptitle(
        f'Predictive Edge-Dilation Spreading Model\n'
        f'Estimated Economic Damages by 2040: '
        f'${total_damage_usd / 1_000_000:,.1f} Million',
        color='#1a1a2e', fontsize=16, fontweight='bold', y=1.01,
    )
    (FIG_DIR / 'predictive').mkdir(parents=True, exist_ok=True)
    fig_out = FIG_DIR / 'predictive' / 'forecast_forest_cover_2040.png'
    fig.savefig(fig_out, dpi=250, bbox_inches='tight', facecolor='#1a1a2e')
    plt.close()
    print(f"  ✅  Figure saved → {fig_out.name}")

    print(f"\n{SEP}")
    print("  ✅  Predictive 2040 Simulation & ESV Complete!")
    print(f"      CSV    : results/economic_damage_2040_forecast.csv")
    print(f"      Figure : figures/predictive/forecast_forest_cover_2040.png")
    print(SEP + "\n")

