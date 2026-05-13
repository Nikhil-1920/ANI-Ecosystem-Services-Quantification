"""
ANI Ecosystem Services — Multi-Scenario Forecasting (2024–2060)
================================================================
Models three socio-political futures based on the historic deforestation
velocity, using the morphological edge-dilation spreading engine:

  1. Conservation (0.2× rate + high-quality core habitat protected)
  2. Business-As-Usual (1.0× rate)
  3. Escalation (2.5× rate — mega-project / uncontrolled development)

Outputs economic damage trajectories at 2030, 2040, 2050, 2060 for
each scenario; renders both a tri-panel spatial map and an economic
cone-of-uncertainty line chart.

Inputs  : data/processed/ANI_ESA_WorldCover_mosaic_clipped.tif
          data/processed/ANI_GFW_Forest_Loss_2001_2023_clipped.tif
          results/habitat_quality_index.tif
Outputs : results/economic_scenarios_2024_2060.csv
          figures/predictive/forecast_economic_damages_2060.png
          figures/predictive/forecast_tri_scenario_2060.png

Run with: venv/bin/python src/predictive_scenarios.py
"""

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import rasterio
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
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
HISTORIC_YEARS   = 24    # Years of observed GFW data (2000–2024)
FORECAST_YEARS   = 36    # Forecast horizon (2024–2060)

# ── Physical Constants ─────────────────────────────────────────────────
PIXEL_AREA_HA    = 0.09  # 30 m × 30 m pixel = 0.09 ha

# ── Economic Valuation Constants (USD) ────────────────────────────────
SCC_PER_TONNE_CO2E  = 51.00      # EPA Social Cost of Carbon (2024)
DREDGE_PER_TONNE    = 5.00       # Coastal sediment dredging unit cost
HABITAT_RESTORE_HA  = 12_000.00  # Tropical habitat replacement cost per ha
# Landscape-level average carbon density for ESV approximation
LANDSCAPE_CO2E_HA   = 150.0 * 1.24 * 0.47 * (44/12)  # ≈ 320 t CO₂e/ha
SOIL_EROSION_SPIKE  = 50.0       # Induced erosion uplift vs. intact forest (t/ha/yr)

# ── Scenario Configuration ─────────────────────────────────────────────
SCENARIO_CONFIG = {
    'Conservation (Best Case)': {'rate_multiplier': 0.2, 'protect_core': True},
    'Business-As-Usual':        {'rate_multiplier': 1.0, 'protect_core': False},
    'Escalation (Worst Case)':  {'rate_multiplier': 2.5, 'protect_core': False},
}

# ── Scenario Timeline Epochs ───────────────────────────────────────────
EPOCH_YEARS         = [2030, 2040, 2050, 2060]
EPOCH_DELTA_YEARS   = [6,    16,   26,   36]   # Years from 2024

# ── Spatial Map Display ────────────────────────────────────────────────
DISPLAY_DOWNSAMPLE  = 10
MAP_PALETTE         = ['#1a1a2e', '#2e7d32', '#ffb300', '#e53935']
# Index:               0=ocean    1=forest   2=hist_loss  3=future_loss

# ── Core Habitat Protection Threshold ─────────────────────────────────
CORE_HABITAT_Q_MIN  = 0.8   # Habitat quality ≥ 0.8 → protected in Conservation

# ── Dilation Engine Safety Cap ─────────────────────────────────────────
DILATION_MAX_ITER   = 1000


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


# ══════════════════════════════════════════════════════════════════════
# SPATIAL EDGE-DILATION ENGINE
# ══════════════════════════════════════════════════════════════════════
def expand_deforestation_scenario(target_pixel_count: int,
                                  intact_forest_mask: np.ndarray,
                                  historic_loss_mask: np.ndarray,
                                  hab_quality_grid:   np.ndarray = None,
                                  protect_core:       bool       = False
                                  ) -> tuple[np.ndarray, list[int]]:
    """Run the morphological edge-dilation deforestation spreading model.

    If protect_core=True and hab_quality_grid is provided, pixels with
    Habitat Quality ≥ CORE_HABITAT_Q_MIN are excluded from spreading.

    Returns:
        future_loss_mask : boolean array of predicted future loss pixels
        epoch_px_counts  : pixel counts achieved at [2030, 2040, 2050, 2060]
    """
    future_loss_mask = np.zeros_like(historic_loss_mask, dtype=bool)
    active_frontier  = historic_loss_mask.copy()
    n_consumed       = 0

    # Apply core-habitat protection if requested
    legal_expansion  = intact_forest_mask.copy()
    if protect_core and hab_quality_grid is not None:
        legal_expansion = legal_expansion & (hab_quality_grid < CORE_HABITAT_Q_MIN)

    # Pre-compute pixel quotas at each epoch proportionally
    quota_fractions  = [6/36, 16/36, 26/36, 1.0]
    epoch_quotas     = [int(target_pixel_count * f) for f in quota_fractions]
    epoch_achieved   = []

    for iteration in range(DILATION_MAX_ITER):
        dilated_mask  = binary_dilation(active_frontier)
        new_frontier  = dilated_mask & legal_expansion & (~future_loss_mask)

        if new_frontier.sum() == 0:
            print(f"    ⚠️  Saturation reached at iteration {iteration}.")
            break

        future_loss_mask = future_loss_mask | new_frontier
        active_frontier  = new_frontier
        n_consumed       = int(future_loss_mask.sum())

        # Log checkpoint when an epoch quota is crossed
        if epoch_quotas and n_consumed >= epoch_quotas[0]:
            epoch_achieved.append(n_consumed)
            epoch_quotas.pop(0)

        if n_consumed >= target_pixel_count:
            break

    # Pad epoch counts if saturation stopped early
    while len(epoch_achieved) < 4:
        epoch_achieved.append(n_consumed)

    return future_loss_mask, epoch_achieved


# ══════════════════════════════════════════════════════════════════════
# CHART RENDERER (shared by the full pipeline and the standalone
# regenerator at src/render_economic_cone.py)
# ══════════════════════════════════════════════════════════════════════
def _render_economic_cone(df: pd.DataFrame, out_path: Path) -> None:
    """Render the light-themed economic-damage trajectory cone.

    `df` must have columns `Year`, `Scenario`, `Economic_Damages_USD`.
    Scenario labels are matched on a `startswith` basis so both
    "Conservation" and "Conservation (Best Case)" work.
    """
    df = df.copy()
    df['millions'] = df['Economic_Damages_USD'] / 1_000_000

    def _series(prefix: str) -> pd.Series:
        mask = df['Scenario'].str.startswith(prefix)
        return df[mask].sort_values('Year').set_index('Year')['millions']

    cons = _series('Conservation')
    bau  = _series('Business')
    esc  = _series('Escalation')
    years = cons.index.values

    final_cons = cons.iloc[-1]
    final_bau  = bau.iloc[-1]
    final_esc  = esc.iloc[-1]
    gap_bau    = final_bau - final_cons
    gap_esc    = final_esc - final_cons
    final_year = int(years[-1])

    plt.style.use('default')
    fig, ax = plt.subplots(figsize=(11, 6.2), facecolor='white')
    ax.set_facecolor('white')

    colour_cons = '#2e7d32'
    colour_bau  = '#ef6c00'
    colour_esc  = '#c62828'

    ax.fill_between(years, cons.values, bau.values,
                    color='#ffb300', alpha=0.13, zorder=1)
    ax.fill_between(years, bau.values, esc.values,
                    color='#c62828', alpha=0.10, zorder=1)

    ax.plot(years, cons.values, color=colour_cons, lw=2.6,
            marker='o', ms=6, label='Conservation (Best Case)', zorder=3)
    ax.plot(years, bau.values,  color=colour_bau,  lw=2.6,
            marker='s', ms=6, label='Business-As-Usual', zorder=3)
    ax.plot(years, esc.values,  color=colour_esc,  lw=2.6,
            marker='^', ms=7, label='Escalation (Worst Case)', zorder=3)

    # Endpoint value labels
    label_x = final_year + 0.4
    ax.text(label_x, final_cons, f'${final_cons:,.0f}M',
            color=colour_cons, fontsize=10, fontweight='bold', va='center')
    ax.text(label_x, final_bau,  f'${final_bau:,.0f}M',
            color=colour_bau,  fontsize=10, fontweight='bold', va='center')
    ax.text(label_x, final_esc,  f'${final_esc:,.0f}M',
            color=colour_esc,  fontsize=10, fontweight='bold', va='center')

    # Two paired callouts — split the previously ambiguous "saves $X" box.
    mid_idx       = len(years) // 2
    mid_year_bau  = years[mid_idx]
    mid_year_esc  = years[min(mid_idx + 1, len(years) - 1)]
    mid_y_bau     = (cons.iloc[mid_idx] + bau.iloc[mid_idx]) / 2
    mid_y_esc     = (bau.iloc[min(mid_idx + 1, len(years) - 1)] +
                      esc.iloc[min(mid_idx + 1, len(years) - 1)]) / 2

    ax.annotate(
        f'BAU → Conservation\navoided damage: ${gap_bau:,.0f}M',
        xy=(mid_year_bau, mid_y_bau),
        ha='center', va='center', fontsize=9, color='#5d4037',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='#fff8e1',
                  edgecolor='#ef6c00', linewidth=0.9),
    )
    ax.annotate(
        f'Escalation → Conservation\navoided damage: ${gap_esc:,.0f}M',
        xy=(mid_year_esc, mid_y_esc),
        ha='center', va='center', fontsize=9, color='#5d2a2a',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='#fdecea',
                  edgecolor='#c62828', linewidth=0.9),
    )

    # Direct in-plot zone labels (clearer than legend-only)
    label_year = years[max(1, len(years) // 4)]
    ax.text(label_year,
            (cons.loc[label_year] + bau.loc[label_year]) / 2,
            'Opportunity gap\n(BAU vs. Conservation)',
            ha='center', va='center', fontsize=8.5,
            color='#7a5a16', style='italic', alpha=0.85)
    ax.text(label_year,
            (bau.loc[label_year] + esc.loc[label_year]) / 2,
            'Escalation risk band\n(Escalation vs. BAU)',
            ha='center', va='center', fontsize=8.5,
            color='#7a2a2a', style='italic', alpha=0.85)

    ax.set_title(
        'Multi-Scenario Projected Economic Damages (2024–2060)\n'
        'Andaman & Nicobar Islands — Ecosystem Services Loss',
        fontsize=13, fontweight='bold', color='#1a1a2e', pad=14,
    )
    ax.set_xlabel('Year', fontsize=11, color='#1a1a2e')
    ax.set_ylabel('Cumulative Ecological Damages (Million USD)',
                  fontsize=11, color='#1a1a2e')
    ax.yaxis.set_major_formatter(
        plt.FuncFormatter(lambda v, _: f'${v:,.0f}M')
    )
    ax.grid(True, color='#dddddd', linewidth=0.7, alpha=0.9)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_color('#bbbbbb')
    ax.set_xlim(years[0] - 0.5, final_year + 4.5)
    ax.tick_params(colors='#1a1a2e')

    leg = ax.legend(loc='upper left', facecolor='white',
                    edgecolor='#cccccc', fontsize=9.5)
    for text in leg.get_texts():
        text.set_color('#1a1a2e')

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, facecolor='white', bbox_inches='tight')
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════
# ENTRY POINT
# ══════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    print(f"\n{SEP}")
    print("  ANI Ecosystem Services — Multi-Scenario Trajectories (2060)")
    print(SEP)

    # 1. Load spatial layers ───────────────────────────────────────────
    print("  Loading spatial layers …")
    landcover_grid, _ = load_raster(PROC_DIR / 'ANI_ESA_WorldCover_mosaic_clipped.tif')
    defor_yr_grid,  _ = load_raster(PROC_DIR / 'ANI_GFW_Forest_Loss_2001_2023_clipped.tif')
    hab_quality,    _ = load_raster(RES_DIR  / 'habitat_quality_index.tif')

    valid_domain   = (~np.isnan(landcover_grid)) & (~np.isnan(defor_yr_grid))
    intact_forest  = (landcover_grid == 10) & (defor_yr_grid == 0) & valid_domain
    historic_loss  = (defor_yr_grid > 0) & valid_domain

    # 2. Historic deforestation velocity ──────────────────────────────
    px_lost_historic = int(np.sum(historic_loss))
    ha_lost_historic = px_lost_historic * PIXEL_AREA_HA
    annual_px_rate   = px_lost_historic / HISTORIC_YEARS
    annual_ha_rate   = ha_lost_historic / HISTORIC_YEARS
    print(f"  Historic Rate    : {annual_ha_rate:,.1f} ha / year")

    # 3. Run all three scenarios ───────────────────────────────────────
    trajectory_rows = []
    scenario_masks  = {}

    # Seed each scenario with a zero-damage 2024 baseline row
    for scenario_name in SCENARIO_CONFIG:
        trajectory_rows.append({
            'year': 2024, 'scenario': scenario_name, 'economic_damage_usd': 0.0
        })

    print(f"\n{SEP}")
    print("  EXECUTING SPATIAL EDGE-DILATION SIMULATIONS")
    print(SEP)

    for scenario_name, cfg in SCENARIO_CONFIG.items():
        print(f"\n  Simulating: {scenario_name}")
        target_px = int(annual_px_rate * FORECAST_YEARS * cfg['rate_multiplier'])
        print(f"    Target Quota : {target_px * PIXEL_AREA_HA:,.1f} ha")

        fut_mask, epoch_px = expand_deforestation_scenario(
            target_pixel_count = target_px,
            intact_forest_mask = intact_forest,
            historic_loss_mask = historic_loss,
            hab_quality_grid   = hab_quality,
            protect_core       = cfg['protect_core'],
        )
        scenario_masks[scenario_name] = fut_mask

        # Compute cumulative economic damage at each epoch
        for i, n_px in enumerate(epoch_px):
            epoch_ha      = n_px * PIXEL_AREA_HA
            carbon_dmg    = epoch_ha * LANDSCAPE_CO2E_HA * SCC_PER_TONNE_CO2E
            dredging_dmg  = epoch_ha * SOIL_EROSION_SPIKE * EPOCH_DELTA_YEARS[i] * DREDGE_PER_TONNE
            habitat_dmg   = epoch_ha * HABITAT_RESTORE_HA
            total_usd     = carbon_dmg + dredging_dmg + habitat_dmg
            trajectory_rows.append({
                'year':                 EPOCH_YEARS[i],
                'scenario':             scenario_name,
                'projected_loss_ha':    round(epoch_ha,   1),
                'carbon_damage_usd':    round(carbon_dmg,   2),
                'dredging_damage_usd':  round(dredging_dmg, 2),
                'habitat_damage_usd':   round(habitat_dmg,  2),
                'economic_damage_usd':  round(total_usd,    2),
            })

    trajectory_df = pd.DataFrame(trajectory_rows)
    csv_out       = RES_DIR / 'economic_scenarios_2024_2060.csv'
    trajectory_df.to_csv(csv_out, index=False)
    print(f"\n  ✅  CSV saved → {csv_out.name}")

    # 4. Render economic trajectory cone chart ─────────────────────────
    print("\n  Generating Economic Trajectory Cone …")
    trajectory_df['damage_millions_usd'] = trajectory_df['economic_damage_usd'] / 1_000_000
    fig_econ = FIG_DIR / 'predictive' / 'forecast_economic_damages_2060.png'
    _render_economic_cone(
        trajectory_df.rename(columns={
            'year': 'Year',
            'scenario': 'Scenario',
            'economic_damage_usd': 'Economic_Damages_USD',
        })[['Year', 'Scenario', 'Economic_Damages_USD']],
        fig_econ,
    )
    print(f"  ✅  Figure saved → {fig_econ.name}")

    # 5. Render tri-panel spatial scenario map ─────────────────────────
    print("  Generating Tri-Panel Spatial Scenario Map …")

    ds           = DISPLAY_DOWNSAMPLE
    domain_ds    = valid_domain[::ds, ::ds]
    forest_ds    = intact_forest[::ds, ::ds]
    hist_loss_ds = historic_loss[::ds, ::ds]

    # Per-scenario ha + $ summary (final epoch = 2060)
    scenario_summary = {}
    for scen in SCENARIO_CONFIG:
        sub = trajectory_df[(trajectory_df['scenario'] == scen) &
                            (trajectory_df['year']     == EPOCH_YEARS[-1])].iloc[0]
        scenario_summary[scen] = (sub['projected_loss_ha'],
                                  sub['economic_damage_usd'] / 1_000_000)

    panel_titles = [
        '(A) Conservation\n(0.2× historic rate, core-protected)',
        '(B) Business-As-Usual\n(1.0× historic rate)',
        '(C) Escalation\n(2.5× historic rate)',
    ]
    scenario_names = list(SCENARIO_CONFIG.keys())

    # Light theme + faded-forest backdrop + dilated red so the predicted
    # spread is actually legible at this zoom (full grid is 25 k × 7 k).
    cmap_delta = ListedColormap(['#dce3ed', '#b8d4ba', '#ffb300', '#e53935'])

    fig2, panel_axes = plt.subplots(1, 3, figsize=(18, 10), facecolor='white')

    for ax_panel, scenario_name, panel_title in zip(
            panel_axes, scenario_names, panel_titles):
        fut_loss_ds = scenario_masks[scenario_name][::ds, ::ds]
        new_loss    = fut_loss_ds & ~hist_loss_ds
        new_loss_v  = binary_dilation(new_loss, iterations=2)
        hist_loss_v = binary_dilation(hist_loss_ds, iterations=1)

        class_map = np.zeros_like(domain_ds, dtype=int)
        class_map[domain_ds]    = 0  # ocean / non-forest
        class_map[forest_ds]    = 1  # surviving forest (faded backdrop)
        class_map[hist_loss_v]  = 2  # historic anthropogenic loss
        class_map[new_loss_v]   = 3  # predicted scenario spread

        ax_panel.set_facecolor('white')
        ax_panel.imshow(class_map, cmap=cmap_delta, vmin=0, vmax=3,
                        interpolation='nearest')
        ax_panel.set_title(panel_title, color='#1a1a2e',
                           fontsize=13, fontweight='bold', pad=10)
        ax_panel.axis('off')

        ha_loss, dmg_m = scenario_summary[scenario_name]
        ax_panel.text(
            0.5, -0.04,
            f'+{ha_loss:,.0f} ha new loss   |   ${dmg_m:,.0f}M damage',
            transform=ax_panel.transAxes,
            ha='center', va='top', fontsize=11, color='#1a1a2e',
            bbox=dict(boxstyle='round,pad=0.45', facecolor='#f5f7fa',
                      edgecolor='#cccccc', linewidth=0.8),
        )

    legend_patches = [
        mpatches.Patch(color='#dce3ed', label='Ocean / Non-Forest'),
        mpatches.Patch(color='#b8d4ba', label='Surviving Forest (faded backdrop)'),
        mpatches.Patch(color='#ffb300', label='Historic Loss 2000–24 (19,502 ha)'),
        mpatches.Patch(color='#e53935', label='New Predicted Loss 2024–60'),
    ]
    fig2.legend(handles=legend_patches, loc='lower center',
                bbox_to_anchor=(0.5, 0.02), ncol=4,
                facecolor='white', edgecolor='#cccccc',
                labelcolor='#1a1a2e', fontsize=11)

    fig2.suptitle(
        'Three-Scenario Predictive Forecast 2024–2060\n'
        'Andaman & Nicobar Islands — Forest-Collapse Scenarios',
        color='#1a1a2e', fontsize=15, fontweight='bold', y=0.97,
    )

    fig2.tight_layout(rect=[0, 0.08, 1, 0.93])
    fig_tri = FIG_DIR / 'predictive' / 'forecast_tri_scenario_2060.png'
    fig2.savefig(fig_tri, dpi=250, facecolor='white', bbox_inches='tight')
    plt.close()
    print(f"  ✅  Figure saved → {fig_tri.name}")

    print(f"\n{SEP}")
    print("  ✅  Multi-Scenario Forecasting Complete!")
    print(f"      CSV     : results/economic_scenarios_2024_2060.csv")
    print(f"      Figures : figures/predictive/")
    print(SEP + "\n")
