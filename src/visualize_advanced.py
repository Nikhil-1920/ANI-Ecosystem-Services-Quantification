"""
ANI Ecosystem Services — Advanced Synthesis Visualizations
===========================================================
Generates three publication-quality multi-dimensional figures:
  1. Trade-off Radar Chart  — Carbon, Habitat, Soil by land cover class
  2. Bivariate KDE Plot     — Habitat quality vs. soil erosion correlation
  3. Hexbin Density Map     — Spatial clustering of ECI triple-collapse zones

Inputs  : data/processed/ANI_ESA_WorldCover_mosaic_clipped.tif
          data/processed/ANI_GEDI_Biomass_Density_clipped.tif
          results/habitat_quality_index.tif
          results/rusle_soil_loss.tif
          results/eci_collapse_hotspots.tif
Outputs : figures/synthesis/tradeoff_radar_chart.png
          figures/synthesis/bivariate_habitat_erosion_kde.png
          figures/synthesis/hotspot_hexbin_density_map.png

Run with: venv/bin/python src/visualize_advanced.py
"""

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import rasterio
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# ── Directory Paths ────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).parent
PROC_DIR   = SCRIPT_DIR.parent / 'data' / 'processed'
RES_DIR    = SCRIPT_DIR.parent / 'results'
FIG_DIR    = SCRIPT_DIR.parent / 'figures'

# ── ESA WorldCover Class Labels ────────────────────────────────────────
ESA_CLASS_LABELS = {
    10: 'Tree Cover',
    20: 'Shrubland',
    40: 'Cropland',
    50: 'Built-up',
    60: 'Bare / Sparse',
}

# ── Carbon Conversion Constants ────────────────────────────────────────
ROOT_SHOOT_R   = 1.24    # Root-to-shoot multiplier applied to AGB (IPCC)
IPCC_CARBON_F  = 0.47    # Biomass to carbon fraction

# ── Random Sample Size for KDE speed ──────────────────────────────────
KDE_SAMPLE_N   = 100_000


# ══════════════════════════════════════════════════════════════════════
# UTILITY
# ══════════════════════════════════════════════════════════════════════
def load_raster_as_float(file_path: Path):
    """Load a GeoTIFF band as float64 with nodata replaced by nan.

    Returns data_array or None if the file does not exist.
    """
    if not file_path.exists():
        return None
    with rasterio.open(file_path) as src:
        data_array = src.read(1).astype(float)
        nodata_val = src.nodata
    if nodata_val is not None:
        data_array = np.where(data_array == nodata_val, np.nan, data_array)
    return data_array


# ══════════════════════════════════════════════════════════════════════
# ENTRY POINT
# ══════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    print("=" * 60)
    print("  ANI Ecosystem Services — Advanced Synthesis Visualizations")
    print("=" * 60)

    # 1. Load all spatial data layers ──────────────────────────────────
    print("  Loading spatial data layers …")
    landcover_grid = load_raster_as_float(PROC_DIR / 'ANI_ESA_WorldCover_mosaic_clipped.tif')
    hab_quality    = load_raster_as_float(RES_DIR  / 'habitat_quality_index.tif')
    soil_erosion   = load_raster_as_float(RES_DIR  / 'rusle_soil_loss.tif')
    agb_density    = load_raster_as_float(PROC_DIR / 'ANI_GEDI_Biomass_Density_clipped.tif')
    eci_hotspots   = load_raster_as_float(RES_DIR  / 'eci_collapse_hotspots.tif')

    if landcover_grid is None or hab_quality is None or soil_erosion is None:
        print("  ❌  Core raster files missing.")
        exit(1)

    # 2. Build aligned sample dataframe ────────────────────────────────
    print("  Extracting valid pixel sample …")
    valid_pixels   = (
        (landcover_grid > 0)
        & (~np.isnan(hab_quality))
        & (~np.isnan(soil_erosion))
        & (~np.isnan(agb_density))
    )

    lc_flat  = landcover_grid[valid_pixels].astype(int)
    hab_flat = hab_quality[valid_pixels]
    rus_flat = soil_erosion[valid_pixels]
    bio_flat = agb_density[valid_pixels] * ROOT_SHOOT_R * IPCC_CARBON_F   # → MgC/ha

    # Stratified random sample for KDE computation speed
    if len(lc_flat) > KDE_SAMPLE_N:
        sample_idx = np.random.choice(len(lc_flat), size=KDE_SAMPLE_N, replace=False)
        lc_flat    = lc_flat[sample_idx]
        hab_flat   = hab_flat[sample_idx]
        rus_flat   = rus_flat[sample_idx]
        bio_flat   = bio_flat[sample_idx]

    class_labels = np.array([ESA_CLASS_LABELS.get(k, 'Other') for k in lc_flat])
    sample_df = pd.DataFrame({
        'Land Cover':    class_labels,
        'Habitat Quality': hab_flat,
        'Soil Erosion':    rus_flat,
        'Carbon Storage':  bio_flat,
    })
    # Keep only visible land cover classes
    sample_df = sample_df[sample_df['Land Cover'].isin(
        ['Tree Cover', 'Cropland', 'Built-up', 'Bare / Sparse']
    )]

    # ── PLOT 1: Trade-off Radar Chart ─────────────────────────────────
    print("\n  Rendering Plot 1: Ecosystem Services Trade-Off Radar …")

    class_means  = sample_df.groupby('Land Cover')[
        ['Habitat Quality', 'Soil Erosion', 'Carbon Storage']
    ].mean()

    # Normalise all metrics to 0–100 % scale
    hab_pct      = class_means['Habitat Quality'] * 100
    carbon_pct   = np.clip((class_means['Carbon Storage'] / 120.0) * 100, 0, 100)
    # Soil Retention = inverse of erosion (0 erosion → 100 %; > 200 → 0 %)
    retention_pct = np.clip(100 - (class_means['Soil Erosion'] / 2.0), 0, 100)

    radar_labels = ['Habitat\nQuality', 'Carbon\nStorage', 'Soil\nRetention']
    n_axes       = len(radar_labels)
    axis_angles  = np.linspace(0, 2 * np.pi, n_axes, endpoint=False).tolist()
    axis_angles  += axis_angles[:1]   # close the loop

    fig1, ax1 = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True),
                             facecolor='#0d0d1a')
    ax1.set_facecolor('#0d0d1a')

    class_colors = {
        'Tree Cover': '#2e7d32',
        'Cropland':   '#ffb300',
        'Built-up':   '#e53935',
    }

    for class_name in ['Tree Cover', 'Cropland', 'Built-up']:
        if class_name not in class_means.index:
            continue
        radar_values  = [hab_pct[class_name], carbon_pct[class_name], retention_pct[class_name]]
        radar_values += radar_values[:1]
        ax1.plot(axis_angles, radar_values,
                 color=class_colors[class_name], linewidth=2, label=class_name)
        ax1.fill(axis_angles, radar_values,
                 color=class_colors[class_name], alpha=0.25)

    ax1.set_theta_offset(np.pi / 2)
    ax1.set_theta_direction(-1)
    ax1.set_thetagrids(np.degrees(axis_angles[:-1]), radar_labels,
                       color='white', fontsize=12, fontweight='bold')
    ax1.set_rlabel_position(180)
    ax1.set_yticks([20, 40, 60, 80, 100])
    ax1.set_yticklabels(['20', '40', '60', '80', '100'], color='#aaa')
    ax1.grid(color='#333')
    ax1.spines['polar'].set_color('#555')

    plt.title('Multidimensional Ecosystem Services Trade-Offs',
              color='white', pad=20, fontsize=14, fontweight='bold')
    legend = plt.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1),
                         facecolor='#14142a')
    for text in legend.get_texts():
        text.set_color('white')

    fig1.tight_layout()
    out_radar = FIG_DIR / 'synthesis' / 'tradeoff_radar_chart.png'
    fig1.savefig(out_radar, dpi=200)
    plt.close()
    print(f"    ✅  Saved → {out_radar.name}")

    # ── PLOT 2: Bivariate KDE — Habitat Quality vs. Soil Erosion ──────
    print("  Rendering Plot 2: Bivariate KDE Correlation …")

    sns.set_theme(style="white", rc={
        "axes.facecolor":   "#0d0d1a",
        "figure.facecolor": "#0d0d1a",
        "text.color":       "white",
        "axes.labelcolor":  "white",
        "xtick.color":      "white",
        "ytick.color":      "white",
    })

    # Restrict to forested and cropped land, trim extreme erosion tail
    kde_df = sample_df[sample_df['Land Cover'].isin(['Tree Cover', 'Cropland'])]
    kde_df = kde_df[kde_df['Soil Erosion'] < 300]

    kde_plot = sns.JointGrid(
        data=kde_df, x='Habitat Quality', y='Soil Erosion',
        hue='Land Cover',
        palette={'Tree Cover': '#2e7d32', 'Cropland': '#ffb300'},
        space=0.1, height=8, ratio=4,
    )
    kde_plot.plot_joint(sns.kdeplot,     fill=False, thresh=0.01, alpha=0.8,
                        linewidths=1.5, levels=5)
    kde_plot.plot_joint(sns.scatterplot, s=5, alpha=0.2)
    kde_plot.plot_marginals(sns.kdeplot, fill=True, linewidth=2, alpha=0.4)

    kde_plot.ax_joint.set_facecolor('#14142a')
    for panel_ax in [kde_plot.ax_joint, kde_plot.ax_marg_x, kde_plot.ax_marg_y]:
        panel_ax.spines['bottom'].set_color('#555')
        panel_ax.spines['left'].set_color('#555')
        panel_ax.spines['right'].set_visible(False)
        panel_ax.spines['top'].set_visible(False)
        panel_ax.tick_params(colors='white')

    kde_plot.figure.suptitle('Bivariate Correlation: Habitat Quality vs. Soil Erosion',
                             y=1.02, fontsize=15, fontweight='bold', color='white')

    out_kde = FIG_DIR / 'synthesis' / 'bivariate_habitat_erosion_kde.png'
    kde_plot.savefig(out_kde, dpi=200, facecolor='#0d0d1a')
    plt.close()
    print(f"    ✅  Saved → {out_kde.name}")

    # ── PLOT 3: Hexbin Spatial Density of Hotspot Pixels ──────────────
    print("  Rendering Plot 3: Hexbin Hotspot Density Map …")

    if eci_hotspots is not None:
        # Extract row/col coordinates of extreme collapse pixels (ECI > 0.8)
        hot_rows, hot_cols = np.where(eci_hotspots > 0.8)
        hot_rows_inv       = -hot_rows      # invert Y so north is up

        fig3, ax3 = plt.subplots(figsize=(10, 8), facecolor='#0d0d1a')
        ax3.set_facecolor('#0d0d1a')
        hb = ax3.hexbin(hot_cols, hot_rows_inv,
                        gridsize=40, cmap='inferno',
                        mincnt=1, bins='log', edgecolors='none')
        ax3.set_aspect('equal')
        ax3.axis('off')

        cb = plt.colorbar(hb, ax=ax3, fraction=0.03, pad=0.04)
        cb.set_label('Hotspot Macro-Density (Log Scale)', color='white')
        cb.ax.yaxis.set_tick_params(color='white')
        plt.setp(cb.ax.yaxis.get_ticklabels(), color='white')

        plt.title('Spatial Hexbin Clustering of Triple-Collapse Hotspots\n'
                  'Andaman & Nicobar Region',
                  color='white', fontsize=14, fontweight='bold', pad=15)

        out_hexbin = FIG_DIR / 'synthesis' / 'hotspot_hexbin_density_map.png'
        fig3.savefig(out_hexbin, dpi=200, bbox_inches='tight')
        plt.close()
        print(f"    ✅  Saved → {out_hexbin.name}")
    else:
        print("    ⚠️  eci_collapse_hotspots.tif not found — skipping hexbin plot.")

    print("\n" + "=" * 60)
    print("  ✅  Advanced Synthesis Visualizations Complete!")
    print("=" * 60 + "\n")
