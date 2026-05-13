"""
ANI Ecosystem Services — Statistical Distribution Visualizations
================================================================
Generates publication-quality violin plots of Habitat Quality and
Soil Erosion distributions stratified by ESA WorldCover land class.

Inputs  : data/processed/ANI_ESA_WorldCover_mosaic_clipped.tif
          results/habitat_quality_index.tif
          results/rusle_soil_loss.tif
Outputs : figures/synthesis/stat_distribution_habitat_quality.png
          figures/synthesis/stat_distribution_soil_erosion.png

Run with: venv/bin/python src/visualize_statistics.py
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

# ── Print Separator ────────────────────────────────────────────────────
SEP = '=' * 60

# ── ESA WorldCover Class Labels ────────────────────────────────────────
ESA_CLASS_LABELS = {
    10: 'Tree Cover',
    20: 'Shrubland',
    30: 'Grassland',
    40: 'Cropland',
    50: 'Built-up',
    60: 'Bare / Sparse',
    80: 'Open water',
    90: 'Herbaceous wetland',
    95: 'Mangroves',
}

# ── Sampling and Display Constants ────────────────────────────────────
SAMPLE_SIZE_N      = 150_000   # Max pixels sampled (prevents KDE OOM)
RUSLE_DISPLAY_CAP  = 350.0     # Clip extreme erosion tail for violin scale
TARGET_CLASSES     = [
    'Tree Cover', 'Mangroves', 'Cropland',
    'Built-up', 'Bare / Sparse', 'Grassland',
]

# ── Per-class Colour Palette ───────────────────────────────────────────
CLASS_COLOR_MAP = {
    'Tree Cover':   '#2e7d32',
    'Mangroves':    '#00b0ff',
    'Cropland':     '#ffb300',
    'Built-up':     '#e53935',
    'Bare / Sparse':'#8d6e63',
    'Grassland':    '#aed581',
}


# ══════════════════════════════════════════════════════════════════════
# UTILITY
# ══════════════════════════════════════════════════════════════════════
def load_raster_band(file_path: Path):
    """Load band 1 of a GeoTIFF as a raw numpy array.

    Returns None if the file does not exist.
    """
    if not file_path.exists():
        return None
    with rasterio.open(file_path) as src:
        return src.read(1)


# ══════════════════════════════════════════════════════════════════════
# ENTRY POINT
# ══════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    print(SEP)
    print("  ANI Ecosystem Services — Statistical Distribution Plots")
    print(SEP)

    # 1. Load raster data ──────────────────────────────────────────────
    print("  Loading raster arrays …")
    landcover_raw = load_raster_band(PROC_DIR / 'ANI_ESA_WorldCover_mosaic_clipped.tif')
    hab_quality   = load_raster_band(RES_DIR  / 'habitat_quality_index.tif')
    soil_erosion  = load_raster_band(RES_DIR  / 'rusle_soil_loss.tif')

    if landcover_raw is None or hab_quality is None or soil_erosion is None:
        print("  ❌  Missing required raster layers — aborting.")
        exit(1)

    # 2. Build valid-pixel mask and extract flat arrays ─────────────────
    print("  Aligning masks and extracting valid pixel pools …")
    valid_mask   = (
        (landcover_raw > 0)
        & (~np.isnan(hab_quality))
        & (~np.isnan(soil_erosion))
        & (soil_erosion >= 0)
    )

    lc_valid  = landcover_raw[valid_mask].astype(int)
    hab_valid = hab_quality[valid_mask]
    rus_valid = soil_erosion[valid_mask]
    n_valid   = len(lc_valid)
    print(f"  Total valid pixels : {n_valid:,}")

    # 3. Stratified random sample ──────────────────────────────────────
    np.random.seed(42)
    if n_valid > SAMPLE_SIZE_N:
        sample_idx = np.random.choice(n_valid, size=SAMPLE_SIZE_N, replace=False)
        lc_sample  = lc_valid[sample_idx]
        hab_sample = hab_valid[sample_idx]
        rus_sample = rus_valid[sample_idx]
    else:
        lc_sample  = lc_valid
        hab_sample = hab_valid
        rus_sample = rus_valid

    # 4. Build labelled dataframe ──────────────────────────────────────
    class_labels = np.array([ESA_CLASS_LABELS.get(k, 'Other') for k in lc_sample])
    class_filter = np.isin(class_labels, TARGET_CLASSES)

    plot_df = pd.DataFrame({
        'Land Cover':           class_labels[class_filter],
        'Habitat Quality':      hab_sample[class_filter],
        'Soil Erosion (t/ha/yr)': rus_sample[class_filter],
    })
    plot_df['Land Cover'] = pd.Categorical(
        plot_df['Land Cover'], categories=TARGET_CLASSES, ordered=True
    )
    print(f"  Plotting dataframe size : {len(plot_df):,} rows")

    # 5. Configure seaborn theme ───────────────────────────────────────
    sns.set_theme(style="whitegrid", rc={
        "axes.facecolor":  "#14142a",
        "figure.facecolor":"#0d0d1a",
        "axes.edgecolor":  "#555",
        "grid.color":      "#2c2c3e",
        "text.color":      "white",
        "axes.labelcolor": "white",
        "xtick.color":     "white",
        "ytick.color":     "white",
    })

    # ── PLOT 1: Habitat Quality Violin ────────────────────────────────
    print(f"\n{SEP}")
    print("  Rendering Habitat Quality Distribution …")
    print(f"{SEP}")

    fig1, ax1 = plt.subplots(figsize=(10, 6))
    sns.violinplot(
        data=plot_df, x='Land Cover', y='Habitat Quality',
        palette=CLASS_COLOR_MAP, ax=ax1,
        inner="box", linewidth=1.5, cut=0,
    )
    ax1.set_title(
        "Statistical Distribution of Habitat Quality by Land Class\n"
        "Andaman & Nicobar Islands (2024)",
        fontsize=14, fontweight='bold', pad=15,
    )
    ax1.set_ylabel("Habitat Quality Index (Q)", fontsize=11)
    ax1.set_xlabel("")
    plt.xticks(rotation=30)
    fig1.tight_layout()

    out_hab = FIG_DIR / 'synthesis' / 'stat_distribution_habitat_quality.png'
    fig1.savefig(out_hab, dpi=200)
    plt.close()
    print(f"  ✅  Figure saved → {out_hab.name}")

    # ── PLOT 2: Soil Erosion Violin ───────────────────────────────────
    print(f"\n{SEP}")
    print("  Rendering Soil Erosion Distribution …")
    print(f"{SEP}")

    # Trim extreme upper tail to prevent scale collapse
    erosion_df = plot_df[plot_df['Soil Erosion (t/ha/yr)'] < RUSLE_DISPLAY_CAP].copy()

    fig2, ax2 = plt.subplots(figsize=(10, 6))
    sns.violinplot(
        data=erosion_df, x='Land Cover', y='Soil Erosion (t/ha/yr)',
        palette=CLASS_COLOR_MAP, ax=ax2,
        inner="box", linewidth=1.5, cut=0,
    )
    ax2.set_title(
        "Statistical Density of Soil Erosion Vulnerability by Land Class\n"
        "Andaman & Nicobar Islands (2024)",
        fontsize=14, fontweight='bold', pad=15,
    )
    ax2.set_ylabel("Annual Soil Loss (tonnes / ha / year)", fontsize=11)
    ax2.set_xlabel("")
    plt.xticks(rotation=30)
    fig2.tight_layout()

    out_soil = FIG_DIR / 'synthesis' / 'stat_distribution_soil_erosion.png'
    fig2.savefig(out_soil, dpi=200)
    plt.close()
    print(f"  ✅  Figure saved → {out_soil.name}")

    print(f"\n{SEP}")
    print("  ✅  Statistical Distribution Visualizations Complete!")
    print(SEP + "\n")
