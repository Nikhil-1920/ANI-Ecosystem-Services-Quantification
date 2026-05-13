"""
ANI Ecosystem Services — Habitat Quality Sensitivity Analysis
=============================================================
Validates the robustness of the InVEST Habitat Quality model by
iterating over a sweep of half-saturation constants (k) and
demonstrating that relative degradation rankings remain stable.

    Q_x = H_x × [1 − D_x^z / (D_x^z + k^z)]

Inputs  : data/processed/ANI_ESA_WorldCover_mosaic_clipped.tif
          (threat rasters derived internally via habitat_quality.py)
Outputs : results/habitat_sensitivity_analysis.csv
          figures/habitat/habitat_sensitivity_analysis.png

Run with: venv/bin/python src/sensitivity_analysis.py
"""

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from pathlib import Path

# Import core InVEST habitat functions from sibling module
sys.path.insert(0, str(Path(__file__).parent))
from habitat_quality import (
    load_raster,
    build_habitat_map,
    build_threat_rasters,
    compute_degradation,
    Z_SCALE,
)

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

# ── Sensitivity Sweep Parameters ──────────────────────────────────────
K_SWEEP_VALUES  = [0.1, 0.25, 0.5, 0.75, 0.9]
TARGET_CLASSES  = ['Tree Cover', 'Mangroves', 'Grassland',
                   'Cropland', 'Built-up', 'Bare / Sparse']

# ── Per-class Colour Palette ───────────────────────────────────────────
CLASS_COLOR_MAP = {
    'Tree Cover':   '#2e7d32',
    'Mangroves':    '#00b0ff',
    'Grassland':    '#aed581',
    'Cropland':     '#ffb300',
    'Built-up':     '#e53935',
    'Bare / Sparse':'#8d6e63',
}


# ══════════════════════════════════════════════════════════════════════
# ENTRY POINT
# ══════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    print(SEP)
    print("  ANI Ecosystem Services — Habitat Quality Sensitivity Analysis")
    print(SEP)

    # 1. Load physical environment─────────────────────────────────────
    landcover_grid, lc_profile = load_raster('ANI_ESA_WorldCover_mosaic_clipped.tif')
    if landcover_grid is None:
        exit(1)

    print("\n  Computing baseline model parameters (≈3 min) …")
    habitat_score  = build_habitat_map(landcover_grid, lc_profile)
    threat_rasters = build_threat_rasters(landcover_grid, lc_profile)
    degradation_d  = compute_degradation(landcover_grid, threat_rasters, lc_profile)

    lc_int         = np.nan_to_num(landcover_grid, nan=0).astype(int)
    valid_pixels   = ~np.isnan(landcover_grid)

    # Pre-compute D^z once; only k changes per iteration
    degradation_dz = np.power(degradation_d, Z_SCALE)

    # 2. Sweep half-saturation constant k ─────────────────────────────
    print(f"\n{SEP}")
    print(f"  SWEEPING k ∈ {K_SWEEP_VALUES}")
    print(f"{SEP}")

    sensitivity_records = []

    for k_val in K_SWEEP_VALUES:
        print(f"  Evaluating k = {k_val} …")
        k_z  = k_val ** Z_SCALE

        # Q_x = H_x × [1 − D_x^z / (D_x^z + k^z)]
        q_grid = habitat_score * (1.0 - degradation_dz / (degradation_dz + k_z))
        q_grid = np.clip(q_grid, 0, 1)

        for class_code, class_label in ESA_CLASS_LABELS.items():
            class_mask = valid_pixels & (lc_int == class_code)
            if class_mask.sum() > 0:
                sensitivity_records.append({
                    'k_value':             k_val,
                    'land_cover':          class_label,
                    'mean_habitat_quality': round(float(np.nanmean(q_grid[class_mask])), 4),
                })

    sensitivity_df = pd.DataFrame(sensitivity_records)
    sensitivity_df = sensitivity_df[sensitivity_df['land_cover'].isin(TARGET_CLASSES)]

    csv_out = RES_DIR / 'habitat_sensitivity_analysis.csv'
    sensitivity_df.to_csv(csv_out, index=False)
    print(f"  ✅  CSV saved → {csv_out.name}")

    # 3. Render sensitivity curve figure ──────────────────────────────
    print(f"\n{SEP}")
    print("  RENDERING SENSITIVITY CURVES")
    print(f"{SEP}")

    sns.set_theme(style="whitegrid")

    # Order classes by final (largest-k) mean Q so legend matches visual rank
    final_k = max(K_SWEEP_VALUES)
    rank_order = (
        sensitivity_df[sensitivity_df['k_value'] == final_k]
        .sort_values('mean_habitat_quality', ascending=False)['land_cover']
        .tolist()
    )

    fig, ax = plt.subplots(figsize=(10, 6), facecolor='white')
    ax.set_facecolor('white')
    sns.lineplot(
        data=sensitivity_df,
        x='k_value', y='mean_habitat_quality',
        hue='land_cover', hue_order=rank_order,
        marker='o', linewidth=2.5, markersize=8,
        palette=CLASS_COLOR_MAP, ax=ax,
    )

    ax.set_title(
        "Habitat Quality Model Sensitivity to Parameter $k$\n"
        "Structural Robustness Validation — Andaman & Nicobar Islands",
        fontsize=14, fontweight='bold', color='#0d1b2a', pad=15,
    )
    ax.set_xlabel("Half-Saturation Constant ($k$)", fontsize=12, color='#0d1b2a')
    ax.set_ylabel("Mean Habitat Quality Index ($Q_x$)", fontsize=12, color='#0d1b2a')
    ax.tick_params(colors='#0d1b2a', which='both')
    for spine_pos in ['top', 'bottom', 'left', 'right']:
        ax.spines[spine_pos].set_color('#888')

    # End-of-curve value labels for quick reading
    for cls in rank_order:
        last = sensitivity_df[
            (sensitivity_df['land_cover'] == cls)
            & (sensitivity_df['k_value'] == final_k)
        ]['mean_habitat_quality'].iloc[0]
        ax.annotate(
            f"{last:.2f}",
            xy=(final_k, last),
            xytext=(6, 0), textcoords='offset points',
            color=CLASS_COLOR_MAP[cls], fontsize=9, fontweight='bold',
            va='center',
        )

    legend = ax.legend(
        loc='center left', bbox_to_anchor=(1.06, 0.5),
        title='Land Cover', frameon=True,
    )
    legend.get_title().set_color('#0d1b2a')
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_edgecolor('#888')
    for text in legend.get_texts():
        text.set_color('#0d1b2a')

    fig.tight_layout()
    fig_out = FIG_DIR / 'habitat' / 'habitat_sensitivity_analysis.png'
    fig.savefig(fig_out, dpi=200, facecolor='white', bbox_inches='tight')
    plt.close()
    print(f"  ✅  Figure saved → {fig_out.name}")

    print(f"\n{SEP}")
    print("  ✅  Habitat Sensitivity Analysis Complete!")
    print(f"      CSV    : results/habitat_sensitivity_analysis.csv")
    print(f"      Figure : figures/habitat/habitat_sensitivity_analysis.png")
    print(SEP + "\n")
