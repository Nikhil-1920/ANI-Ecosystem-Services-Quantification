"""
ANI Ecosystem Services — Week 7: Habitat Quality Assessment
============================================================
Implements an InVEST-equivalent Habitat Quality model using
pure NumPy (no InVEST installation required).

Model:
    Threat rasters → per-pixel degradation D_x
    D_x + habitat sensitivity H_x → quality Q_x ∈ [0, 1]

Threat layers  : OSM Roads (rasterized) + ESA built-up (class 50)
                 + ESA cropland (class 40)
Habitat scores : derived from ESA WorldCover land cover classes
Half-saturation: k = 0.5 (standard InVEST default)

Inputs  : data/processed/  (clipped GeoTIFFs)
          data/processed/ANI_OSM_Roads_32646.shp
Outputs : results/habitat_quality_index.tif
          results/habitat_quality_delta.tif
          results/habitat_quality_by_landcover.csv
          figures/habitat/habitat_quality_index_map.png
          figures/habitat/habitat_quality_by_landcover.png
          figures/habitat/habitat_quality_delta_map.png

Run with: venv/bin/python src/habitat_quality.py
"""

import os
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import rasterio
from rasterio.features import rasterize
from rasterio.transform import from_bounds
import geopandas as gpd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import pandas as pd
from pathlib import Path
from scipy.ndimage import (distance_transform_edt, binary_erosion,
                            binary_dilation, label)

# ── Directory Paths ────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).parent
PROC_DIR   = SCRIPT_DIR.parent / 'data' / 'processed'
RES_DIR    = SCRIPT_DIR.parent / 'results'
FIG_DIR    = SCRIPT_DIR.parent / 'figures'
RES_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ── Print Separator ────────────────────────────────────────────────────
SEP = '=' * 60

# ── ESA WorldCover classes → habitat sensitivity ───────────────────────
# Scale 0 (non-habitat) to 1.0 (pristine primary habitat)
HABITAT_SENSITIVITY = {
    10:  1.00,   # Tree cover (tropical evergreen)
    20:  0.50,   # Shrubland
    30:  0.40,   # Grassland
    40:  0.10,   # Cropland — low habitat value
    50:  0.00,   # Built-up — no habitat
    60:  0.15,   # Bare / sparse vegetation
    70:  0.00,   # Snow/ice
    80:  0.20,   # Open water
    90:  0.45,   # Herbaceous wetland
    95:  0.90,   # Mangroves — high value
    100: 0.30,   # Moss/Lichen
}

# ── Threat parameters ──────────────────────────────────────────────────
# Each threat: weight, max effect distance (pixels @30m), decay type
# Distance units: number of 30m pixels
# 500m ≈ 17 pixels;  2km ≈ 67 pixels;  5km ≈ 167 pixels

THREATS = {
    'roads': {
        'weight':    0.7,
        'max_dist':  100,     # 3 km effect radius
        'decay':     'linear',
        'sensitivity': {      # habitat sensitivity to THIS threat (by ESA class)
            10: 0.85,  20: 0.60,  30: 0.50,  40: 0.10,
            50: 0.00,  60: 0.20,  80: 0.30,  90: 0.55, 95: 0.80,
        },
    },
    'builtup': {
        'weight':    0.9,
        'max_dist':  50,      # 1.5 km effect radius
        'decay':     'exponential',
        'sensitivity': {
            10: 0.90,  20: 0.65,  30: 0.55,  40: 0.15,
            50: 0.00,  60: 0.25,  80: 0.35,  90: 0.60, 95: 0.85,
        },
    },
    'cropland': {
        'weight':    0.5,
        'max_dist':  25,      # 750 m edge effect
        'decay':     'linear',
        'sensitivity': {
            10: 0.60,  20: 0.40,  30: 0.30,  40: 0.05,
            50: 0.00,  60: 0.15,  80: 0.25,  90: 0.35, 95: 0.55,
        },
    },
}

# ── InVEST constants ───────────────────────────────────────────────────
HALF_SAT_K = 0.5    # k in Q = H*(1 - D^z/(D^z + k^z))
Z_SCALE    = 2.5    # z (scaling parameter, default InVEST value)


# ── Split-map row boundaries (Ten Degree Channel) ────────────────────
# Andaman group: rows 0 → ANDAMAN_ROW_END
# Nicobar group: rows NICOBAR_ROW_START → end
ANDAMAN_ROW_END   = 11750
NICOBAR_ROW_START = 16200
_PAD = 200   # pixel padding around each extent


def _split_extent(arr):
    """Return (andaman_crop, nicobar_crop) slices of a 2-D masked/ndarray."""
    h = arr.shape[0]
    r0_a = 0
    r1_a = min(ANDAMAN_ROW_END + _PAD, h)
    r0_n = max(0, NICOBAR_ROW_START - _PAD)
    r1_n = h
    return arr[r0_a:r1_a, :], arr[r0_n:r1_n, :]


def _tight_crop(arr):
    """Crop 2-D array to bounding box of non-NaN values + small pad."""
    if np.ma.is_masked(arr):
        valid = ~arr.mask
    else:
        valid = ~np.isnan(arr)
    if not valid.any():
        return arr
    rows = np.any(valid, axis=1)
    cols = np.any(valid, axis=0)
    r0, r1 = np.where(rows)[0][[0, -1]]
    c0, c1 = np.where(cols)[0][[0, -1]]
    pad = 60
    r0 = max(0, r0 - pad); r1 = min(arr.shape[0], r1 + pad)
    c0 = max(0, c0 - pad); c1 = min(arr.shape[1], c1 + pad)
    return arr[r0:r1, c0:c1]


def split_imshow(fig, gs_left, gs_right, data, cmap, vmin, vmax,
                 label_left='Andaman Islands', label_right='Nicobar Islands',
                 interpolation='nearest', aspect='equal'):
    """Plot a raster split into Andaman (left) and Nicobar (right) sub-axes."""
    and_raw, nic_raw = _split_extent(data)
    and_crop = _tight_crop(and_raw)
    nic_crop = _tight_crop(nic_raw)

    ax_a = fig.add_subplot(gs_left)
    ax_n = fig.add_subplot(gs_right)
    for ax, crop, title in [(ax_a, and_crop, label_left),
                             (ax_n, nic_crop, label_right)]:
        ax.set_facecolor('white')
        ax.imshow(crop, cmap=cmap, vmin=vmin, vmax=vmax,
                  interpolation=interpolation, aspect=aspect)
        ax.set_title(title, color='#1a1a2e', fontsize=11, fontweight='bold', pad=8)
        ax.axis('off')
    return ax_a, ax_n


# ── Helper: load raster ────────────────────────────────────────────────
def load_raster(name):
    path = PROC_DIR / name
    if not path.exists():
        print(f"  ⚠️  Not found: {name}")
        return None, None
    with rasterio.open(path) as src:
        data    = src.read(1).astype(float)
        profile = src.profile.copy()
        nodata  = src.nodata
    if nodata is not None:
        data = np.where(data == nodata, np.nan, data)
    data = np.where(data < -1e9, np.nan, data)
    return data, profile


# ══════════════════════════════════════════════════════════════════════
# PART 1 — Build Habitat Sensitivity Map (H_x)
# ══════════════════════════════════════════════════════════════════════
def build_habitat_map(esa, profile):
    print(f"\n{SEP}")
    print("PART 1 — Building Habitat Sensitivity Map (H_x)")
    print(SEP)

    H = np.zeros_like(esa, dtype=float)
    valid = ~np.isnan(esa)
    for cls_code, h_val in HABITAT_SENSITIVITY.items():
        mask = valid & (np.nan_to_num(esa, nan=-1).astype(int) == cls_code)
        H[mask] = h_val
    # Pixels with no ESA class = 0 sensitivity
    H[~valid] = np.nan

    unique, counts = np.unique(np.nan_to_num(esa[valid], nan=-1).astype(int), return_counts=True)
    print(f"  ESA classes present: {list(zip(unique, counts))}")
    print(f"  Mean habitat sensitivity: {np.nanmean(H):.4f}")
    print(f"  High-quality pixels (H≥0.8): {np.sum(H >= 0.8):,}")
    return H


# ══════════════════════════════════════════════════════════════════════
# PART 2 — Build Threat Rasters
# ══════════════════════════════════════════════════════════════════════
def build_threat_rasters(esa, profile):
    print(f"\n{SEP}")
    print("PART 2 — Building Threat Rasters")
    print(SEP)

    valid = ~np.isnan(esa)
    esa_int = np.nan_to_num(esa, nan=0).astype(int)
    threats_out = {}

    # --- Threat 1: Built-up (ESA class 50) ---
    builtup_mask = (esa_int == 50) & valid
    print(f"  Built-up pixels : {builtup_mask.sum():,}")
    threats_out['builtup'] = builtup_mask.astype(float)

    # --- Threat 2: Cropland (ESA class 40) ---
    cropland_mask = (esa_int == 40) & valid
    print(f"  Cropland pixels : {cropland_mask.sum():,}")
    threats_out['cropland'] = cropland_mask.astype(float)

    # --- Threat 3: Roads (rasterize OSM shapefile) ---
    roads_path = PROC_DIR / 'ANI_OSM_Roads_32646.shp'
    if roads_path.exists():
        print(f"  Rasterizing OSM roads from {roads_path.name}...")
        gdf = gpd.read_file(roads_path)
        # Rasterize road geometries onto the ESA grid
        transform = profile['transform']
        height    = profile['height']
        width     = profile['width']
        road_shapes = [(geom, 1) for geom in gdf.geometry if geom is not None]
        if road_shapes:
            road_raster = rasterize(
                road_shapes,
                out_shape=(height, width),
                transform=transform,
                fill=0,
                dtype=np.uint8,
                all_touched=True
            ).astype(float)
        else:
            road_raster = np.zeros((height, width), dtype=float)
        print(f"  Road pixels rasterized: {(road_raster > 0).sum():,}")
        threats_out['roads'] = road_raster
    else:
        print("  ⚠️  OSM Roads shapefile not found — using GFW loss as roads proxy")
        # Fallback: GFW loss pixels as proxy threat source
        gfw, _ = load_raster('ANI_GFW_Forest_Loss_2001_2023_clipped.tif')
        threats_out['roads'] = np.where(~np.isnan(gfw) & (gfw > 0), 1.0, 0.0)

    return threats_out


# ══════════════════════════════════════════════════════════════════════
# PART 3 — Compute Degradation Score (D_x) with Distance Decay
# ══════════════════════════════════════════════════════════════════════
def compute_degradation(esa, threat_rasters, profile):
    print(f"\n{SEP}")
    print("PART 3 — Computing Degradation Index D_x")
    print(SEP)

    esa_int  = np.nan_to_num(esa, nan=0).astype(int)
    valid    = ~np.isnan(esa)

    total_weight = sum(T['weight'] for T in THREATS.values())
    D_total      = np.zeros_like(esa, dtype=float)

    for threat_name, T in THREATS.items():
        threat_arr = threat_rasters.get(threat_name)
        if threat_arr is None:
            print(f"  ⚠️  No raster for threat '{threat_name}' — skipping.")
            continue

        print(f"  Processing threat: {threat_name} (weight={T['weight']}, "
              f"max_dist={T['max_dist']} px, decay={T['decay']})")

        # Distance-from-threat (in pixels) via EDT on the inverted mask
        no_threat = (threat_arr == 0).astype(np.uint8)
        dist_px   = distance_transform_edt(no_threat).astype(float)

        # Decay function: linear or exponential
        max_d = T['max_dist']
        if T['decay'] == 'linear':
            decay_factor = np.clip(1.0 - dist_px / max_d, 0, 1)
        else:  # exponential
            decay_factor = np.exp(-2.99 * dist_px / max_d)
            decay_factor = np.where(dist_px > max_d, 0.0, decay_factor)

        # Per-pixel habitat sensitivity to this specific threat
        sens_map = np.zeros_like(esa, dtype=float)
        for cls_code, s_val in T['sensitivity'].items():
            sens_map[esa_int == cls_code] = s_val

        # Weighted degradation contribution
        D_j = (T['weight'] / total_weight) * sens_map * decay_factor
        D_total += D_j
        print(f"    Mean D contribution: {np.nanmean(D_j[valid]):.4f}")

    D_total[~valid] = np.nan
    # Clip to [0, 1]
    D_total = np.clip(D_total, 0, 1)
    print(f"\n  D_x summary: mean={np.nanmean(D_total):.4f}  "
          f"max={np.nanmax(D_total):.4f}")
    return D_total


# ══════════════════════════════════════════════════════════════════════
# PART 4 — Compute Habitat Quality (Q_x)
# ══════════════════════════════════════════════════════════════════════
def compute_quality(H, D):
    print(f"\n{SEP}")
    print("PART 4 — Computing Habitat Quality Index Q_x")
    print(SEP)

    # InVEST formula: Q_x = H_x × [1 − D_x^z / (D_x^z + k^z)]
    D_z = np.power(D, Z_SCALE)
    k_z = HALF_SAT_K ** Z_SCALE
    Q   = H * (1.0 - D_z / (D_z + k_z))
    Q   = np.clip(Q, 0, 1)

    valid = ~np.isnan(Q)
    print(f"  Q_x summary: mean={np.nanmean(Q):.4f}  "
          f"min={np.nanmin(Q):.4f}  max={np.nanmax(Q):.4f}")
    print(f"  High quality (Q≥0.7) : {np.sum(Q[valid] >= 0.7):,} pixels "
          f"({100*np.sum(Q[valid] >= 0.7)/valid.sum():.1f}%)")
    print(f"  Degraded    (Q≤0.3)  : {np.sum(Q[valid] <= 0.3):,} pixels "
          f"({100*np.sum(Q[valid] <= 0.3)/valid.sum():.1f}%)")

    # ── Summary statistics by ESA class ───────────────────────────
    class_labels = {
        10: 'Tree cover', 20: 'Shrubland', 30: 'Grassland',
        40: 'Cropland',  50: 'Built-up',  80: 'Open water',
        90: 'Wetland',   95: 'Mangroves', 100: 'Moss',
    }
    print(f"\n  Mean Habitat Quality by Land Cover:")
    print(f"  {'Class':<22} {'Mean Q':>7}  {'Pixels':>10}")
    print(f"  {'-'*45}")
    return Q


# ══════════════════════════════════════════════════════════════════════
# PART 5 — Save Outputs & Figures
# ══════════════════════════════════════════════════════════════════════

# ── Shared palette constants ───────────────────────────────────────────
OCEAN_COLOR  = '#b3d9f2'   # light blue ocean
LAND_COLOR   = '#DDECC5'   # sage green for intact land
BG_COLOR     = 'white'     # all figure backgrounds
TEXT_COLOR   = '#1a1a2e'   # dark navy text
SPINE_COLOR  = '#cccccc'


def _draw_land_ocean_base(ax, land_binary, ocean_rgba, land_rgba):
    """Render ocean + land base layers before overlaying data."""
    h, w = land_binary.shape
    base = np.ones((h, w, 4), dtype=float)
    # Ocean
    oc = mcolors.to_rgba(ocean_rgba)
    lc = mcolors.to_rgba(land_rgba)
    ocean_px = (land_binary == 0)
    base[ocean_px] = oc
    land_px  = (land_binary == 1)
    base[land_px] = lc
    ax.imshow(base, aspect='equal', interpolation='nearest')


def save_outputs(Q, D, H, esa, profile, Q_2000=None):
    print(f"\n{SEP}")
    print("PART 5 — Saving Outputs & Figures")
    print(SEP)

    from matplotlib.colors import LinearSegmentedColormap as LSC
    import matplotlib.gridspec as gridspec

    esa_int   = np.nan_to_num(esa, nan=0).astype(int)
    land_mask = ~np.isnan(esa) & (esa_int != 0)   # True = land pixel

    # ── Save GeoTIFF ───────────────────────────────────────────────
    out_profile = profile.copy()
    out_profile.update({'dtype': 'float32', 'nodata': -9999,
                        'count': 1, 'compress': 'lzw'})
    q_export = np.where(np.isnan(Q), -9999, Q).astype('float32')
    tif_path = RES_DIR / 'habitat_quality_map.tif'
    with rasterio.open(tif_path, 'w', **out_profile) as dst:
        dst.write(q_export, 1)
    print(f"  ✅  GeoTIFF saved → {tif_path}")

    delta_Q = None
    if Q_2000 is not None:
        delta_Q = Q - Q_2000
        delta_export = np.where(np.isnan(delta_Q), -9999, delta_Q).astype('float32')
        delta_tif_path = RES_DIR / 'habitat_quality_delta_map.tif'
        with rasterio.open(delta_tif_path, 'w', **out_profile) as dst:
            dst.write(delta_export, 1)
        print(f"  ✅  GeoTIFF saved → {delta_tif_path}")

    # ── CSV statistics by ESA class ────────────────────────────────
    class_labels = {
        10: 'Tree cover', 20: 'Shrubland', 30: 'Grassland',
        40: 'Cropland',  50: 'Built-up',  60: 'Bare',
        80: 'Open water', 90: 'Wetland',  95: 'Mangroves', 100: 'Moss',
    }
    records = []
    for code, lbl in class_labels.items():
        mask = (esa_int == code) & ~np.isnan(Q)
        if mask.sum() == 0:
            continue
        q_vals = Q[mask]
        record = {
            'esa_class': code, 'label': lbl,
            'pixel_count': int(mask.sum()),
            'area_ha': round(mask.sum() * 0.09, 1),
            'mean_Q': round(float(np.mean(q_vals)), 4),
            'std_Q':  round(float(np.std(q_vals)), 4),
            'min_Q':  round(float(np.min(q_vals)), 4),
            'max_Q':  round(float(np.max(q_vals)), 4),
        }
        if delta_Q is not None:
            record['mean_delta_Q'] = round(float(np.mean(delta_Q[mask])), 4)
        records.append(record)
        print(f"  {lbl:<22} mean Q = {np.mean(q_vals):.4f}  n={mask.sum():,}")

    csv_path = RES_DIR / 'habitat_quality_stats.csv'
    pd.DataFrame(records).to_csv(csv_path, index=False)
    print(f"\n  ✅  CSV saved → {csv_path}")

    (FIG_DIR / 'habitat').mkdir(parents=True, exist_ok=True)

    # ── Shared colormaps ───────────────────────────────────────────
    cmap_q = LSC.from_list(
        'habitat',
        ['#7f0000', '#d73027', '#fdae61', '#ffffbf', '#a6d96a', '#1a9641', '#003d00'], N=256
    )
    cmap_q.set_bad(color=OCEAN_COLOR)

    # Threat: yellow→orange→red (no black = land visible)
    cmap_threat = LSC.from_list(
        'threat',
        ['#ffffcc', '#fecc5c', '#fd8d3c', '#f03b20', '#bd0026', '#67000d'], N=256
    )
    cmap_threat.set_bad(color=OCEAN_COLOR)

    # ── Helper: render split map (Andaman left, Nicobar right) ───────
    def _split_map(fig, gs_l, gs_r, arr, cmap, vmin, vmax,
                   title_l, title_r, stat_text_l=None, stat_text_r=None):
        """White background split map with ocean coloring."""
        ocean_px = ~land_mask
        arr_m    = np.ma.masked_where(ocean_px, arr)
        and_arr, nic_arr = _split_extent(arr_m)
        and_crop = _tight_crop(and_arr)
        nic_crop = _tight_crop(nic_arr)

        ax_a = fig.add_subplot(gs_l)
        ax_n = fig.add_subplot(gs_r)

        for ax, crop, ttl, stat in [
            (ax_a, and_crop, title_l, stat_text_l),
            (ax_n, nic_crop, title_r, stat_text_r)
        ]:
            ax.set_facecolor(OCEAN_COLOR)
            ax.imshow(crop, cmap=cmap, vmin=vmin, vmax=vmax,
                      interpolation='nearest', aspect='equal')
            ax.set_title(ttl, color=TEXT_COLOR, fontsize=11,
                         fontweight='bold', pad=6)
            ax.axis('off')
            if stat:
                ax.text(0.03, 0.04, stat, transform=ax.transAxes,
                        color=TEXT_COLOR, fontsize=8.5,
                        bbox=dict(boxstyle='round,pad=0.35',
                                  facecolor='white', edgecolor=SPINE_COLOR,
                                  alpha=0.92))
        return ax_a, ax_n

    def _add_colorbar(fig, gs_slot, cmap, vmin, vmax, label,
                      ticks=None, ticklabels=None):
        cax  = fig.add_subplot(gs_slot)
        sm   = plt.cm.ScalarMappable(cmap=cmap,
                                      norm=plt.Normalize(vmin=vmin, vmax=vmax))
        sm.set_array([])
        cbar = fig.colorbar(sm, cax=cax)
        cbar.set_label(label, color=TEXT_COLOR, fontsize=9.5)
        cbar.ax.tick_params(colors=TEXT_COLOR, labelsize=8)
        if ticks is not None:
            cbar.set_ticks(ticks)
        if ticklabels is not None:
            cbar.set_ticklabels(ticklabels)
        plt.setp(cbar.ax.yaxis.get_ticklabels(), color=TEXT_COLOR)
        cbar.outline.set_edgecolor(SPINE_COLOR)
        return cbar

    # ══════════════════════════════════════════════════════════════
    # Figure — habitat_quality_index_map.png
    # Split view: rows = [Q, D]; cols = [Andaman, Nicobar, colorbar]
    # ══════════════════════════════════════════════════════════════
    fig2 = plt.figure(figsize=(16, 20), facecolor=BG_COLOR)
    gs2  = gridspec.GridSpec(2, 3, figure=fig2,
                              width_ratios=[2, 2, 0.09],
                              hspace=0.14, wspace=0.06)
    fig2.patch.set_facecolor(BG_COLOR)

    d_vmax = float(np.nanmax(D[land_mask]))

    _split_map(fig2, gs2[0, 0], gs2[0, 1], Q,
               cmap_q, 0, 1,
               'Andaman Islands — Habitat Quality (Q_x)',
               'Nicobar Islands — Habitat Quality (Q_x)',
               stat_text_l=f'Mean Q = {np.nanmean(Q[land_mask]):.3f}')
    _add_colorbar(fig2, gs2[0, 2], cmap_q, 0, 1,
                  'Habitat Quality Index (0–1)',
                  ticks=[0, 0.2, 0.4, 0.6, 0.8, 1.0],
                  ticklabels=['0.0\n(Degraded)', '0.2', '0.4',
                              '0.6', '0.8', '1.0\n(Pristine)'])

    _split_map(fig2, gs2[1, 0], gs2[1, 1], D,
               cmap_threat, 0, d_vmax,
               'Andaman Islands — Cumulative Threat (D_x)',
               'Nicobar Islands — Cumulative Threat (D_x)',
               stat_text_l=f'Max D = {d_vmax:.3f}')
    _add_colorbar(fig2, gs2[1, 2], cmap_threat, 0, d_vmax,
                  'Cumulative Threat Level (D_x)\nRoads + Built-up + Cropland')

    fig2.suptitle(
        'InVEST Habitat Quality Model — Andaman & Nicobar Islands (2024)',
        color=TEXT_COLOR, fontsize=14, fontweight='bold', y=0.995)
    fig2.savefig(FIG_DIR / 'habitat' / 'habitat_quality_index_map.png',
                 dpi=200, bbox_inches='tight', facecolor=BG_COLOR)
    plt.close()
    print("  ✅  Figure saved → figures/habitat/habitat_quality_index_map.png")

    # ══════════════════════════════════════════════════════════════
    # Figure 3 — habitat_quality_by_landcover.png
    # Lollipop chart — white background, proper labels both files
    # ══════════════════════════════════════════════════════════════
    if records:
        df_rec = pd.DataFrame(records).sort_values('mean_Q', ascending=True)
        fig3, ax3 = plt.subplots(figsize=(12, 7), facecolor=BG_COLOR)
        ax3.set_facecolor('#f5f8fc')
        bar_colors = [cmap_q(v) for v in df_rec['mean_Q']]
        y_pos = list(range(len(df_rec)))

        # Threshold zone shading — turns the dashed reference lines into
        # readable regions. Drawn at zorder=0 so they sit behind every
        # lollipop and dot.
        ax3.axvspan(0.00, 0.30, color='#e53935', alpha=0.07, zorder=0)
        ax3.axvspan(0.70, 1.00, color='#2e7d32', alpha=0.07, zorder=0)
        ax3.text(0.15, len(df_rec) - 0.4, 'Degraded zone',
                 ha='center', va='center', fontsize=8.5,
                 color='#7a2a2a', style='italic', alpha=0.85, zorder=2)
        ax3.text(0.85, len(df_rec) - 0.4, 'High-quality zone',
                 ha='center', va='center', fontsize=8.5,
                 color='#1b5e20', style='italic', alpha=0.85, zorder=2)

        for i, (_, row) in enumerate(df_rec.iterrows()):
            ax3.plot([0, row['mean_Q']], [i, i],
                     color='#aaaacc', linewidth=1.8, zorder=1)
        ax3.scatter(df_rec['mean_Q'], y_pos, color=bar_colors, s=240,
                    zorder=3, edgecolors=TEXT_COLOR, linewidths=0.9)

        # Per-row label: Q value + landscape area (ha). Truly-zero Q
        # values are tagged as such instead of being printed "0.000",
        # which previously hid whether the result was rounded.
        has_area = 'area_ha' in df_rec.columns
        for i, (_, row) in enumerate(df_rec.iterrows()):
            val = row['mean_Q']
            if val < 0.0005:
                q_text = '< 0.001'
            else:
                q_text = f'{val:.3f}'
            area_text = (f'   ({row["area_ha"]:,.0f} ha)'
                         if has_area else '')
            ax3.text(val + 0.018, i, f'{q_text}{area_text}',
                     va='center', color=TEXT_COLOR,
                     fontsize=9.5, fontweight='bold')

        ax3.set_yticks(y_pos)
        ax3.set_yticklabels(df_rec['label'], color=TEXT_COLOR, fontsize=10.5)
        ax3.set_xlabel('Mean Habitat Quality Index (Q_x)',
                       color=TEXT_COLOR, fontsize=11)
        ax3.set_xlim(-0.02, 1.35)
        ax3.tick_params(axis='x', colors=TEXT_COLOR)
        ax3.tick_params(axis='y', colors=TEXT_COLOR)
        for spine in ['top', 'right']:
            ax3.spines[spine].set_visible(False)
        for spine in ['bottom', 'left']:
            ax3.spines[spine].set_color(SPINE_COLOR)

        ax3.axvline(0.7, color='#2e7d32', linestyle='--', alpha=0.65,
                    linewidth=1.6, label='High Quality (Q ≥ 0.7)')
        ax3.axvline(0.3, color='#e53935', linestyle='--', alpha=0.65,
                    linewidth=1.6, label='Degraded (Q ≤ 0.3)')
        ax3.legend(facecolor='white', edgecolor=SPINE_COLOR,
                   labelcolor=TEXT_COLOR, fontsize=9, loc='lower right')
        ax3.set_title('Mean Habitat Quality by Land Cover Class\n'
                      'Andaman & Nicobar Islands (2024)',
                      color=TEXT_COLOR, fontsize=13, fontweight='bold', pad=14)

        # Open-water footnote (clarifies the Q score for that class).
        if 'Open water' in df_rec['label'].values or \
           'open water' in [str(v).lower() for v in df_rec['label']]:
            fig3.text(
                0.5, -0.02,
                '★ Open-water Q reflects inland water bodies (lakes, '
                'lagoons) — not the surrounding ocean, which is excluded '
                'from the InVEST land domain.',
                ha='center', va='top', fontsize=8.5, color='#4a5568',
                bbox=dict(boxstyle='round,pad=0.35', facecolor='#fffde7',
                          edgecolor='#f0c000', alpha=0.9),
            )

        fig3.tight_layout()
        fig3.savefig(FIG_DIR / 'habitat' / 'habitat_quality_by_landcover.png',
                     dpi=180, bbox_inches='tight', facecolor=BG_COLOR)
        plt.close()
        print("  ✅  Figure saved → figures/habitat/habitat_quality_by_landcover.png")

    # ══════════════════════════════════════════════════════════════
    # Figure 4 — habitat_quality_delta_map.png
    # Temporal change 2000→2024
    # ══════════════════════════════════════════════════════════════
    if delta_Q is not None:
        delta_land_vals = delta_Q[land_mask]
        changed_mask    = np.abs(delta_land_vals) > 0.001

        # Crop the colormap to the 2nd–98th percentile of *changed*
        # pixels so the gradient actually pops. Floor at 0.05 so the
        # range is always meaningful even if changes are tiny.
        if changed_mask.sum() > 0:
            abs_max = max(
                abs(np.percentile(delta_land_vals[changed_mask], 2)),
                abs(np.percentile(delta_land_vals[changed_mask], 98)),
                0.05,
            )
        else:
            abs_max = 0.1
        vmin_d, vmax_d = -abs_max, abs_max

        cmap_delta = LSC.from_list(
            'delta_hab',
            ['#b2182b', '#ef8a62', '#fddbc7',
             '#f7f7f7', '#d1e5f0', '#4393c3', '#053061'], N=256,
        )
        cmap_delta.set_bad(alpha=0)

        # Per-class hectare tallies (real, not dilated).
        lost_total   = float(np.sum(delta_land_vals < -0.01) * 0.09)
        gained_total = float(np.sum(delta_land_vals >  0.01) * 0.09)

        # Dilated copy of delta_Q for *rendering* only — single-pixel
        # changes are otherwise invisible at full-island zoom. We
        # propagate the change *value* outward, not just a binary mask,
        # via grey_dilation so dilated pixels keep their ΔQ colour.
        changed_2d = (np.abs(delta_Q) > 0.001) & land_mask
        changed_2d = binary_dilation(changed_2d, iterations=4)
        delta_for_vis = np.where(np.isnan(delta_Q), 0.0, delta_Q)
        # Use the max-abs neighbourhood value so each dilated pixel
        # takes the colour of the strongest nearby real change.
        from scipy.ndimage import grey_dilation as _grey_dil
        pos = np.where(delta_for_vis > 0,  delta_for_vis, 0)
        neg = np.where(delta_for_vis < 0, -delta_for_vis, 0)
        pos = _grey_dil(pos, size=(9, 9))
        neg = _grey_dil(neg, size=(9, 9))
        combined = np.where(neg > pos, -neg, pos)
        delta_vis  = np.where(changed_2d, combined, np.nan)
        delta_vis_m = np.ma.masked_where(np.isnan(delta_vis), delta_vis)

        fig4 = plt.figure(figsize=(20, 12), facecolor=BG_COLOR)
        gs4  = gridspec.GridSpec(
            1, 4, figure=fig4,
            width_ratios=[2.0, 1.4, 0.06, 1.6], wspace=0.10,
            left=0.04, right=0.97, top=0.92, bottom=0.07,
        )

        # Per-region splits. land_nan = land mask with ocean=NaN so
        # _tight_crop can find the land bbox.
        land_nan = np.where(land_mask, 1.0, np.nan)
        and_delta, nic_delta = _split_extent(delta_vis_m)
        and_lm,    nic_lm    = _split_extent(land_nan)
        and_lm_b,  nic_lm_b  = _split_extent(land_mask)

        # Per-region hectare tallies for panel titles (computed on
        # raw delta_Q, restricted to land).
        and_drow_end = and_lm_b.shape[0]
        nic_drow_beg = delta_Q.shape[0] - nic_lm_b.shape[0]
        and_dvals = delta_Q[:and_drow_end][and_lm_b]
        nic_dvals = delta_Q[nic_drow_beg:][nic_lm_b]
        and_lost   = float(np.sum(and_dvals < -0.01) * 0.09)
        and_gained = float(np.sum(and_dvals >  0.01) * 0.09)
        nic_lost   = float(np.sum(nic_dvals < -0.01) * 0.09)
        nic_gained = float(np.sum(nic_dvals >  0.01) * 0.09)

        panels = [
            (gs4[0, 0], and_delta, and_lm,
             'Andaman Islands',
             f'Degraded: {and_lost:,.0f} ha   •   Improved: {and_gained:,.0f} ha'),
            (gs4[0, 1], nic_delta, nic_lm,
             'Nicobar Islands',
             f'Degraded: {nic_lost:,.0f} ha   •   Improved: {nic_gained:,.0f} ha'),
        ]
        for gs_pos, arr_rgn, lm_rgn, region, stats_line in panels:
            ax = fig4.add_subplot(gs_pos)
            ax.set_facecolor('#e6eef5')          # pale blue-grey ocean

            # Crop both layers to the same land bbox using lm_rgn (NaN
            # over ocean) so _tight_crop returns a tight land window.
            # Bbox derived from lm_rgn validity, then sliced identically.
            valid_rows = np.any(~np.isnan(lm_rgn), axis=1)
            valid_cols = np.any(~np.isnan(lm_rgn), axis=0)
            if valid_rows.any() and valid_cols.any():
                r0, r1 = np.where(valid_rows)[0][[0, -1]]
                c0, c1 = np.where(valid_cols)[0][[0, -1]]
                pad = 30
                r0 = max(0, r0 - pad); r1 = min(lm_rgn.shape[0], r1 + pad)
                c0 = max(0, c0 - pad); c1 = min(lm_rgn.shape[1], c1 + pad)
                lm_crop  = lm_rgn[r0:r1, c0:c1]
                arr_crop = arr_rgn[r0:r1, c0:c1]
            else:
                lm_crop  = lm_rgn
                arr_crop = arr_rgn

            # Pale-green land underlay (so unchanged land reads as land,
            # not ocean) drawn beneath the ΔQ overlay.
            land_only = np.where(np.isnan(lm_crop), np.nan, 1.0)
            ax.imshow(
                land_only,
                cmap=LSC.from_list('flatland', ['#dfeadd', '#dfeadd']),
                vmin=0, vmax=1, interpolation='nearest',
            )

            # ΔQ overlay (dilated for visibility).
            ax.imshow(arr_crop, cmap=cmap_delta,
                      vmin=vmin_d, vmax=vmax_d,
                      interpolation='nearest', aspect='equal')

            # Bold coastline (white halo + dark line) — drawn from a
            # 0/1 land grid (use np.where because contour can't take NaN).
            lm_for_contour = np.where(np.isnan(lm_crop), 0.0, 1.0)
            ax.contour(lm_for_contour, levels=[0.5], colors=['#ffffff'],
                       linewidths=2.4, alpha=0.95)
            ax.contour(lm_for_contour, levels=[0.5], colors=['#0d1b2a'],
                       linewidths=1.4, alpha=1.0)

            ax.set_title(f'{region}\n{stats_line}',
                         color=TEXT_COLOR, fontsize=12,
                         fontweight='bold', pad=10)
            ax.axis('off')

        _add_colorbar(
            fig4, gs4[0, 2], cmap_delta, vmin_d, vmax_d,
            'ΔQ  (2000 → 2024)\n← Degraded   |   Improved →',
        )

        # Histogram panel — its own clean axes (no longer wedged
        # alongside the colorbar).
        ax_h = fig4.add_subplot(gs4[0, 3])
        ax_h.set_facecolor('#f5f8fc')
        hist_vals = delta_land_vals[np.abs(delta_land_vals) > 1e-5]
        if len(hist_vals) > 0:
            n, bins, patches = ax_h.hist(
                hist_vals, bins=60,
                edgecolor='#cccccc', linewidth=0.4,
            )
            # Colour each bar by its bin-centre ΔQ value.
            bin_centres = 0.5 * (bins[:-1] + bins[1:])
            for patch, c in zip(patches, bin_centres):
                frac = (c - vmin_d) / (vmax_d - vmin_d)
                frac = max(0.0, min(1.0, frac))
                patch.set_facecolor(cmap_delta(frac))
        ax_h.axvline(0, color='#555', linewidth=1.2,
                     linestyle='--', alpha=0.8)
        ax_h.set_xlabel('ΔQ Value', color=TEXT_COLOR, fontsize=10)
        ax_h.set_ylabel('Pixel Count (log)', color=TEXT_COLOR, fontsize=10)
        ax_h.set_yscale('log')
        ax_h.tick_params(colors=TEXT_COLOR)
        for s in ['top', 'right']:
            ax_h.spines[s].set_visible(False)
        for s in ['bottom', 'left']:
            ax_h.spines[s].set_color(SPINE_COLOR)
        ax_h.set_title(
            f'ΔQ Distribution\n'
            f'Net: −{lost_total:,.0f} ha degraded, +{gained_total:,.0f} ha improved',
            color=TEXT_COLOR, fontsize=10.5, fontweight='bold',
        )

        fig4.suptitle(
            'Temporal Habitat Quality Change — '
            'Andaman & Nicobar Islands (2000 → 2024)',
            color=TEXT_COLOR, fontsize=14, fontweight='bold', y=0.97,
        )
        fig4.savefig(FIG_DIR / 'habitat' / 'habitat_quality_delta_map.png',
                     dpi=200, bbox_inches='tight', facecolor=BG_COLOR)
        plt.close()
        print("  ✅  Figure saved → figures/habitat/habitat_quality_delta_map.png")



if __name__ == '__main__':
    print(SEP)
    print("  ANI Ecosystem Services — Week 7: Habitat Quality")
    print(f"  Input : {PROC_DIR}")
    print(f"  Output: {RES_DIR}  |  {FIG_DIR}")
    print(SEP)

    esa, profile = load_raster('ANI_ESA_WorldCover_mosaic_clipped.tif')
    if esa is None:
        print("❌  ESA WorldCover layer not found. Aborting.")
        exit(1)

    loss_yr, _ = load_raster('ANI_GFW_Forest_Loss_2001_2023_clipped.tif')
    
    print(f"\n{SEP}")
    print("  PHASE A: RECONSTRUCTING YEAR 2000 BASELINE")
    print(f"{SEP}")
    
    esa_2000 = esa.copy()
    if loss_yr is not None:
        # If forest was lost after 2000, it MUST have been forest in 2000!
        loss_mask = (np.nan_to_num(loss_yr, nan=0) > 0)
        esa_2000[loss_mask] = 10  # Set to class 10 (Tree cover / Primary habitat)
        print(f"  Re-forested {loss_mask.sum():,} deforested pixels to restore Year 2000 baseline.")

    H_2000 = build_habitat_map(esa_2000, profile)
    threat_rasts_2000 = build_threat_rasters(esa_2000, profile)
    
    # Remove threats on reforested pixels (they were forest, not roads/farms back then)
    if loss_yr is not None:
        for t_name in threat_rasts_2000:
            threat_rasts_2000[t_name][loss_mask] = 0.0
            
    D_2000 = compute_degradation(esa_2000, threat_rasts_2000, profile)
    Q_2000 = compute_quality(H_2000, D_2000)

    print(f"\n{SEP}")
    print("  PHASE A2: STRUCTURAL FRAGMENTATION (YEAR 2000)")
    print(f"{SEP}")
    
    def compute_structural_fragmentation(esa_array, year_label):
        esa_int = np.nan_to_num(esa_array, nan=0).astype(int)
        # Tree cover (10) and Mangroves (95)
        forest_mask = (esa_int == 10) | (esa_int == 95)
        
        # Shrink forest by ~90 meters (3 pixels) to find CORE forest
        core_mask = binary_erosion(forest_mask, iterations=3)
        edge_mask = forest_mask & ~core_mask
        
        core_ha = core_mask.sum() * 0.09
        edge_ha = edge_mask.sum() * 0.09
        total_ha = forest_mask.sum() * 0.09
        
        # FRAGSTATS Metrics: Patch Number (NP) and Mean Patch Size (MPS)
        # 8-connectivity labeling (a patch is contiguous if pixels touch diagonally)
        structure = np.ones((3, 3), dtype=int)
        labeled_cores, num_patches = label(core_mask, structure=structure)
        
        mps_ha = (core_ha / num_patches) if num_patches > 0 else 0
        
        print(f"  Structural Fragmentation Analysis ({year_label}):")
        if total_ha > 0:
            print(f"    Total Forest Area : {total_ha:,.1f} ha")
            print(f"    Core Forest Area  : {core_ha:,.1f} ha ({(core_ha/total_ha)*100:.1f}%)")
            print(f"    Edge Forest Area  : {edge_ha:,.1f} ha ({(edge_ha/total_ha)*100:.1f}%)")
            print(f"    FRAGSTATS (NP)    : {num_patches:,} Core Patches")
            print(f"    FRAGSTATS (MPS)   : {mps_ha:,.2f} ha / patch")
            
        return core_ha, edge_ha, num_patches, mps_ha

    core_2000, edge_2000, np_2000, mps_2000 = compute_structural_fragmentation(esa_2000, "Year 2000 Baseline")

    print(f"\n{SEP}")
    print("  PHASE B: CURRENT STATE (2024)")
    print(f"{SEP}")

    H            = build_habitat_map(esa, profile)
    threat_rasts = build_threat_rasters(esa, profile)
    D            = compute_degradation(esa, threat_rasts, profile)
    Q            = compute_quality(H, D)
    
    print(f"\n{SEP}")
    print("  PHASE B2: STRUCTURAL FRAGMENTATION (CURRENT 2024)")
    print(f"{SEP}")
    
    core_2024, edge_2024, np_2024, mps_2024 = compute_structural_fragmentation(esa, "Year 2024 Current")
    
    print(f"\n  FRAGMENTATION DELTA (2000 -> 2024):")
    print(f"    Lost Core Forest : {core_2000 - core_2024:,.1f} ha")
    print(f"    Change in Edge   : {edge_2024 - edge_2000:+,.1f} ha")
    print(f"    New Core Patches : {np_2024 - np_2000:+,} islands created (NP)")
    print(f"    Change in Area/P : {mps_2024 - mps_2000:+,.2f} ha (MPS)")

    Q            = save_outputs(Q, D, H, esa, profile, Q_2000=Q_2000)

    print(f"\n{SEP}")
    print("  ✅  Week 7 Habitat Quality Assessment Complete!")
    print(f"      GeoTIFF : results/habitat_quality_index.tif")
    print(f"      CSV     : results/habitat_quality_by_landcover.csv")
    print(f"      Figures : figures/habitat/")
    print(f"  Next → run: venv/bin/python src/soil_retention.py")
    print(SEP)
