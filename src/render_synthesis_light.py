"""
Render all synthesis figures on light backgrounds — publication-ready.

Fixes applied:
 1. bivariate_habitat_erosion_kde   — seaborn-style JointGrid, light bg, r annotation
 2. eci_collapse_hotspots_map       — land-masked ECI + proper histogram (not flat)
 3. hotspot_hexbin_density_map      — hexbin on white + island boundary overlay
 4. stat_distribution_habitat_quality — light violin + annotated medians
 5. stat_distribution_soil_erosion   — same
 6. tradeoff_radar_chart             — properly normalised (high=good), light bg

Run: ./venv/venv/bin/python src/render_synthesis_light.py
"""

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import rasterio
import geopandas as gpd
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats
from scipy.stats import gaussian_kde
from pathlib import Path

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'axes.spines.top':   False,
    'axes.spines.right': False,
})

PROC = Path('data/processed')
RAW  = Path('data/raw')
RES  = Path('results')
FIG  = Path('figures/synthesis')
FIG.mkdir(parents=True, exist_ok=True)

BG = 'white'
TC = '#1a1a2e'

# ── Load rasters ──────────────────────────────────────────────────────
def load(path):
    with rasterio.open(path) as src:
        arr = src.read(1).astype(float)
        nd  = src.nodata
        pro = src.profile
    if nd is not None:
        arr[arr == nd] = np.nan
    arr[arr < -1e9] = np.nan
    return arr, pro

print("Loading rasters …")
lc,  lc_p = load(PROC / 'ANI_ESA_WorldCover_mosaic_clipped.tif')
agb, _    = load(PROC / 'ANI_GEDI_Biomass_Density_clipped.tif')
agb       = np.where((agb > 0) & (agb <= 800), agb, np.nan)
hq, hq_p  = load(RES  / 'habitat_quality_index.tif')
eci, eci_p = load(RES / 'eci_collapse_hotspots.tif')
dy, _     = load(PROC / 'ANI_GFW_Forest_Loss_2001_2023_clipped.tif')

land_mask = ~np.isnan(lc)

print(f"  LC  : {lc.shape}   HQ: {hq.shape}   ECI: {eci.shape}")
print(f"  HQ range   : {np.nanmin(hq):.3f} – {np.nanmax(hq):.3f}")
print(f"  ECI range  : {np.nanmin(eci):.3f} – {np.nanmax(eci):.3f}")

# ── Align HQ/ECI to LC grid ───────────────────────────────────────────
def align_to(arr, ref):
    """Crop/pad arr to match ref shape."""
    if arr.shape == ref.shape:
        return arr
    out = np.full(ref.shape, np.nan, dtype=float)
    r  = min(arr.shape[0], ref.shape[0])
    c  = min(arr.shape[1], ref.shape[1])
    out[:r, :c] = arr[:r, :c]
    return out

hq  = align_to(hq,  lc)
eci = align_to(eci, lc)

# Clip valid values
hq  = np.where(land_mask & (hq  >= 0) & (hq  <= 1),   hq,  np.nan)
eci = np.where(land_mask & (eci >= 0) & (eci <= 1),    eci, np.nan)

# ── Real per-pixel RUSLE soil-loss raster ─────────────────────────────
#    results/rusle_soil_loss.tif: 25500×7503, EPSG:32646, t/ha/yr,
#    NoData = -9999. Aligns 1:1 with the LC grid, so no resampling needed.
se, _ = load(RES / 'rusle_soil_loss.tif')
se = align_to(se, lc)
se = np.where(land_mask & (se >= 0), se, np.nan)
np.random.seed(42)   # kept for reproducibility of downstream samplers

# ── Boundary overlay ──────────────────────────────────────────────────
# Store boundary in BOTH pixel space (for imshow plots) and world CRS
# coordinates (for plots whose axes are in geographic units like the hexbin).
bnd = gpd.read_file(RAW / 'ANI_Administrative_Boundary.shp').to_crs(lc_p['crs'])
bnd_pixel = []   # list of (cols, rows)  — pixel indices into the LC raster
bnd_world = []   # list of (xs,   ys)    — same units as lc_p['crs']
for geom in bnd.geometry:
    if geom is None:
        continue
    polys = list(geom.geoms) if geom.geom_type == 'MultiPolygon' else [geom]
    for p in polys:
        xs, ys = p.exterior.xy
        xs, ys = list(xs), list(ys)
        rows, cols = rasterio.transform.rowcol(lc_p['transform'], xs, ys)
        bnd_pixel.append((list(cols), list(rows)))
        bnd_world.append((xs, ys))

def draw_bnd(ax, DS=1, col='#444444', lw=1.0):
    """Pixel-space boundary — for imshow plots."""
    for c, r in bnd_pixel:
        ax.plot([x/DS for x in c], [x/DS for x in r],
                color=col, lw=lw, alpha=0.75, zorder=10)

def draw_bnd_world(ax, col='#444444', lw=1.0):
    """CRS-world-space boundary — for plots whose axes are in CRS units."""
    for xs, ys in bnd_world:
        ax.plot(xs, ys, color=col, lw=lw, alpha=0.75, zorder=10)

# ── Sample pixels per land-cover class ───────────────────────────────
def sample_class(arr, cls, max_n=50_000):
    mask = (lc == cls) & ~np.isnan(arr)
    vals = arr[mask].ravel()
    if len(vals) > max_n:
        vals = np.random.choice(vals, max_n, replace=False)
    return vals

LC_DEFS = [
    (10, 'Tree Cover',  '#2E7D32'),
    (95, 'Mangroves',   '#00838F'),
    (40, 'Cropland',    '#FF8F00'),
    (50, 'Built-up',    '#D32F2F'),
    (60, 'Bare/Sparse', '#8D6E63'),
    (30, 'Grassland',   '#8BC34A'),
    (20, 'Shrubland',   '#A5D6A7'),
    (90, 'Wetland',     '#26A69A'),
]

# ══════════════════════════════════════════════════════════════════════
# FIG 1 — bivariate_habitat_erosion_kde
# ══════════════════════════════════════════════════════════════════════
print("\nRendering: bivariate_habitat_erosion_kde.png …")

valid_both = ~np.isnan(hq) & ~np.isnan(se) & land_mask
tree_px    = valid_both & (lc == 10)
crop_px    = valid_both & (lc == 40)

def rand_sample_xy(mask, q_arr, s_arr, n=15_000):
    idx = np.where(mask.ravel())[0]
    if len(idx) > n:
        idx = np.random.choice(idx, n, replace=False)
    r, c = np.unravel_index(idx, mask.shape)
    return q_arr[r, c], s_arr[r, c]

np.random.seed(42)
hq_t, se_t = rand_sample_xy(tree_px, hq, se, 12_000)
hq_c, se_c = rand_sample_xy(crop_px, hq, se,  4_000)

hq_all = np.concatenate([hq_t, hq_c])
se_all = np.concatenate([se_t, se_c])
r_val, p_val = stats.pearsonr(hq_all, np.log1p(se_all))

fig = plt.figure(figsize=(10, 9), facecolor=BG)
gs  = fig.add_gridspec(4, 4, hspace=0.05, wspace=0.05)
ax_main  = fig.add_subplot(gs[1:, :3])
ax_top   = fig.add_subplot(gs[0,  :3], sharex=ax_main)
ax_right = fig.add_subplot(gs[1:,  3], sharey=ax_main)
for a in [ax_main, ax_top, ax_right]:
    a.set_facecolor('white')   # ← pure white background (not #fafafa grey)

ax_main.scatter(hq_t, se_t, s=4, alpha=0.20, color='#2E7D32', label='Tree Cover', rasterized=True)
ax_main.scatter(hq_c, se_c, s=8, alpha=0.45, color='#FF8F00', label='Cropland',   rasterized=True)

x_line = np.linspace(hq_all.min(), hq_all.max(), 100)
m, b   = np.polyfit(hq_all, np.log1p(se_all), 1)
ax_main.plot(x_line, np.expm1(m*x_line + b), '--', color='#B71C1C', lw=1.8, label='Log-linear fit', zorder=5)
ax_main.set_xlabel('Habitat Quality Index (Q)', fontsize=12, color=TC)
ax_main.set_ylabel('Annual Soil Erosion (t/ha/yr)', fontsize=12, color=TC)
ylim_top_main = np.nanpercentile(se_all, 99) * 1.1
ax_main.set_ylim(0, ylim_top_main)
ax_main.set_xlim(0.05, 1.02)
ax_main.legend(fontsize=10, framealpha=0.9)
p_str = 'p < 0.001' if p_val < 0.001 else f'p = {p_val:.3f}'
ax_main.text(0.97, 0.97,
             f'r = {r_val:.3f}\n{p_str}\nn = {len(hq_all):,}',
             transform=ax_main.transAxes, ha='right', va='top', fontsize=11,
             bbox=dict(facecolor='white', edgecolor='#ccc', boxstyle='round'))

for vals, col in [(hq_t, '#2E7D32'), (hq_c, '#FF8F00')]:
    kde = gaussian_kde(vals, bw_method=0.12)
    xs  = np.linspace(0.05, 1.0, 200)
    ax_top.fill_between(xs, kde(xs), alpha=0.35, color=col)
    ax_top.plot(xs, kde(xs), color=col, lw=1.5)
ax_top.axis('off')

for vals, col in [(se_t, '#2E7D32'), (se_c, '#FF8F00')]:
    # Clip to visible y range to prevent outlier spikes
    ylim_top = ylim_top_main
    vals_clipped = vals[vals <= ylim_top]
    if len(vals_clipped) < 10:
        continue
    ys = np.linspace(0, ylim_top, 300)
    kde = gaussian_kde(vals_clipped, bw_method=0.35)
    density = kde(ys)
    density = density / density.max() * 0.8   # normalise width
    ax_right.fill_betweenx(ys, density, alpha=0.35, color=col)
    ax_right.plot(density, ys, color=col, lw=1.5)
ax_right.set_ylim(0, ylim_top_main)
ax_right.axis('off')

plt.setp(ax_top.get_xticklabels(),   visible=False)
plt.setp(ax_right.get_yticklabels(), visible=False)
fig.suptitle('Bivariate Correlation: Habitat Quality vs. Soil Erosion\nAndaman & Nicobar Islands (2024)',
             fontsize=13, fontweight='bold', color=TC, y=0.97)
fig.savefig(FIG / 'bivariate_habitat_erosion_kde.png', dpi=180, bbox_inches='tight', facecolor=BG)
plt.close()
print("✅  bivariate_habitat_erosion_kde.png")

# ══════════════════════════════════════════════════════════════════════
# FIG 2 — eci_collapse_hotspots_map
# ══════════════════════════════════════════════════════════════════════
print("\nRendering: eci_collapse_hotspots_map.png …")

eci_land = eci[~np.isnan(eci)]
thresh   = np.nanpercentile(eci_land, 90)
hotspot  = (eci >= thresh) & land_mask
print(f"  ECI: {eci_land.min():.3f}–{eci_land.max():.3f}  threshold:{thresh:.3f}  "
      f"hotspots: {hotspot.sum()*0.09:,.0f} ha")

# The ECI raster is effectively binary here (>99% of valid pixels sit
# at ECI=1.0), so a continuous colormap and a density histogram both
# read as empty. Render categorically: pale-sage land base + dilated
# red hotspots, and switch the histogram to log pixel-counts so the
# rare intermediate-ECI pixels actually become visible. The Andaman and
# Nicobar groups are plotted in separate panels so the dilation does
# not visually merge them across the 10-degree channel.
from scipy.ndimage import binary_dilation
DS              = 8
hot_full        = binary_dilation(hotspot, iterations=3)

# Row-split: the ANI raster has Andamans in the upper rows, Nicobars in
# the lower rows. We use the same split convention as the hexbin map.
H = land_mask.shape[0]
ANDAMAN_ROW_END   = min(int(H * 0.59), H)   # ~rows 0  → 15000 at full res
NICOBAR_ROW_START = min(int(H * 0.66), H)   # ~rows 16800 → H
and_slice = slice(0, ANDAMAN_ROW_END)
nic_slice = slice(NICOBAR_ROW_START, H)

land_a = land_mask[and_slice][::DS, ::DS]
land_n = land_mask[nic_slice][::DS, ::DS]
hot_a  = hot_full[and_slice][::DS, ::DS]
hot_n  = hot_full[nic_slice][::DS, ::DS]

and_ha = hotspot[and_slice].sum() * 0.09
nic_ha = hotspot[nic_slice].sum() * 0.09

fig = plt.figure(figsize=(15, 10), facecolor=BG)
gs_eci = fig.add_gridspec(1, 3, width_ratios=[1.1, 1.0, 1.1], wspace=0.18)
ax_a    = fig.add_subplot(gs_eci[0, 0])
ax_n    = fig.add_subplot(gs_eci[0, 1])
ax_hist = fig.add_subplot(gs_eci[0, 2])

def _render_island_panel(ax, land_ds, hot_ds, title, hot_ha_val,
                          bnd_row_offset=0):
    ax.set_facecolor('#b3d9f2')
    # Land base
    land_rgba = np.zeros((*land_ds.shape, 4))
    land_rgba[land_ds] = [0.91, 0.94, 0.85, 1.0]
    ax.imshow(land_rgba, aspect='equal', interpolation='nearest')
    # Hotspot overlay
    hot_rgba = np.zeros((*hot_ds.shape, 4))
    hot_rgba[hot_ds] = [0.72, 0.10, 0.10, 0.95]
    ax.imshow(hot_rgba, aspect='equal', interpolation='nearest')
    # Boundary (only the polygons whose pixel rows fall inside this slice)
    for c, r in bnd_pixel:
        r_arr = np.asarray(r); c_arr = np.asarray(c)
        in_slice = (r_arr >= bnd_row_offset) & (r_arr < bnd_row_offset + land_ds.shape[0] * DS)
        if in_slice.any():
            ax.plot((c_arr[in_slice]) / DS,
                    (r_arr[in_slice] - bnd_row_offset) / DS,
                    color='#333333', lw=0.9, alpha=0.75, zorder=10)
    ax.set_title(title, fontsize=12, fontweight='bold', color=TC, pad=8)
    ax.axis('off')
    ax.text(0.02, 0.02,
            f'Triple-collapse hotspots:\n{hot_ha_val:,.0f} ha at ECI $\\geq$ {thresh:.2f}',
            transform=ax.transAxes, fontsize=9, color='#B71C1C',
            bbox=dict(facecolor='white', alpha=0.9, boxstyle='round',
                      edgecolor='#ccc'))

_render_island_panel(ax_a, land_a, hot_a, 'Andaman Islands', and_ha,
                     bnd_row_offset=0)
_render_island_panel(ax_n, land_n, hot_n, 'Nicobar Islands', nic_ha,
                     bnd_row_offset=NICOBAR_ROW_START)

# One shared legend below the two map panels
legend_handles = [
    mpatches.Patch(facecolor='#e8f0d8', edgecolor='#888',
                    label='Land (not in hotspot bin)'),
    mpatches.Patch(facecolor='#B71C1C',
                    label=f'Triple-collapse hotspot (ECI $\\geq$ {thresh:.2f})'),
]
fig.legend(handles=legend_handles, loc='lower center', ncol=2,
           facecolor='white', edgecolor='#ccc', fontsize=10,
           framealpha=0.95, bbox_to_anchor=(0.36, 0.02))

# Histogram of ECI — log pixel count so the rare intermediate values
# (only 15 pixels between 0.667 and 0.99 in this raster) are visible
# alongside the >26000 pixels at 1.0. Y axis cropped to the data range.
y_min = max(0.6, float(eci_land.min()) - 0.02)
y_max = 1.02
ax_hist.hist(eci_land, bins=40, color='#EF9A9A', edgecolor='white',
             linewidth=0.4, orientation='horizontal', log=True)
ax_hist.axhline(thresh, color='#B71C1C', lw=2, ls='--',
                label=f'Top-10\\% threshold = {thresh:.2f}')
ax_hist.fill_betweenx([thresh, y_max],
                      0.5, 1e6,
                      color='#B71C1C', alpha=0.10, label='Hotspot zone (top 10\\%)')
ax_hist.set_xlabel('Pixel count (log scale)', fontsize=10)
ax_hist.set_ylabel('ECI Score', fontsize=10)
ax_hist.set_ylim(y_min, y_max)
ax_hist.set_title('ECI Score Distribution', fontsize=11, fontweight='bold', color=TC)
ax_hist.legend(fontsize=9, framealpha=0.9, loc='upper right')
ax_hist.text(0.02, 0.02,
             f'{(eci_land >= 0.999).sum():,} pixels at ECI = 1.00\n'
             f'{((eci_land > 0.667) & (eci_land < 0.999)).sum():,} pixels at intermediate values\n'
             f'$n_{{\\mathrm{{total}}}} = {eci_land.size:,}$',
             transform=ax_hist.transAxes, fontsize=8.5, color='#444',
             bbox=dict(facecolor='white', alpha=0.9, boxstyle='round', edgecolor='#ccc'))

fig.suptitle('Triple-Collapse Ecosystem Risk — ANI (2000–2024)',
             fontsize=14, fontweight='bold', color=TC, y=1.01)
fig.tight_layout()
fig.savefig(FIG / 'eci_collapse_hotspots_map.png', dpi=180, bbox_inches='tight', facecolor=BG)
plt.close()
print("✅  eci_collapse_hotspots_map.png")

# ══════════════════════════════════════════════════════════════════════
# FIG 3 — hotspot_hexbin_density_map
# ══════════════════════════════════════════════════════════════════════
print("\nRendering: hotspot_hexbin_density_map.png …")

rows_h, cols_h = np.where(hotspot)
if len(rows_h) > 40_000:
    idx = np.random.choice(len(rows_h), 40_000, replace=False)
    rows_h, cols_h = rows_h[idx], cols_h[idx]

# Convert pixel row/col to real geographic coordinates
tf      = lc_p['transform']
# xs_h and ys_h are already in the raster's CRS units
# tf.c = west edge (lon or easting), tf.f = north edge (lat or northing)
# tf.a = pixel width (positive), tf.e = pixel height (negative)
xs_h    = tf.c + cols_h * tf.a   # longitude / easting
ys_h    = tf.f + rows_h * tf.e   # latitude  / northing (negative step flips correctly)

# Auto-detect degree coords vs projected (metres)
is_geo  = abs(tf.a) < 1.0   # degrees have |step| << 1
x_label = 'Longitude (°E)' if is_geo else 'Easting (m)'
y_label = 'Latitude (°N)'  if is_geo else 'Northing (m)'

from matplotlib.colors import LogNorm

# Geographic split: the 10°N channel between Little Andaman and Car
# Nicobar lands at ~1.10e6 m northing in this raster's CRS (or 10°
# latitude if the raster is in geographic degrees).
SPLIT_Y = 10.0 if is_geo else 1_100_000.0

and_mask = ys_h >= SPLIT_Y
nic_mask = ys_h <  SPLIT_Y

def _bnds_for(y_min, y_max):
    out = []
    for xs, ys in bnd_world:
        ys_arr = np.asarray(ys)
        if ys_arr.mean() >= y_min and ys_arr.mean() < y_max:
            out.append((xs, ys))
    return out

bnd_and = _bnds_for(SPLIT_Y,       float('inf'))
bnd_nic = _bnds_for(float('-inf'), SPLIT_Y)

# Shared log colour scale so both panels are directly comparable.
def _approx_vmax(x, y, n=55):
    if len(x) == 0:
        return 1.0
    x_span = float(np.ptp(x))
    y_span = float(np.ptp(y))
    ny = max(1, int(n * (y_span / max(x_span, 1.0))))
    H, _, _ = np.histogram2d(x, y, bins=[n, ny])
    return max(float(H.max()), 1.0)

vmax = max(
    _approx_vmax(xs_h[and_mask], ys_h[and_mask]),
    _approx_vmax(xs_h[nic_mask], ys_h[nic_mask]),
)
norm = LogNorm(vmin=1, vmax=vmax)

fig, (ax_a, ax_n) = plt.subplots(
    1, 2, figsize=(14, 11), facecolor=BG,
    gridspec_kw={'wspace': 0.18},
)

def _draw_panel(ax, xm, bnds, title):
    ax.set_facecolor('#b3d9f2')
    hx = ax.hexbin(
        xs_h[xm], ys_h[xm], gridsize=45, cmap='YlOrRd',
        norm=norm, mincnt=1, linewidths=0.3, edgecolors='white',
    )
    for xs, ys in bnds:
        ax.plot(xs, ys, color='#333333', lw=1.0, alpha=0.8, zorder=10)
    ax.set_aspect('equal')
    if bnds:
        bx = np.concatenate([np.asarray(xs) for xs, _ in bnds])
        by = np.concatenate([np.asarray(ys) for _, ys in bnds])
    else:
        bx, by = xs_h[xm], ys_h[xm]
    fx = np.concatenate([bx, xs_h[xm]])
    fy = np.concatenate([by, ys_h[xm]])
    px = (fx.max() - fx.min()) * 0.05
    py = (fy.max() - fy.min()) * 0.03
    ax.set_xlim(fx.min() - px, fx.max() + px)
    ax.set_ylim(fy.min() - py, fy.max() + py)
    if is_geo:
        ax.set_xlabel(x_label, fontsize=11, color=TC)
        ax.set_ylabel(y_label, fontsize=11, color=TC)
    else:
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f'{v/1000:.0f} km'))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f'{v/1000:.0f} km'))
        ax.set_xlabel('Easting (km)',  fontsize=11, color=TC)
        ax.set_ylabel('Northing (km)', fontsize=11, color=TC)
    ax.set_title(title, fontsize=12, fontweight='bold', color=TC, pad=8)
    return hx

hx_a = _draw_panel(ax_a, and_mask, bnd_and, 'Andaman Islands')
_     = _draw_panel(ax_n, nic_mask, bnd_nic, 'Nicobar Islands')

cb = fig.colorbar(hx_a, ax=[ax_a, ax_n], fraction=0.025, pad=0.02)
cb.set_label('Hotspot Pixel Density (log scale)', fontsize=10)

fig.suptitle(
    'Spatial Hexbin Clustering of Triple-Collapse Hotspots\n'
    'Andaman vs. Nicobar Islands',
    fontsize=14, fontweight='bold', color=TC, y=0.98,
)
fig.savefig(FIG / 'hotspot_hexbin_density_map.png',
            dpi=180, bbox_inches='tight', facecolor=BG)
plt.close()
print("✅  hotspot_hexbin_density_map.png")

# ══════════════════════════════════════════════════════════════════════
# FIGS 4 & 5 — stat_distribution_*
# ══════════════════════════════════════════════════════════════════════
def violin_figure(entries, ylabel, title, outname, ylim_top=None,
                  log_y=False, log_y_min=0.3, log_y_max=600):
    fig, ax = plt.subplots(figsize=(12, 7), facecolor=BG)
    ax.set_facecolor('#f7f7f7')
    vals_list  = [e[2] for e in entries]
    colors     = [e[1] for e in entries]
    labels     = [e[0] for e in entries]
    positions  = list(range(len(entries)))

    # For log axes we plot log10(values) so the KDE bandwidth lives in log
    # space — that's what gives a properly shaped violin on a log y-axis.
    # Tick labels are then formatted back to the original (un-logged) values.
    if log_y:
        plot_vals = [np.log10(np.clip(v, 1e-3, None)) for v in vals_list]
    else:
        plot_vals = vals_list

    parts = ax.violinplot(plot_vals, positions=positions,
                          showmedians=True, showextrema=False)
    # Mark small-sample classes (n < 1000) with a dashed outline so the
    # violin shape isn't read with the same weight as a 50k-pixel class.
    LOW_N = 1000
    for pc, col, vals in zip(parts['bodies'], colors, vals_list):
        pc.set_facecolor(col); pc.set_alpha(0.72); pc.set_edgecolor('#333')
        if len(vals) < LOW_N:
            pc.set_linestyle('--')
            pc.set_linewidth(1.4)
            pc.set_alpha(0.45)
    parts['cmedians'].set_color('#111'); parts['cmedians'].set_linewidth(2.5)

    for i, vals in enumerate(vals_list):
        med = float(np.median(vals))
        y_med = np.log10(max(med, 1e-3)) if log_y else med
        ax.text(i, y_med, f'{med:.2f}', ha='center', va='bottom', fontsize=9,
                fontweight='bold', color='#111',
                bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.15', edgecolor='none'))
        n_px = len(vals)
        if n_px >= 1000:
            n_label = f'n={n_px/1000:.0f}k' if n_px >= 10_000 else f'n={n_px/1000:.1f}k'
        else:
            n_label = f'n={n_px:,}'
        if log_y:
            span = np.log10(log_y_max) - np.log10(log_y_min)
            y_n  = np.log10(log_y_min) - 0.04 * span
        else:
            y_n = -0.04*ylim_top if ylim_top else 0
        ax.text(i, y_n, n_label, ha='center', va='top', fontsize=8, color='#666')

    ax.set_xticks(positions); ax.set_xticklabels(labels, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=12, color=TC)
    ax.set_title(title, fontsize=13, fontweight='bold', color=TC)
    ax.grid(axis='y', alpha=0.35, linestyle='--')
    if log_y:
        ax.set_ylim(np.log10(log_y_min), np.log10(log_y_max))
        ticks = [t for t in [0.1, 1, 10, 100, 1000] if log_y_min <= t <= log_y_max]
        ax.set_yticks([np.log10(t) for t in ticks])
        ax.set_yticklabels([f'{t:g}' for t in ticks])
    elif ylim_top:
        ax.set_ylim(0, ylim_top)
    fig.tight_layout()
    fig.savefig(FIG / outname, dpi=180, bbox_inches='tight', facecolor=BG)
    plt.close()

print("\nRendering: stat_distribution_habitat_quality.png …")
hq_entries = []
for cls, lbl, col in LC_DEFS:
    vals = sample_class(hq, cls)
    # Keep Q==0 pixels (Built-up has Q≈0 by InVEST construction — dropping
    # them silently removes the class from the chart, which hides the very
    # contrast the figure is meant to communicate).
    vals = vals[~np.isnan(vals) & (vals >= 0) & (vals <= 1)]
    if len(vals) >= 50:
        hq_entries.append((lbl, col, vals))
# Sort ascending by median for consistency with the soil-erosion chart.
hq_entries.sort(key=lambda e: float(np.median(e[2])))
violin_figure(hq_entries,
              'Habitat Quality Index (Q)',
              'Statistical Distribution of Habitat Quality by Land Cover\nAndaman & Nicobar Islands (2024)',
              'stat_distribution_habitat_quality.png', ylim_top=1.10)
print("✅  stat_distribution_habitat_quality.png")

print("\nRendering: stat_distribution_soil_erosion.png …")
se_entries = []
for cls, lbl, col in LC_DEFS:
    vals = sample_class(se, cls)
    # log axis requires strictly positive values; this also drops Built-up,
    # which is uniformly 0 under RUSLE (impervious → no sediment yield) and
    # contributes nothing informative to a per-class erosion distribution.
    vals = vals[~np.isnan(vals) & (vals > 0)]
    if len(vals) >= 50:
        se_entries.append((lbl, col, vals))
# Sort classes by median erosion so the chart reads naturally left→right.
se_entries.sort(key=lambda e: float(np.median(e[2])))
violin_figure(se_entries,
              'Annual Soil Loss (t / ha / yr, log scale)',
              'Statistical Density of Soil Erosion Vulnerability by Land Cover\nAndaman & Nicobar Islands (2024)',
              'stat_distribution_soil_erosion.png',
              log_y=True, log_y_min=0.3, log_y_max=600)
print("✅  stat_distribution_soil_erosion.png")

# ══════════════════════════════════════════════════════════════════════
# FIG 6 — tradeoff_radar_chart
# ══════════════════════════════════════════════════════════════════════
print("\nRendering: tradeoff_radar_chart.png …")

# All three axes are computed from the same rasters the other figures use,
# so the radar can never drift away from the per-pixel data shown elsewhere.
RADAR_CLS = [10, 95, 40, 50, 30, 60]

# Habitat Quality (Q × 100) — per-class mean of habitat_quality_index.tif
hq_vals = {}
for cls in RADAR_CLS:
    v = sample_class(hq, cls)
    v = v[~np.isnan(v) & (v >= 0) & (v <= 1)]
    hq_vals[cls] = float(v.mean() * 100) if len(v) > 0 else 0.0

# Carbon Storage — per-class mean AGB, normalised to the empirical class max
# (not a hand-set 200 Mg/ha reference, which previously pushed Tree past 100).
agb_raw = {}
for cls in RADAR_CLS:
    v = sample_class(agb, cls)
    v = v[~np.isnan(v)]
    agb_raw[cls] = float(v.mean()) if len(v) > 0 else 0.0
agb_max = max(agb_raw.values()) or 1.0
agb_vals = {cls: (v / agb_max * 100) for cls, v in agb_raw.items()}

# Soil Retention — derived from the real RUSLE raster (the same one the
# violin chart uses), normalised against the worst class (Bare/Sparse).
# Built-up: no override. With Built-up impervious in RUSLE, mean erosion ≈ 0
# → retention ≈ 100. That reads as "no sediment yield", which is honest about
# the model output even if a richer ecosystem-service framing would discount
# it (there is no soil there to retain in the first place).
erosion_raw = {}
for cls in RADAR_CLS:
    v = sample_class(se, cls)
    v = v[~np.isnan(v) & (v >= 0)]
    erosion_raw[cls] = float(v.mean()) if len(v) > 0 else 0.0
erosion_max = max(erosion_raw.values()) or 1.0
soil_ret = {cls: 100 - (v / erosion_max * 100) for cls, v in erosion_raw.items()}

radar_classes = [
    (10, 'Tree Cover', '#2E7D32', 0.25),
    (95, 'Mangroves',  '#00838F', 0.30),
    (40, 'Cropland',   '#FF8F00', 0.35),
    (50, 'Built-up',   '#D32F2F', 0.35),
]
categories = ['Habitat\nQuality', 'Carbon\nStorage', 'Soil\nRetention']
N      = len(categories)
angles = [n / N * 2 * np.pi for n in range(N)] + [0]   # close

fig, ax = plt.subplots(figsize=(10, 9), subplot_kw=dict(polar=True), facecolor=BG)
ax.set_facecolor('white')   # ← pure white (was #f5f5f5)

# Per-class radial nudge offsets so labels don't pile on top of each other
# Each entry: [nudge for HQ axis, nudge for Carbon axis, nudge for Soil axis]
label_nudge = {
    10: [ 6,  6,  6],   # Tree Cover  — push outward
    95: [-8, -8, -8],   # Mangroves   — pull slightly inward
    40: [ 5,  5,  5],   # Cropland
    50: [-6, -6, -6],   # Built-up
}

for cls, label, color, alpha in radar_classes:
    vals = [hq_vals[cls], agb_vals[cls], soil_ret[cls]] + [hq_vals[cls]]
    ax.plot(angles, vals, color=color, lw=2.5, label=label)
    ax.fill(angles, vals, color=color, alpha=alpha)
    # Place value labels with per-class radial nudge to avoid overlap
    nudges = label_nudge.get(cls, [5, 5, 5])
    for ang, val, nudge in zip(angles[:-1], vals[:-1], nudges):
        r_pos = max(val + nudge, 2)   # never go below centre
        ax.text(ang, r_pos, f'{val:.0f}', fontsize=9, ha='center', va='center',
                color=color, fontweight='bold',
                bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.1',
                          edgecolor='none'))

ax.set_xticks(angles[:-1])
ax.set_xticklabels(categories, fontsize=12, fontweight='bold', color=TC)
ax.set_ylim(0, 115)
ax.set_yticks([20, 40, 60, 80, 100])
ax.set_yticklabels(['20', '40', '60', '80', '100'], fontsize=8, color='#777')
ax.grid(color='#cccccc', linewidth=0.8)
ax.spines['polar'].set_color('#cccccc')
# Legend placed inside the figure bbox so it doesn't clip
ax.legend(loc='upper left', bbox_to_anchor=(-0.18, 1.18), fontsize=11,
          framealpha=0.95, edgecolor='#ccc')
ax.set_title('Multidimensional Ecosystem Services Trade-Offs\nby Land-Cover Class — ANI (2024)',
             fontsize=13, fontweight='bold', color=TC, pad=28)
fig.savefig(FIG / 'tradeoff_radar_chart.png', dpi=180, bbox_inches='tight', facecolor=BG)
plt.close()
print("✅  tradeoff_radar_chart.png")

print("\n✅  All 6 synthesis figures regenerated on light background.")
