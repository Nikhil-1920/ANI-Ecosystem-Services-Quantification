"""
ANI Ecosystem Services — Split-Region Maps (Andaman / Nicobar Separately)
===========================================================================
Renders each major spatial raster as TWO large, high-resolution maps:
  • Left panel  → Andaman Islands (north group)
  • Right panel → Nicobar Islands (south group)

The Ten Degree Channel (a 138-km ocean gap between rows 11662 and 16288)
provides the natural split boundary.

Maps generated (all saved to figures/split_maps/):
  1.  andaman_nicobar_esa_landcover.png       — ESA WorldCover land cover
  2.  andaman_nicobar_agb_baseline.png        — GEDI AGB biomass density
  3.  andaman_nicobar_habitat_quality.png     — Habitat Quality Index (Q)
  4.  andaman_nicobar_habitat_delta.png       — Δ Habitat Quality 2000→2024
  5.  andaman_nicobar_carbon_loss.png         — Carbon loss hotspots
  6.  andaman_nicobar_soil_loss.png           — RUSLE soil erosion risk
  7.  andaman_nicobar_eci_hotspots.png        — Ecosystem Collapse Index
  8.  andaman_nicobar_forest_loss.png         — GFW deforestation year
  9.  andaman_nicobar_rainfall.png            — CHIRPS annual rainfall

Run:  ./venv/venv/bin/python src/render_split_maps.py
"""

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import rasterio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm, ListedColormap
from pathlib import Path

# ─────────────────────────────────────────────────────────────────────
# Paths
# ─────────────────────────────────────────────────────────────────────
ROOT    = Path(__file__).parent.parent
PROC    = ROOT / 'data' / 'processed'
RES     = ROOT / 'results'
FIG_OUT = ROOT / 'figures' / 'split_maps'
FIG_OUT.mkdir(parents=True, exist_ok=True)

SEP = '=' * 65

# ─────────────────────────────────────────────────────────────────────
# GEOGRAPHIC SPLIT  (EPSG:32646 / UTM 46N)
# Full extent: rows 0-25499, northing 1511840 → 746870
# Largest ocean gap: rows 11662–16288  (Ten Degree Channel, ~138 km)
# ─────────────────────────────────────────────────────────────────────
ANDAMAN_ROW_END   = 11750   # Include a small ocean buffer below last Andaman pixel
NICOBAR_ROW_START = 16200   # Include a small ocean buffer above first Nicobar pixel

# Extra padding (rows) added around each island group bounding box
PAD_ROWS = 80    # ~2.4 km top / bottom
PAD_COLS = 80    # ~2.4 km left / right


# ─────────────────────────────────────────────────────────────────────
# UTILITY FUNCTIONS
# ─────────────────────────────────────────────────────────────────────
def load(path: Path, band: int = 1) -> tuple:
    """Load a GeoTIFF as a float array. Returns (array, profile)."""
    if not path.exists():
        print(f"  ⚠  Not found: {path.name}")
        return None, None
    with rasterio.open(path) as src:
        arr     = src.read(band).astype(float)
        profile = src.profile.copy()
        nd      = src.nodata
    if nd is not None:
        arr[arr == nd] = np.nan
    arr[arr < -1e9] = np.nan
    return arr, profile


def split_regions(arr: np.ndarray):
    """Split a full-extent array into Andaman (north) and Nicobar (south) crops.

    Returns (andaman_crop, nicobar_crop) – both with land-tight bounding boxes
    and PAD_ROWS/PAD_COLS padding.
    """
    if arr is None:
        return None, None

    H, W = arr.shape

    # ── Andaman crop (rows 0 → ANDAMAN_ROW_END) ────────────────────────
    ani_slice = arr[:ANDAMAN_ROW_END, :]
    land_a    = ~np.isnan(ani_slice)
    if land_a.any():
        r_a, c_a  = np.where(land_a)
        r0a = max(0, r_a.min() - PAD_ROWS)
        r1a = min(ani_slice.shape[0], r_a.max() + PAD_ROWS + 1)
        c0a = max(0, c_a.min() - PAD_COLS)
        c1a = min(W, c_a.max() + PAD_COLS + 1)
        and_crop  = ani_slice[r0a:r1a, c0a:c1a]
    else:
        and_crop  = ani_slice

    # ── Nicobar crop (rows NICOBAR_ROW_START → end) ────────────────────
    nic_slice = arr[NICOBAR_ROW_START:, :]
    land_n    = ~np.isnan(nic_slice)
    if land_n.any():
        r_n, c_n  = np.where(land_n)
        r0n = max(0, r_n.min() - PAD_ROWS)
        r1n = min(nic_slice.shape[0], r_n.max() + PAD_ROWS + 1)
        c0n = max(0, c_n.min() - PAD_COLS)
        c1n = min(W, c_n.max() + PAD_COLS + 1)
        nic_crop  = nic_slice[r0n:r1n, c0n:c1n]
    else:
        nic_crop  = nic_slice

    return and_crop, nic_crop


def side_by_side(
    and_arr, nic_arr,
    cmap, vmin, vmax,
    title_and, title_nic,
    suptitle,
    outname,
    cbar_label    = '',
    ocean_color   = '#0a1929',
    fig_bg        = '#0d0d1a',
    title_color   = 'white',
    dpi           = 250,
    cbar_ticks    = None,
    cbar_ticklabs = None,
    extra_legend  = None,
    norm          = None,
):
    """
    Render two large maps side-by-side (Andaman | Nicobar).

    Parameters
    ----------
    and_arr, nic_arr   : 2-D numpy arrays (already cropped)
    cmap               : matplotlib colormap
    vmin, vmax         : colour scale limits
    title_and/nic      : per-panel subtitles
    suptitle           : overall figure title
    outname            : output filename (stem only, no extension)
    cbar_label         : colour-bar axis label
    ocean_color        : background fill for NaN (ocean)
    fig_bg             : figure background colour
    dpi                : output resolution
    cbar_ticks/labs    : optional tick overrides for colour bar
    extra_legend       : list of matplotlib Patch objects for a legend
    norm               : optional BoundaryNorm (overrides vmin/vmax)
    """
    if cmap is None:
        return

    # Compute figure aspect from the two panels
    h_a, w_a = and_arr.shape if and_arr is not None else (1, 1)
    h_n, w_n = nic_arr.shape if nic_arr is not None else (1, 1)
    max_h    = max(h_a, h_n)

    # Each panel gets width proportional to its column count; min 5 in wide
    panel_width = 8.5   # inches per panel
    fig_h       = panel_width * (max_h / max(w_a, w_n, 1)) * 1.05
    fig_h       = max(fig_h, 7.0)

    fig, axes = plt.subplots(
        1, 2,
        figsize=(panel_width * 2 + 1.2, fig_h),
        facecolor=fig_bg,
        gridspec_kw={'wspace': 0.08},
    )

    ims = []
    for ax, arr, subtitle in zip(axes, [and_arr, nic_arr], [title_and, title_nic]):
        ax.set_facecolor(ocean_color)
        if arr is not None:
            masked = np.ma.masked_where(np.isnan(arr), arr)
            kw = dict(interpolation='nearest', aspect='auto')
            if norm is not None:
                im = ax.imshow(masked, cmap=cmap, norm=norm, **kw)
            else:
                im = ax.imshow(masked, cmap=cmap, vmin=vmin, vmax=vmax, **kw)
            ims.append(im)
        ax.set_title(subtitle, color=title_color, fontsize=12,
                     fontweight='bold', pad=9)
        ax.axis('off')

    # Shared colorbar (attached to right panel)
    if ims:
        cbar = fig.colorbar(
            ims[-1], ax=axes, fraction=0.025, pad=0.02,
            shrink=0.80, aspect=28,
        )
        cbar.set_label(cbar_label, color=title_color, fontsize=10.5)
        cbar.ax.tick_params(colors=title_color, labelsize=8.5)
        plt.setp(cbar.ax.yaxis.get_ticklabels(), color=title_color)
        if cbar_ticks is not None:
            cbar.set_ticks(cbar_ticks)
        if cbar_ticklabs is not None:
            cbar.set_ticklabels(cbar_ticklabs)

    if extra_legend:
        axes[1].legend(
            handles=extra_legend, loc='lower right', fontsize=9,
            framealpha=0.85, facecolor='#14142a', edgecolor='#555',
            labelcolor=title_color,
        )

    fig.suptitle(suptitle, color=title_color,
                 fontsize=14, fontweight='bold', y=1.01)
    out_path = FIG_OUT / f'{outname}.png'
    fig.savefig(out_path, dpi=dpi, bbox_inches='tight', facecolor=fig_bg)
    plt.close(fig)
    print(f'  ✅  Saved → figures/split_maps/{outname}.png')


# ─────────────────────────────────────────────────────────────────────
# LOAD ALL RASTERS
# ─────────────────────────────────────────────────────────────────────
print(SEP)
print('  ANI Split-Region Map Renderer')
print(f'  Output → {FIG_OUT}')
print(SEP)

print('\nLoading rasters …')
esa,    esa_p    = load(PROC / 'ANI_ESA_WorldCover_mosaic_clipped.tif')
agb,    _        = load(PROC / 'ANI_GEDI_Biomass_Density_clipped.tif')
gfw,    _        = load(PROC / 'ANI_GFW_Forest_Loss_2001_2023_clipped.tif')
chirps, _        = load(PROC / 'ANI_CHIRPS_Annual_Total_Precip_clipped.tif')
hq,     _        = load(RES  / 'habitat_quality_map.tif')
hq_delta, _      = load(RES  / 'habitat_quality_delta_map.tif')
soil,   _        = load(RES  / 'rusle_soil_loss.tif')
eci,    _        = load(RES  / 'eci_collapse_hotspots.tif')

# Clean AGB
if agb is not None:
    agb = np.where((agb > 0) & (agb <= 800), agb, np.nan)

# Align all to ESA grid (all should match, but be safe)
def align(arr, ref):
    if arr is None or ref is None:
        return arr
    if arr.shape == ref.shape:
        return arr
    out = np.full(ref.shape, np.nan)
    r, c = min(arr.shape[0], ref.shape[0]), min(arr.shape[1], ref.shape[1])
    out[:r, :c] = arr[:r, :c]
    return out

for name in ['agb', 'gfw', 'chirps', 'hq', 'hq_delta', 'soil', 'eci']:
    globals()[name] = align(globals()[name], esa)

print('  All rasters loaded and aligned.\n')


# ─────────────────────────────────────────────────────────────────────
# HELPER: split once and reuse
# ─────────────────────────────────────────────────────────────────────
# Pre-split the ESA land mask so we can use it to build ocean-masked
# versions of any raster quickly.
def masked_by_esa(arr):
    """Return arr with ocean pixels (NaN in ESA) set to NaN."""
    if arr is None or esa is None:
        return arr
    return np.where(np.isnan(esa), np.nan, arr)


# ═════════════════════════════════════════════════════════════════════
# MAP 1 — ESA WorldCover Land Cover
# ═════════════════════════════════════════════════════════════════════
print('─' * 40)
print('MAP 1 — ESA WorldCover Land Cover')

ESA_CLASSES = {
    10: ('#1b7837', 'Tree Cover'),
    20: ('#a6d96a', 'Shrubland'),
    30: ('#ffffbf', 'Grassland'),
    40: ('#fdae61', 'Cropland'),
    50: ('#d7191c', 'Built-up'),
    60: ('#cab2d6', 'Bare / Sparse'),
    80: ('#4393c3', 'Open Water'),
    90: ('#74c476', 'Wetland'),
    95: ('#00441b', 'Mangroves'),
   100: ('#f7f7f7', 'Moss / Lichen'),
}
esa_classes_sorted = sorted(ESA_CLASSES.keys())
esa_colors  = [ESA_CLASSES[k][0] for k in esa_classes_sorted]
esa_cmap    = ListedColormap(esa_colors)
esa_bounds  = [k - 5 for k in esa_classes_sorted] + [esa_classes_sorted[-1] + 5]
esa_norm    = BoundaryNorm(esa_bounds, len(esa_colors))
esa_legend  = [
    mpatches.Patch(facecolor=ESA_CLASSES[k][0], label=ESA_CLASSES[k][1])
    for k in esa_classes_sorted if k in ESA_CLASSES
]

esa_and, esa_nic = split_regions(esa)
side_by_side(
    esa_and, esa_nic,
    cmap=esa_cmap, vmin=None, vmax=None, norm=esa_norm,
    title_and='Andaman Islands',
    title_nic='Nicobar Islands',
    suptitle='ESA WorldCover Land Cover — Andaman & Nicobar Islands (2021)',
    outname='andaman_nicobar_esa_landcover',
    cbar_label='',
    extra_legend=esa_legend,
    ocean_color='#0a1929',
)


# ═════════════════════════════════════════════════════════════════════
# MAP 2 — GEDI AGB Biomass Density
# ═════════════════════════════════════════════════════════════════════
print('─' * 40)
print('MAP 2 — GEDI L4B Aboveground Biomass')

agb_masked = masked_by_esa(agb)
agb_and, agb_nic = split_regions(agb_masked)

cmap_agb = plt.cm.YlGn.copy()
cmap_agb.set_bad(color='#0a1929')
side_by_side(
    agb_and, agb_nic,
    cmap=cmap_agb, vmin=0, vmax=350,
    title_and='Andaman Islands',
    title_nic='Nicobar Islands',
    suptitle='GEDI L4B Aboveground Biomass Density — Andaman & Nicobar Islands',
    outname='andaman_nicobar_agb_baseline',
    cbar_label='Aboveground Biomass (Mg/ha)',
    ocean_color='#0a1929',
)


# ═════════════════════════════════════════════════════════════════════
# MAP 3 — Habitat Quality Index
# ═════════════════════════════════════════════════════════════════════
print('─' * 40)
print('MAP 3 — Habitat Quality Index (Q)')

cmap_hq = LinearSegmentedColormap.from_list(
    'habitat',
    ['#7f0000', '#d73027', '#fdae61', '#ffffbf', '#a6d96a', '#1a9641', '#003d00'],
    N=256,
)
cmap_hq.set_bad(color='#0a1929')

hq_masked = masked_by_esa(hq)
hq_and, hq_nic = split_regions(hq_masked)

hq_ticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
hq_labs  = ['0.0\nCritical', '0.2', '0.4', '0.6', '0.8', '1.0\nPristine']
side_by_side(
    hq_and, hq_nic,
    cmap=cmap_hq, vmin=0, vmax=1,
    title_and='Andaman Islands',
    title_nic='Nicobar Islands',
    suptitle='InVEST Habitat Quality Index (Q) — Andaman & Nicobar Islands (2024)',
    outname='andaman_nicobar_habitat_quality',
    cbar_label='Habitat Quality Index (0 = Critical  →  1 = Pristine)',
    cbar_ticks=hq_ticks, cbar_ticklabs=hq_labs,
    ocean_color='#0a1929',
)


# ═════════════════════════════════════════════════════════════════════
# MAP 4 — Delta Habitat Quality (2000 → 2024)
# ═════════════════════════════════════════════════════════════════════
print('─' * 40)
print('MAP 4 — Δ Habitat Quality 2000 → 2024')

if hq_delta is not None:
    hq_delta_masked = masked_by_esa(hq_delta)
    land_vals = hq_delta_masked[~np.isnan(hq_delta_masked)]
    if len(land_vals) > 0:
        abs_max = max(abs(np.percentile(land_vals, 1)),
                      abs(np.percentile(land_vals, 99)))
        abs_max = max(abs_max, 0.01)
    else:
        abs_max = 0.1

    cmap_delta = LinearSegmentedColormap.from_list(
        'delta',
        ['#b2182b', '#ef8a62', '#fddbc7', '#f7f7f7', '#d1e5f0', '#4393c3', '#053061'],
        N=256,
    )
    cmap_delta.set_bad(color='#0a1929')
    hd_and, hd_nic = split_regions(hq_delta_masked)
    side_by_side(
        hd_and, hd_nic,
        cmap=cmap_delta, vmin=-abs_max, vmax=abs_max,
        title_and='Andaman Islands',
        title_nic='Nicobar Islands',
        suptitle='Δ Habitat Quality Change (2000 → 2024) — Red = Degraded  |  Blue = Improved',
        outname='andaman_nicobar_habitat_delta',
        cbar_label='ΔQ (negative = habitat loss)',
        ocean_color='#0a1929',
    )
else:
    print('  ⚠  habitat_quality_delta_map.tif not found — skipping MAP 4.')


# ═════════════════════════════════════════════════════════════════════
# MAP 5 — Carbon Loss Hotspots (GFW × GEDI)
# ═════════════════════════════════════════════════════════════════════
print('─' * 40)
print('MAP 5 — Carbon Loss Hotspots')

if gfw is not None and agb is not None:
    IPCC_C = 0.47
    loss_mask   = (gfw >= 1) & (gfw <= 24)
    carbon_map  = np.where(loss_mask, agb * IPCC_C, np.nan)
    carbon_map  = masked_by_esa(carbon_map)

    cmap_c = plt.cm.hot_r.copy()
    cmap_c.set_bad(color='#0a1929')
    cm_and, cm_nic = split_regions(carbon_map)
    side_by_side(
        cm_and, cm_nic,
        cmap=cmap_c, vmin=0, vmax=200,
        title_and='Andaman Islands',
        title_nic='Nicobar Islands',
        suptitle='Carbon Loss Hotspots 2001–2024 — Andaman & Nicobar Islands',
        outname='andaman_nicobar_carbon_loss',
        cbar_label='Carbon Lost (MgC/ha)',
        ocean_color='#0a1929',
    )
else:
    print('  ⚠  GFW or AGB layer missing — skipping MAP 5.')


# ═════════════════════════════════════════════════════════════════════
# MAP 6 — RUSLE Soil Erosion Risk
# ═════════════════════════════════════════════════════════════════════
print('─' * 40)
print('MAP 6 — RUSLE Soil Erosion Risk')

if soil is not None:
    soil_masked = masked_by_esa(soil)
    soil_log    = np.log1p(soil_masked)
    vmax_s      = np.nanpercentile(soil_log, 98)

    cmap_s = plt.cm.YlOrRd.copy()
    cmap_s.set_bad(color='#0a1929')
    soil_and, soil_nic = split_regions(soil_log)
    side_by_side(
        soil_and, soil_nic,
        cmap=cmap_s, vmin=0, vmax=vmax_s,
        title_and='Andaman Islands',
        title_nic='Nicobar Islands',
        suptitle='RUSLE Annual Soil Erosion Risk — Andaman & Nicobar Islands (2024)',
        outname='andaman_nicobar_soil_loss',
        cbar_label='Soil Loss  log(1 + A)  [A in t/ha/yr]',
        ocean_color='#0a1929',
    )
else:
    print('  ⚠  rusle_soil_loss.tif not found — skipping MAP 6.')


# ═════════════════════════════════════════════════════════════════════
# MAP 7 — Ecosystem Collapse Index (ECI)
# ═════════════════════════════════════════════════════════════════════
print('─' * 40)
print('MAP 7 — Ecosystem Collapse Index (ECI)')

if eci is not None:
    eci_masked = masked_by_esa(eci)
    eci_and, eci_nic = split_regions(eci_masked)

    fire_colors = [
        '#0d0d1a', '#2c1e45', '#7b2a59',
        '#cd3e45', '#f97e20', '#ffe600', '#ffffff',
    ]
    cmap_eci = LinearSegmentedColormap.from_list('eci_fire', fire_colors, N=256)
    cmap_eci.set_bad(color='#0a1929')

    side_by_side(
        eci_and, eci_nic,
        cmap=cmap_eci, vmin=0, vmax=1,
        title_and='Andaman Islands',
        title_nic='Nicobar Islands',
        suptitle='Ecosystem Collapse Index (ECI) — Triple Risk Convergence (2000–2024)',
        outname='andaman_nicobar_eci_hotspots',
        cbar_label='ECI  (0 = Stable  →  1 = Collapse)',
        ocean_color='#0a1929',
    )
else:
    print('  ⚠  eci_collapse_hotspots.tif not found — skipping MAP 7.')


# ═════════════════════════════════════════════════════════════════════
# MAP 8 — GFW Deforestation Year
# ═════════════════════════════════════════════════════════════════════
print('─' * 40)
print('MAP 8 — GFW Forest Loss Year')

if gfw is not None:
    # Show only pixels where deforestation occurred (year code 1–23 = 2001–2023)
    gfw_loss = np.where((gfw >= 1) & (gfw <= 23), gfw + 2000, np.nan)
    gfw_loss  = masked_by_esa(gfw_loss)

    cmap_gfw = plt.cm.plasma.copy()
    cmap_gfw.set_bad(color='#0a1929')
    gfw_and, gfw_nic = split_regions(gfw_loss)
    side_by_side(
        gfw_and, gfw_nic,
        cmap=cmap_gfw, vmin=2001, vmax=2023,
        title_and='Andaman Islands',
        title_nic='Nicobar Islands',
        suptitle='GFW Forest Loss Year (2001–2023) — Andaman & Nicobar Islands',
        outname='andaman_nicobar_forest_loss',
        cbar_label='Year of Forest Loss',
        ocean_color='#0a1929',
    )
else:
    print('  ⚠  GFW layer missing — skipping MAP 8.')


# ═════════════════════════════════════════════════════════════════════
# MAP 9 — CHIRPS Annual Rainfall
# ═════════════════════════════════════════════════════════════════════
print('─' * 40)
print('MAP 9 — CHIRPS Annual Rainfall')

if chirps is not None:
    chirps_masked = masked_by_esa(chirps)
    vmin_r = np.nanpercentile(chirps_masked[~np.isnan(chirps_masked)], 2)
    vmax_r = np.nanpercentile(chirps_masked[~np.isnan(chirps_masked)], 98)

    cmap_r = plt.cm.Blues.copy()
    cmap_r.set_bad(color='#0a1929')
    rain_and, rain_nic = split_regions(chirps_masked)
    side_by_side(
        rain_and, rain_nic,
        cmap=cmap_r, vmin=vmin_r, vmax=vmax_r,
        title_and='Andaman Islands',
        title_nic='Nicobar Islands',
        suptitle='CHIRPS Annual Rainfall — Andaman & Nicobar Islands (Mean 2000–2023)',
        outname='andaman_nicobar_rainfall',
        cbar_label='Annual Precipitation (mm/yr)',
        ocean_color='#0a1929',
    )
else:
    print('  ⚠  CHIRPS layer missing — skipping MAP 9.')


# ─────────────────────────────────────────────────────────────────────
print()
print(SEP)
print('  ✅  All split-region maps saved to:')
print(f'      {FIG_OUT}')
print(SEP)
