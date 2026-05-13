"""
render_carbon_v4.py  —  Fixes:
  1. Maps: vmax = 95th percentile (124 MgCha) so colors spread properly
  2. Timeseries: gain shown correctly — only after 2012 as a step-up on cumulative
  3. Net balance panel: starts from 0, gain applied only when it was detected (post-2012)
"""
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import rasterio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import geopandas as gpd
from pathlib import Path
from matplotlib.patches import Patch

PROC_DIR      = Path('data/processed')
RAW_DIR       = Path('data/raw')
FIG_DIR       = Path('figures/carbon')
RES_DIR       = Path('results')
FIG_DIR.mkdir(parents=True, exist_ok=True)

PIXEL_AREA_HA = 0.09
IPCC_CARBON_F = 0.47
CO2E_MW_RATIO = 44 / 12
ROOT_SHOOT_R  = 0.24

BG     = 'white'
AX_BG  = '#f7f7f7'
TC     = '#1a1a2e'
LOSS_C = '#C62828'
GAIN_C = '#1B5E20'
NET_C  = '#0D47A1'
CUM_C  = '#880E4F'

plt.rcParams.update({
    'font.family'     : 'DejaVu Sans',
    'axes.facecolor'  : AX_BG,
    'figure.facecolor': BG,
    'axes.labelcolor' : TC,
    'xtick.color'     : TC,
    'ytick.color'     : TC,
    'text.color'      : TC,
    'grid.color'      : '#cccccc',
    'grid.alpha'      : 0.5,
})


def load(fname):
    p = PROC_DIR / fname
    if not p.exists():
        print(f'  ⚠  Not found: {fname}')
        return None, None
    with rasterio.open(p) as src:
        arr     = src.read(1).astype(float)
        nd      = src.nodata
        profile = src.profile
    if nd is not None:
        arr[arr == nd] = np.nan
    arr[arr < -1e9] = np.nan
    return arr, profile


def load_boundary(profile):
    shp = gpd.read_file(RAW_DIR / 'ANI_Administrative_Boundary.shp')
    return shp.to_crs(profile['crs'])


def boundary_coords(gdf, profile):
    coords = []
    for geom in gdf.geometry:
        if geom is None:
            continue
        polys = list(geom.geoms) if geom.geom_type == 'MultiPolygon' else [geom]
        for poly in polys:
            xs, ys = poly.exterior.xy
            rows, cols = rasterio.transform.rowcol(profile['transform'], xs, ys)
            coords.append((list(cols), list(rows)))
    return coords


def draw_border(ax, coords, lw=1.4, color='#222222', alpha=0.85):
    for c, r in coords:
        ax.plot(c, r, color=color, lw=lw, alpha=alpha, zorder=10)


# ── Load ──────────────────────────────────────────────────────────────
print('Loading rasters …')
agb_grid,   agb_p   = load('ANI_GEDI_Biomass_Density_clipped.tif')
defor_grid, defor_p = load('ANI_GFW_Forest_Loss_2001_2023_clipped.tif')
gain_grid,  gain_p  = load('ANI_GFW_Forest_Gain_clipped.tif')

if agb_grid is not None:
    agb_grid = np.where((agb_grid >= 0) & (agb_grid <= 800), agb_grid, np.nan)

bnd_gdf  = load_boundary(agb_p)
bnd_coords = boundary_coords(bnd_gdf, agb_p)
print(f'Boundary: {len(bnd_coords)} polygons')

annual_df = pd.read_csv(RES_DIR / 'carbon_annual_loss_by_year.csv')

# ── Gain stats ────────────────────────────────────────────────────────
gain_px_mask = None
gain_stats   = None
if gain_grid is not None and agb_grid is not None:
    gain_px_mask = (gain_grid == 1)
    g_agb        = agb_grid[gain_px_mask]
    g_agb        = g_agb[~np.isnan(g_agb)]
    gain_area_ha = gain_px_mask.sum() * PIXEL_AREA_HA
    gain_bio_mg  = g_agb.sum() * PIXEL_AREA_HA * (1 + ROOT_SHOOT_R)
    gain_carbon  = (gain_bio_mg * IPCC_CARBON_F) / 1e3
    gain_co2e    = gain_carbon * CO2E_MW_RATIO
    gain_stats   = dict(area_ha=gain_area_ha, carbon_ggc=gain_carbon, co2e_ggco2e=gain_co2e)
    print(f'Gain: {gain_area_ha:,.0f} ha  |  {gain_co2e:.2f} GgCO₂e')

# ── Carbon loss map ───────────────────────────────────────────────────
defor_int   = np.nan_to_num(defor_grid, nan=0).astype(int)
loss_pixels = (defor_int >= 1) & (defor_int <= 24)
carbon_map  = np.where(loss_pixels, agb_grid * IPCC_CARBON_F, np.nan)

# Percentile-based vmax so colors SPREAD across actual data range
carbon_vals = carbon_map[~np.isnan(carbon_map)]
VMAX = float(np.percentile(carbon_vals, 95))   # = ~124 MgC/ha
print(f'Carbon vmax (p95): {VMAX:.1f} MgC/ha  (mean: {carbon_vals.mean():.1f})')


# ════════════════════════════════════════════════════════════════════════
# FIG 1 — AGB Baseline + Orange Gain + Boundary
# ════════════════════════════════════════════════════════════════════════
print('\nFIG 1: AGB Baseline + Gain …')
fig, ax = plt.subplots(figsize=(9, 12), facecolor=BG)
ax.set_facecolor('#dce9dc')

# Use percentile vmax for AGB too so greens spread
agb_vals = agb_grid[~np.isnan(agb_grid)]
agb_vmax = float(np.percentile(agb_vals, 97))
print(f'AGB vmax (p97): {agb_vmax:.1f} Mg/ha')

im   = ax.imshow(agb_grid, cmap='YlGn', vmin=0, vmax=agb_vmax, alpha=0.92)
cbar = plt.colorbar(im, ax=ax, fraction=0.028, pad=0.02, shrink=0.85)
cbar.set_label(f'AGB (Mg/ha, capped at p97={agb_vmax:.0f})', fontsize=9, color=TC)
plt.setp(cbar.ax.yaxis.get_ticklabels(), color=TC)

if gain_px_mask is not None:
    h, w = agb_grid.shape
    gan  = gain_px_mask[:h, :w]
    ov   = np.zeros((*agb_grid.shape, 4), dtype=float)
    ov[gan, 0] = 1.00   # bright orange — contrasts green AGB
    ov[gan, 1] = 0.45
    ov[gan, 2] = 0.00
    ov[gan, 3] = 0.90
    ax.imshow(ov, zorder=5)

draw_border(ax, bnd_coords)
ax.set_title('GEDI L4B AGB Baseline + Forest Gain (orange)\nAndaman & Nicobar Islands',
             fontsize=12, fontweight='bold', color=TC, pad=10)
ax.legend(handles=[
    Patch(facecolor='#FF7300', alpha=0.90,
          label=f'GFW Forest Gain 2000–2012\n'
                f'{gain_area_ha:,.0f} ha  ·  {gain_co2e:.1f} GgCO₂e sequestered'),
], loc='lower left', framealpha=0.9, fontsize=9)
ax.axis('off')
fig.tight_layout()
fig.savefig(FIG_DIR / 'agb_gedi_baseline_map.png', dpi=180, bbox_inches='tight', facecolor=BG)
plt.close()
print('✅  agb_gedi_baseline_map.png')


# ════════════════════════════════════════════════════════════════════════
# FIG 2 — Carbon Loss Hotspots (standalone, correct vmax, boundary)
# ════════════════════════════════════════════════════════════════════════
print('\nFIG 2: Loss hotspot map …')
fig, ax = plt.subplots(figsize=(8, 12), facecolor=BG)
ax.set_facecolor('#dce9dc')
im   = ax.imshow(carbon_map, cmap='YlOrRd', vmin=0, vmax=VMAX, alpha=0.93)
cbar = plt.colorbar(im, ax=ax, fraction=0.028, pad=0.02, shrink=0.82)
cbar.set_label(f'Carbon Lost (MgC/ha, p95={VMAX:.0f})', fontsize=9)
plt.setp(cbar.ax.yaxis.get_ticklabels(), color=TC)
draw_border(ax, bnd_coords)
ax.set_title('Carbon Loss Hotspots (2001–2024)\nAndaman & Nicobar Islands',
             fontsize=12, fontweight='bold', color=TC, pad=10)
ax.axis('off')
fig.tight_layout()
fig.savefig(FIG_DIR / 'carbon_loss_hotspots_map.png', dpi=200, bbox_inches='tight', facecolor=BG)
plt.close()
print('✅  carbon_loss_hotspots_map.png')


# ════════════════════════════════════════════════════════════════════════
# FIG 3 — Spatial: Loss (YlOrRd) + Gain (ORANGE overlay) + Boundary
# ════════════════════════════════════════════════════════════════════════
print('\nFIG 3: Spatial loss + gain overlay …')
if gain_px_mask is not None:
    h, w = carbon_map.shape
    gan  = gain_px_mask[:h, :w]

    ov_gain = np.zeros((*carbon_map.shape, 4), dtype=float)
    ov_gain[gan, 0] = 1.00
    ov_gain[gan, 1] = 0.45
    ov_gain[gan, 2] = 0.00
    ov_gain[gan, 3] = 0.90

    fig, ax = plt.subplots(figsize=(9, 13), facecolor=BG)
    ax.set_facecolor('#dce9dc')
    im = ax.imshow(carbon_map, cmap='YlOrRd', vmin=0, vmax=VMAX, alpha=0.90, zorder=2)
    cbar = plt.colorbar(im, ax=ax, fraction=0.028, pad=0.02, shrink=0.80)
    cbar.set_label(f'Carbon Lost (MgC/ha, p95={VMAX:.0f})', fontsize=9)
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color=TC)
    ax.imshow(ov_gain, zorder=3)
    draw_border(ax, bnd_coords)
    ax.legend(handles=[
        Patch(facecolor='#C62828', alpha=0.85,
              label=f'Carbon Loss  {loss_pixels.sum()*PIXEL_AREA_HA:,.0f} ha  ·  {annual_df["co2e_ggco2e"].sum():.0f} GgCO₂e'),
        Patch(facecolor='#FF7300', alpha=0.90,
              label=f'Forest Gain  {gan.sum()*PIXEL_AREA_HA:,.0f} ha  ·  {gain_co2e:.0f} GgCO₂e seq.')
    ], loc='lower left', framealpha=0.92, fontsize=9)
    ax.set_title(
        'Carbon Loss + Forest Gain Overlay\nANI  (Loss 2001–2024  |  Gain 2000–2012)',
        fontsize=12, fontweight='bold', color=TC, pad=10)
    ax.axis('off')
    fig.tight_layout()
    fig.savefig(FIG_DIR / 'carbon_loss_vs_gain_spatial.png', dpi=200, bbox_inches='tight', facecolor=BG)
    plt.close()
print('✅  carbon_loss_vs_gain_spatial.png')


# ════════════════════════════════════════════════════════════════════════
# FIG 4 — Timeseries + Net Balance (CORRECTED gain placement)
#
# GFW gain layer = cumulative binary mask, 2000–2012.
# We don't know the year-breakdown, so we show gain correctly as:
#   • Top panel: annual loss bars + cumulative loss curve
#               + green dashed line at total gain = 508 on SAME axis as cumulative
#               + annotation "by 2012, ANI forests regrew 508 GgCO₂e"
#   • Bottom: net = cumulative loss curve.
#             We add gain as a STEP UP at year 2012 (one-time credit).
#             Before 2012: cumulative loss only.
#             At 2012: subtract gain credit → net drops.
#             After 2012: net = cumulative loss - gain (flat credit).
# ════════════════════════════════════════════════════════════════════════
print('\nFIG 4: Timeseries …')
if annual_df is not None and gain_stats is not None:
    df         = annual_df.copy()
    df_active  = df[df['area_ha'] > 0].copy()
    total_loss = df['co2e_ggco2e'].sum()
    total_gain = gain_stats['co2e_ggco2e']

    all_years  = df['year'].values
    all_co2e   = df['co2e_ggco2e'].values
    cum_loss   = np.cumsum(all_co2e)

    # Net cumulative: before 2013, no gain credit. From 2013 onwards, subtract total gain.
    net_cum = np.where(all_years >= 2013, cum_loss - total_gain, cum_loss)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10), facecolor=BG, sharex=True)

    # ─ TOP: Annual bars + cumulative ────────────────────────────────────
    ax1.set_facecolor(AX_BG)
    ax1.bar(all_years, all_co2e, color=LOSS_C, alpha=0.80,
            label='Annual Deforestation Emissions',
            edgecolor='white', linewidth=0.4, zorder=3)
    ax1.set_ylabel('Annual CO₂e Loss (GgCO₂e)', fontsize=11, color=LOSS_C)
    ax1.tick_params(axis='y', colors=LOSS_C)
    ax1.spines['left'].set_color(LOSS_C)
    ax1.grid(axis='y', zorder=0)

    ax1b = ax1.twinx()
    ax1b.plot(all_years, cum_loss, color=CUM_C, lw=2.5, marker='o', ms=3.5,
              label='Cumulative Loss', zorder=4)
    # Show gain as a vertical drop at 2013 (when GFW gain credit kicks in)
    ax1b.annotate('',
                  xy     =(2013, cum_loss[list(all_years).index(2013)] - total_gain),
                  xytext =(2013, cum_loss[list(all_years).index(2013)]),
                  arrowprops=dict(arrowstyle='->', color=GAIN_C, lw=2.0))
    ax1b.annotate(
        f'Forest Gain credit applied\n(2000–2012 = {total_gain:.0f} GgCO₂e seq.)',
        xy=(2013, cum_loss[list(all_years).index(2013)] - total_gain/2),
        xytext=(2015.5, 750),
        fontsize=9, color=GAIN_C, fontweight='bold',
        arrowprops=dict(arrowstyle='->', color=GAIN_C, lw=1.2),
        bbox=dict(boxstyle='round,pad=0.35', facecolor='#e8f5e9', alpha=0.9, edgecolor=GAIN_C)
    )
    ax1b.set_ylabel('Cumulative CO₂e (GgCO₂e)', fontsize=11, color=CUM_C)
    ax1b.tick_params(axis='y', colors=CUM_C)
    ax1b.spines['right'].set_color(CUM_C)

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax1b.get_legend_handles_labels()
    ax1.legend(h1+h2, l1+l2, fontsize=9, framealpha=0.92, loc='upper left')
    ax1.set_title('Annual Deforestation Emissions + Cumulative Loss\n'
                  '(GFW Forest Gain 2000–2012 credit shown at 2013)',
                  fontsize=11, fontweight='bold')

    # ─ BOTTOM: Net cumulative (0 before gain, drops at 2013) ────────────
    ax2.set_facecolor(AX_BG)
    ax2.fill_between(all_years, cum_loss, alpha=0.12, color=LOSS_C)
    ax2.plot(all_years, cum_loss, color=LOSS_C, lw=2.5,
             label=f'Gross Cumulative Loss = {total_loss:.0f} GgCO₂e')

    ax2.fill_between(all_years, net_cum, alpha=0.12, color=NET_C)
    ax2.plot(all_years, net_cum, color=NET_C, lw=2.8, marker='s', ms=3.5,
             label=f'Net Flux (after gain credit at 2013) = {net_cum[-1]:.0f} GgCO₂e')

    # Shade the gap = gain credit zone
    ax2.fill_between(all_years, cum_loss, net_cum,
                     where=(all_years >= 2013),
                     alpha=0.18, color=GAIN_C,
                     label=f'Gain credit  = {total_gain:.0f} GgCO₂e  ({total_gain/total_loss*100:.1f}% offset)')

    ax2.axhline(0, color='black', lw=0.8, ls=':', alpha=0.4)
    ax2.set_xlabel('Year', fontsize=11)
    ax2.set_ylabel('Cumulative GgCO₂e', fontsize=11)
    ax2.tick_params(axis='x', rotation=45)
    ax2.grid(axis='y')
    ax2.legend(fontsize=9, framealpha=0.92, loc='upper left')
    ax2.set_title(
        f'Net Cumulative Carbon Flux — Gain credit applied from 2013\n'
        f'Gain offsets {total_gain/total_loss*100:.1f}% of total gross deforestation emissions',
        fontsize=11, fontweight='bold')

    fig.suptitle('Carbon Emission & Sequestration Balance — Andaman & Nicobar Islands (2001–2024)',
                 fontsize=13, fontweight='bold', y=1.005)
    plt.tight_layout()
    fig.savefig(FIG_DIR / 'carbon_annual_loss_timeseries.png',
                dpi=200, bbox_inches='tight', facecolor=BG)
    plt.close()
    print('✅  carbon_annual_loss_timeseries.png')


# ════════════════════════════════════════════════════════════════════════
# FIG 5 — Net Balance Summary Bars
# ════════════════════════════════════════════════════════════════════════
print('\nFIG 5: Net balance summary …')
if annual_df is not None and gain_stats is not None:
    total_loss = annual_df['co2e_ggco2e'].sum()
    total_gain = gain_stats['co2e_ggco2e']
    net_flux   = total_loss - total_gain

    fig, ax = plt.subplots(figsize=(9, 6), facecolor=BG)
    ax.set_facecolor(AX_BG)
    cats   = ['Gross\nDeforestation\nEmissions', 'Forest Gain\nSequestration\n(2000–2012)', 'Net Carbon\nFlux']
    vals   = [total_loss, total_gain, net_flux]
    colors = [LOSS_C, GAIN_C, NET_C]
    bars   = ax.bar(cats, vals, color=colors, alpha=0.85,
                    edgecolor='white', linewidth=0.8, width=0.5)
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 12,
                f'{val:.1f}\nGgCO₂e', ha='center', fontsize=11.5,
                fontweight='bold', color=TC)
    ax.set_ylabel('GgCO₂e', fontsize=12)
    ax.set_title(
        'Carbon Balance Summary — Andaman & Nicobar Islands\n'
        f'Forest Gain offsets {total_gain/total_loss*100:.1f}% of gross emissions',
        fontsize=11, fontweight='bold')
    ax.grid(axis='y', alpha=0.5)
    ax.set_ylim(0, max(vals) * 1.25)
    fig.tight_layout()
    fig.savefig(FIG_DIR / 'carbon_net_balance_figure.png', dpi=200, bbox_inches='tight', facecolor=BG)
    plt.close()
    print('✅  carbon_net_balance_figure.png')

print('\n✅  All figures complete.')
