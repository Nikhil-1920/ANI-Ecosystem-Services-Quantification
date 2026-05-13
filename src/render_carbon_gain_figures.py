"""
Quick-render: Carbon Gain vs Loss Figures
==========================================
Loads pre-computed rasters and renders the 3 carbon gain figures
without re-running the full pipeline. Takes ~60 seconds.
"""
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import rasterio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
from matplotlib.patches import Patch

PROC_DIR       = Path('data/processed')
FIG_DIR        = Path('figures/carbon')
RES_DIR        = Path('results')
FIG_DIR.mkdir(parents=True, exist_ok=True)

PIXEL_AREA_HA  = 0.09
IPCC_CARBON_F  = 0.47
CO2E_MW_RATIO  = 44 / 12
ROOT_SHOOT_R   = 0.24


def load(fname):
    p = PROC_DIR / fname
    if not p.exists():
        print(f'  ⚠  Not found: {fname}')
        return None
    with rasterio.open(p) as src:
        arr = src.read(1).astype(float)
        nd  = src.nodata
    if nd is not None:
        arr[arr == nd] = np.nan
    arr[arr < -1e9] = np.nan
    return arr


print('Loading rasters …')
agb_grid      = load('ANI_GEDI_Biomass_Density_clipped.tif')
defor_grid    = load('ANI_GFW_Forest_Loss_2001_2023_clipped.tif')
gain_grid     = load('ANI_GFW_Forest_Gain_clipped.tif')

# Sanitize
agb_grid  = np.where((agb_grid  is not None) & (agb_grid  >= 0) & (agb_grid  <= 800),
                     agb_grid,  np.nan) if agb_grid  is not None else None
print('Rasters loaded.')

# ── Annual loss CSV (already computed) ────────────────────────────────
csv_path = RES_DIR / 'carbon_annual_loss_by_year.csv'
annual_df = pd.read_csv(csv_path) if csv_path.exists() else None

# ── Gain stats ────────────────────────────────────────────────────────
gain_stats   = None
gain_px_mask = None

if gain_grid is not None and agb_grid is not None:
    gain_px_mask = (gain_grid == 1)
    # Only count gain pixels that overlap with valid GEDI data
    gain_agb     = agb_grid[gain_px_mask]
    gain_agb     = gain_agb[~np.isnan(gain_agb)]
    gain_area_ha = gain_px_mask.sum() * PIXEL_AREA_HA
    gain_agb_mg  = gain_agb.sum() * PIXEL_AREA_HA
    gain_bio_mg  = gain_agb_mg * (1 + ROOT_SHOOT_R)
    gain_carbon  = (gain_bio_mg * IPCC_CARBON_F) / 1e3   # GgC
    gain_co2e    = gain_carbon * CO2E_MW_RATIO            # GgCO2e
    gain_stats = dict(area_ha=gain_area_ha, carbon_ggc=gain_carbon, co2e_ggco2e=gain_co2e)
    print(f'Gain area : {gain_area_ha:,.0f} ha')
    print(f'Gain CO2e : {gain_co2e:.3f} GgCO2e')


# ════════════════════════════════════════════════════════
# FIGURE 1 — Net Carbon Balance (bars + cumulative lines)
# ════════════════════════════════════════════════════════
if annual_df is not None and gain_stats is not None:
    df = annual_df[annual_df['area_ha'] > 0].copy()
    total_loss  = df['co2e_ggco2e'].sum()
    total_gain  = gain_stats['co2e_ggco2e']
    net_flux    = total_loss - total_gain
    n_yrs       = len(df)
    # GFW gain is cumulative 2000-2012 — distribute proportionally
    # to loss years as a visual proxy (makes bars visible for comparison)
    gain_per_yr = total_gain / n_yrs

    fig, axes = plt.subplots(1, 2, figsize=(16, 6), facecolor='#0d0d1a')
    bg = '#14142a'

    # Left — grouped annual bars
    ax = axes[0]
    ax.set_facecolor(bg)
    x     = df['year'].values
    width = 0.38
    ax.bar(x - width/2, df['co2e_ggco2e'],
           width=width, color='#ef5350', alpha=0.88,
           label='Gross Loss (deforestation)')
    ax.bar(x + width/2, [gain_per_yr] * n_yrs,
           width=width, color='#66BB6A', alpha=0.88,
           label=f'Gain (cumul. total = {total_gain:.2f} GgCO₂e,\nshown as equal annual share)')
    ax.set_xlabel('Year', color='white', fontsize=10)
    ax.set_ylabel('GgCO₂e / yr', color='white', fontsize=10)
    ax.tick_params(colors='white', axis='both')
    ax.tick_params(axis='x', rotation=45)
    ax.spines[:].set_color('#444')
    ax.legend(framealpha=0.25, labelcolor='white', fontsize=7.5, loc='upper right')
    ax.set_title('Annual Gross Loss  vs  Forest Gain\n'
                 '(Gain = GFW cumulative 2000–2012, distributed equally for visualisation)',
                 color='white', fontsize=10, fontweight='bold')

    # Right — cumulative lines
    ax2 = axes[1]
    ax2.set_facecolor(bg)
    cum_loss = df['co2e_ggco2e'].cumsum().values
    cum_gain = np.array([gain_per_yr * (i + 1) for i in range(n_yrs)])
    cum_net  = cum_loss - cum_gain
    ax2.fill_between(x, cum_loss, alpha=0.22, color='#ef5350')
    ax2.plot(x, cum_loss, color='#ef5350', lw=2.2, label='Cumulative Loss')
    ax2.fill_between(x, cum_gain, alpha=0.22, color='#66BB6A')
    ax2.plot(x, cum_gain, color='#66BB6A', lw=2.2, label='Cumulative Gain')
    ax2.plot(x, cum_net,  color='#4FC3F7', lw=2.5, marker='o', ms=3.5,
             label='Net Flux (Loss − Gain)')
    ax2.axhline(0, color='white', lw=0.7, ls='--', alpha=0.35)
    ax2.set_xlabel('Year', color='white', fontsize=10)
    ax2.set_ylabel('Cumulative GgCO₂e', color='white', fontsize=10)
    ax2.tick_params(colors='white', axis='both')
    ax2.tick_params(axis='x', rotation=45)
    ax2.spines[:].set_color('#444')
    ax2.legend(framealpha=0.25, labelcolor='white', fontsize=8)
    ax2.set_title(f'Cumulative Net Flux\nNet = {net_flux:.2f} GgCO₂e  |  Gain offsets {total_gain/total_loss*100:.1f}% of Loss',
                  color='white', fontsize=10, fontweight='bold')

    fig.tight_layout()
    out = FIG_DIR / 'carbon_net_balance_figure.png'
    fig.savefig(out, dpi=200, bbox_inches='tight', facecolor='#0d0d1a')
    plt.close()
    print(f'✅  Saved → {out}')


# ════════════════════════════════════════════════════════
# FIGURE 2 — Spatial: Loss vs Gain side-by-side
# ════════════════════════════════════════════════════════
if defor_grid is not None and agb_grid is not None and gain_px_mask is not None:
    defor_int   = np.nan_to_num(defor_grid, nan=0).astype(int)
    loss_pixels = (defor_int >= 1) & (defor_int <= 24)
    carbon_map  = np.where(loss_pixels, agb_grid * IPCC_CARBON_F, np.nan)

    h, w = carbon_map.shape
    gain_aligned = gain_px_mask[:h, :w]

    # Binary RGBA overlays (no AGB lookup for gain — young regrowth has near-NaN GEDI)
    loss_rgba = np.zeros((*carbon_map.shape, 4), dtype=float)
    loss_rgba[loss_pixels, 0] = 0.90
    loss_rgba[loss_pixels, 1] = 0.20
    loss_rgba[loss_pixels, 2] = 0.20
    loss_rgba[loss_pixels, 3] = 0.85

    gain_rgba = np.zeros((*carbon_map.shape, 4), dtype=float)
    gain_rgba[gain_aligned, 0] = 0.10
    gain_rgba[gain_aligned, 1] = 0.88
    gain_rgba[gain_aligned, 2] = 0.35
    gain_rgba[gain_aligned, 3] = 0.85

    fig, axes = plt.subplots(1, 2, figsize=(15, 10), facecolor='#0d0d1a')

    # Left — loss colormap
    im0 = axes[0].imshow(carbon_map, cmap='hot_r', vmin=0, vmax=200)
    axes[0].set_title('Carbon Loss Hotspots\n(Deforestation 2001–2024)',
                      color='white', fontsize=11, fontweight='bold')
    axes[0].axis('off')
    cbar0 = plt.colorbar(im0, ax=axes[0], fraction=0.035, pad=0.03)
    cbar0.set_label('Carbon Lost (MgC/ha)', color='white', fontsize=8)
    cbar0.ax.yaxis.set_tick_params(color='white', labelsize=7)
    plt.setp(cbar0.ax.yaxis.get_ticklabels(), color='white')

    # Right — binary overlay
    axes[1].set_facecolor('#0d0d1a')
    axes[1].imshow(loss_rgba)
    axes[1].imshow(gain_rgba)
    axes[1].set_title(
        'Forest Loss (red) vs Forest Gain (green)\n'
        'GFW: Loss 2001–2024  |  Gain 2000–2012',
        color='white', fontsize=11, fontweight='bold')
    axes[1].axis('off')
    axes[1].legend(handles=[
        Patch(facecolor='#16E05A', alpha=0.9,
              label=f'Forest Gain  ({int(gain_aligned.sum()):,} px = {gain_aligned.sum()*PIXEL_AREA_HA:,.0f} ha)'),
        Patch(facecolor='#E63333', alpha=0.85,
              label=f'Forest Loss  ({int(loss_pixels.sum()):,} px = {loss_pixels.sum()*PIXEL_AREA_HA:,.0f} ha)'),
    ], loc='lower right', framealpha=0.35, labelcolor='white', fontsize=8)

    fig.suptitle('Carbon Forest Loss vs Forest Gain — Andaman & Nicobar Islands',
                 color='white', fontsize=12, fontweight='bold')
    fig.tight_layout()
    out = FIG_DIR / 'carbon_loss_vs_gain_spatial.png'
    fig.savefig(out, dpi=200, bbox_inches='tight', facecolor='#0d0d1a')
    plt.close()
    print(f'✅  Saved → {out}')


# ════════════════════════════════════════════════════════
# FIGURE 3 — Baseline AGB map with gain overlay
# ════════════════════════════════════════════════════════
if agb_grid is not None and gain_px_mask is not None:
    h, w = agb_grid.shape
    gain_aligned = gain_px_mask[:h, :w]

    fig, ax = plt.subplots(figsize=(6, 9), facecolor='#0d0d1a')
    ax.set_facecolor('#0d0d1a')
    im   = ax.imshow(agb_grid, cmap='YlGn', vmin=0, vmax=400)
    cbar = plt.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label('Aboveground Biomass (Mg/ha)', color='white', fontsize=9)
    cbar.ax.yaxis.set_tick_params(color='white', labelsize=8)
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color='white')

    # Cyan overlay for gain pixels
    gain_rgba = np.zeros((*agb_grid.shape, 4), dtype=float)
    gain_rgba[gain_aligned, 0] = 0.00
    gain_rgba[gain_aligned, 1] = 0.90
    gain_rgba[gain_aligned, 2] = 0.90
    gain_rgba[gain_aligned, 3] = 0.75
    ax.imshow(gain_rgba)

    ax.set_title('GEDI L4B AGB Baseline\n+ Forest Gain pixels (cyan) — ANI',
                 color='white', fontsize=10, fontweight='bold', pad=8)
    ax.legend(handles=[
        Patch(facecolor='#00E5E5', alpha=0.85,
              label=f'GFW Forest Gain 2000–2012\n({gain_aligned.sum()*PIXEL_AREA_HA:,.0f} ha)'),
    ], loc='lower right', framealpha=0.3, labelcolor='white', fontsize=8)
    ax.axis('off')
    fig.tight_layout()
    out = FIG_DIR / 'agb_gedi_baseline_map.png'
    fig.savefig(out, dpi=150, bbox_inches='tight', facecolor='#0d0d1a')
    plt.close()
    print(f'✅  Saved → {out}')


# ════════════════════════════════════════════════════════
# Update CSV with correct gain columns
# ════════════════════════════════════════════════════════
if annual_df is not None and gain_stats is not None:
    annual_df['cumulative_gain_co2e_ggco2e'] = round(gain_stats['co2e_ggco2e'], 6)
    annual_df['cumulative_gain_area_ha']     = round(gain_stats['area_ha'], 2)
    net_total = annual_df['co2e_ggco2e'].sum() - gain_stats['co2e_ggco2e']
    annual_df['net_co2e_total_ggco2e']       = round(net_total, 6)
    csv_out = RES_DIR / 'carbon_annual_loss_by_year.csv'
    annual_df.to_csv(csv_out, index=False)
    print(f'✅  CSV updated → {csv_out}')
    print(f'\n   Gross Loss : {annual_df["co2e_ggco2e"].sum():.2f} GgCO₂e')
    print(f'   Gain       : {gain_stats["co2e_ggco2e"]:.2f} GgCO₂e')
    print(f'   Net Flux   : {net_total:.2f} GgCO₂e')
    print(f'   Gain offsets {gain_stats["co2e_ggco2e"]/annual_df["co2e_ggco2e"].sum()*100:.1f}% of gross loss')

print('\nAll done!')
