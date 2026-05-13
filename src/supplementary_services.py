"""
ANI Ecosystem Services — Supplementary Modules (Lightweight)
============================================================
Quantifies four additional ecosystem services that the main report
acknowledged as out-of-scope, using only datasets already in
data/processed/.  None require external downloads.

Modules:
  1. Coastal protection      (mangrove × elevation × distance-from-coast)
  2. Freshwater yield        (Budyko-style P − ET, per WorldCover class)
  3. Pollination suitability (Lonsdorf-style index from WorldCover)
  4. Soil organic carbon     (SoilGrids clay+silt × GFW loss × Don 2011 25%)

Each module produces a single number for the report's "lower-bound"
caveat and a one-panel figure.  Combined results upgrade the BAU 2040
economic-damage forecast from $502M to an upper-bound estimate.

Inputs : data/processed/*.tif
         results/carbon_annual_loss_by_year.csv
         results/economic_damage_2040_forecast.csv

Outputs: results/supplementary_services.csv
         results/supplementary_services_summary.json
         results/economic_damage_2040_upperbound.csv
         figures/supplementary/supplementary_services_panel.png
"""

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import ndimage as ndi

warnings.filterwarnings("ignore")

ROOT = Path(__file__).parent.parent
PROC = ROOT / "data" / "processed"
RES = ROOT / "results"
FIG = ROOT / "figures" / "supplementary"
FIG.mkdir(parents=True, exist_ok=True)

PIXEL_HA = 0.09          # 30 m × 30 m
PIXEL_KM2 = 0.0009
USD_PER_TC = 51.0        # Social Cost of Carbon (US EPA 2023 SC-CO2 → SC-C)
USD_PER_TCO2 = 51.0      # commonly cited 2023 SCC value (lower bound)

# ESA WorldCover class codes
TREE = 10; SHRUB = 20; GRASS = 30; CROP = 40; BUILT = 50
BARE = 60; WATER = 80; WETLAND = 90; MANGROVE = 95


def _read(path: Path, band: int = 1):
    with rasterio.open(path) as s:
        a = s.read(band).astype(np.float32)
        if s.nodata is not None:
            a = np.where(a == s.nodata, np.nan, a)
        prof = s.profile
    return a, prof


def _save_geotiff(arr: np.ndarray, ref_profile: dict, out_path: Path):
    p = ref_profile.copy()
    p.update(dtype="float32", count=1, compress="deflate", nodata=-9999.0)
    arr = np.where(np.isnan(arr), -9999.0, arr).astype(np.float32)
    with rasterio.open(out_path, "w", **p) as dst:
        dst.write(arr, 1)


# ──────────────────────────────────────────────────────────────────────
# 1.  COASTAL PROTECTION
# ──────────────────────────────────────────────────────────────────────
def module_coastal_protection():
    """Mangrove-based storm-surge / tsunami buffer index.

    CPI per pixel = mangrove_presence × elev_factor × coast_factor
      elev_factor   = 1 if elev <= 10 m, decays linearly to 0 at 30 m
      coast_factor  = 1 if dist_to_water <= 1 km, decays to 0 at 5 km

    Captures the idea that mangroves close to the shore on low-lying
    land provide the most surge buffering (Spalding et al. 2014;
    Menendez et al. 2020).
    """
    print("\n[1] Coastal protection (mangrove buffer index)")
    lc, prof = _read(PROC / "ANI_ESA_WorldCover_mosaic_clipped.tif")
    dem, _ = _read(PROC / "ANI_SRTM_DEM_30m_clipped.tif")
    loss, _ = _read(PROC / "ANI_GFW_Forest_Loss_2001_2023_clipped.tif")

    # Crop to common shape
    h = min(lc.shape[0], dem.shape[0], loss.shape[0])
    w = min(lc.shape[1], dem.shape[1], loss.shape[1])
    lc, dem, loss = lc[:h, :w], dem[:h, :w], loss[:h, :w]

    mangrove = (lc == MANGROVE)
    water = np.isin(lc, [WATER]) | (lc == 0) | np.isnan(lc)

    # Distance to nearest water pixel (km).  Each pixel = 30 m.
    print("    computing distance-to-water (Euclidean)...")
    dist_px = ndi.distance_transform_edt(~water)
    dist_km = dist_px * 0.030

    # Coast factor: 1 inside 1 km, linear ramp to 0 at 5 km
    coast_factor = np.clip(1.0 - (dist_km - 1.0) / 4.0, 0.0, 1.0)
    coast_factor = np.where(dist_km <= 1.0, 1.0, coast_factor)

    # Elevation factor: full effect ≤10 m, linear ramp to 0 at 30 m
    elev_factor = np.clip(1.0 - (dem - 10.0) / 20.0, 0.0, 1.0)
    elev_factor = np.where(dem <= 10.0, 1.0, elev_factor)

    cpi = mangrove.astype(np.float32) * elev_factor * coast_factor

    mangrove_area_ha = mangrove.sum() * PIXEL_HA
    protective_area_ha = (cpi > 0.5).sum() * PIXEL_HA
    high_value_area_ha = (cpi > 0.8).sum() * PIXEL_HA

    # Mangrove lost (any year, GFW)
    mangrove_lost_px = mangrove & (loss > 0)
    mangrove_lost_ha = mangrove_lost_px.sum() * PIXEL_HA
    high_value_lost_ha = (mangrove_lost_px & (cpi > 0.5)).sum() * PIXEL_HA

    # Per Menendez et al. 2020 (Sci Rep): mangroves reduce expected annual
    # damages from coastal flooding by ~US$65 B globally over 2.9M ha
    # = ~US$22,400/ha/yr.  Conservative ANI scaling: 1/10th this for
    # storm-only (no tsunami exposure modeling).
    USD_PER_HA_YR = 2240.0
    annual_protection_value_usd = protective_area_ha * USD_PER_HA_YR

    # 15-yr NPV: previously this was an *undiscounted* 15× multiplier,
    # which mislabelled a gross-value figure as NPV. Apply a proper
    # annuity-factor discount at 4% (standard public-sector real rate).
    DISCOUNT_RATE = 0.04
    NPV_YEARS = 15
    annuity_factor = (1 - (1 + DISCOUNT_RATE) ** (-NPV_YEARS)) / DISCOUNT_RATE
    value_lost_usd = high_value_lost_ha * USD_PER_HA_YR * annuity_factor

    print(f"    Total mangrove area   : {mangrove_area_ha:>9,.0f} ha")
    print(f"    Protective (CPI>0.5)  : {protective_area_ha:>9,.0f} ha")
    print(f"    High-value (CPI>0.8)  : {high_value_area_ha:>9,.0f} ha")
    print(f"    Mangrove lost (any)   : {mangrove_lost_ha:>9,.0f} ha")
    print(f"    High-value lost       : {high_value_lost_ha:>9,.0f} ha")
    print(f"    Protection value/yr   : ${annual_protection_value_usd/1e6:>6.1f} M")
    print(f"    NPV value lost (15yr) : ${value_lost_usd/1e6:>6.1f} M")

    _save_geotiff(cpi, prof, RES / "coastal_protection_index.tif")

    return dict(
        mangrove_area_ha=float(mangrove_area_ha),
        protective_area_ha=float(protective_area_ha),
        high_value_area_ha=float(high_value_area_ha),
        mangrove_lost_ha=float(mangrove_lost_ha),
        high_value_lost_ha=float(high_value_lost_ha),
        annual_protection_value_usd=float(annual_protection_value_usd),
        npv_value_lost_usd=float(value_lost_usd),
        cpi_map=cpi,
        landcover=lc,
    )


# ──────────────────────────────────────────────────────────────────────
# 2.  FRESHWATER YIELD  (Budyko-style P − ET)
# ──────────────────────────────────────────────────────────────────────
def module_freshwater_yield():
    """Annual water yield = P − AET, where AET is class-specific.

    AET coefficients (fraction of P consumed by ET) follow tropical
    values from Zhang et al. 2001, Sun et al. 2006:
        Tree cover   0.85    Mangrove 0.90
        Cropland     0.65    Grassland 0.55
        Built-up     0.30    Bare      0.25
    """
    print("\n[2] Freshwater yield (Budyko-style P − ET)")
    lc, prof = _read(PROC / "ANI_ESA_WorldCover_mosaic_clipped.tif")
    chirps, _ = _read(PROC / "ANI_CHIRPS_Annual_Total_Precip_clipped.tif")

    h = min(lc.shape[0], chirps.shape[0])
    w = min(lc.shape[1], chirps.shape[1])
    lc, chirps = lc[:h, :w], chirps[:h, :w]

    et_coef = np.full_like(lc, np.nan, dtype=np.float32)
    coef_map = {TREE: 0.85, MANGROVE: 0.90, SHRUB: 0.70,
                 GRASS: 0.55, CROP: 0.65, BARE: 0.25,
                 WETLAND: 0.95, BUILT: 0.30}
    for cls, c in coef_map.items():
        et_coef[lc == cls] = c

    land_mask = np.isfinite(et_coef) & np.isfinite(chirps)
    aet_mm = chirps * et_coef
    yield_mm = chirps - aet_mm                  # = (1 − coef) × P
    # mm × 30 m × 30 m × 1e−3 m / 1e3 = m³ per pixel; / 1000 = ML; /1e6 = GL
    pixel_area_m2 = 900.0
    yield_m3 = yield_mm * 1e-3 * pixel_area_m2

    total_yield_GL = np.nansum(yield_m3) / 1e6  # GL = 1e9 L = 1e6 m³

    # Per-class yield summary (Shrubland previously omitted from display
    # even though its ET coefficient was applied — adds back here).
    rows = []
    for cls, label in [(TREE,"Tree cover"),(MANGROVE,"Mangrove"),
                       (SHRUB,"Shrubland"),
                       (GRASS,"Grassland"),(CROP,"Cropland"),
                       (BUILT,"Built-up"),(BARE,"Bare"),(WETLAND,"Wetland")]:
        m = (lc == cls) & land_mask
        if m.sum() == 0:
            continue
        rows.append(dict(
            esa_class=cls, label=label,
            area_ha=float(m.sum() * PIXEL_HA),
            mean_P_mm=float(np.nanmean(chirps[m])),
            mean_yield_mm=float(np.nanmean(yield_mm[m])),
            total_yield_GL=float(np.nansum(yield_m3[m]) / 1e6),
        ))
    per_class = pd.DataFrame(rows)
    print(per_class.to_string(index=False))

    # Loss of yield if forested area converted to cropland (counter-factual)
    # Use ESA tree-cover area × yield-per-ha difference (tree→cropland)
    tree_m = (lc == TREE) & land_mask
    yield_under_tree = float(np.nansum(yield_m3[tree_m])) / max(tree_m.sum(), 1) * pixel_area_m2 * 0  # placeholder
    # Compute the delta-yield directly: replace 0.85 with 0.65 on tree pixels
    delta_yield_m3 = chirps[tree_m] * (0.85 - 0.65) * 1e-3 * pixel_area_m2
    incremental_yield_if_converted_GL = float(np.nansum(delta_yield_m3)) / 1e6
    print(f"\n    Total water yield      : {total_yield_GL:,.1f} GL/yr")
    print(f"    Incremental yield      : {incremental_yield_if_converted_GL:,.1f} GL/yr")
    print(f"    (if all forest → crop, more runoff but ↑ erosion/flood risk)")

    _save_geotiff(yield_mm, prof, RES / "freshwater_yield_mm.tif")

    return dict(
        total_yield_GL=float(total_yield_GL),
        incremental_yield_if_converted_GL=float(incremental_yield_if_converted_GL),
        per_class=per_class,
        yield_map=yield_mm,
        landcover=lc,
    )


# ──────────────────────────────────────────────────────────────────────
# 3.  POLLINATION SUITABILITY (Lonsdorf-style)
# ──────────────────────────────────────────────────────────────────────
def module_pollination():
    """Pollinator habitat-suitability index per pixel.

    Lonsdorf et al. 2009: S = sqrt(F × N), where
       F = floral resource (0-1) by class
       N = nesting resource (0-1) by class
    Native bee guild (no foraging-kernel convolution applied — local
    proxy only).  Values from Koh et al. 2016 tropical defaults.
    """
    print("\n[3] Pollination habitat suitability")
    lc, prof = _read(PROC / "ANI_ESA_WorldCover_mosaic_clipped.tif")

    # Lonsdorf-style floral and nesting scores (tropical defaults)
    floral = {TREE: 0.65, MANGROVE: 0.55, SHRUB: 0.80, GRASS: 0.75,
               CROP: 0.40, WETLAND: 0.60, BUILT: 0.10, BARE: 0.05,
               WATER: 0.00}
    nesting = {TREE: 0.90, MANGROVE: 0.60, SHRUB: 0.70, GRASS: 0.50,
                CROP: 0.10, WETLAND: 0.30, BUILT: 0.00, BARE: 0.20,
                WATER: 0.00}

    F = np.zeros_like(lc, dtype=np.float32)
    N = np.zeros_like(lc, dtype=np.float32)
    for cls, v in floral.items():
        F[lc == cls] = v
    for cls, v in nesting.items():
        N[lc == cls] = v

    S = np.sqrt(F * N).astype(np.float32)
    land = np.isin(lc, [TREE, SHRUB, GRASS, CROP, BUILT, BARE,
                          WETLAND, MANGROVE])
    S[~land] = np.nan

    mean_S = float(np.nanmean(S))
    rows = []
    for cls, label in [(TREE,"Tree cover"),(MANGROVE,"Mangrove"),
                       (SHRUB,"Shrubland"),(GRASS,"Grassland"),
                       (CROP,"Cropland"),(BUILT,"Built-up"),
                       (BARE,"Bare"),(WETLAND,"Wetland")]:
        m = (lc == cls)
        if m.sum() == 0:
            continue
        rows.append(dict(
            esa_class=cls, label=label,
            area_ha=float(m.sum() * PIXEL_HA),
            floral=floral.get(cls, 0),
            nesting=nesting.get(cls, 0),
            mean_S=float(np.nanmean(S[m])) if np.isfinite(S[m]).any() else 0.0,
        ))
    per_class = pd.DataFrame(rows)
    print(per_class.to_string(index=False))

    # Counter-factual: pollinator support loss from forest → cropland
    forest_S = np.sqrt(0.65 * 0.90)        # = 0.764
    crop_S = np.sqrt(0.40 * 0.10)          # = 0.200
    delta_per_ha = forest_S - crop_S       # 0.564 unit drop
    # Apply to the documented 19,502 ha deforested
    df = pd.read_csv(RES / "carbon_annual_loss_by_year.csv")
    deforested_ha = float(df["area_ha"].sum())
    pollination_units_lost = deforested_ha * delta_per_ha
    print(f"\n    Mean S (whole ANI)        : {mean_S:.3f}")
    print(f"    Forest→Cropland ΔS        : {delta_per_ha:.3f} per ha")
    print(f"    Pollination units lost    : {pollination_units_lost:,.0f} unit-ha")

    _save_geotiff(S, prof, RES / "pollination_suitability.tif")

    return dict(
        mean_S=float(mean_S),
        forest_S=float(forest_S),
        crop_S=float(crop_S),
        delta_per_ha=float(delta_per_ha),
        deforested_ha=float(deforested_ha),
        pollination_units_lost=float(pollination_units_lost),
        per_class=per_class,
        S_map=S,
        landcover=lc,
    )


# ──────────────────────────────────────────────────────────────────────
# 4.  SOIL ORGANIC CARBON  (Don et al. 2011)
# ──────────────────────────────────────────────────────────────────────
def module_soc_loss():
    """SOC loss from forest conversion.

    Baseline tropical surface (0–30 cm) SOC stock estimated as:
        SOC_Mg_per_ha = 40 + 0.5 × (clay% + silt%)
    (Tropical empirical fit; mean ~65 Mg C/ha, range 40–95).

    Conversion loss factor: 25% (Don et al. 2011 global meta-analysis,
    tropical forest → cropland).
    """
    print("\n[4] Soil organic carbon (SOC) loss")
    lc, prof = _read(PROC / "ANI_ESA_WorldCover_mosaic_clipped.tif")
    loss, _ = _read(PROC / "ANI_GFW_Forest_Loss_2001_2023_clipped.tif")
    clay, _ = _read(PROC / "ANI_SoilGrids_Clay_Silt_GEE_clipped.tif", band=1)
    silt, _ = _read(PROC / "ANI_SoilGrids_Clay_Silt_GEE_clipped.tif", band=2)

    h = min(lc.shape[0], loss.shape[0], clay.shape[0])
    w = min(lc.shape[1], loss.shape[1], clay.shape[1])
    lc, loss, clay, silt = lc[:h, :w], loss[:h, :w], clay[:h, :w], silt[:h, :w]

    # Baseline SOC stock (Mg C/ha) per pixel.
    # The SoilGrids clay+silt export from GEE is mostly zero-filled at
    # 30 m (most pixels have no soil data — the underlying SoilGrids
    # grid is at ~250 m and the resample left gaps). Using the formula
    # 40 + 0.5·(clay+silt) on those pixels collapses to a uniform
    # 40 Mg C/ha floor — implausibly low for tropical forest soils.
    # Hybrid baseline:
    #   • Where SoilGrids reports valid data (clay+silt > 0): keep the
    #     empirical formula — it captures real local soil-texture variance.
    #   • Where it is zero-filled: fall back to a literature class-mean
    #     SOC stock (0–30 cm) by land-cover type.
    # Class means follow Powers et al. 2018, Donato et al. 2011
    # (mangroves), FAO 2017 tropical soil-carbon synthesis.
    CLASS_SOC_BASELINE = {
        TREE:     75.0,   # closed tropical forest
        MANGROVE: 100.0,  # mangrove soils are SOC-rich
        SHRUB:    60.0,
        GRASS:    55.0,
        CROP:     50.0,
        BUILT:    35.0,
        BARE:     30.0,
        WETLAND:  80.0,
    }
    soc_empirical = 40.0 + 0.5 * (clay + silt)
    soc_fallback = np.full_like(lc, np.nan, dtype=np.float32)
    for cls, val in CLASS_SOC_BASELINE.items():
        soc_fallback[lc == cls] = val
    has_soil_data = (clay + silt) > 0
    soc_stock = np.where(has_soil_data, soc_empirical, soc_fallback)

    n_valid_pixels = int(np.nansum(has_soil_data))
    n_total_pixels = int(np.isfinite(soc_stock).sum())
    print(f"    SoilGrids coverage    : {n_valid_pixels:,} / {n_total_pixels:,} "
          f"pixels ({100*n_valid_pixels/max(n_total_pixels,1):.1f}% empirical, "
          f"rest from class-mean fallback)")

    # Pixels lost (any GFW year)
    lost_mask = (loss > 0)
    DON_FACTOR = 0.25         # 25% SOC loss after tropical conversion
    soc_lost_Mg_per_ha = soc_stock * DON_FACTOR
    soc_lost_Mg_total = float(np.nansum(soc_lost_Mg_per_ha[lost_mask] * PIXEL_HA))

    # CO2e equivalence (44/12 ratio)
    co2e_lost_Mg = soc_lost_Mg_total * (44/12)
    co2e_lost_Gg = co2e_lost_Mg / 1e3

    # Economic value @ social cost of carbon
    soc_damage_usd = co2e_lost_Mg * USD_PER_TCO2

    mean_baseline = float(np.nanmean(soc_stock))
    print(f"    Mean baseline SOC     : {mean_baseline:.1f} Mg C/ha")
    print(f"    Deforested pixels     : {lost_mask.sum():,}")
    print(f"    Total SOC lost        : {soc_lost_Mg_total:,.0f} Mg C")
    print(f"    Equivalent CO2e       : {co2e_lost_Gg:,.1f} Gg CO2e")
    print(f"    Economic damage @ SCC : ${soc_damage_usd/1e6:,.1f} M USD")

    _save_geotiff(soc_stock, prof, RES / "soc_baseline_Mg_per_ha.tif")

    return dict(
        mean_baseline_Mg_per_ha=float(mean_baseline),
        deforested_px=int(lost_mask.sum()),
        soc_lost_Mg_C=float(soc_lost_Mg_total),
        co2e_lost_Gg=float(co2e_lost_Gg),
        soc_damage_usd=float(soc_damage_usd),
        soc_map=soc_stock,
        landcover=lc,
    )


# ──────────────────────────────────────────────────────────────────────
# DRIVER
# ──────────────────────────────────────────────────────────────────────
def render_summary_figure(coastal, water, pollin, soc):
    """4-panel synthesis: each module gets one panel."""
    fig = plt.figure(figsize=(15, 11), facecolor="white")
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.22,
                            left=0.06, right=0.97, top=0.88, bottom=0.06)

    TEXT = "#1a1a2e"

    # Panel 1: Coastal protection — bar of area buckets
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_facecolor("#fafbfc")
    cats = ["Total\nmangrove", "Protective\n(CPI > 0.5)",
             "High-value\n(CPI > 0.8)", "Mangrove\nlost", "High-value\nlost"]
    vals = [coastal["mangrove_area_ha"], coastal["protective_area_ha"],
             coastal["high_value_area_ha"], coastal["mangrove_lost_ha"],
             coastal["high_value_lost_ha"]]
    colors = ["#2e7d32", "#388e3c", "#1b5e20", "#c62828", "#7f0000"]
    bars = ax1.bar(cats, vals, color=colors, alpha=0.88, edgecolor="white")
    for b, v in zip(bars, vals):
        ax1.text(b.get_x() + b.get_width() / 2, b.get_height(),
                  f"{v:,.0f}", ha="center", va="bottom", fontsize=9,
                  color=TEXT)
    ax1.set_ylabel("Area (ha)", color=TEXT)
    ax1.set_title(f"Coastal protection — mangrove buffer index\n"
                   f"Annual value: ${coastal['annual_protection_value_usd']/1e6:.1f}M | "
                   f"15-yr NPV lost: ${coastal['npv_value_lost_usd']/1e6:.1f}M",
                   fontweight="bold", color=TEXT, fontsize=11)
    ax1.tick_params(colors=TEXT, labelsize=8)
    ax1.spines[:].set_color("#cccccc")

    # Panel 2: Freshwater yield — bar per land-cover class
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.set_facecolor("#fafbfc")
    pc = water["per_class"].sort_values("total_yield_GL", ascending=True)
    bars = ax2.barh(pc["label"], pc["total_yield_GL"],
                     color="#1976d2", alpha=0.88, edgecolor="white")
    for b, v in zip(bars, pc["total_yield_GL"]):
        ax2.text(b.get_width(), b.get_y() + b.get_height() / 2,
                  f" {v:,.0f} GL", va="center", fontsize=8, color=TEXT)
    ax2.set_xlabel("Annual water yield (GL = 10⁶ m³)", color=TEXT)
    # Total = sum of displayed bars, so the headline always matches what the reader sees.
    panel_total_GL = float(pc["total_yield_GL"].sum())
    ax2.set_title(f"Freshwater yield (Budyko-style P − ET)\n"
                   f"Total (shown classes): {panel_total_GL:,.0f} GL/yr | "
                   f"Forest→Crop: +{water['incremental_yield_if_converted_GL']:,.0f} GL "
                   f"(↑ flood/erosion risk)",
                   fontweight="bold", color=TEXT, fontsize=11)
    ax2.tick_params(colors=TEXT, labelsize=9)
    ax2.spines[:].set_color("#cccccc")

    # Panel 3: Pollination — bar of S by class
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.set_facecolor("#fafbfc")
    pp = pollin["per_class"].sort_values("mean_S", ascending=True)
    bars = ax3.barh(pp["label"], pp["mean_S"],
                     color="#f9a825", alpha=0.92, edgecolor="white")
    for b, v in zip(bars, pp["mean_S"]):
        ax3.text(b.get_width() + 0.01, b.get_y() + b.get_height() / 2,
                  f"{v:.2f}", va="center", fontsize=9, color=TEXT)
    ax3.set_xlim(0, 1.0)
    ax3.set_xlabel("Pollination suitability S = √(F·N)", color=TEXT)
    ax3.set_title(f"Pollination suitability (Lonsdorf 2009)\n"
                   f"Forest→Cropland ΔS: -{pollin['delta_per_ha']:.2f}/ha | "
                   f"Total lost: {pollin['pollination_units_lost']:,.0f} unit-ha",
                   fontweight="bold", color=TEXT, fontsize=11)
    ax3.tick_params(colors=TEXT, labelsize=9)
    ax3.spines[:].set_color("#cccccc")

    # Panel 4: SOC loss — big-number card matching the other three chart panels
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.set_facecolor("#fafbfc")
    ax4.set_xlim(0, 1); ax4.set_ylim(0, 1)
    ax4.set_xticks([]); ax4.set_yticks([])
    ax4.spines[:].set_color("#cccccc")

    pct_of_main = soc['co2e_lost_Gg'] / 1390 * 100

    # Headline number
    ax4.text(0.5, 0.78, f"{soc['co2e_lost_Gg']:,.0f}",
             ha="center", va="center", transform=ax4.transAxes,
             fontsize=44, fontweight="bold", color="#6a1b1b")
    ax4.text(0.5, 0.63, "Gg CO₂e released from soil",
             ha="center", va="center", transform=ax4.transAxes,
             fontsize=11, color=TEXT)
    ax4.text(0.5, 0.56,
             f"(+{pct_of_main:.0f}% on top of the {1390:,} Gg CO₂e "
             f"above-ground loss in the main report)",
             ha="center", va="center", transform=ax4.transAxes,
             fontsize=9, color="#555", style="italic")

    # Supporting figures, three columns
    col_y = 0.34
    for x, (label, value) in zip(
        [0.18, 0.50, 0.82],
        [("Baseline SOC",
          f"{soc['mean_baseline_Mg_per_ha']:.1f} Mg C/ha"),
         ("Total SOC lost",
          f"{soc['soc_lost_Mg_C']/1000:,.0f} kt C"),
         ("Damage @ SCC",
          f"${soc['soc_damage_usd']/1e6:,.1f} M")]
    ):
        ax4.text(x, col_y + 0.07, value, ha="center", va="center",
                 transform=ax4.transAxes,
                 fontsize=14, fontweight="bold", color=TEXT)
        ax4.text(x, col_y - 0.02, label, ha="center", va="center",
                 transform=ax4.transAxes,
                 fontsize=9, color="#555")

    ax4.text(0.5, 0.07,
             "Method: 40 + 0.5(clay% + silt%) baseline · Don et al. 2011 25% loss factor",
             ha="center", va="center", transform=ax4.transAxes,
             fontsize=8, color="#888", style="italic")

    ax4.set_title("Soil organic carbon (SOC) loss from deforestation",
                   fontweight="bold", color=TEXT, fontsize=11)

    fig.suptitle("Supplementary Ecosystem Services — Andaman & Nicobar Islands\n"
                  "Three additional services plus a soil-carbon supplement to the main carbon budget",
                  fontsize=13, fontweight="bold", color=TEXT, y=0.96)
    fig.savefig(FIG / "supplementary_services_panel.png", dpi=200,
                 bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"\n  Saved -> {FIG / 'supplementary_services_panel.png'}")


def update_economic_upper_bound(coastal, soc):
    """Combine main BAU damage + coastal NPV + SOC damage → upper bound."""
    print("\n[5] Updated economic damages (upper bound)")
    base = pd.read_csv(RES / "economic_damage_2040_forecast.csv").iloc[0]
    main_usd = float(base["total_economic_damage_usd"])
    coastal_usd = float(coastal["npv_value_lost_usd"])
    soc_usd = float(soc["soc_damage_usd"])

    upper_usd = main_usd + coastal_usd + soc_usd
    pct_increase = (upper_usd - main_usd) / main_usd * 100.0

    print(f"    Main report (lower bound)  : ${main_usd/1e6:,.1f} M")
    print(f"    + Coastal protection NPV   : ${coastal_usd/1e6:,.1f} M")
    print(f"    + SOC loss @ SCC           : ${soc_usd/1e6:,.1f} M")
    print(f"    UPPER-BOUND TOTAL          : ${upper_usd/1e6:,.1f} M  "
          f"({pct_increase:+.0f}% vs lower bound)")

    pd.DataFrame([dict(
        component="Carbon + sediment + habitat (lower bound, main report)",
        usd=main_usd),
        dict(component="Coastal protection NPV (15-yr)",
              usd=coastal_usd),
        dict(component="SOC loss @ social cost of carbon",
              usd=soc_usd),
        dict(component="UPPER BOUND TOTAL", usd=upper_usd),
    ]).to_csv(RES / "economic_damage_2040_upperbound.csv", index=False)
    print(f"    Saved -> results/economic_damage_2040_upperbound.csv")
    return dict(main_usd=main_usd, coastal_usd=coastal_usd,
                 soc_usd=soc_usd, upper_usd=upper_usd,
                 pct_increase=pct_increase)


def main():
    coastal = module_coastal_protection()
    water = module_freshwater_yield()
    pollin = module_pollination()
    soc = module_soc_loss()
    econ = update_economic_upper_bound(coastal, soc)

    # Summary CSV
    rows = [
        ("Coastal protection", "Total mangrove area (ha)",
         f"{coastal['mangrove_area_ha']:,.0f}"),
        ("Coastal protection", "Protective mangrove (CPI > 0.5) (ha)",
         f"{coastal['protective_area_ha']:,.0f}"),
        ("Coastal protection", "Mangrove lost (ha)",
         f"{coastal['mangrove_lost_ha']:,.0f}"),
        ("Coastal protection", "Annual protection value (USD)",
         f"{coastal['annual_protection_value_usd']:,.0f}"),
        ("Coastal protection", "15-yr NPV of value lost (USD)",
         f"{coastal['npv_value_lost_usd']:,.0f}"),
        ("Freshwater yield", "Total water yield (GL/yr)",
         f"{water['total_yield_GL']:,.1f}"),
        ("Freshwater yield", "Forest→Cropland incremental yield (GL/yr)",
         f"{water['incremental_yield_if_converted_GL']:,.1f}"),
        ("Pollination", "Mean ANI suitability S",
         f"{pollin['mean_S']:.3f}"),
        ("Pollination", "Forest→Cropland ΔS per ha",
         f"-{pollin['delta_per_ha']:.3f}"),
        ("Pollination", "Total suitability lost (unit-ha)",
         f"{pollin['pollination_units_lost']:,.0f}"),
        ("SOC", "Mean baseline SOC (Mg C/ha)",
         f"{soc['mean_baseline_Mg_per_ha']:.1f}"),
        ("SOC", "Total SOC lost (Mg C)",
         f"{soc['soc_lost_Mg_C']:,.0f}"),
        ("SOC", "Equivalent CO2e (Gg)",
         f"{soc['co2e_lost_Gg']:,.1f}"),
        ("SOC", "Economic damage @ SCC (USD)",
         f"{soc['soc_damage_usd']:,.0f}"),
        ("Economic", "Main report damage (lower bound, USD)",
         f"{econ['main_usd']:,.0f}"),
        ("Economic", "Upper-bound damage (with all services, USD)",
         f"{econ['upper_usd']:,.0f}"),
        ("Economic", "Uplift over lower bound (%)",
         f"+{econ['pct_increase']:.0f}%"),
    ]
    pd.DataFrame(rows, columns=["module", "metric", "value"]).to_csv(
        RES / "supplementary_services.csv", index=False)
    print(f"\n  Saved -> results/supplementary_services.csv")

    # JSON dump (numbers only)
    json_safe = {
        "coastal_protection": {k: v for k, v in coastal.items()
                                 if not isinstance(v, np.ndarray)},
        "freshwater_yield": {k: v for k, v in water.items()
                              if not isinstance(v, (np.ndarray, pd.DataFrame))},
        "pollination": {k: v for k, v in pollin.items()
                          if not isinstance(v, (np.ndarray, pd.DataFrame))},
        "soc_loss": {k: v for k, v in soc.items()
                      if not isinstance(v, np.ndarray)},
        "economic": econ,
    }
    (RES / "supplementary_services_summary.json").write_text(
        json.dumps(json_safe, indent=2))

    render_summary_figure(coastal, water, pollin, soc)


if __name__ == "__main__":
    main()
