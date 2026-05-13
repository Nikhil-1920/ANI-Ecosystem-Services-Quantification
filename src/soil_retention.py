"""
ANI Ecosystem Services — Week 8: Soil Retention Modeling (RUSLE)
=================================================================
Implements the Revised Universal Soil Loss Equation:
    A = R × K × LS × C × P   (tonnes / ha / yr)

Factor derivations:
  R  (Rainfall Erosivity)     — CHIRPS: R = 0.04830 × P^1.610
                                [Renard & Freimund 1994]
  K  (Soil Erodibility)       — SoilGrids clay+silt fractions
                                [Williams 1995 EPIC simplified]
  LS (Slope Length+Steepness) — SRTM 30 m DEM central-difference gradient
  C  (Cover Management)       — ESA WorldCover class lookup table
  P  (Support Practice)       — 1.0 (natural tropical, no practices)

Inputs  : data/processed/  (clipped GeoTIFFs)
Outputs : results/rusle_soil_loss.tif
          results/rusle_soil_loss_delta.tif
          results/rusle_erosion_by_landcover.csv
          figures/soil/rusle_factor_components_map.png
          figures/soil/rusle_soil_loss_map.png
          figures/soil/rusle_erosion_by_landcover.png
          figures/soil/rusle_soil_loss_delta_map.png

Run with: venv/bin/python src/soil_retention.py
"""

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import rasterio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import pandas as pd
from pathlib import Path

# ── Directory Paths ────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).parent
PROC_DIR   = SCRIPT_DIR.parent / 'data' / 'processed'
RES_DIR    = SCRIPT_DIR.parent / 'results'
FIG_DIR    = SCRIPT_DIR.parent / 'figures'
RES_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ── Print Separator ────────────────────────────────────────────────────
SEP = '=' * 60

# ── Split-map row boundaries (Ten Degree Channel) ─────────────────────
ANDAMAN_ROW_END   = 11750
NICOBAR_ROW_START = 16200
_PAD = 200


def _split_extent(arr):
    h = arr.shape[0]
    return arr[:min(ANDAMAN_ROW_END + _PAD, h), :], \
           arr[max(0, NICOBAR_ROW_START - _PAD):, :]


def _tight_crop(arr):
    if np.ma.is_masked(arr):
        valid = ~arr.mask
    else:
        valid = ~np.isnan(arr)
    if not valid.any():
        return arr
    rows = np.any(valid, axis=1); cols = np.any(valid, axis=0)
    r0, r1 = np.where(rows)[0][[0, -1]]
    c0, c1 = np.where(cols)[0][[0, -1]]
    pad = 60
    r0 = max(0, r0-pad); r1 = min(arr.shape[0], r1+pad)
    c0 = max(0, c0-pad); c1 = min(arr.shape[1], c1+pad)
    return arr[r0:r1, c0:c1]


# ── Shared design tokens ──────────────────────────────────────────────
OCEAN_COLOR = '#b3d9f2'   # light sky-blue ocean fill
LAND_COLOR  = '#DDECC5'   # soft sage-green land base
BG_COLOR    = 'white'
TEXT_COLOR  = '#1a1a2e'   # dark navy for all text


def _draw_land_ocean_base(ax, land_mask_crop, ocean_color=OCEAN_COLOR, land_color=LAND_COLOR):
    """Fill ocean (NaN / masked) with ocean_color and land with land_color as base layers."""
    h, w = land_mask_crop.shape
    # Ocean fill
    ocean_rgba = np.zeros((h, w, 4), dtype=np.float32)
    import matplotlib.colors as _mc
    oc = _mc.to_rgba(ocean_color)
    lc = _mc.to_rgba(land_color)
    ocean_rgba[~land_mask_crop] = oc
    ocean_rgba[land_mask_crop]  = lc
    ax.imshow(ocean_rgba, interpolation='nearest')


def _tight_crop_bbox(arr):
    """Return the (r0, r1, c0, c1) bounding box for valid (True/nonzero) pixels in arr."""
    if arr.dtype == bool:
        valid = arr
    elif np.ma.is_masked(arr):
        valid = ~arr.mask
    else:
        valid = ~np.isnan(arr.astype(float)) & (arr != 0)
    if not valid.any():
        return 0, arr.shape[0], 0, arr.shape[1]
    rows = np.any(valid, axis=1); cols = np.any(valid, axis=0)
    r0, r1 = np.where(rows)[0][[0, -1]]
    c0, c1 = np.where(cols)[0][[0, -1]]
    pad = 60
    r0 = max(0, r0 - pad); r1 = min(arr.shape[0], r1 + pad)
    c0 = max(0, c0 - pad); c1 = min(arr.shape[1], c1 + pad)
    return r0, r1, c0, c1


def split_imshow(fig, gs_l, gs_r, data, cmap, vmin, vmax,
                 land_mask=None,
                 label_left='Andaman Islands', label_right='Nicobar Islands',
                 interpolation='nearest', norm=None):
    """Render data as Andaman (left) + Nicobar (right) sub-axes with land/ocean base.

    Pass `norm` (e.g. LogNorm) to override the default linear vmin/vmax scaling.
    """
    and_r, nic_r = _split_extent(data)
    ax_a = fig.add_subplot(gs_l)
    ax_n = fig.add_subplot(gs_r)

    if land_mask is not None:
        lm_a, lm_n = _split_extent(land_mask)
    else:
        lm_a = lm_n = None

    for ax, raw_crop, lm_raw, title in [
        (ax_a, and_r, lm_a, label_left),
        (ax_n, nic_r, lm_n, label_right),
    ]:
        ax.set_facecolor(BG_COLOR)
        ax.set_aspect('equal')

        if lm_raw is not None:
            # Derive crop bbox from land mask so both arrays get the same shape
            r0, r1, c0, c1 = _tight_crop_bbox(lm_raw.astype(float))
            lm_crop   = lm_raw[r0:r1, c0:c1]
            data_crop = raw_crop[r0:r1, c0:c1]
            _draw_land_ocean_base(ax, lm_crop)
            masked = np.ma.masked_where(~lm_crop, data_crop)
        else:
            masked = _tight_crop(raw_crop)

        if norm is not None:
            ax.imshow(masked, cmap=cmap, norm=norm,
                      interpolation=interpolation, alpha=0.92)
        else:
            ax.imshow(masked, cmap=cmap, vmin=vmin, vmax=vmax,
                      interpolation=interpolation, alpha=0.92)
        ax.set_title(title, color=TEXT_COLOR, fontsize=11, fontweight='bold', pad=8)
        ax.axis('off')
    return ax_a, ax_n

# ── Model Constants ────────────────────────────────────────────────────
PIXEL_SIZE_M    = 30.0    # DEM / ESA resolution in metres
PIXEL_AREA_HA   = 0.09    # 30 m × 30 m = 0.09 ha
P_FACTOR        = 1.0     # Support-practice factor (natural land = 1.0)
LS_CAP          = 50.0    # Maximum LS factor (avoids cliff-edge artefacts)
RUSLE_A_CAP     = 500.0   # Maximum soil loss A (t/ha/yr) — physical ceiling
K_MIN           = 0.005   # Minimum plausible K (t·ha·h)/(ha·MJ·mm)
K_MAX           = 0.060   # Maximum plausible K (t·ha·h)/(ha·MJ·mm)

# ── RUSLE C-Factor Lookup: ESA WorldCover class → C value ─────────────
# Sources: Wischmeier & Smith (1978); Panagos et al. (2015); Borrelli (2020)
C_FACTOR_BY_CLASS = {
    10:  0.001,   # Tree cover       — closed canopy, minimal erosion
    20:  0.035,   # Shrubland
    30:  0.060,   # Grassland
    40:  0.280,   # Cropland
    50:  0.000,   # Built-up         — impervious, RUSLE not applicable
    60:  0.450,   # Bare / sparse    — high erosion risk
    70:  0.000,   # Snow / Ice
    80:  0.000,   # Open water
    90:  0.020,   # Herbaceous wetland
    95:  0.010,   # Mangroves
    100: 0.025,   # Moss / Lichen
}


# ══════════════════════════════════════════════════════════════════════
# UTILITY
# ══════════════════════════════════════════════════════════════════════
def load_clipped_raster(filename: str, band: int = 1):
    """Load a single band from a clipped GeoTIFF in PROC_DIR.

    Returns (data_array, raster_profile).
    Returns (None, None) when the file is absent.
    """
    file_path = PROC_DIR / filename
    if not file_path.exists():
        print(f"  ⚠️  Not found: {filename}")
        return None, None
    with rasterio.open(file_path) as src:
        data_array   = src.read(band).astype(float)
        rast_profile = src.profile.copy()
        nodata_val   = src.nodata
    if nodata_val is not None:
        data_array = np.where(data_array == nodata_val, np.nan, data_array)
    data_array = np.where(np.abs(data_array) > 1e10, np.nan, data_array)
    return data_array, rast_profile


# ══════════════════════════════════════════════════════════════════════
# FACTOR R — Rainfall Erosivity
# ══════════════════════════════════════════════════════════════════════
def compute_rainfall_erosivity():
    """Derive RUSLE R-factor from CHIRPS precipitation via Renard & Freimund (1994)."""
    print(f"\n{SEP}")
    print("FACTOR R — Rainfall Erosivity  (CHIRPS → Renard & Freimund 1994)")
    print(SEP)

    # Prefer annual total; fall back to daily mean × 365.25
    precip_annual, r_profile = load_clipped_raster('ANI_CHIRPS_Annual_Total_Precip_clipped.tif')
    if precip_annual is None:
        precip_mean, r_profile = load_clipped_raster('ANI_CHIRPS_Mean_Precip_2000_2023_clipped.tif')
        if precip_mean is None:
            print("  ❌  No CHIRPS data found. Using default R = 1200 MJ·mm/(ha·h·yr).")
            return np.full((1, 1), 1200.0), None
        precip_annual = precip_mean * 365.25
        print("  Using CHIRPS daily mean × 365.25 to approximate annual total.")

    precip_annual = np.where(precip_annual <= 0, np.nan, precip_annual)

    # R = 0.04830 × P^1.610  [MJ·mm/(ha·h·yr)]
    r_factor = 0.04830 * (precip_annual ** 1.610)

    print(f"  CHIRPS P_annual  : mean={np.nanmean(precip_annual):.1f} mm/yr  "
          f"range=[{np.nanmin(precip_annual):.0f}, {np.nanmax(precip_annual):.0f}]")
    print(f"  R-Factor         : mean={np.nanmean(r_factor):.1f}  "
          f"range=[{np.nanmin(r_factor):.1f}, {np.nanmax(r_factor):.1f}] MJ·mm/(ha·h·yr)")
    return r_factor, r_profile


# ══════════════════════════════════════════════════════════════════════
# FACTOR K — Soil Erodibility
# ══════════════════════════════════════════════════════════════════════
def compute_soil_erodibility():
    """Derive RUSLE K-factor from SoilGrids clay+silt fractions (Williams 1995)."""
    print(f"\n{SEP}")
    print("FACTOR K — Soil Erodibility  (SoilGrids clay+silt)")
    print(SEP)

    soil_b1, _ = load_clipped_raster('ANI_SoilGrids_Clay_Silt_GEE_clipped.tif', band=1)
    if soil_b1 is None:
        print("  ⚠️  SoilGrids not found — using default K = 0.025.")
        return np.full((1, 1), 0.025), None

    # Read both bands at once to detect clay + sand
    soil_path = PROC_DIR / 'ANI_SoilGrids_Clay_Silt_GEE_clipped.tif'
    with rasterio.open(soil_path) as src:
        n_bands      = src.count
        clay_raw     = src.read(1).astype(float)
        nodata_val   = src.nodata
        k_profile    = src.profile.copy()
        sand_raw     = src.read(2).astype(float) if n_bands >= 2 else clay_raw.copy()

    if nodata_val is not None:
        clay_raw = np.where(clay_raw == nodata_val, np.nan, clay_raw)
        sand_raw = np.where(sand_raw == nodata_val, np.nan, sand_raw)

    clay_raw = np.where(clay_raw < 0, np.nan, clay_raw)
    sand_raw = np.where(sand_raw < 0, np.nan, sand_raw)

    # Auto-detect unit scale: g/kg (>100), % (>1), fraction (<1)
    clay_max = np.nanmax(clay_raw) if not np.all(np.isnan(clay_raw)) else 0.0
    if   clay_max > 100: unit_scale = 1000
    elif clay_max > 1.1: unit_scale = 100
    else:                unit_scale = 1.0

    clay_frac = np.clip(clay_raw / unit_scale, 0, 1)
    sand_frac = np.clip(sand_raw / unit_scale, 0, 1)
    silt_frac = np.clip(1.0 - clay_frac - sand_frac, 0, 1)

    # Williams (1995) EPIC simplified K equation
    sand_safe = np.where(sand_frac > 0.01, sand_frac, 0.01)
    clay_safe = np.where(clay_frac + 0.3 > 0, clay_frac + 0.3, 1e-6)
    denom_cs  = np.where((clay_frac + silt_frac) > 0, clay_frac + silt_frac, 1e-6)

    k_factor = (
        (0.2 + 0.3 * np.exp(-0.0256 * sand_frac * 100 * (1 - silt_frac / sand_safe)))
        * ((silt_frac / denom_cs) ** 0.3)
        * (1 - 0.25 * clay_frac / clay_safe)
    )
    k_factor = np.clip(k_factor, K_MIN, K_MAX)
    k_factor = np.where(np.isnan(clay_frac), np.nan, k_factor)

    print(f"  Clay fraction    : mean={np.nanmean(clay_frac):.3f}")
    print(f"  Sand fraction    : mean={np.nanmean(sand_frac):.3f}")
    print(f"  K-Factor         : mean={np.nanmean(k_factor):.4f}  "
          f"range=[{np.nanmin(k_factor):.4f}, {np.nanmax(k_factor):.4f}]")
    return k_factor, k_profile


# ══════════════════════════════════════════════════════════════════════
# FACTOR LS — Slope Length & Steepness
# ══════════════════════════════════════════════════════════════════════
def compute_slope_steepness(dem_profile):
    """Compute RUSLE LS-factor from SRTM 30 m DEM via McCool et al. (1987)."""
    print(f"\n{SEP}")
    print("FACTOR LS — Slope Length & Steepness  (SRTM 30m DEM)")
    print(SEP)

    dem_grid, ls_profile = load_clipped_raster('ANI_SRTM_DEM_30m_clipped.tif')
    if dem_grid is None:
        print("  ⚠️  DEM not found — using default LS = 1.0.")
        return np.ones((1, 1)), None

    # Central-difference gradient (rise/run = tan slope)
    dem_filled   = np.nan_to_num(dem_grid, nan=0)
    dz_dx        = np.gradient(dem_filled, PIXEL_SIZE_M, axis=1)
    dz_dy        = np.gradient(dem_filled, PIXEL_SIZE_M, axis=0)
    slope_rad    = np.arctan(np.sqrt(dz_dx**2 + dz_dy**2))
    slope_deg    = np.degrees(slope_rad)
    slope_pct    = np.tan(slope_rad) * 100

    # S-factor: McCool et al. (1987) / Renard et al. (1997)
    s_factor = np.where(slope_deg < 9.0,
                        10.8 * np.sin(slope_rad) + 0.03,
                        16.8 * np.sin(slope_rad) - 0.50)

    # L-factor: unit plot length proxy per slope class
    m_exp    = np.where(slope_pct >= 9, 0.5,
               np.where(slope_pct >= 3, 0.4,
               np.where(slope_pct >= 1, 0.3, 0.2)))
    l_factor = (PIXEL_SIZE_M / 22.13) ** m_exp

    ls_factor = np.clip(l_factor * s_factor, 0, LS_CAP)
    ls_factor = np.where(np.isnan(dem_grid), np.nan, ls_factor)

    print(f"  Slope range      : {np.nanmin(slope_deg):.1f}° – {np.nanmax(slope_deg):.1f}°")
    print(f"  Mean slope       : {np.nanmean(slope_deg):.2f}°")
    print(f"  LS-Factor        : mean={np.nanmean(ls_factor):.3f}  "
          f"range=[{np.nanmin(ls_factor):.3f}, {np.nanmax(ls_factor):.3f}]")
    return ls_factor, ls_profile


# ══════════════════════════════════════════════════════════════════════
# FACTOR C — Cover Management
# ══════════════════════════════════════════════════════════════════════
def compute_cover_management(landcover_grid):
    """Assign RUSLE C-factor from ESA WorldCover class lookup table."""
    print(f"\n{SEP}")
    print("FACTOR C — Cover Management  (ESA WorldCover class lookup)")
    print(SEP)

    c_factor     = np.zeros_like(landcover_grid, dtype=float)
    valid_pixels = ~np.isnan(landcover_grid)
    lc_int       = np.nan_to_num(landcover_grid, nan=0).astype(int)

    for class_code, c_val in C_FACTOR_BY_CLASS.items():
        class_mask   = valid_pixels & (lc_int == class_code)
        c_factor[class_mask] = c_val
    c_factor[~valid_pixels] = np.nan

    print("  C-Factor by ESA class:")
    for class_code, c_val in sorted(C_FACTOR_BY_CLASS.items()):
        class_mask = valid_pixels & (lc_int == class_code)
        if class_mask.sum() > 0:
            print(f"    Class {class_code:>3} : C={c_val:.3f}  ({class_mask.sum():>9,} pixels)")

    print(f"\n  Overall mean C   : {np.nanmean(c_factor):.4f}")
    return c_factor


# ══════════════════════════════════════════════════════════════════════
# RUSLE INTEGRATION: A = R × K × LS × C × P
# ══════════════════════════════════════════════════════════════════════
def integrate_rusle_equation(r_factor, k_factor, ls_factor,
                             c_factor, reference_profile):
    """Combine RUSLE factors into annual soil loss A (t/ha/yr)."""
    print(f"\n{SEP}")
    print("INTEGRATING RUSLE: A = R × K × LS × C × P")
    print(SEP)

    def _align_to_target(arr, target_shape):
        """Resize an array to target_shape by crop or constant-fill."""
        if arr is None or arr.shape == (1, 1):
            fill_val = np.nanmean(arr) if arr is not None else 1.0
            return np.full(target_shape, fill_val)
        h_tar, w_tar = target_shape
        h_arr, w_arr = arr.shape
        if (h_arr, w_arr) == (h_tar, w_tar):
            return arr
        h_crop = min(h_arr, h_tar)
        w_crop = min(w_arr, w_tar)
        result = np.full(target_shape, np.nan)
        result[:h_crop, :w_crop] = arr[:h_crop, :w_crop]
        return result

    target_shape  = (reference_profile['height'], reference_profile['width'])
    r_aligned     = _align_to_target(r_factor,  target_shape)
    k_aligned     = _align_to_target(k_factor,  target_shape)
    ls_aligned    = _align_to_target(ls_factor, target_shape)
    c_aligned     = _align_to_target(c_factor,  target_shape)
    p_grid        = np.full(target_shape, P_FACTOR)

    soil_loss_a   = r_aligned * k_aligned * ls_aligned * c_aligned * p_grid
    soil_loss_a   = np.clip(soil_loss_a, 0, RUSLE_A_CAP)
    soil_loss_a   = np.where(np.isnan(c_aligned), np.nan, soil_loss_a)

    valid_land    = ~np.isnan(soil_loss_a)
    total_loss_mt = np.nansum(soil_loss_a) * PIXEL_AREA_HA / 1e6   # megalitres

    print(f"  RUSLE A summary:")
    print(f"    Mean A           : {np.nanmean(soil_loss_a):.2f} t/ha/yr")
    print(f"    Median A         : {np.nanmedian(soil_loss_a):.2f} t/ha/yr")
    print(f"    Max A            : {np.nanmax(soil_loss_a):.2f} t/ha/yr")
    print(f"    Valid pixels     : {valid_land.sum():,}")
    print(f"    Total soil loss  : {total_loss_mt:.4f} Mt/yr")
    return soil_loss_a, r_aligned, k_aligned, ls_aligned


# ══════════════════════════════════════════════════════════════════════
# OUTPUTS — GeoTIFFs, CSVs, Figures
# ══════════════════════════════════════════════════════════════════════
def save_rusle_outputs(soil_loss_a, r_aligned, k_aligned, ls_aligned,
                       c_factor, landcover_grid, rast_profile,
                       soil_loss_baseline=None):
    """Save all RUSLE GeoTIFFs, statistics CSV, and publication figures."""
    print(f"\n{SEP}")
    print("SAVING OUTPUTS & FIGURES")
    print(SEP)

    # ── GeoTIFF: current soil loss A ──────────────────────────────────
    export_profile = rast_profile.copy()
    export_profile.update({'dtype': 'float32', 'nodata': -9999,
                           'count': 1, 'compress': 'lzw'})
    a_export  = np.where(np.isnan(soil_loss_a), -9999, soil_loss_a).astype('float32')
    tif_path  = RES_DIR / 'rusle_soil_loss.tif'
    with rasterio.open(tif_path, 'w', **export_profile) as dst:
        dst.write(a_export, 1)
    print(f"  ✅  GeoTIFF saved → {tif_path.name}")

    # ── GeoTIFF: delta soil loss (2000 baseline vs 2024) ──────────────
    delta_soil_loss = None
    if soil_loss_baseline is not None:
        delta_soil_loss = soil_loss_a - soil_loss_baseline
        delta_export    = np.where(np.isnan(delta_soil_loss), -9999, delta_soil_loss).astype('float32')
        delta_path      = RES_DIR / 'rusle_soil_loss_delta.tif'
        with rasterio.open(delta_path, 'w', **export_profile) as dst:
            dst.write(delta_export, 1)
        print(f"  ✅  GeoTIFF saved → {delta_path.name}")

    # ── CSV: erosion statistics by ESA land cover class ───────────────
    class_labels = {
        10: 'Tree cover (Forest)', 20: 'Shrubland',    30: 'Grassland',
        40: 'Cropland',           50: 'Built-up',      60: 'Bare / Sparse',
        90: 'Wetland',            95: 'Mangroves',
    }
    lc_int      = np.nan_to_num(landcover_grid, nan=0).astype(int)
    class_stats = []

    for class_code, class_label in class_labels.items():
        class_mask = (lc_int == class_code) & ~np.isnan(soil_loss_a)
        if class_mask.sum() == 0:
            continue
        a_vals      = soil_loss_a[class_mask]
        class_area  = class_mask.sum() * PIXEL_AREA_HA
        total_loss  = np.sum(a_vals) * PIXEL_AREA_HA
        record = {
            'esa_class':            class_code,
            'land_cover_label':     class_label,
            'pixel_count':          int(class_mask.sum()),
            'area_ha':              round(class_area,            1),
            'mean_soil_loss_t_ha':  round(float(np.mean(a_vals)), 3),
            'median_soil_loss_t_ha':round(float(np.median(a_vals)), 3),
            'max_soil_loss_t_ha':   round(float(np.max(a_vals)),   3),
            'total_soil_loss_t_yr': round(float(total_loss),       1),
            'c_factor':             C_FACTOR_BY_CLASS.get(class_code, 'N/A'),
        }
        if delta_soil_loss is not None:
            delta_vals = delta_soil_loss[class_mask]
            record['mean_delta_soil_loss']   = round(float(np.mean(delta_vals)), 3)
            record['total_new_erosion_t_yr'] = round(float(np.sum(delta_vals) * PIXEL_AREA_HA), 1)
        class_stats.append(record)
        print(f"  {class_label:<28} mean={np.mean(a_vals):>7.2f}  total={total_loss:>12,.0f} t/yr")

    csv_path = RES_DIR / 'rusle_erosion_by_landcover.csv'
    pd.DataFrame(class_stats).to_csv(csv_path, index=False)
    print(f"\n  ✅  CSV saved → {csv_path.name}")

    # ── Build land mask (valid, non-ocean pixels) ─────────────────────
    lc_int_lm = np.nan_to_num(landcover_grid, nan=0).astype(int)
    land_mask  = (~np.isnan(landcover_grid)) & (lc_int_lm != 0)

    # ── Figure 1: RUSLE four-factor component maps ─────────────────────
    # Panel-by-panel scaling choices:
    #   R: linear, Blues, p2-p98 (rainfall erosivity varies smoothly).
    #   K: cividis (high-contrast, colorblind-safe) — the previous
    #      YlOrBr ramp flattened over the [0.005, 0.060] clip range and
    #      produced an almost-uniform brown wash. Use p2-p98 of actual
    #      land pixels so the gradient maps to real variation.
    #   LS: linear with empirical p2-p98 cap (most slopes on ANI are
    #       gentle; a fixed 0-15 range left subtle gradients invisible).
    #   C: log scale — C is a per-LC-class lookup with values spanning
    #      three orders of magnitude (Tree 0.001 → Bare 0.45). A linear
    #      0-0.5 scale collapses the 98%-forest landscape to a single
    #      green tone. Log spreads the class values so cropland, bare,
    #      and forest read as distinct shades.
    from matplotlib.colors import LogNorm
    factor_panels = [
        (r_aligned,  'R-Factor\n(Rainfall Erosivity)', 'Blues',     None,  None, 'MJ·mm/(ha·h·yr)', 'linear'),
        (k_aligned,  'K-Factor\n(Soil Erodibility)',   'cividis',   None,  None, 't·ha·h/(ha·MJ·mm)', 'linear'),
        (ls_aligned, 'LS-Factor\n(Slope × Steepness)', 'copper_r',  None,  None, 'Dimensionless',     'linear'),
        (c_factor,   'C-Factor\n(Cover Management)',   'magma_r',   0.001, 0.5,  'Dimensionless (log)', 'log'),
    ]
    (FIG_DIR / 'soil').mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(14, 24), facecolor=BG_COLOR)
    gs_f = gridspec.GridSpec(4, 3, figure=fig,
                             width_ratios=[2, 2, 0.08],
                             hspace=0.28, wspace=0.04)
    fig.suptitle('RUSLE Factor Component Maps — Andaman & Nicobar Islands',
                 color=TEXT_COLOR, fontsize=15, fontweight='bold', y=0.995)
    for row_i, (factor_arr, panel_title, cmap_name, vmin, vmax, unit_label, scale) in enumerate(factor_panels):
        land_vals = factor_arr[land_mask] if factor_arr.shape == land_mask.shape else factor_arr.ravel()
        if vmin is None: vmin = float(np.nanpercentile(land_vals, 2))
        if vmax is None: vmax = float(np.nanpercentile(land_vals, 98))
        if scale == 'log':
            norm = LogNorm(vmin=max(vmin, 1e-4), vmax=vmax)
        else:
            norm = plt.Normalize(vmin=vmin, vmax=vmax)
        ax_a, ax_n = split_imshow(fig, gs_f[row_i, 0], gs_f[row_i, 1],
                                   factor_arr, cmap_name, vmin, vmax,
                                   land_mask=land_mask, norm=norm)
        # Push titles further from the panel so they don't crowd the
        # island shapes underneath (previous pad=7 left titles touching
        # the image edge in tall thin panels).
        ax_a.set_title(f'Andaman — {panel_title}',
                       color=TEXT_COLOR, fontsize=11, fontweight='bold',
                       pad=16, y=1.01)
        ax_n.set_title(f'Nicobar — {panel_title}',
                       color=TEXT_COLOR, fontsize=11, fontweight='bold',
                       pad=16, y=1.01)
        cax = fig.add_subplot(gs_f[row_i, 2])
        sm  = plt.cm.ScalarMappable(cmap=cmap_name, norm=norm)
        sm.set_array([])
        cb = fig.colorbar(sm, cax=cax)
        cb.set_label(unit_label, color=TEXT_COLOR, fontsize=8)
        cb.ax.yaxis.set_tick_params(color=TEXT_COLOR, labelsize=7)
        plt.setp(cb.ax.yaxis.get_ticklabels(), color=TEXT_COLOR)
    fig.savefig(FIG_DIR / 'soil' / 'rusle_factor_components_map.png',
                dpi=160, bbox_inches='tight', facecolor=BG_COLOR)
    plt.close()
    print("  ✅  Figure saved → figures/soil/rusle_factor_components_map.png")

    # ── Figure 2: Annual Soil Erosion Risk Map ─────────────────────────
    a_land_vals = soil_loss_a[land_mask]
    a_log       = np.where(land_mask, np.log1p(soil_loss_a), np.nan)
    a_vmax      = float(np.nanpercentile(np.log1p(a_land_vals), 98))
    tick_raw    = np.linspace(0, float(np.nanpercentile(a_land_vals, 98)), 6)
    tick_locs   = np.log1p(tick_raw)
    cmap_erosion = plt.cm.YlOrRd.copy()
    cmap_erosion.set_bad(color=OCEAN_COLOR)

    a_max_real  = float(np.nanmax(a_land_vals))
    cap_hit     = a_max_real >= RUSLE_A_CAP - 1e-6
    max_label   = f'≥{RUSLE_A_CAP:.0f} (capped)' if cap_hit else f'{a_max_real:.0f}'

    fig = plt.figure(figsize=(14, 13), facecolor=BG_COLOR)
    gs_s = gridspec.GridSpec(1, 3, figure=fig,
                             width_ratios=[2, 2, 0.08], wspace=0.04)
    ax_sa, ax_sn = split_imshow(fig, gs_s[0, 0], gs_s[0, 1],
                                 a_log, cmap_erosion, 0, a_vmax,
                                 land_mask=land_mask)
    ax_sa.set_title('Andaman Islands\nRUSLE Annual Soil Erosion Risk',
                    color=TEXT_COLOR, fontsize=12, fontweight='bold', pad=8)
    ax_sn.set_title('Nicobar Islands\nRUSLE Annual Soil Erosion Risk',
                    color=TEXT_COLOR, fontsize=12, fontweight='bold', pad=8)
    ax_sa.text(0.03, 0.03,
               f'Mean: {np.nanmean(a_land_vals):.1f} t/ha/yr\n'
               f'Max:  {max_label} t/ha/yr',
               transform=ax_sa.transAxes, color=TEXT_COLOR, fontsize=9,
               bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.9, edgecolor='#cccccc'))
    cax_s = fig.add_subplot(gs_s[0, 2])
    sm_s  = plt.cm.ScalarMappable(cmap=cmap_erosion,
                                   norm=plt.Normalize(vmin=0, vmax=a_vmax))
    sm_s.set_array([])
    cbar_s = fig.colorbar(sm_s, cax=cax_s, extend='max')
    cbar_s.set_ticks(tick_locs)
    cbar_s.set_ticklabels([f'{v:.0f}' for v in tick_raw])
    cbar_s.set_label('Soil Loss  A  (t / ha / yr)  — log scale (tick values are raw)',
                     color=TEXT_COLOR, fontsize=9)
    cbar_s.ax.yaxis.set_tick_params(color=TEXT_COLOR, labelsize=8)
    plt.setp(cbar_s.ax.yaxis.get_ticklabels(), color=TEXT_COLOR)
    fig.suptitle('RUSLE Annual Soil Erosion Risk  (A = R · K · LS · C · P)\n'
                 'Andaman & Nicobar Islands — 2024',
                 color=TEXT_COLOR, fontsize=14, fontweight='bold', y=0.995)
    fig.savefig(FIG_DIR / 'soil' / 'rusle_soil_loss_map.png',
                dpi=200, bbox_inches='tight', facecolor=BG_COLOR)
    plt.close()
    print("  ✅  Figure saved → figures/soil/rusle_soil_loss_map.png")

    # ── Figure 3: Lollipop chart — mean erosion by land cover ─────────
    if class_stats:
        stats_df  = pd.DataFrame(class_stats).sort_values('mean_soil_loss_t_ha', ascending=True)
        n_rows    = len(stats_df)
        y_pos     = list(range(n_rows))
        norm_vals = stats_df['mean_soil_loss_t_ha'].values
        norm_c    = plt.Normalize(vmin=0, vmax=max(float(norm_vals.max()), 1))

        fig, ax = plt.subplots(figsize=(12, max(5, n_rows * 0.85 + 1.5)),
                               facecolor=BG_COLOR)
        ax.set_facecolor('#f5f8fc')
        ax.grid(axis='x', color='#dde3ec', linewidth=0.8, zorder=0)
        for i, (_, row) in enumerate(stats_df.iterrows()):
            ax.plot([0, row['mean_soil_loss_t_ha']], [i, i],
                    color='#9fb3c8', linewidth=1.8, zorder=1)
        sc = ax.scatter(norm_vals, y_pos, c=norm_vals, cmap='YlOrRd',
                        norm=norm_c, s=260, zorder=3, edgecolors='white', linewidths=1.2)
        for i, val in enumerate(norm_vals):
            ax.text(val + max(float(norm_vals.max()) * 0.015, 0.5), i,
                    f'{val:.2f}', va='center', ha='left',
                    color=TEXT_COLOR, fontsize=9, fontweight='bold')
        ax.set_yticks(y_pos)
        ax.set_yticklabels(stats_df['land_cover_label'], color=TEXT_COLOR, fontsize=10)
        ax.set_xlabel('Mean Soil Loss  A  (t / ha / yr)', color=TEXT_COLOR, fontsize=11)
        ax.set_xlim(left=-float(norm_vals.max()) * 0.04)
        ax.tick_params(axis='x', colors=TEXT_COLOR, labelsize=9)
        ax.tick_params(axis='y', length=0)
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)
        for spine in ['bottom', 'left']:
            ax.spines[spine].set_color('#cccccc')
        ax.axvline(50,  color='#f57c00', linestyle='--', alpha=0.7, linewidth=1.4, label='High Risk (50 t/ha/yr)')
        ax.axvline(200, color='#c62828', linestyle='--', alpha=0.7, linewidth=1.4, label='Severe Risk (200 t/ha/yr)')
        ax.legend(facecolor='white', edgecolor='#cccccc', labelcolor=TEXT_COLOR, fontsize=9, loc='lower right')
        cbar_ax = fig.add_axes([0.92, 0.12, 0.015, 0.76])
        cb3 = plt.colorbar(sc, cax=cbar_ax)
        cb3.set_label('t / ha / yr', color=TEXT_COLOR, fontsize=9)
        cb3.ax.yaxis.set_tick_params(color=TEXT_COLOR, labelsize=8)
        plt.setp(cb3.ax.yaxis.get_ticklabels(), color=TEXT_COLOR)
        ax.set_title('Mean Annual Soil Erosion by Land Cover Class (RUSLE)\nAndaman & Nicobar Islands — 2024',
                     color=TEXT_COLOR, fontsize=13, fontweight='bold', pad=14)
        fig.tight_layout(rect=[0, 0, 0.91, 1])
        fig.savefig(FIG_DIR / 'soil' / 'rusle_erosion_by_landcover.png',
                    dpi=180, bbox_inches='tight', facecolor=BG_COLOR)
        plt.close()
        print("  ✅  Figure saved → figures/soil/rusle_erosion_by_landcover.png")

    # ── Figure 4: Delta soil loss map — split Andaman / Nicobar ───────
    if delta_soil_loss is not None:
        from scipy.ndimage import maximum_filter
        from matplotlib.colors import LinearSegmentedColormap as LSC

        delta_land_vals = delta_soil_loss[land_mask]
        pos_delta       = delta_land_vals[delta_land_vals > 0]
        # Tight colorbar range so small but real deltas show clearly
        # (p99 was washing everything to near-white; p90 is far more readable)
        vmax_delta      = float(np.nanpercentile(pos_delta, 90)) if len(pos_delta) > 0 else 50.0
        DELTA_THRESHOLD = 0.5   # t/ha/yr — below this is noise

        # Change mask, then dilate so sparse single-pixel changes are visible
        # at print/render resolution (3x3 max-filter widens dots ~90m → 270m).
        change_mask    = land_mask & (delta_soil_loss > DELTA_THRESHOLD)
        delta_for_disp = np.where(change_mask,
                                   np.clip(delta_soil_loss, DELTA_THRESHOLD, vmax_delta),
                                   0.0)
        delta_dilated  = maximum_filter(delta_for_disp, size=3)
        delta_display  = np.where((delta_dilated > 0) & land_mask,
                                   delta_dilated, np.nan)

        total_new_t    = float(np.nansum(pos_delta)) * PIXEL_AREA_HA
        new_area_ha    = int(np.sum(delta_land_vals > 1.0)) * PIXEL_AREA_HA

        # Bolder palette: start at saturated orange, not pale peach, so even
        # threshold-level pixels register against the sage land base.
        cmap_delta = LSC.from_list('erosion_delta',
                                   ['#fd8d3c', '#e6550d', '#a63603', '#7f2704'], N=256)
        cmap_delta.set_bad(color='none')

        fig = plt.figure(figsize=(14, 13), facecolor=BG_COLOR)
        gs_d = gridspec.GridSpec(1, 3, figure=fig,
                                 width_ratios=[2, 2, 0.08], wspace=0.04)
        ax_da, ax_dn = split_imshow(fig, gs_d[0, 0], gs_d[0, 1],
                                     delta_display, cmap_delta, DELTA_THRESHOLD, vmax_delta,
                                     land_mask=land_mask)
        ax_da.set_title('Andaman Islands\nΔ Soil Erosion Risk (2000 → 2024)',
                        color=TEXT_COLOR, fontsize=12, fontweight='bold', pad=8)
        ax_dn.set_title('Nicobar Islands\nΔ Soil Erosion Risk (2000 → 2024)',
                        color=TEXT_COLOR, fontsize=12, fontweight='bold', pad=8)
        ax_da.text(0.03, 0.03,
                   f'New erosion area: {new_area_ha:,.0f} ha\n'
                   f'Added loss: {total_new_t:,.0f} t/yr\n'
                   f'(colour capped at p90 = {vmax_delta:.0f} t/ha/yr)',
                   transform=ax_da.transAxes, color=TEXT_COLOR, fontsize=9,
                   bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.9, edgecolor='#cccccc'))
        cax_d = fig.add_subplot(gs_d[0, 2])
        sm_d  = plt.cm.ScalarMappable(cmap=cmap_delta,
                                       norm=plt.Normalize(vmin=0, vmax=vmax_delta))
        sm_d.set_array([])
        cbar_d = fig.colorbar(sm_d, cax=cax_d, extend='max')
        cbar_d.set_label('Increase in Soil Loss (t/ha/yr)\nSince Year 2000',
                         color=TEXT_COLOR, fontsize=9)
        cbar_d.ax.yaxis.set_tick_params(color=TEXT_COLOR, labelsize=8)
        plt.setp(cbar_d.ax.yaxis.get_ticklabels(), color=TEXT_COLOR)
        fig.suptitle('Δ Soil Erosion Risk (2000 vs 2024) — Anthropogenic Impact\n'
                     'Andaman & Nicobar Islands',
                     color=TEXT_COLOR, fontsize=14, fontweight='bold', y=0.995)
        fig.savefig(FIG_DIR / 'soil' / 'rusle_soil_loss_delta_map.png',
                    dpi=200, bbox_inches='tight', facecolor=BG_COLOR)
        plt.close()
        print("  ✅  Figure saved → figures/soil/rusle_soil_loss_delta_map.png")


# ══════════════════════════════════════════════════════════════════════
# ENTRY POINT
# ══════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    print(SEP)
    print("  ANI Ecosystem Services — Week 8: RUSLE Soil Retention")
    print(f"  Input : {PROC_DIR}")
    print(f"  Output: {RES_DIR}  |  {FIG_DIR}")
    print(SEP)

    # Load ESA WorldCover (defines valid land domain)
    landcover_grid, esa_profile = load_clipped_raster('ANI_ESA_WorldCover_mosaic_clipped.tif')
    if landcover_grid is None:
        print("❌  ESA WorldCover not found — aborting.")
        exit(1)

    # Load GFW deforestation year layer (needed for 2000 baseline)
    defor_yr_grid, _ = load_clipped_raster('ANI_GFW_Forest_Loss_2001_2023_clipped.tif')

    # ── Compute all RUSLE factors ──────────────────────────────────────
    r_factor,  _ = compute_rainfall_erosivity()
    k_factor,  _ = compute_soil_erodibility()
    ls_factor, _ = compute_slope_steepness(esa_profile)

    # ── Phase A: Year 2000 baseline (pre-deforestation) ───────────────
    print(f"\n{SEP}")
    print("  PHASE A: RECONSTRUCTING YEAR 2000 BASELINE")
    print(f"{SEP}")

    lc_year2000 = landcover_grid.copy()
    if defor_yr_grid is not None:
        deforested_pixels = (np.nan_to_num(defor_yr_grid, nan=0) > 0)
        lc_year2000[deforested_pixels] = 10   # Restore to tree cover / primary habitat
        print(f"  Re-forested {deforested_pixels.sum():,} pixels to reconstruct Year 2000 baseline.")

    c_year2000       = compute_cover_management(lc_year2000)
    soil_loss_2000, _, _, _ = integrate_rusle_equation(
        r_factor, k_factor, ls_factor, c_year2000, esa_profile
    )

    # ── Phase B: Current state (2024) ─────────────────────────────────
    print(f"\n{SEP}")
    print("  PHASE B: CURRENT STATE (2024)")
    print(f"{SEP}")

    c_factor_2024    = compute_cover_management(landcover_grid)
    soil_loss_2024, r_aligned, k_aligned, ls_aligned = integrate_rusle_equation(
        r_factor, k_factor, ls_factor, c_factor_2024, esa_profile
    )

    save_rusle_outputs(
        soil_loss_2024, r_aligned, k_aligned, ls_aligned,
        c_factor_2024, landcover_grid, esa_profile,
        soil_loss_baseline=soil_loss_2000
    )

    print(f"\n{SEP}")
    print("  ✅  Week 8 RUSLE Soil Retention Complete!")
    print(f"      GeoTIFF : results/rusle_soil_loss.tif")
    print(f"      CSV     : results/rusle_erosion_by_landcover.csv")
    print(f"      Figures : figures/soil/")
    print(f"  Next → run: venv/bin/python src/synthesis_hotspots.py")
    print(SEP)
