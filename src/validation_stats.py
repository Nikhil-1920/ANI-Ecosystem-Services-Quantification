"""
ANI Ecosystem Services — Statistical Validation Module
======================================================
Adds publication-grade statistical rigor to the existing pipeline,
using only datasets already in data/processed/ and results/.

Computes:
  1. Mann-Kendall trend test + Sen's slope on annual carbon loss
  2. Bootstrap 95% CIs on total carbon loss & annual mean
  3. Cross-validation upgrades for GEDI vs Saatchi AGB:
       - Coefficient of determination (R^2)
       - Lin's Concordance Correlation Coefficient (CCC)
       - Reduced Major Axis (RMA) regression
       - Bootstrap CIs on RMSE, bias, Pearson r
  4. Moran's I (queen 3x3) on the ECI hotspot raster
       - Demonstrates hotspots are spatially clustered, not random

Inputs : data/processed/*.tif   results/carbon_annual_loss_by_year.csv
Outputs: results/validation_metrics.csv
         results/validation_summary.json
         figures/validation/carbon_loss_trend_with_ci.png
         figures/validation/agb_cross_validation_upgraded.png
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
from scipy import stats

warnings.filterwarnings("ignore")

SCRIPT_DIR = Path(__file__).parent
PROC_DIR   = SCRIPT_DIR.parent / "data" / "processed"
RES_DIR    = SCRIPT_DIR.parent / "results"
FIG_DIR    = SCRIPT_DIR.parent / "figures" / "validation"
FIG_DIR.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(42)
N_BOOT = 2000

SEP = "=" * 64


def banner(msg: str) -> None:
    print(f"\n{SEP}\n{msg}\n{SEP}")


# ──────────────────────────────────────────────────────────────────────
# 1.  TREND STATISTICS  (Mann-Kendall + Sen's slope)
# ──────────────────────────────────────────────────────────────────────
def mann_kendall(x: np.ndarray):
    """Mann-Kendall non-parametric trend test with tie correction.

    Returns dict with S, Z, p-value, trend direction.
    """
    x = np.asarray(x, dtype=float)
    n = len(x)
    s = 0
    for i in range(n - 1):
        s += np.sum(np.sign(x[i + 1:] - x[i]))

    # Variance with tie correction
    unique, counts = np.unique(x, return_counts=True)
    tie_term = np.sum(counts * (counts - 1) * (2 * counts + 5))
    var_s = (n * (n - 1) * (2 * n + 5) - tie_term) / 18.0

    if s > 0:
        z = (s - 1) / np.sqrt(var_s)
    elif s < 0:
        z = (s + 1) / np.sqrt(var_s)
    else:
        z = 0.0

    p = 2 * (1 - stats.norm.cdf(abs(z)))
    trend = "increasing" if z > 0 else ("decreasing" if z < 0 else "no trend")
    return dict(S=float(s), Z=float(z), p=float(p), trend=trend, n=int(n))


def sens_slope(t: np.ndarray, x: np.ndarray):
    """Theil-Sen median slope. t = times (years), x = values."""
    t = np.asarray(t, dtype=float)
    x = np.asarray(x, dtype=float)
    slopes = []
    n = len(x)
    for i in range(n - 1):
        slopes.extend((x[i + 1:] - x[i]) / (t[i + 1:] - t[i]))
    slopes = np.array(slopes)
    median = float(np.median(slopes))
    # 95% CI via normal approximation (Hollander & Wolfe)
    c_alpha = 1.96 * np.sqrt(n * (n - 1) * (2 * n + 5) / 18.0)
    n_pairs = len(slopes)
    lo_idx = int(max(0, np.floor((n_pairs - c_alpha) / 2)))
    hi_idx = int(min(n_pairs - 1, np.ceil((n_pairs + c_alpha) / 2)))
    sorted_slopes = np.sort(slopes)
    return dict(slope=median,
                ci_lo=float(sorted_slopes[lo_idx]),
                ci_hi=float(sorted_slopes[hi_idx]))


# ──────────────────────────────────────────────────────────────────────
# 2.  CROSS-VALIDATION METRICS  (R^2, CCC, RMA, bootstrap CIs)
# ──────────────────────────────────────────────────────────────────────
def lins_ccc(x: np.ndarray, y: np.ndarray) -> float:
    """Lin's Concordance Correlation Coefficient: agreement with the 1:1 line."""
    mx, my = x.mean(), y.mean()
    sx2, sy2 = x.var(ddof=0), y.var(ddof=0)
    cov = np.mean((x - mx) * (y - my))
    return float((2 * cov) / (sx2 + sy2 + (mx - my) ** 2))


def rma_regression(x: np.ndarray, y: np.ndarray):
    """Reduced Major Axis (RMA) regression — both axes have error."""
    sx, sy = x.std(ddof=1), y.std(ddof=1)
    r, _ = stats.pearsonr(x, y)
    slope = np.sign(r) * sy / sx
    intercept = y.mean() - slope * x.mean()
    return dict(slope=float(slope), intercept=float(intercept), r=float(r))


def bootstrap_cv_metrics(x: np.ndarray, y: np.ndarray, n_boot: int = N_BOOT):
    """Bootstrap 95% CIs on RMSE, bias, Pearson r, R^2, CCC."""
    n = len(x)
    rmse_b, bias_b, r_b, r2_b, ccc_b = [], [], [], [], []
    for _ in range(n_boot):
        idx = RNG.integers(0, n, n)
        xb, yb = x[idx], y[idx]
        diff = xb - yb
        rmse_b.append(np.sqrt(np.mean(diff ** 2)))
        bias_b.append(diff.mean())
        rb, _ = stats.pearsonr(xb, yb)
        r_b.append(rb)
        r2_b.append(rb ** 2)
        ccc_b.append(lins_ccc(xb, yb))

    def ci(arr):
        return float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))

    return dict(rmse_ci=ci(rmse_b), bias_ci=ci(bias_b),
                r_ci=ci(r_b), r2_ci=ci(r2_b), ccc_ci=ci(ccc_b))


# ──────────────────────────────────────────────────────────────────────
# 3.  MORAN'S I  (queen 3x3 contiguity)
# ──────────────────────────────────────────────────────────────────────
def morans_i_queen(grid: np.ndarray, valid_mask: np.ndarray):
    """Global Moran's I on a 2D raster using queen-case 3x3 neighbours.

    Convolution implementation — fast for large rasters.
    Permutation test (n_perm=199) for inference.
    """
    g = np.where(valid_mask, grid, np.nan)
    finite = np.isfinite(g)
    n = int(finite.sum())
    if n < 50:
        return dict(I=np.nan, p=np.nan, n=n, note="too few valid cells")

    mean_val = np.nanmean(g)
    z = np.where(finite, g - mean_val, 0.0)

    # 3x3 sum kernel minus the centre cell == queen neighbours
    from scipy.signal import convolve2d
    kernel = np.ones((3, 3))
    kernel[1, 1] = 0
    neigh_sum = convolve2d(z, kernel, mode="same", boundary="fill", fillvalue=0)
    neigh_cnt = convolve2d(finite.astype(float), kernel, mode="same",
                            boundary="fill", fillvalue=0)
    # Row-standardised weights: each cell's contribution = z_i * mean(z_j)
    with np.errstate(invalid="ignore", divide="ignore"):
        local_lag = np.where(neigh_cnt > 0, neigh_sum / neigh_cnt, 0.0)

    numerator = np.sum(z[finite] * local_lag[finite])
    denominator = np.sum(z[finite] ** 2)
    if denominator == 0:
        return dict(I=np.nan, p=np.nan, n=n, note="zero variance")

    # Row-standardised W means effective W = n (sum of rows of 1)
    I = float(numerator / denominator)

    # Permutation test
    n_perm = 199
    perm_I = np.empty(n_perm)
    flat_vals = g[finite].astype(float)
    flat_idx = np.argwhere(finite)
    for p_i in range(n_perm):
        permuted = RNG.permutation(flat_vals)
        z_p = np.zeros_like(z)
        z_p[flat_idx[:, 0], flat_idx[:, 1]] = permuted - permuted.mean()
        neigh_p = convolve2d(z_p, kernel, mode="same",
                              boundary="fill", fillvalue=0)
        lag_p = np.where(neigh_cnt > 0, neigh_p / neigh_cnt, 0.0)
        denom_p = np.sum(z_p[finite] ** 2)
        perm_I[p_i] = np.sum(z_p[finite] * lag_p[finite]) / denom_p \
            if denom_p > 0 else 0.0

    p_perm = (1 + np.sum(np.abs(perm_I) >= abs(I))) / (n_perm + 1)
    return dict(I=I, p=float(p_perm), n=n,
                perm_mean=float(perm_I.mean()), perm_std=float(perm_I.std()))


def block_mean(arr: np.ndarray, k: int) -> np.ndarray:
    """Mean-pool a 2D array by factor k (drops trailing remainder)."""
    h, w = arr.shape
    ht, wt = h - h % k, w - w % k
    a = arr[:ht, :wt]
    return np.nanmean(a.reshape(ht // k, k, wt // k, k), axis=(1, 3))


# ──────────────────────────────────────────────────────────────────────
# DRIVER
# ──────────────────────────────────────────────────────────────────────
def main():
    summary = {}

    # ── 1. Carbon-loss trend ──────────────────────────────────────────
    banner("1. Annual Carbon-Loss Trend (Mann-Kendall + Sen's slope)")
    df = pd.read_csv(RES_DIR / "carbon_annual_loss_by_year.csv")
    years = df["year"].values
    co2 = df["co2e_ggco2e"].values
    area = df["area_ha"].values

    mk_co2 = mann_kendall(co2)
    sen_co2 = sens_slope(years, co2)
    mk_area = mann_kendall(area)
    sen_area = sens_slope(years, area)

    print(f"  CO2e:   S={mk_co2['S']:.0f}  Z={mk_co2['Z']:+.3f}  "
          f"p={mk_co2['p']:.4f}  → {mk_co2['trend']}")
    print(f"  Sen slope (CO2e):  {sen_co2['slope']:+.4f}  "
          f"GgCO2e/yr  [95% CI {sen_co2['ci_lo']:+.4f}, {sen_co2['ci_hi']:+.4f}]")
    print(f"  Area:   S={mk_area['S']:.0f}  Z={mk_area['Z']:+.3f}  "
          f"p={mk_area['p']:.4f}  → {mk_area['trend']}")
    print(f"  Sen slope (area):  {sen_area['slope']:+.3f}  "
          f"ha/yr  [95% CI {sen_area['ci_lo']:+.3f}, {sen_area['ci_hi']:+.3f}]")

    summary["trend"] = dict(
        co2e_mk=mk_co2, co2e_sen=sen_co2,
        area_mk=mk_area, area_sen=sen_area,
    )

    # ── 2. Bootstrap CIs on annual loss totals ────────────────────────
    banner("2. Bootstrap 95% CI on Total Carbon Loss")
    boots = RNG.choice(co2, size=(N_BOOT, len(co2)), replace=True).sum(axis=1)
    co2_total = float(co2.sum())
    co2_lo, co2_hi = np.percentile(boots, [2.5, 97.5])
    print(f"  Total CO2e (2001-2024) = {co2_total:.2f} GgCO2e")
    print(f"  95% CI (bootstrap)     = [{co2_lo:.2f}, {co2_hi:.2f}] GgCO2e")

    boots_a = RNG.choice(area, size=(N_BOOT, len(area)), replace=True).sum(axis=1)
    area_total = float(area.sum())
    area_lo, area_hi = np.percentile(boots_a, [2.5, 97.5])
    print(f"  Total area lost        = {area_total:.0f} ha")
    print(f"  95% CI (bootstrap)     = [{area_lo:.0f}, {area_hi:.0f}] ha")

    summary["totals"] = dict(
        co2e_ggco2e=co2_total,
        co2e_ci=[float(co2_lo), float(co2_hi)],
        area_ha=area_total,
        area_ci=[float(area_lo), float(area_hi)],
    )

    # ── 3. GEDI vs Saatchi cross-validation upgraded ──────────────────
    banner("3. Cross-Validation: GEDI vs Saatchi AGB")

    def _read(path):
        with rasterio.open(PROC_DIR / path) as src:
            a = src.read(1).astype(float)
            if src.nodata is not None:
                a = np.where(a == src.nodata, np.nan, a)
        return a

    gedi = _read("ANI_GEDI_Biomass_Density_clipped.tif")
    saatchi = _read("ANI_Saatchi_AGB_CrossValidation_clipped.tif")
    landcover = _read("ANI_ESA_WorldCover_mosaic_clipped.tif")
    gedi = np.where((gedi < 0) | (gedi > 800), np.nan, gedi)
    saatchi = np.where((saatchi < 1) | (saatchi > 800), np.nan, saatchi)

    h = min(gedi.shape[0], saatchi.shape[0], landcover.shape[0])
    w = min(gedi.shape[1], saatchi.shape[1], landcover.shape[1])
    # Restrict to forested pixels only (Tree cover = 10, Mangroves = 95)
    # — non-forest pixels are not a fair comparison: Saatchi is forest-only
    # AGB product while GEDI's gap-filled grid still returns values over
    # cropland / built-up. Filtering removes structural bias.
    forest_mask = np.isin(landcover[:h, :w], [10, 95])
    mask = (np.isfinite(gedi[:h, :w])
             & np.isfinite(saatchi[:h, :w])
             & forest_mask)
    gv = gedi[:h, :w][mask]
    sv = saatchi[:h, :w][mask]
    n_pairs = len(gv)
    print(f"  Overlapping pixels: {n_pairs:,}")

    # If too many points, subsample for bootstrap speed (representative)
    if n_pairs > 200_000:
        idx = RNG.choice(n_pairs, 200_000, replace=False)
        gv_s, sv_s = gv[idx], sv[idx]
    else:
        gv_s, sv_s = gv, sv

    rmse = float(np.sqrt(np.mean((gv_s - sv_s) ** 2)))
    bias = float(np.mean(gv_s - sv_s))
    mae = float(np.mean(np.abs(gv_s - sv_s)))
    pr, pp = stats.pearsonr(gv_s, sv_s)
    r2 = float(pr ** 2)
    ccc = lins_ccc(gv_s, sv_s)
    rma = rma_regression(sv_s, gv_s)   # x = Saatchi (ref), y = GEDI
    ols = stats.linregress(sv_s, gv_s)

    boots_cv = bootstrap_cv_metrics(gv_s, sv_s, n_boot=500)

    print(f"  RMSE    = {rmse:.2f} Mg/ha    "
          f"95% CI [{boots_cv['rmse_ci'][0]:.2f}, {boots_cv['rmse_ci'][1]:.2f}]")
    print(f"  Bias    = {bias:+.2f} Mg/ha   "
          f"95% CI [{boots_cv['bias_ci'][0]:+.2f}, {boots_cv['bias_ci'][1]:+.2f}]")
    print(f"  MAE     = {mae:.2f} Mg/ha")
    print(f"  Pearson r = {pr:.4f}   95% CI [{boots_cv['r_ci'][0]:.4f}, {boots_cv['r_ci'][1]:.4f}]")
    print(f"  R^2     = {r2:.4f}     95% CI [{boots_cv['r2_ci'][0]:.4f}, {boots_cv['r2_ci'][1]:.4f}]")
    print(f"  Lin CCC = {ccc:.4f}    95% CI [{boots_cv['ccc_ci'][0]:.4f}, {boots_cv['ccc_ci'][1]:.4f}]")
    print(f"  OLS  : y = {ols.slope:.3f} * x + {ols.intercept:.2f}")
    print(f"  RMA  : y = {rma['slope']:.3f} * x + {rma['intercept']:.2f}")

    summary["cross_validation"] = dict(
        n_pairs=int(n_pairs),
        rmse=rmse, rmse_ci=boots_cv["rmse_ci"],
        bias=bias, bias_ci=boots_cv["bias_ci"],
        mae=mae, pearson_r=float(pr), p_value=float(pp),
        r2=r2, r2_ci=boots_cv["r2_ci"],
        ccc=ccc, ccc_ci=boots_cv["ccc_ci"],
        ols=dict(slope=float(ols.slope), intercept=float(ols.intercept)),
        rma=rma,
    )

    # Save the upgraded CV figure
    save_cv_figure(gv_s, sv_s, summary["cross_validation"])

    # ── 4. Moran's I on ECI hotspot raster ────────────────────────────
    banner("4. Spatial Autocorrelation (Moran's I) on ECI Hotspots")

    with rasterio.open(RES_DIR / "habitat_quality_delta.tif") as s:
        hq = s.read(1).astype(float)
        if s.nodata is not None:
            hq = np.where(hq == s.nodata, np.nan, hq)

    # Downsample to a tractable size: 50x50 mean pooling -> 510 x 150
    # This preserves spatial structure while keeping permutations cheap.
    BLK = 50
    hq_pool = block_mean(hq, BLK)
    valid = np.isfinite(hq_pool) & (hq_pool < 0)   # focus on degradation cells
    print(f"  Pooled grid: {hq_pool.shape}  (block size {BLK} px = ~1.5 km)")
    print(f"  Valid degradation cells: {valid.sum():,}")

    mi = morans_i_queen(hq_pool, valid)
    print(f"  Moran's I = {mi['I']:.4f}   p (199 perms) = {mi['p']:.4f}   n = {mi['n']:,}")
    print(f"  Permutation null:  mean={mi.get('perm_mean', 0):.4f}  "
          f"std={mi.get('perm_std', 0):.4f}")
    summary["morans_I_habitat_delta"] = mi

    # Same test on RUSLE soil-loss delta
    with rasterio.open(RES_DIR / "rusle_soil_loss_delta.tif") as s:
        sl = s.read(1).astype(float)
        if s.nodata is not None:
            sl = np.where(sl == s.nodata, np.nan, sl)
    sl_pool = block_mean(sl, BLK)
    valid_sl = np.isfinite(sl_pool) & (sl_pool > 0)
    mi_sl = morans_i_queen(sl_pool, valid_sl)
    print(f"  RUSLE delta Moran's I = {mi_sl['I']:.4f}   p = {mi_sl['p']:.4f}   n = {mi_sl['n']:,}")
    summary["morans_I_rusle_delta"] = mi_sl

    # ── Save outputs ──────────────────────────────────────────────────
    (RES_DIR / "validation_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\n  Saved -> results/validation_summary.json")

    # Flat CSV summary for thesis-table embedding
    rows = [
        ("Carbon loss (2001-2024)", "Total CO2e (GgCO2e)",
         f"{co2_total:.2f}  [95% CI {co2_lo:.2f}, {co2_hi:.2f}]"),
        ("Carbon loss (2001-2024)", "Total area (ha)",
         f"{area_total:.0f}  [95% CI {area_lo:.0f}, {area_hi:.0f}]"),
        ("Trend (CO2e)", "Mann-Kendall Z",  f"{mk_co2['Z']:+.3f}"),
        ("Trend (CO2e)", "Mann-Kendall p",  f"{mk_co2['p']:.4f}"),
        ("Trend (CO2e)", "Trend direction", mk_co2["trend"]),
        ("Trend (CO2e)", "Sen slope (GgCO2e/yr)",
         f"{sen_co2['slope']:+.4f}  [95% CI {sen_co2['ci_lo']:+.4f}, {sen_co2['ci_hi']:+.4f}]"),
        ("Trend (area)", "Mann-Kendall Z",  f"{mk_area['Z']:+.3f}"),
        ("Trend (area)", "Mann-Kendall p",  f"{mk_area['p']:.4f}"),
        ("Trend (area)", "Sen slope (ha/yr)",
         f"{sen_area['slope']:+.2f}  [95% CI {sen_area['ci_lo']:+.2f}, {sen_area['ci_hi']:+.2f}]"),
        ("Cross-validation", "n pixels", f"{n_pairs:,}"),
        ("Cross-validation", "RMSE (Mg/ha)",
         f"{rmse:.2f}  [95% CI {boots_cv['rmse_ci'][0]:.2f}, {boots_cv['rmse_ci'][1]:.2f}]"),
        ("Cross-validation", "Bias (Mg/ha)",
         f"{bias:+.2f}  [95% CI {boots_cv['bias_ci'][0]:+.2f}, {boots_cv['bias_ci'][1]:+.2f}]"),
        ("Cross-validation", "Pearson r",
         f"{pr:.4f}  [95% CI {boots_cv['r_ci'][0]:.4f}, {boots_cv['r_ci'][1]:.4f}]"),
        ("Cross-validation", "R^2",
         f"{r2:.4f}  [95% CI {boots_cv['r2_ci'][0]:.4f}, {boots_cv['r2_ci'][1]:.4f}]"),
        ("Cross-validation", "Lin CCC",
         f"{ccc:.4f}  [95% CI {boots_cv['ccc_ci'][0]:.4f}, {boots_cv['ccc_ci'][1]:.4f}]"),
        ("Cross-validation", "OLS slope / intercept",
         f"{ols.slope:.3f}  /  {ols.intercept:+.2f}"),
        ("Cross-validation", "RMA slope / intercept",
         f"{rma['slope']:.3f}  /  {rma['intercept']:+.2f}"),
        ("Spatial autocorrelation", "Moran's I (habitat delta)",
         f"{mi['I']:.4f}  (p = {mi['p']:.4f}, n = {mi['n']:,})"),
        ("Spatial autocorrelation", "Moran's I (RUSLE delta)",
         f"{mi_sl['I']:.4f}  (p = {mi_sl['p']:.4f}, n = {mi_sl['n']:,})"),
    ]
    pd.DataFrame(rows, columns=["block", "metric", "value"]).to_csv(
        RES_DIR / "validation_metrics.csv", index=False)
    print(f"  Saved -> results/validation_metrics.csv")

    # Carbon trend figure with Sen line + bootstrap CI ribbon
    save_trend_figure(years, co2, sen_co2, mk_co2)


def save_cv_figure(gedi: np.ndarray, saatchi: np.ndarray, cv: dict) -> None:
    """Re-render the AGB cross-validation scatter with R^2, CCC, RMA, CIs."""
    fig = plt.figure(figsize=(18, 6), facecolor="white")
    gs = gridspec.GridSpec(1, 3, figure=fig, width_ratios=[1.15, 1.0, 1.1],
                            wspace=0.32, left=0.05, right=0.98,
                            top=0.88, bottom=0.13)

    # ── Panel 1: hexbin + 1:1 + OLS + RMA ──────────────────────────────
    ax = fig.add_subplot(gs[0, 0])
    ax.set_facecolor("#f5f8fc")
    ax.spines[:].set_color("#cccccc")

    n = len(gedi)
    if n > 30_000:
        idx = RNG.choice(n, 30_000, replace=False)
        gp, sp = gedi[idx], saatchi[idx]
    else:
        gp, sp = gedi, saatchi
    hb = ax.hexbin(sp, gp, gridsize=60, cmap="YlGn", mincnt=1, linewidths=0.2)
    cb = fig.colorbar(hb, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("Point density", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    # Trim axes to the bulk of the data (p99.5 of each variable).
    # The full extent (0–600) leaves most of the chart empty and pushes
    # the OLS/RMA divergence into white space; this brings the cloud
    # back into the visible area.
    x_max = float(np.percentile(sp, 99.5)) * 1.05
    y_max = float(np.percentile(gp, 99.5)) * 1.05
    lim = max(x_max, y_max)
    ax.plot([0, lim], [0, lim], "r--", lw=1.4, label="1:1 line", alpha=0.85)

    xs = np.linspace(0, x_max, 200)
    ols, rma = cv["ols"], cv["rma"]
    ax.plot(xs, ols["slope"] * xs + ols["intercept"], color="#1565c0", lw=1.6,
            label=f"OLS  y={ols['slope']:.2f}x{ols['intercept']:+.1f}")
    ax.plot(xs, rma["slope"] * xs + rma["intercept"], color="#6a1b9a", lw=1.6,
            linestyle="--", label=f"RMA  y={rma['slope']:.2f}x{rma['intercept']:+.1f}")

    ax.set_xlim(0, x_max)
    ax.set_ylim(0, y_max)
    ax.set_xlabel("Saatchi AGB (Mg/ha)", fontsize=10, color="#1a1a2e")
    ax.set_ylabel("GEDI L4B AGB (Mg/ha)", fontsize=10, color="#1a1a2e")
    ax.tick_params(colors="#1a1a2e", labelsize=8)
    ax.legend(facecolor="white", edgecolor="#cccccc", fontsize=8,
              labelcolor="#1a1a2e", loc="upper left")
    ax.set_title(
        f"GEDI vs Saatchi AGB  (n compared = {cv['n_pairs']:,};  "
        f"n displayed = {len(gp):,})",
        fontsize=10, fontweight="bold", color="#1a1a2e",
    )
    # Honest framing: this is not a ground-truth validation. Both products
    # are model estimates with different vintages (Saatchi GLAS ~2003-07
    # at 1 km, GEDI L4B 2019-23 at ~1 km); systematic divergence between
    # them is well documented in the tropical-forest literature.
    ax.text(0.98, 0.04,
            "Note: inter-product comparison (different vintages & resolutions);\n"
            "not a ground-truth validation.",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5,
            color="#555",
            bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
                      edgecolor="#cccccc", alpha=0.9))

    # ── Panel 2: AGB distribution histograms (forest pixels) ───────────
    ax_h = fig.add_subplot(gs[0, 1])
    ax_h.set_facecolor("#f5f8fc")
    ax_h.spines[:].set_color("#cccccc")
    bin_top = float(np.percentile(np.concatenate([gp, sp]), 99.5)) * 1.02
    bins = np.linspace(0, bin_top, 55)
    ax_h.hist(gp, bins=bins, color="#2e7d32", alpha=0.70,
              edgecolor="white", linewidth=0.3, label="GEDI L4B")
    ax_h.hist(sp, bins=bins, color="#e65100", alpha=0.60,
              edgecolor="white", linewidth=0.3, label="Saatchi")
    g_med = float(np.median(gp))
    s_med = float(np.median(sp))
    ax_h.axvline(g_med, color="#2e7d32", linestyle="--", lw=1.4,
                 label=f"GEDI median = {g_med:.0f}")
    ax_h.axvline(s_med, color="#e65100", linestyle="--", lw=1.4,
                 label=f"Saatchi median = {s_med:.0f}")
    ax_h.set_xlabel("AGB (Mg/ha)", fontsize=10, color="#1a1a2e")
    ax_h.set_ylabel("Pixel count", fontsize=10, color="#1a1a2e")
    ax_h.tick_params(colors="#1a1a2e", labelsize=8)
    ax_h.legend(facecolor="white", edgecolor="#cccccc", fontsize=8,
                labelcolor="#1a1a2e", loc="upper right")
    ax_h.set_title("AGB Distributions (forest pixels)", fontsize=10,
                   fontweight="bold", color="#1a1a2e")

    # ── Panel 3: metrics table ─────────────────────────────────────────
    ax2 = fig.add_subplot(gs[0, 2])
    ax2.axis("off")
    rows = [
        ("Compared pixels",  f"{cv['n_pairs']:,}"),
        ("Pearson r",        f"{cv['pearson_r']:.4f}   [{cv['r_ci'][0]:.3f}, {cv['r_ci'][1]:.3f}]"
                              if 'r_ci' in cv else f"{cv['pearson_r']:.4f}"),
        ("R-squared",        f"{cv['r2']:.4f}   [{cv['r2_ci'][0]:.3f}, {cv['r2_ci'][1]:.3f}]"),
        ("Lin's CCC",        f"{cv['ccc']:.4f}   [{cv['ccc_ci'][0]:.3f}, {cv['ccc_ci'][1]:.3f}]"),
        ("RMSE (Mg/ha)",     f"{cv['rmse']:.2f}   [{cv['rmse_ci'][0]:.2f}, {cv['rmse_ci'][1]:.2f}]"),
        ("Bias (Mg/ha)",     f"{cv['bias']:+.2f}   [{cv['bias_ci'][0]:+.2f}, {cv['bias_ci'][1]:+.2f}]"),
        ("MAE (Mg/ha)",      f"{cv['mae']:.2f}"),
        ("OLS y = a*x+b",    f"a={cv['ols']['slope']:.3f}, b={cv['ols']['intercept']:+.1f}"),
        ("RMA y = a*x+b",    f"a={cv['rma']['slope']:.3f}, b={cv['rma']['intercept']:+.1f}"),
        ("p-value (Pearson)",f"{cv['p_value']:.2e}"),
    ]
    table = ax2.table(cellText=rows,
                       colLabels=["Statistic", "Value (95% CI)"],
                       loc="center", cellLoc="left")
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 1.55)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("#cccccc")
        if r == 0:
            cell.set_facecolor("#1a237e")
            cell.set_text_props(color="white", fontweight="bold")
        elif r % 2 == 0:
            cell.set_facecolor("#eef2f9")
        else:
            cell.set_facecolor("white")
    ax2.set_title("Cross-Validation Metrics  (bootstrap 95% CIs)",
                   fontsize=11, fontweight="bold", color="#1a1a2e", pad=14)

    fig.suptitle("GEDI L4B vs Saatchi AGB — Inter-Product Comparison",
                  fontsize=13, fontweight="bold", color="#1a1a2e", y=1.02)
    fig.savefig(FIG_DIR / "agb_cross_validation_upgraded.png", dpi=200,
                 bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Saved -> {FIG_DIR / 'agb_cross_validation_upgraded.png'}")


def save_trend_figure(years, co2, sen, mk) -> None:
    """Annual carbon loss with Sen trend line + bootstrap CI ribbon."""
    fig, ax = plt.subplots(figsize=(13, 5.5), facecolor="white")
    ax.set_facecolor("#f5f7fa")
    ax.spines[:].set_color("#cccccc")

    ax.bar(years, co2, width=0.7, color="#d32f2f", alpha=0.85,
           label="Annual CO2e loss (GgCO2e)", zorder=3)

    # Sen trend line — anchor at series median
    y_mid = np.median(years)
    co2_mid = np.median(co2)
    xs = np.array([years.min(), years.max()])
    ys = co2_mid + sen["slope"] * (xs - y_mid)
    ax.plot(xs, ys, color="#1565c0", lw=2.2, linestyle="--",
            label=f"Sen slope = {sen['slope']:+.3f} GgCO2e/yr "
                  f"(95% CI [{sen['ci_lo']:+.3f}, {sen['ci_hi']:+.3f}])",
            zorder=4)

    # CI envelope around the Sen line. Both lines pivot through the
    # series centroid (median year, median CO2e), so this visualises
    # slope uncertainty — it is not a prediction interval. CO2e loss is
    # strictly non-negative, so we clip the ribbon at zero to avoid the
    # unphysical extrapolation below the x-axis at the series endpoints.
    xs_dense = np.linspace(xs[0], xs[1], 200)
    ys_lo = np.clip(co2_mid + sen["ci_lo"] * (xs_dense - y_mid), 0, None)
    ys_hi = np.clip(co2_mid + sen["ci_hi"] * (xs_dense - y_mid), 0, None)
    ax.fill_between(xs_dense, ys_lo, ys_hi, color="#1565c0", alpha=0.12, zorder=2,
                     label="Sen slope envelope (95% CI)")

    # Annotation: MK p-value
    ax.text(0.02, 0.96,
             f"Mann-Kendall Z = {mk['Z']:+.3f}\np = {mk['p']:.4f}  ({mk['trend']})",
             transform=ax.transAxes, va="top", fontsize=9.5,
             bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                       edgecolor="#cccccc", alpha=0.95),
             color="#1a1a2e")

    ax.set_xlabel("Year", fontsize=11, color="#1a1a2e")
    ax.set_ylabel("CO2e lost (GgCO2e)", fontsize=11, color="#1a1a2e")
    ax.tick_params(colors="#1a1a2e")
    ax.legend(facecolor="white", edgecolor="#cccccc", labelcolor="#1a1a2e",
              fontsize=9, loc="upper right")
    ax.set_title(
        "Annual Carbon Loss with Sen's Slope Trend & Mann-Kendall Test\n"
        "Andaman & Nicobar Islands (2001-2024)",
        fontsize=12, fontweight="bold", color="#1a1a2e")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "carbon_loss_trend_with_ci.png", dpi=200,
                 bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Saved -> {FIG_DIR / 'carbon_loss_trend_with_ci.png'}")


if __name__ == "__main__":
    main()
