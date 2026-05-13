# 07 · Validation — Trend, Cross-Comparison, Spatial Statistics

This addendum documents the additional statistical rigour applied to
the ANI Ecosystem Services pipeline. All tests were run on data
already present in `data/processed/` and `results/`; no new external
dataset was introduced. Companion to `04_methods_carbon.md`,
`05_methods_habitat.md`, and `06_methods_soil.md`.

**Source code.** `src/validation_stats.py`.
**Outputs.** `results/validation_summary.json`,
`results/validation_metrics.csv`, `figures/validation/`.

---

## 1. Trend Analysis of Annual Deforestation (2001–2024)

### 1.1 Mann-Kendall Test

We applied the non-parametric Mann-Kendall (MK) test to the annual
CO₂-equivalent loss and the annual deforested-area series. MK is
robust to non-normality and seasonality and is the standard test for
monotonic trend in environmental time-series (Mann 1945; Kendall 1975).
Tie-corrected variance was used:

$$ S = \sum_{i<j} \mathrm{sgn}(x_j - x_i),\qquad
   \mathrm{Var}(S) = \frac{n(n-1)(2n+5) - \sum_t t(t-1)(2t+5)}{18} $$

### 1.2 Sen's Slope Estimator

The Theil-Sen median slope was used as a robust trend-magnitude
estimator; the 95 % confidence interval was derived from the normal
approximation of Hollander & Wolfe (1973).

**Note on the slope-envelope band shown in `carbon_loss_trend_with_ci.png`.**
The upper and lower bounding lines both pivot through the data
centroid `(median(year), median(CO₂e))`, so the visible band is a
**slope envelope** rather than a predictive confidence interval. The
band is also **clipped at zero** because annual CO₂e loss is strictly
non-negative — an un-clipped extrapolation of the lower CI line would
dip below the x-axis at the right edge of the chart, which is
unphysical. The legend label was updated from "95 % CI" to
"Sen envelope (95 % CI)" to set expectations correctly.

### 1.3 Results

| Series | MK Z | p-value | Trend direction | Sen slope (95% CI) |
|---|---|---|---|---|
| Annual CO₂e (GgCO₂e/yr) | −2.505 | 0.0122 | decreasing | −2.668 [−5.254, −0.656] |
| Annual area (ha/yr)     | −3.349 | 0.0008 | decreasing | −54.1 [−83.8, −20.4]    |

Interpretation: deforestation in ANI has slowed significantly over the
study window. The decrease is most plausibly attributable to (a) earlier
loss of the most accessible forest, (b) post-2008 implementation of the
A&N Islands (Protection of Aboriginal Tribes) Regulation tightening land
conversion, and (c) reduced agricultural-expansion pressure after the
2004 tsunami displaced settlements. This trend should **not** be
interpreted as ecosystem-services recovery — accumulated losses persist
even as new losses slow.

---

## 2. Bootstrap Confidence Intervals on Carbon-Loss Totals

To quantify sampling uncertainty in the 24-year total, we drew 2,000
non-parametric bootstrap resamples (with replacement) of the annual
series and computed the sum each time.

| Metric                            | Point estimate | 95% bootstrap CI |
|---|---|---|
| Total CO₂e lost, 2001–2024 (GgCO₂e) | 1,390.34       | [936.13, 1,887.49] |
| Total area lost, 2001–2024 (ha)     | 19,502         | [12,351, 27,826]   |

These intervals capture variability in the annual loss series itself,
*not* measurement error in GEDI biomass — the latter is reported
separately in the cross-validation section.

---

## 3. GEDI L4B vs Saatchi AGB Cross-Validation (Upgraded)

The original pipeline reported Pearson r, RMSE, bias, and MAE between
the GEDI L4B aboveground biomass density product and the Saatchi
pan-tropical AGB raster. We extended this with four additional
agreement statistics, restricted the comparison to forested land
cover only, and bootstrapped 95% CIs on every metric.

### 3.1 Pixel-selection rule

Comparison was restricted to ESA WorldCover classes **10 (Tree cover)**
and **95 (Mangroves)**. Cropland, built-up, water, and bare pixels were
excluded because Saatchi's AGB raster is conditioned on forest cover,
whereas GEDI's gap-filled grid returns non-zero values over any class
with measured returns — comparing across all land would inject a
structural bias unrelated to instrument accuracy.

### 3.2 Statistics computed

| Statistic | Definition | Use |
|---|---|---|
| Pearson r | linear correlation | classical fit, sensitive to outliers |
| R² | r-squared | variance explained |
| Lin's CCC | concordance correlation | agreement with the 1:1 line (not just association) |
| RMSE | √mean[(x−y)²] | mean disagreement magnitude |
| Bias | mean(GEDI − Saatchi) | systematic offset |
| MAE | mean |GEDI − Saatchi| | robust mean error |
| OLS | y = ax + b (least squares) | standard regression |
| RMA | y = ±(σ_y/σ_x)x + b | reduced major axis — error on both axes |

Bootstrap CIs (500 resamples) were computed on RMSE, bias, r, R², and
CCC.

### 3.3 Results (n = 2,314,513 forest pixels)

| Statistic | Value (95 % CI) |
|---|---|
| Pearson r (linear)  | 0.254 [0.250, 0.259]   |
| Pearson r (log–log) | **0.358** (single-pass; bootstrap CI not computed) |
| R²              | 0.064 [0.062, 0.067]   |
| Lin CCC         | 0.021 [0.021, 0.022]   |
| RMSE (Mg/ha)    | 168.76 [168.45, 169.11] |
| Bias (Mg/ha)    | +152.59 [+152.27, +152.93] (GEDI > Saatchi) |
| MAE (Mg/ha)     | 152.76                  |
| OLS  | y = 1.158·x + 142.74    |
| RMA  | y = 4.559·x − 68.89     |

**Linear vs log–log r.** The jump from linear r = 0.254 to log–log
r = 0.358 is informative: the two products **rank pixels in a similar
relative order** (the log–log signal) even though their absolute
scales differ by a factor of ~3.5× (the bias-dominated linear signal).
A diagnostic test that rescaled Saatchi by 2.13 (the Mg C/ha → Mg AGB
conversion) closed only **half** of the bias gap, ruling out a simple
unit-conversion error and confirming that the residual disagreement
is a true product-level divergence.

### 3.4 Interpretation

The two products **disagree systematically**, with GEDI showing on
average +153 Mg/ha higher biomass than Saatchi over forested ANI. This
is not an instrument failure: it is consistent with the published
literature on Saatchi 2011 (Mitchard et al. 2014; Avitabile et al. 2016),
which reports that the original 1-km pan-tropical AGB product
**saturates near 200–250 Mg/ha** in dense tropical evergreen forests,
systematically underestimating the high-biomass tail captured by GEDI's
laser waveforms. The GEDI L4B product is the more recent and physically
more sensitive estimator for tropical evergreen canopies.

**Practical implication for this study.** We use the Saatchi RMSE as a
*conservative upper bound* on GEDI absolute uncertainty when reporting
carbon-loss ranges, while acknowledging that the true GEDI uncertainty
(per Dubayah et al. 2022) is substantially lower (~30–40 Mg/ha for
30-m gridded means). The Saatchi-based ±value should therefore be read
as worst-case, not best-estimate.

---

## 4. Spatial Autocorrelation of Ecosystem Collapse (Moran's I)

To test whether the spatial clustering of ecosystem degradation
hotspots is statistically significant (vs. arising by chance from
random pixel-level damage), we computed global Moran's I with
queen-case (3 × 3) contiguity on two of the delta rasters:

- `habitat_quality_delta.tif` — change in habitat-quality index
- `rusle_soil_loss_delta.tif` — change in annual soil-loss rate

Both rasters were mean-pooled to a 50× coarser grid (~1.5 km
effective resolution, 510 × 150 cells) to make the permutation test
tractable while preserving spatial structure. Inference was performed
via 199 random label permutations.

### 4.1 Results

| Layer | Moran's I | p (perm) | n cells |
|---|---|---|---|
| Habitat-quality Δ (degradation cells, <0)  | **0.333** | **0.005** | 2,260 |
| RUSLE soil-loss Δ (increase cells, >0)     | **0.202** | **0.005** | 1,374 |
| Permutation null distribution (habitat)    | mean ≈ 0, std ≈ 0.013 | — | — |

Both layers show **highly significant positive spatial autocorrelation**.
Translation: ecosystem-services collapse in ANI does not occur as
isolated random pixels — it **clusters**. Degraded pixels are
neighboured by other degraded pixels far more often than chance allows.
This justifies the synthesis-level hotspot mapping (`synthesis_hotspots.py`,
ECI raster) and the prioritisation of compact contiguous areas
for intervention, rather than scattered remediation across the
archipelago.

---

## 5. Outputs

| File | Contents |
|---|---|
| `results/validation_summary.json` | Full machine-readable results |
| `results/validation_metrics.csv` | Flat thesis-table-ready metrics |
| `figures/validation/agb_cross_validation_upgraded.png` | Upgraded scatter + metric table |
| `figures/validation/carbon_loss_trend_with_ci.png` | Annual loss with Sen line & MK p-value |

---

## References

- Mann, H. B. (1945). *Nonparametric tests against trend.* Econometrica.
- Kendall, M. G. (1975). *Rank Correlation Methods.* Griffin, London.
- Hollander, M. & Wolfe, D. A. (1973). *Nonparametric Statistical Methods.* Wiley.
- Lin, L. I.-K. (1989). *A concordance correlation coefficient to evaluate reproducibility.* Biometrics.
- Moran, P. A. P. (1950). *Notes on continuous stochastic phenomena.* Biometrika.
- Mitchard, E. T. A., et al. (2014). *Markedly divergent estimates of Amazon forest carbon density from ground plots and satellites.* Global Ecology and Biogeography.
- Avitabile, V., et al. (2016). *An integrated pan-tropical biomass map using multiple reference datasets.* Global Change Biology.
- Dubayah, R., et al. (2022). *GEDI L4A footprint level aboveground biomass density.* Earth System Science Data.
- Harris, N. L., et al. (2021). *Global maps of twenty-first century forest carbon fluxes.* Nature Climate Change.
