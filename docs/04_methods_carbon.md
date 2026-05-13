# 04 · Methods — Carbon Storage & Deforestation Emissions

## 4.1 Concept

Tropical forests are the world's most powerful land carbon sinks.
Quantifying carbon loss from ANI deforestation lets us measure the
islands' contribution to (and exposure to) global climate change.

The analytical framework follows **IPCC Tier 1 allometric accounting**:

- *Aboveground Biomass* (AGB) — measured in Mg ha⁻¹, sourced from GEDI.
- *Carbon stock* — 1 Mg of dry biomass ≈ 0.47 Mg of carbon (IPCC default).
- *CO₂ equivalent* — 1 unit of carbon ≈ 3.667 units of CO₂ (the
  molecular-mass ratio 44/12).

Source code: `src/carbon_analysis.py`. Outputs in `results/`:
`carbon_annual_loss_by_year.csv`, `carbon_loss_co2e_tonnes.tif`.
Figures in `figures/carbon/` and `figures/validation/`.

---

## 4.2 Mathematical foundation

For every pixel `(x, y)` and every year `t ∈ {2001, …, 2023}` encoded
by the GFW lossyear raster:

$$
\Delta C_t(x, y) \;=\; B(x, y) \cdot M_t(x, y) \cdot 0.09 \cdot 0.47
$$

where *B(x, y)* is the GEDI baseline AGB in Mg ha⁻¹, *M_t(x, y)* is the
binary mask `lossyear == t − 2000`, the constant 0.09 converts
per-hectare biomass to per-pixel biomass for a 30-m pixel, and 0.47 is
the IPCC biomass-to-carbon fraction. Annual CO₂-equivalent emissions
follow from the molecular-mass ratio:

$$
\text{CO}_2\text{e}_t \;=\; \Delta C_t \cdot \frac{44}{12} \;\approx\; 3.667 \cdot \Delta C_t .
$$

Aggregated to a national time-series, this yields the
`carbon_annual_loss_by_year.csv` table consumed by the trend tests
(see **07_validation.md §1**) and by every figure in
`figures/carbon/`.

---

## 4.3 Part 1 — Establishing the GEDI baseline

The GEDI L4B raster `ANI_GEDI_Biomass_Density_clipped.tif` is loaded
and any obviously invalid values (`AGB < 0` or `AGB > 800`) are masked
to NaN. The Andaman mean AGB is **208 Mg ha⁻¹** with a
95th-percentile of 336 Mg ha⁻¹; Nicobar is **184 Mg ha⁻¹** with a
95th-percentile of 322 Mg ha⁻¹. The densest biomass is on the central
spine of North and Middle Andaman (Saddle Peak block) and on the deep
interior of Great Nicobar.

The script also stratifies AGB by ESA WorldCover class:

| ESA class | Description | Use in this study |
|---|---|---|
| 10 | Tree cover | Tropical evergreen / semi-evergreen forest |
| 20 | Shrubland | Edge / regrowth scrub |
| 95 | Mangroves | Coastal mangrove fringe (separately reported) |

These two forest classes are also the pixel set used for the GEDI-vs-
Saatchi inter-product comparison.

The full baseline surface is rendered in
`figures/carbon/agb_gedi_baseline_map.png` with the Andaman and
Nicobar groups in separate panels and the GFW gain pixels overlaid in
magenta. Section 3.1 of the LaTeX report discusses the spatial pattern
in detail.

---

## 4.4 Part 2 — Inter-product comparison against Saatchi 2011

The independent Saatchi 2011 pan-tropical AGB raster is used to
quantify product-level divergence in the GEDI baseline. Both products
are reprojected to the same UTM 46N grid; only forest-only pixels
(ESA 10 + 95, n = 2,314,513) are retained.

Eight agreement metrics are computed (Pearson r, R², Lin's CCC, RMSE,
bias, MAE, OLS slope/intercept, RMA slope/intercept), each with a
500-resample bootstrap 95 % CI on metrics that admit one.

**Why "inter-product comparison" and not "validation."** Both Saatchi
and GEDI are model estimates with different vintages and methodologies
(Saatchi 2011 calibrated against ~2003–07 GLAS waveforms at 1 km; GEDI
L4B against 2019–23 full-waveform shots). Disagreement between them
quantifies *product-level divergence* rather than the error of either
against true ground biomass. Full methodology and the result table are
in **07_validation.md §3**.

The figure for this comparison is
`figures/validation/agb_cross_validation_upgraded.png`.

---

## 4.5 Part 3 — Annual carbon-loss calculation

The Hansen GFW `lossyear` band encodes the year of loss as integers
1–23 (= 2001–2023). For each year code, the script extracts the binary
mask, multiplies pixel-wise by the GEDI baseline, converts to carbon
and to CO₂e, and writes one row to the annual CSV. Headline outputs:

| Quantity | Value (95 % bootstrap CI) |
|---|---|
| Total area lost 2001–2023 | **19,502 ha** (12,351 – 27,826) |
| Total CO₂e lost 2001–2023 | **1,390 Gg CO₂e** (936 – 1,887) |
| Largest single-year spike | 2005, ~180 Gg CO₂e |
| Forest *gain* 2000–2012 (offset) | 5,806 ha → ~508 Gg CO₂e sequestered |
| Cumulative **net** balance to end-2024 | ~ 882 Gg CO₂e of standing debt |

The annual time-series and the cumulative curve are plotted in
`figures/carbon/carbon_annual_loss_timeseries.png` and
`figures/validation/carbon_loss_trend_with_ci.png`. The spatial
distribution of cumulative emissions is in
`figures/carbon/carbon_loss_hotspots_map.png`.

---

## 4.6 Part 4 — Loss vs gain and the net balance

The GFW gain layer covers 2000–2012 only — there is no equivalent
gain product for 2013–2024 yet — so the net-balance calculation has
an inherent forward bias. The figure
`figures/carbon/carbon_loss_vs_gain_spatial.png` overlays the cumulative
loss (red) and cumulative gain (yellow) on the intact-land basemap.
Three observations follow:

1. **Loss and gain are spatially anticorrelated.** Loss clusters along
   roads and settlement edges; gain occurs in the deep interior.
2. **Loss exceeds gain by ~3.4×** (19,502 ha vs 5,806 ha).
3. **The gain ledger is incomplete** for post-2012 regrowth.

Despite the offset, the net-balance time-series remains negative
(loss-dominated) in every year of the record — the archipelago is a
net carbon source even after granting full credit to the GFW gain
product.

---

## 4.7 Uncertainty quantification

The bootstrap 95 % CI on the 24-year total (936 – 1,887 Gg CO₂e) is
dominated by the spatial variance of the GEDI baseline, not by
temporal uncertainty in the GFW loss count. A 10 % systematic
overestimate in mean baseline AGB would inflate the total by the same
10 %. The Saatchi-vs-GEDI bias of +152.6 Mg ha⁻¹ (Section 4.4)
gives a conservative upper bound on absolute GEDI uncertainty, though
the documented GEDI L4B accuracy (Dubayah et al. 2022) is much
tighter (~30–40 Mg ha⁻¹ for gridded means).

---

## 4.8 Key equations summary (Methods section of the paper)

| Step | Equation | Source |
|---|---|---|
| AGB baseline | GEDI L4B 1 km product, bilinear-resampled to 30 m | Dubayah et al. (2022) |
| Biomass → Carbon | `C = AGB × 0.47` | IPCC Tier 1 (2006) |
| Carbon → CO₂e | `CO₂e = C × (44/12)` | IPCC (2006) |
| Annual loss | `ΔC_t = Σ (AGB_px × 0.47 × 0.09 ha)` for `lossyear == t` | Hansen et al. (2013) |
| Bootstrap CI | 2,000 resamples of the annual series, sum each | Efron & Tibshirani (1993) |
| Inter-product check | Pearson r, CCC, RMA on forest-only pixels | Lin (1989); Warton (2006) |

---

## 4.9 Action checklist (reproducibility)

- [ ] Confirm `data/processed/ANI_GEDI_Biomass_Density_clipped.tif`,
      `ANI_GFW_Forest_Loss_2001_2023_clipped.tif`,
      `ANI_GFW_Forest_Gain_clipped.tif`,
      `ANI_ESA_WorldCover_mosaic_clipped.tif`,
      `ANI_Saatchi_AGB_CrossValidation_clipped.tif` exist.
- [ ] Run `python src/carbon_analysis.py` — confirm
      `results/carbon_annual_loss_by_year.csv` and the eight
      `figures/carbon/*.png` are produced.
- [ ] Run `python src/validation_stats.py` — confirm
      `results/validation_summary.json` contains the GEDI–Saatchi
      block and the bootstrap CIs on the 24-year totals.
