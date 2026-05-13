# 01 · Research Plan & Theoretical Framework

> **Project:** *Quantifying Changes in Ecosystem Services in the Andaman
> & Nicobar Islands due to Forest Conversion using Open Global Carbon
> and Biomass Datasets.*
> **Type:** Independent Study Report (12-week schedule, completed).

---

## 1.1 Overarching Goal

Estimate carbon-storage, habitat-quality and soil-retention losses
arising from forest conversion in the Andaman & Nicobar Islands (ANI)
between 2001 and 2024, using only **open global remote-sensing
products** plus published parameter values, and deliver a publication-
ready synthesis with figures, tables and economic interpretation.

---

## 1.2 Research Questions

A high-impact paper needs 3–4 specific, quantifiable research questions
that guide the entire methodology:

- **RQ1 — Spatio-temporal patterns.** What are the dominant spatial
  and temporal patterns of forest conversion across ANI from 2001 to
  2024, and how do these patterns differ between primary evergreen
  forest and insular mangrove ecosystems?
- **RQ2 — Carbon emissions.** How much aboveground biomass (AGB) and
  corresponding CO₂-equivalent has been lost due to deforestation, and
  how do baseline GEDI L4B estimates — cross-compared with the regional
  Saatchi 2011 product — constrain the uncertainty of these emissions?
- **RQ3 — Habitat degradation.** Using an InVEST-style Habitat Quality
  model, how has proximity to anthropogenic land-cover transitions
  (quantified via ESA WorldCover) degraded the habitat quality of
  remaining contiguous forest blocks, and how sensitive is this
  degradation to varying threat decay distances?
- **RQ4 — Ecosystem-service synergies.** Where are the critical
  spatial hotspots in which high carbon loss, severe habitat
  degradation and elevated soil-erosion risk (RUSLE) intersect, and
  what are the implications for targeted conservation policy in the
  archipelago?

---

## 1.3 Literature-Survey Matrix

To write a defensible Introduction and Methods section, papers were
organised into three buckets.

### Bucket A — Ecosystem services in tropical islands (ANI-specific)
**Goal:** establish the ecological importance of ANI and the specific
threats (logging, agriculture, infrastructure) it faces; identify the
research gap that this work fills (no prior paper has combined 3-D
GEDI lidar + GFW + InVEST + RUSLE on this archipelago).
**Search terms:** `("Andaman and Nicobar" OR "Andaman Islands") AND
("deforestation" OR "forest cover change" OR "ecosystem services")`.

### Bucket B — Methodological precedents (GEDI + GFW)
**Goal:** justify the use of GEDI L4B and Global Forest Watch to
reviewers; show that this combination is state-of-the-art for tropical
carbon accounting. **Search terms:** `("GEDI" OR "Global Ecosystem
Dynamics Investigation") AND ("Hansen" OR "Global Forest Watch") AND
"carbon loss"`. **Key focus:** papers that document GEDI's tropical-
canopy uncertainty and how authors cross-validate with other products.

### Bucket C — Spatial modelling (InVEST Habitat + RUSLE Erosion)
**Goal:** find peer-reviewed precedent for the exact parameters used
in the InVEST and RUSLE models — particularly threat weights, decay
distances, and tropical C-factor values.
**Search terms:** `("InVEST Habitat Quality" OR "RUSLE") AND
("tropical forest" OR "island ecosystem") AND "GIS"`.

---

## 1.4 Study Area Definition

- **Geographic scope.** The entire Andaman and Nicobar archipelago
  (Union Territory of India), comprising approximately 572 islands
  scattered between ~6.5°N and ~13.5°N in the eastern Bay of Bengal,
  totalling ~8,249 km² of land area. The Andaman group lies to the
  north and the Nicobar group to the south, separated by the
  Ten-Degree Channel. The FAO GAUL Level 1 boundary is used for
  spatial clipping.
- **Ecological scope.** Both major forest biomes present in the
  archipelago are analysed separately: tropical evergreen / semi-
  evergreen forest (ESA WorldCover class 10) and mangrove forest
  (ESA WorldCover class 95). ESA's discrimination of these two types
  permits stratified carbon and habitat reporting.
- **Coordinate Reference System.** All data are reprojected to
  **WGS 84 / UTM Zone 46N (EPSG:32646)** for planar-area and
  distance-based calculations (RUSLE LS, InVEST distance-to-threat,
  coastal-protection distance-to-water).

---

## 1.5 12-Week Schedule (as executed)

| Month | Window | Outputs |
|---|---|---|
| 1 — Foundation & data | Wk 1–2 | Research questions, lit-review buckets, study area definition |
|   | Wk 3 | GEE asset inventory; Avitabile/Saatchi cross-validation layer procured |
|   | Wk 4–5 | All 9 GeoTIFFs reprojected to EPSG:32646 at 30 m, clipped to ANI boundary, masked to land (`processed/` folder) |
| 2 — Core analysis | Wk 6 | GEDI baseline + GFW × GEDI carbon-loss accounting; CO₂e time-series |
|   | Wk 7 | InVEST-style habitat quality map + sensitivity sweep |
|   | Wk 8 | RUSLE soil-loss surface + counterfactual delta |
| 3 — Synthesis & write-up | Wk 9 | Ecosystem Collapse Index (ECI) raster + bivariate / hexbin syntheses |
|   | Wk 10 | Publication-quality figures (light backgrounds, side-by-side Andaman/Nicobar panels, log-scale violin charts) |
|   | Wk 11–12 | Discussion, limitations, validation addenda, economic forecast, LaTeX report |

Implementation order is `src/preprocess.py` → `src/carbon_analysis.py`
→ `src/habitat_quality.py` → `src/soil_retention.py` →
`src/synthesis_hotspots.py` → `src/supplementary_services.py` →
`src/validation_stats.py` → `src/render_synthesis_light.py`.

---

## 1.6 What Each Subsequent Doc Covers

| Doc | Topic |
|---|---|
| **02_data_inventory.md** | Per-dataset reference card (source, resolution, role, processing notes) |
| **03_preprocessing.md** | The reprojection / clipping / masking pipeline and the integrity check |
| **04_methods_carbon.md** | GEDI baseline, GFW × GEDI loss accounting, Saatchi inter-product comparison |
| **05_methods_habitat.md** | InVEST-style habitat-quality model with threat decay and sensitivity sweep |
| **06_methods_soil.md** | RUSLE factor build, A = R·K·LS·C·P, per-class statistics |
| **07_validation.md** | Mann–Kendall + Sen's slope, bootstrap CIs, Moran's *I*, log-log r |
| **08_results_synthesis.md** | All headline numbers (carbon, habitat, soil, ECI, economic) |
| **09_discussion_and_limitations.md** | Supplementary services, parameter defensibility, scope caveats |
