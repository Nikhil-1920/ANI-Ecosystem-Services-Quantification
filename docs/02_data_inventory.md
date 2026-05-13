# 02 · Data Inventory & Reference Cards

This document details the eight primary spatial datasets used in the
ANI Ecosystem Services project, plus one ancillary vector layer.
Google Earth Engine (GEE) acts as the retrieval engine for most of
these layers; the underlying datasets are produced by different
scientific consortia using different satellite constellations.

---

## 2.1 ESA WorldCover 2021 — "What is where?"

| Field | Value |
|---|---|
| Source | European Space Agency (ESA CCI) |
| Sensors | Sentinel-1 (radar) + Sentinel-2 (multispectral optical) |
| Native resolution | 10 m |
| Native CRS | EPSG:4326 |
| Project copy | `data/processed/ANI_ESA_WorldCover_mosaic_clipped.tif` |

**Role in project.**
- Stratifies biomass analysis (Week 6): separates evergreen-forest
  carbon from mangrove carbon.
- Defines habitat sensitivity *H* and threat sources (Week 7).
- Drives the RUSLE C-factor lookup (Week 8).

ESA WorldCover is the **central tabulation key** for the project — every
per-class statistic in the synthesis is keyed on its class codes.

---

## 2.2 Global Forest Watch (Hansen et al.) — "When was it cut down?"

| Field | Value |
|---|---|
| Source | University of Maryland & WRI / Google Cloud |
| Sensors | Landsat 5 / 7 / 8 (optical) |
| Native resolution | 30 m |
| Native CRS | EPSG:4326 |
| Bands used | `lossyear` (0–23 for 2001–2023), `treecover2000` baseline, `gain` (2000–2012) |
| Project copies | `ANI_GFW_Forest_Loss_2001_2023_clipped.tif`, `ANI_GFW_TreeCover2000_Baseline_clipped.tif`, `ANI_GFW_Forest_Gain_clipped.tif`, `ANI_GFW_DataMask_Land_Water_clipped.tif` |

**Role in project (Week 6).** Per-year forest-loss masks are overlaid on
the GEDI baseline to compute annual carbon emissions. The gain layer
(2000–2012) drives the offset side of the net-balance calculation; the
DataMask is used to construct the binary land/water mask in
preprocessing.

---

## 2.3 GEDI L4B Aboveground Biomass Density — "How much carbon is stored?"

| Field | Value |
|---|---|
| Source | NASA Global Ecosystem Dynamics Investigation (mission v2.1) |
| Instrument | Full-waveform lidar mounted on the International Space Station, 242 laser pulses / sec, $\pm 51.6^\circ$ lat coverage |
| Native resolution | 1 km gridded (calibrated against 25-m footprints) |
| Project copy | `ANI_GEDI_Biomass_Density_clipped.tif` (bilinearly resampled to 30 m) |

**Role in project (Week 6).** Forms the baseline AGB raster
*B(x, y)* — the multiplicand for every per-pixel carbon-loss
calculation. The product reports mean AGBD in Mg ha⁻¹ via a
random-forest model calibrated on field plots.

Optical sensors can only see the *top* of a forest canopy. They
cannot easily differentiate a 10 m and a 40 m canopy. GEDI's lasers
penetrate the canopy and estimate the actual physical weight of the
trees — this is why it is the right baseline for this study.

---

## 2.4 Saatchi 2011 Pan-Tropical Biomass — "Independent comparison product"

| Field | Value |
|---|---|
| Source | Saatchi et al. 2011 (NASA JPL / Wageningen mirror) |
| Methodology | Field plots + GLAS lidar + optical, 1 km gridded |
| Vintage | ~2003–2007 GLAS waveforms |
| Project copy | `ANI_Saatchi_AGB_CrossValidation_clipped.tif` |

**Role in project (Week 6).** Used as an **inter-product comparison**
for GEDI L4B — not as a ground-truth validation. Both products are
model estimates with different vintages and methodologies, so
disagreement between them quantifies *product-level divergence* rather
than the error of either against truth.

The comparison is restricted to forest-only pixels (ESA 10 Tree cover
and 95 Mangroves; n = 2,314,513) and reports Pearson r, R², Lin's CCC,
RMSE, bias, OLS regression and Reduced-Major-Axis regression with
bootstrap CIs.

**Headline results.** Pearson r = 0.254, log–log r = 0.358, CCC = 0.021,
RMSE = 168.8 Mg/ha, bias = +152.6 Mg/ha (GEDI higher than Saatchi),
RMA slope = 4.56. A diagnostic test that rescaled Saatchi by ×2.13
(Mg C/ha → Mg AGB/ha) closed only half of the bias gap, confirming
this is a product-level divergence and not a unit-conversion error.
The pattern is consistent with the documented conservatism of Saatchi
2011 in tall closed-canopy tropical forest (Mitchard et al. 2014).

Full statistical detail is in **07_validation.md §3**.

---

## 2.5 SRTM 30 m DEM — "How steep is the terrain?"

| Field | Value |
|---|---|
| Source | NASA + NGA, Shuttle Radar Topography Mission (Feb 2000) |
| Native resolution | 30 m |
| Project copy | `ANI_SRTM_DEM_30m_clipped.tif` |

**Role in project (Week 8).** Provides slope and flow-accumulation
inputs to the RUSLE LS factor (slope-length-and-steepness). ANI's
mountainous interior (Saddle Peak in north Andaman, central Great
Nicobar ridge) produces high LS values that combine with deforestation
to drive the bare-pixel hotspots seen in the RUSLE delta map.

---

## 2.6 CHIRPS Annual Precipitation — "How hard is it raining?"

| Field | Value |
|---|---|
| Source | UC Santa Barbara Climate Hazards Center |
| Methodology | IR satellite cloud-top temperatures + rain-gauge station calibration |
| Native resolution | ~5.5 km (0.05°) |
| Project copies | `ANI_CHIRPS_Annual_Total_Precip_clipped.tif`, `ANI_CHIRPS_Mean_Precip_2000_2023_clipped.tif` |

**Role in project (Week 8).** Drives the RUSLE R-factor (rainfall
erosivity) and the freshwater-yield supplementary module
(Section 4 of `08_results_synthesis.md`). ANI receives 2,500–4,000
mm yr⁻¹, producing R-factor values among the highest globally and
explaining the very high erosion potential on the few deforested
parcels.

---

## 2.7 SoilGrids Clay+Silt Fractions — "What type of dirt is this?"

| Field | Value |
|---|---|
| Source | ISRIC SoilGrids v2 (via Google Earth Engine asset) |
| Methodology | Machine-learning prediction from 250 k+ soil profiles and environmental covariates |
| Native resolution | 250 m |
| Project copy | `ANI_SoilGrids_Clay_Silt_GEE_clipped.tif` (two-band: clay %, silt %) |

**Role in project (Week 8).** Plugged into the simplified K-factor
relation `K = 0.005 + 0.0006·(clay+silt)`, then clipped to the
plausible physical range [0.005, 0.060]. The Soil Organic Carbon
supplementary module (08 §4) also uses this layer in a **hybrid
baseline**: pedotransfer formula `SOC = 40 + 0.5·(clay+silt)` where
SoilGrids reports valid data, and literature class-mean SOC values
where the SoilGrids raster is zero-filled (a coarse-to-fine resample
artefact). Mean baseline SOC under the hybrid scheme is ~75.3 Mg C/ha.

---

## 2.8 ANI Administrative Boundary — "Cookie-cutter polygon"

| Field | Value |
|---|---|
| Source | Survey of India / GADM Level-1 (filtered to Andaman & Nicobar) |
| Format | Vector shapefile |
| Project copy | `data/raw/ANI_Administrative_Boundary.shp` |

**Role.** Defines the region-of-interest polygon used for clipping
every raster in preprocessing. Once reprojected to EPSG:32646 it is
shared by all downstream scripts.

---

## 2.9 Summary — Dataset × Service matrix

| Dataset | Carbon (W6) | Habitat (W7) | Soil (W8) |
|---|:---:|:---:|:---:|
| ESA WorldCover | ✅ | ✅ | ✅ |
| GFW Forest Loss / Gain | ✅ |   |   |
| GEDI L4B AGBD | ✅ |   |   |
| Saatchi 2011 AGB | (validation) |   |   |
| SRTM DEM |   |   | ✅ |
| CHIRPS |   |   | ✅ |
| SoilGrids Clay+Silt |   |   | ✅ |
| ANI boundary | ✅ | ✅ | ✅ |

**ESA WorldCover** is the most critical dataset: it provides the
common land-class framework for all three primary models.

---

## 2.10 Glossary of Acronyms

| Acronym | Full form |
|---|---|
| ESA | European Space Agency |
| GFW | Global Forest Watch (Hansen et al.) |
| GEDI | Global Ecosystem Dynamics Investigation (NASA lidar) |
| SRTM | Shuttle Radar Topography Mission |
| CHIRPS | Climate Hazards Group InfraRed Precipitation with Station data |
| SoilGrids | ISRIC global digital soil mapping system |
| OSM | OpenStreetMap |
| IPCC | Intergovernmental Panel on Climate Change |
| RUSLE | Revised Universal Soil Loss Equation |
| InVEST | Integrated Valuation of Ecosystem Services and Tradeoffs |
| AGB / AGBD | Aboveground Biomass / Aboveground Biomass Density |
| UTM | Universal Transverse Mercator |
