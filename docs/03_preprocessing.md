# 03 · Spatial Preprocessing & Data Harmonisation

## 3.1 Why this step exists

Raw inputs come from multiple sensors (NASA, ESA, UMD, ISRIC, UCSB)
with different projections, resolutions and nodata conventions, so
they cannot be analysed together in their native state. This phase
transforms **eight heterogeneous global datasets** plus one vector
boundary into a **single, unified 30 m research database** specifically
optimised for ANI.

The preprocessing script is `src/preprocess.py`. Its outputs land in
`data/processed/` and are the entry point for every Week-6 / 7 / 8 /
synthesis script.

---

## 3.2 The "golden standard" target grid

To ensure pixel-perfect alignment for mathematical operations (e.g.
multiplying biomass by forest-loss masks), every raster is forced into
the following grid:

| Parameter | Standard value | Rationale |
|---|---|---|
| Coordinate Reference System | **EPSG:32646 (UTM Zone 46N)** | Metric projection that minimises area distortion for ANI |
| Spatial resolution | **30 m** | Matches the native resolution of GFW and SRTM; finest practical resolution for the whole stack |
| File format | **Cloud-Optimised GeoTIFF** | Standard for fast raster I/O via rasterio |
| Region of interest | **ANI administrative boundary** | Survey-of-India / GADM polygon, reprojected to EPSG:32646 |
| Grid dimensions | **25,500 rows × 7,503 columns** | Single shared transform across the stack |
| Pixel area | **0.09 ha (900 m²)** | Used in every per-pixel mass / area conversion |

---

## 3.3 Six-step pipeline

### Step 1 — ESA WorldCover mosaicking
The raw ESA tiles are merged into a single seamless raster via
`rasterio.merge` so there are no seams between adjacent tiles spanning
North and South Andaman or the Nicobars.

### Step 2 — Reprojection
Raw data from GEE arrives in WGS 84 (EPSG:4326). Each raster is warped
into the UTM 46N metric grid using
`rasterio.warp.calculate_default_transform` + `reproject` so that, for
example, a road in the OSM layer falls geographically on top of the
forest pixels it might be degrading.

### Step 3 — Resampling (pixel size alignment)
Two methods are used, chosen by data type:

- **Nearest-neighbour** for *categorical* layers — ESA WorldCover land
  cover, GFW `lossyear` (years 0–23). Preserves integer codes.
- **Bilinear interpolation** for *continuous* layers — GEDI biomass,
  Saatchi biomass, SRTM elevation, CHIRPS precipitation, SoilGrids
  clay+silt. Produces smooth gradients.

### Step 4 — Land / water masking
The GFW DataMask band is extracted and converted to a binary land mask
(1 = land, 0 = water). This is used to blank out ocean pixels — every
downstream statistic ignores them so that no analysis "averages over
ocean".

### Step 5 — Exact spatial clipping
Using `ANI_Administrative_Boundary.shp` as a cookie-cutter, all rasters
are clipped to the precise island polygon via `rasterio.mask.mask`.
This eliminates ~70 % of the bounding box and shrinks file sizes by an
order of magnitude.

### Step 6 — Vector harmonisation (roads, if used)
Any OSM road shapefile or settlement polygon is reprojected to
EPSG:32646 and clipped to the boundary. (In the final pipeline, road
threats for the habitat-quality module are derived from the built-up
ESA class itself, so OSM rasterisation is only used in sensitivity
runs.)

---

## 3.4 Quality control

After the pipeline runs, a `verify_alignment.py`-style check confirms
every output GeoTIFF has the same width, height, transform and CRS as
the reference ESA raster. The expected output stack:

```
data/processed/
├── ANI_ESA_WorldCover_mosaic_clipped.tif
├── ANI_GFW_Forest_Loss_2001_2023_clipped.tif
├── ANI_GFW_TreeCover2000_Baseline_clipped.tif
├── ANI_GFW_Forest_Gain_clipped.tif
├── ANI_GFW_DataMask_Land_Water_clipped.tif
├── ANI_GEDI_Biomass_Density_clipped.tif
├── ANI_Saatchi_AGB_CrossValidation_clipped.tif
├── ANI_SRTM_DEM_30m_clipped.tif
├── ANI_CHIRPS_Annual_Total_Precip_clipped.tif
├── ANI_CHIRPS_Mean_Precip_2000_2023_clipped.tif
├── ANI_SoilGrids_Clay_Silt_GEE_clipped.tif
└── ANI_land_mask.tif
```

All files share an identical 25 500 × 7 503 grid, with NaN propagated
through `np.nan_to_num` only at the very end of each analysis (so that
ocean pixels never leak into the statistics).

---

## 3.5 Output summary

The preprocessing pipeline produces a deterministic, fully-aligned 30 m
analysis stack. From this point onward a simple Python expression like
`biomass[lossyear == 22]` works instantly across the entire 800 km
island chain — which is the whole reason the preprocessing phase
exists.

## 3.6 Technical stack used

- **Rasterio** — core raster I/O, warping and masking.
- **GeoPandas** — shapefile handling for the ANI boundary.
- **PyProj** — CRS transformations to EPSG:32646.
- **NumPy** — array-based masking and NoData handling.
- **scikit-image** — block-reduce helpers used by the cross-validation
  block-average to 1 km grid (Week 6).
