# Quantifying Changes in Ecosystem Services in the Andaman & Nicobar Islands due to Forest Conversion (2001 – 2024)

An Independent Study project using open global remote-sensing data to
quantify changes in three biophysical ecosystem services — **carbon
storage**, **soil retention** and **habitat quality** — across the
Andaman & Nicobar archipelago between 2001 and 2024, and to synthesise
the three losses into a single Ecosystem Collapse Index (ECI) raster.

**Authors:** Nikhil Singh, Mihir Chaudhary
**Guide:** Prof. Ramachandra Prasad
**Project type:** Independent Study, 12-week schedule
**Repository:** <https://github.com/Nikhil-1920/ANI-Ecosystem-Services-Quantification>

---

## What is in this repository

```
ANI-Ecosystem-Services-Quantification/
├── data/
│   ├── raw/                  # Raw GeoTIFFs (ESA, GFW, GEDI, Saatchi, CHIRPS,
│   │                         #   SoilGrids, SRTM) + ANI admin boundary shapefile
│   └── processed/            # Reprojected to UTM 46N, clipped to ANI, masked
│
├── src/                      # Python pipeline (see "Pipeline order" below)
│   ├── preprocess.py         #   not all scripts are listed here; see folder
│   ├── carbon_analysis.py
│   ├── habitat_quality.py
│   ├── soil_retention.py
│   ├── synthesis_hotspots.py
│   ├── supplementary_services.py
│   ├── validation_stats.py
│   ├── render_synthesis_light.py    # final figure renderer
│   └── ...                   # additional rendering / sensitivity helpers
│
├── results/                  # CSV outputs and intermediate rasters
│   ├── carbon_annual_loss_by_year.csv
│   ├── habitat_quality_by_landcover.csv
│   ├── rusle_erosion_by_landcover.csv
│   ├── eci_collapse_hotspots_summary.csv
│   ├── supplementary_services.csv
│   ├── validation_summary.json
│   └── *.tif intermediate rasters
│
├── figures/                  # PNGs grouped by analytical theme
│   ├── carbon/         5 figures (baseline, time-series, hotspots, loss-vs-gain, net balance)
│   ├── habitat/        4 figures (index map, delta, by-landcover, sensitivity)
│   ├── soil/           4 figures (factor components, soil-loss map, delta, by-landcover)
│   ├── synthesis/      6 figures (violins, bivariate KDE, radar, ECI, hexbin)
│   ├── validation/     2 figures (Mann-Kendall trend, inter-product comparison)
│   ├── supplementary/  1 figure  (4-panel coastal + water + pollination + SOC)
│   └── predictive/     3 figures (forest-cover 2040, tri-scenario 2060, economic 2060)
│
├── docs/                     # Project documentation (9 numbered chapters + index)
│   ├── README.md
│   ├── 01_research_plan.md
│   ├── 02_data_inventory.md
│   ├── 03_preprocessing.md
│   ├── 04_methods_carbon.md
│   ├── 05_methods_habitat.md
│   ├── 06_methods_soil.md
│   ├── 07_validation.md
│   ├── 08_results_synthesis.md
│   └── 09_discussion_and_limitations.md
│
├── report/                   # LaTeX project report
│   ├── main.tex              # 12 pt, 1.5 spacing, Times-like, environmental
│   │                         #   colour theme; ~45 pages with 25 figures
│   └── figures/              # Local copy of the figures used in the report
│
├── requirements.txt
└── README.md                 # this file
```

---

## Pipeline order

The Python scripts in `src/` are deterministic given the same inputs
(NumPy seed 42 is set explicitly). Run them in this order from the
repo root after activating the virtual environment:

```bash
source venv/bin/activate

# 1. Preprocessing (one-time): reproject + clip + mask to a 30 m UTM 46N grid
python src/preprocess.py

# 2. Core service models
python src/carbon_analysis.py        # GEDI baseline + GFW × GEDI carbon-loss
python src/habitat_quality.py        # InVEST-style habitat quality + sensitivity
python src/soil_retention.py         # RUSLE A = R × K × LS × C × P

# 3. Synthesis & supplementary
python src/synthesis_hotspots.py     # Ecosystem Collapse Index (ECI) raster
python src/supplementary_services.py # Coastal / freshwater / pollination / SOC

# 4. Validation (trend + spatial autocorrelation + cross-comparison)
python src/validation_stats.py

# 5. Render the figures used in the report
python src/render_synthesis_light.py
```

Outputs land in `results/` and `figures/`. Re-running a single
downstream stage (e.g. just the figures) does not require re-running
the earlier stages — the intermediate `.tif` and `.csv` files are
cached on disk.

---

## Documentation index

The `docs/` folder is the long-form companion to the LaTeX report. It
reads as a linear sequence of nine numbered chapters:

| # | File | Topic |
|---|---|---|
| 01 | `01_research_plan.md` | Objectives, research questions (RQ1–4), study area, 12-week schedule |
| 02 | `02_data_inventory.md` | Per-dataset reference cards (source, resolution, role) |
| 03 | `03_preprocessing.md` | The reprojection / clipping / masking pipeline |
| 04 | `04_methods_carbon.md` | GEDI baseline + GFW × GEDI accounting, Saatchi inter-product comparison |
| 05 | `05_methods_habitat.md` | InVEST-style habitat quality with threat decay + sensitivity sweep |
| 06 | `06_methods_soil.md` | RUSLE factor build, real-raster per-class medians |
| 07 | `07_validation.md` | Mann–Kendall + Sen's slope (envelope clipped at zero), bootstrap CIs, Moran's I, log–log Pearson r |
| 08 | `08_results_synthesis.md` | Headline numbers for carbon, habitat, soil, ECI, 2024–2060 economic scenarios |
| 09 | `09_discussion_and_limitations.md` | Supplementary services, parameter defensibility, scope caveats |

The previous folder layout (13 overlapping files with weekly diaries
and a separate methodology summary) was consolidated into this linear
9-chapter structure in the final pass.

---

## Headline findings (2001 – 2024)

| Quantity | Value (95 % CI where given) |
|---|---|
| Total forest loss 2001–2023 | 19,502 ha (12,351 – 27,826) |
| Total CO₂e loss | 1,390 Gg CO₂e (936 – 1,887) |
| Net carbon debt (loss − gain) | ~ 882 Gg CO₂e |
| Mann–Kendall trend Z (CO₂e) | −2.51 (p = 0.012, decreasing) |
| Sen slope (CO₂e) | −2.67 Gg / yr (95 % CI −5.25 to −0.66) |
| RUSLE median, Tree Cover | 2.83 t / ha / yr |
| RUSLE median, Bare/Sparse | 424.81 t / ha / yr |
| Habitat Quality, Tree Cover (median) | 1.00 |
| Habitat Quality, Cropland (median) | 0.10 |
| ECI critical hotspots (Andamans + Nicobars) | 2,404 ha |
| Moran's I, habitat-quality Δ | 0.33 (p = 0.005) |
| GEDI vs Saatchi, Pearson r (linear / log–log) | 0.254 / 0.358 |
| Economic damage 2040 (lower bound, BAU) | US$501.8 M |
| Economic damage 2040 (upper bound, with supplementary services) | US$594.4 M (+ 18 %) |

The full headline table with confidence intervals is reproduced in
the LaTeX report (`report/main.tex`, Table 4) and in
`docs/08_results_synthesis.md`.

---

## Reproducibility

- Python 3.11.
- Core libraries: `numpy`, `scipy`, `pandas`, `matplotlib`,
  `rasterio`, `geopandas`, `pyproj`, `scikit-image`, `pillow`.
- See `requirements.txt` for exact pins.
- All scripts respect a single NumPy random seed (42) so that the
  same input rasters produce byte-identical CSV outputs.
- `results/validation_summary.json` is the machine-readable record of
  every numeric quantity reported in the LaTeX document.

---

## Building the report

The LaTeX source under `report/` compiles with any standard pdflatex
distribution (tested on TeX Live 2024). On Overleaf, upload the
`report/` folder and compile twice for the TOC, List of Figures and
List of Tables to populate. The figures referenced by the report are
included locally under `report/figures/` so the project compiles
without requiring the upstream `figures/` tree.

```bash
cd report
pdflatex main.tex && pdflatex main.tex
```

---

## Acknowledgements

This work uses freely-available remote-sensing data: ESA WorldCover
2021 (ESA CCI), Global Forest Watch annual loss layers (Hansen et al.
/ WRI), the NASA GEDI L4B aboveground biomass density product, the
Saatchi 2011 pan-tropical biomass map, CHIRPS annual precipitation
(UCSB Climate Hazards Center), SoilGrids v2 (ISRIC), and the SRTM 30 m
DEM (NASA / USGS). The Survey of India administrative boundary
shapefile was used to delimit the study area. Python, NumPy, SciPy,
matplotlib, pandas, rasterio and GeoPandas are gratefully
acknowledged.
