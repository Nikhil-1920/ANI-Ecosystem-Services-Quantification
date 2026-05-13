# ANI Ecosystem Services — Documentation Index

Linear reading order for the project documentation. Each file's number
matches its position in the natural research-paper flow (plan → data
→ preprocessing → methods → validation → results → discussion).

| # | File | Topic |
|---|---|---|
| 01 | [`01_research_plan.md`](01_research_plan.md) | Project goal, research questions (RQ1–4), literature-survey buckets, study area, 12-week schedule |
| 02 | [`02_data_inventory.md`](02_data_inventory.md) | The eight remote-sensing datasets used (ESA, GFW, GEDI, Saatchi, SRTM, CHIRPS, SoilGrids, ANI boundary) — sources, resolutions, roles |
| 03 | [`03_preprocessing.md`](03_preprocessing.md) | Reprojection / clipping / masking pipeline to a unified 30 m UTM 46N grid (`src/preprocess.py`) |
| 04 | [`04_methods_carbon.md`](04_methods_carbon.md) | GEDI baseline + GFW × GEDI annual carbon-loss accounting, Saatchi inter-product comparison (`src/carbon_analysis.py`) |
| 05 | [`05_methods_habitat.md`](05_methods_habitat.md) | InVEST-style Habitat Quality model with threat decay + sensitivity sweep (`src/habitat_quality.py`) |
| 06 | [`06_methods_soil.md`](06_methods_soil.md) | RUSLE soil-loss surface, factor build, counterfactual delta, real-raster per-class medians (`src/soil_retention.py`) |
| 07 | [`07_validation.md`](07_validation.md) | Mann–Kendall + Sen's slope (with envelope clipping), bootstrap CIs, log-log r, Moran's *I* (`src/validation_stats.py`) |
| 08 | [`08_results_synthesis.md`](08_results_synthesis.md) | Headline carbon / habitat / soil / ECI numbers, 2024–2060 economic-damage scenarios |
| 09 | [`09_discussion_and_limitations.md`](09_discussion_and_limitations.md) | Supplementary services (coastal protection, freshwater, pollination, SOC), parameter defensibility, scope caveats |

## What was consolidated in this restructure

The previous folder contained 13 markdown files with overlap between
`quantification_methodology.md` and the week-by-week diary files
(`week_1_2_foundation.md` … `week_8_soil_retention.md`). The restructure
merged them as follows:

| Original files (deleted) | Merged into |
|---|---|
| `research_plan.md` + `week_1_2_foundation.md` + `week_3_advanced_data.md` | **01_research_plan.md** |
| `dataset_reference.md` | **02_data_inventory.md** (expanded) |
| `preprocessing_report.md` + `week_4_5_preprocessing.md` | **03_preprocessing.md** |
| `quantification_methodology.md` (carbon section) + `week_6_carbon_analysis.md` | **04_methods_carbon.md** |
| `quantification_methodology.md` (habitat section) + `week_7_habitat_quality.md` | **05_methods_habitat.md** |
| `quantification_methodology.md` (soil section) + `week_8_soil_retention.md` | **06_methods_soil.md** |
| `validation_methods.md` | **07_validation.md** (renamed) |
| `final_results_synthesis.md` | **08_results_synthesis.md** (renamed, economic numbers updated) |
| `discussion_and_limitations.md` | **09_discussion_and_limitations.md** (renamed, $603 M → $594 M correction) |

13 docs → 9 docs. All cross-references in the new files use the new
numbered filenames. Outdated values (coastal NPV, SOC baseline, upper-
bound total, hero-figure references, lookup-table RUSLE values) have
been corrected to match the current pipeline outputs.

## Pipeline order (in case you want to re-run from scratch)

```
src/preprocess.py            # → 03_preprocessing.md
src/carbon_analysis.py       # → 04_methods_carbon.md
src/habitat_quality.py       # → 05_methods_habitat.md
src/soil_retention.py        # → 06_methods_soil.md
src/synthesis_hotspots.py    # → 08_results_synthesis.md §5 (ECI)
src/supplementary_services.py# → 09_discussion_and_limitations.md §1.1
src/validation_stats.py      # → 07_validation.md
src/render_synthesis_light.py# → publication figures in figures/synthesis/
```
