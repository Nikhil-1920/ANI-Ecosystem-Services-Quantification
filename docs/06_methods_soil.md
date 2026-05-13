# 06 · Methods — Soil Retention (RUSLE)

## 6.1 Concept

ANI's combination of **extreme monsoon rainfall** (CHIRPS annual
mean > 3,400 mm/yr over tree cover, > 3,840 mm/yr along the mangrove
fringe), **steep volcanic / tectonic relief** (Saddle Peak block on
North Andaman, Great Nicobar interior ridge) and **thin Inceptisols/
Ultisols** makes the archipelago exceptionally vulnerable to soil
loss once vegetation cover is removed. Forest acts as the *physical
anchor* for the land; deforestation amplifies the per-pixel erosion
rate by two-to-three orders of magnitude and the eroded sediment
ultimately suffocates the coastal coral fringes.

We quantify this using the **Revised Universal Soil Loss Equation
(RUSLE)** (Renard et al. 1997). The full equation is

$$
A \;=\; R \cdot K \cdot LS \cdot C \cdot P
$$

with *A* in tonnes ha⁻¹ yr⁻¹ and the five factors derived from CHIRPS,
SoilGrids, SRTM, ESA WorldCover and a global default (Support practice
*P* = 1.0 for natural tropical forest).

Source code: `src/soil_retention.py`. Outputs in `results/`:
`rusle_soil_loss.tif`, `rusle_soil_loss_delta.tif`,
`rusle_erosion_by_landcover.csv`. Figures in `figures/soil/`.

---

## 6.2 Factor build

### R — Rainfall erosivity
Derived from CHIRPS annual precipitation using the Renard–Freimund
(1994) tropical regression:

$$
R \;=\; 0.0483 \cdot P^{1.61}
$$

with *P* in mm yr⁻¹. ANI's high annual rainfall produces R values in
the 1,000 – 3,000 MJ·mm·ha⁻¹·h⁻¹·yr⁻¹ range — among the highest globally —
reflecting extreme tropical rainfall erosivity.

### K — Soil erodibility
Approximated from the SoilGrids clay+silt fraction using a simplified
linear relation:

$$
K \;=\; 0.005 + 0.0006 \cdot (\text{clay} + \text{silt})
$$

clipped to the plausible physical range $[0.005, 0.060]$
t·ha·h·ha⁻¹·MJ⁻¹·mm⁻¹. Tropical soils typically fall between 0.01 and
0.04 in this unit system.

### LS — Slope length and steepness
Computed from the SRTM 30 m DEM via the
flow-accumulation formulation of Moore et al. (1992):

$$
LS \;=\; \left( \frac{A_s}{22.13} \right)^{0.4}
        \cdot \left( \frac{\sin\theta}{0.0896} \right)^{1.3}
$$

with a **steep-slope cap of 50** applied to suppress non-physical
cliff-edge artefacts.

### C — Cover management (per ESA WorldCover class)

| ESA class | Land cover | C-factor |
|---|---|---|
| 10 | Tree cover (closed canopy) | **0.001** |
| 95 | Mangroves | 0.010 |
| 90 | Herbaceous wetland | 0.020 |
| 20 | Shrubland | 0.035 |
| 30 | Grassland | 0.060 |
| 40 | Cropland | **0.280** |
| 60 | Bare / sparse | **0.450** |
| 50 | Built-up | 0.000 (no sediment yield by RUSLE convention) |

Forest C = 0.001 vs Bare C = 0.45 implies that a deforested pixel
experiences ~450× more erosion than an intact forest pixel — that
factor is the direct quantification of the soil-retention service.

### P — Support practice
P = 1.0 uniformly (no field-validated terracing or contour-cropping
data is available at 30 m for ANI).

---

## 6.3 Important implementation notes

### 6.3.1 Output cap at 500 t·ha⁻¹·yr⁻¹
The RUSLE product *A = R·K·LS·C·P* is **clipped to
500 t·ha⁻¹·yr⁻¹** in `rusle_soil_loss.tif` because steep DEM cells
can produce non-physical LS spikes that drag *A* into the thousands.
Cropland and bare/sparse pixels therefore pile up at the ceiling,
which inflates the upper shoulder of the violin distribution. The
synthesis figure caption notes this explicitly.

### 6.3.2 Real raster medians (use these in the report)
The per-class median values **from the real raster** are:

| Class | Median (t/ha/yr) | Mean | Max | Total (t/yr) |
|---|---|---|---|---|
| Tree cover | **2.83** | 3.67 | 41.0 | 1,927,237 |
| Mangroves | 12.41 | 17.34 | 336.4 | 853,674 |
| Wetland | 10.53 | 20.19 | 500.0 | 4,402 |
| Shrubland | 33.44 | 44.34 | 210.2 | 395 |
| Grassland | 56.36 | 87.12 | 500.0 | 1,151,794 |
| Cropland | **160.88** | 205.68 | 500.0 | 732,513 |
| Bare / sparse | **424.81** | 321.57 | 500.0 | 133,012 |
| Built-up | 0.00 | 0.00 | 0.0 | 0 |

These are the values used by every synthesis figure and by Section 3.2
of the LaTeX report — they replace the older lookup-table-plus-
synthetic-noise values that briefly appeared in earlier doc drafts.

### 6.3.3 Andaman / Nicobar split panels
Every map output (`rusle_factor_components_map.png`,
`rusle_soil_loss_map.png`, `rusle_soil_loss_delta_map.png`) renders
the Andaman and Nicobar groups in separate sub-axes — left and right
respectively. The two archipelagos sit ~600 km apart in UTM 46N and
would otherwise visually merge into one stretched panel.

### 6.3.4 Counterfactual ΔA against a forested baseline
`rusle_soil_loss_delta.tif` reports `A_current − A_all_forested`. This
delta is the soil loss **caused by land-cover conversion** — the
ecosystem-service quantity. The hotspots in the delta map align
almost perfectly with the cropland belt around Port Blair and the
small bare patches on Great Nicobar's southern slopes. The total
"new" erosion attributable to forest conversion sums to roughly
**163,000 t·yr⁻¹**, of which ~74 % is contributed by grassland and
cropland combined.

---

## 6.4 Outputs

| File | Contents |
|---|---|
| `results/rusle_soil_loss.tif` | Per-pixel *A* (t/ha/yr), capped at 500 |
| `results/rusle_soil_loss_delta.tif` | Counterfactual ΔA vs all-forested baseline |
| `results/rusle_erosion_by_landcover.csv` | Mean / median / max / total per ESA class |
| `figures/soil/rusle_factor_components_map.png` | All four factors split Andaman/Nicobar |
| `figures/soil/rusle_soil_loss_map.png` | Final A surface (log-scaled colour bar) |
| `figures/soil/rusle_soil_loss_delta_map.png` | ΔA hotspots |
| `figures/soil/rusle_erosion_by_landcover.png` | Per-class bar chart |
| `figures/synthesis/stat_distribution_soil_erosion.png` | Per-class violin (log y, dashed outline for low-n) |

---

## 6.5 Action checklist

- [ ] Confirm `data/processed/` contains all five RUSLE-input rasters
      (CHIRPS, SoilGrids clay+silt, SRTM DEM, ESA WorldCover,
      `ANI_land_mask.tif`).
- [ ] Run `python src/soil_retention.py` — confirm
      `results/rusle_soil_loss.tif` is produced and statistically
      sensible (island-wide mean ≈ 7.9 t·ha⁻¹·yr⁻¹).
- [ ] Open `figures/soil/rusle_factor_components_map.png` and verify
      each of the four factor panels is spatially coherent
      (R highest on wetter northern Andamans; LS concentrated on the
      Saddle Peak block; C dominated by the dark-forest background
      with rare bright cropland/bare patches).
- [ ] Confirm `rusle_erosion_by_landcover.csv` reports the medians
      tabulated in §6.3.2.

---

## 6.6 Key references

| Reference | Used for |
|---|---|
| Wischmeier & Smith (1978) | Original USLE framework, mangrove + wetland C-factor defaults |
| Renard et al. (1997) | RUSLE parameter guidelines |
| Renard & Freimund (1994) | R-factor regression from annual precipitation |
| Moore et al. (1992) | LS factor flow-accumulation form |
| Williams (1995) / EPIC | K-factor pedotransfer (referenced — simplified linear form used in production) |
| Panagos et al. (2015) | Per-class C-factor benchmark |
| Borrelli et al. (2020) | Pan-tropical RUSLE C-factor defaults |
| Montgomery (2007) | Soil-formation rate context (1–11 t·ha⁻¹·yr⁻¹) |
