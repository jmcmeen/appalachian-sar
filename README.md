# Species-Area Relationships of Amphibians in the Southern Appalachians

Artifacts and updated analyses from:

> **Stout, J.B., Jessee, L.D., & McMeen, J.N. (2025).** Nested and island models for determining the species-area relationship of southern Appalachian amphibians. *Journal of North American Herpetology*, 2025(1), 1–7. [https://doi.org/10.17161/jnah.v2025i1.21557](https://doi.org/10.17161/jnah.v2025i1.21557)

## About the Paper

The southern Appalachian ecoregion is the global biodiversity hotspot for salamander genera, with high rates of endemism among amphibians driven by elevational gradients and limited dispersal ability. Despite a rich history of naturalist investigations, no quantitative species–area relationship (SAR) for amphibians had been established for the region until this study.

Stout et al. (2025) built two SAR models using the classical power function S = CA^z fitted via log-linear OLS regression. The **nested model** uses seven concentric land areas ranging from a single wetland pond (0.0018 km²) to all of southern Appalachia (178,758 km²). The **island model** uses 26 independent, well-surveyed sites (parks, counties, and other localities) within the study area. Both models demonstrate that land area is a significant predictor of amphibian species richness, with the nested model producing stronger correlations (R² = 0.93–0.99) than the island model (R² = 0.56–0.76). The paper recommends averaging the two models for practical prediction at unstudied localities.

## Repository Contents

```text
data/
  nested.csv                    Nested area-richness data (7 concentric areas, Appendix 1)
  island.csv                    Island area-richness data (26 independent sites, Appendix 2)
publication/
  Stout_et_al_2025.pdf          Published paper
  figures/                      Figures from the publication
analysis.ipynb                  Reproduction of the original paper's SAR analysis
analysis_updates.ipynb          NLS reanalysis using the sars library
comparative_analysis.md         Full methodological comparison write-up
README.md                       This file
requirements.txt                Python dependencies (includes sars)
```

## Data

### Nested Model (Appendix 1)

Seven concentric areas of increasing size, each encompassing all smaller areas:

| Area ID | Locality | Area (km²) | Amphibians | Frogs | Salamanders |
| --- | --- | ---: | ---: | ---: | ---: |
| 1 | Wetland pond at Steele Creek Park | 0.0018 | 10 | 7 | 3 |
| 2 | Steele Creek Park | 9.3 | 22 | 10 | 12 |
| 3 | Sullivan County, TN | 1,114 | 38 | 13 | 25 |
| 4 | Northeastern Tennessee | 4,137 | 44 | 13 | 31 |
| 5 | Eastern Tennessee | 37,438 | 69 | 18 | 52 |
| 6 | + W North Carolina + W Virginia | 101,271 | 88 | 20 | 68 |
| 7 | Southern Appalachia | 178,758 | 98 | 22 | 76 |

### Island Model (Appendix 2)

Twenty-six independent sites within the southern Appalachian ecoregion, ranging from John's Bog (0.0081 km²) to Great Smoky Mountains NP (2,114 km²). Sites include national and state parks, national recreation areas, counties, and other well-surveyed localities across Tennessee, North Carolina, West Virginia, Virginia, and Georgia. See `data/island.csv` for the full dataset.

## Analysis

### Original SAR Analysis

The main notebook (`analysis.ipynb`) reproduces the paper's log-linear OLS analysis for both datasets:

1. **Log-log regression** of species richness on area for amphibians, frogs, and salamanders — separately for nested and island models
2. **Power model derivation** (S = CA^z) from the regression parameters
3. **Nested vs. island comparison** reproducing the paper's Figure 2, with both regression lines overlaid for each taxonomic group
4. **Predictive estimates** at reference areas, including the averaged model the paper recommends

The nested model produces the paper's reported parameters exactly: amphibians C = 18.87, z = 0.1222, R² = 0.97; salamanders C = 8.36, z = 0.1734, R² = 0.99; frogs C = 9.34, z = 0.0598, R² = 0.93.

### Updated Analysis with `sars`

The updated notebook (`analysis_updates.ipynb`) reanalyzes both datasets using the [`sars`](https://pypi.org/project/sars/) Python library, which fits SAR models via nonlinear least squares (NLS) in arithmetic space rather than the log-linear OLS used in the original calculations. The NLS approach is generally preferred in modern SAR literature because it does not impose the variance structure assumed by log-transformation (Tjørve & Tjørve 2021). The notebook covers:

1. **Power-law fits** for all three taxonomic groups on both datasets, with NLS parameter estimates and comparison to original values
2. **Cross-group and cross-dataset comparison** — parameter tables and overlay plots showing nested vs. island NLS fits
3. **Alternative model comparison** for weak fits — power, linear, logarithmic, and powerR models compared for island groups (all R² < 0.80), with residual tables for island salamanders where the linear model outperforms the power model
4. **Residual diagnostics** — 2×3 panel grid across both datasets
5. **Multi-model comparison** across all 20 SAR models supported by `sars`, for both datasets
6. **Threshold analysis** testing for breakpoints / small-island effects
7. **Model averaging** and **bootstrap confidence intervals** (fully available for the island dataset at n = 26)
8. **Predictive comparison** of NLS vs. original log-linear predictions at reference areas

The full methodological discussion — parameter shifts between fitting methods, the salamander anomaly in the island model, nested vs. island structural comparison, z-value interpretation, and ecological implications — is in [`comparative_analysis.md`](comparative_analysis.md).

## Setup

```bash
pip install -r requirements.txt
jupyter notebook analysis.ipynb
```

The updated analysis requires the [`sars`](https://pypi.org/project/sars/) library (included in `requirements.txt`).
