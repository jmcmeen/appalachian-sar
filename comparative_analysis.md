# Comparative Analysis: NLS vs. Log-Linear SAR Fitting

A companion to [`analysis_updates.ipynb`](analysis_updates.ipynb), comparing the nonlinear least squares (NLS) results from the [`sars`](https://pypi.org/project/sars/) Python library against the original log-linear OLS calculations published in Stout, Jessee & McMeen (2025).

---

## Fitting Method Differences

The original calculations used log₁₀-transformed OLS regression:

> log₁₀(S) = z · log₁₀(A) + log₁₀(c)

This minimizes squared residuals in log-space and implicitly assumes multiplicative error — i.e., error proportional to the predicted value. The c and z parameters are recovered as c = 10^(intercept) and z = slope.

The `sars` library fits the power model directly in arithmetic space via nonlinear least squares (`scipy.optimize.least_squares`):

> S = cA^z

This minimizes squared residuals in the original species-count scale and assumes additive error. The NLS approach is generally preferred in modern SAR literature because it does not impose the variance structure that log-transformation assumes (Tjørve & Tjørve 2021).

## Parameter Estimates

| Group | Method | c | z | R² |
|---|---|---|---|---|
| Total Amphibians | Original (log-linear) | 13.92 | 0.1079 | 0.757 (log-space) |
| Total Amphibians | `sars` NLS | 12.88 | 0.1272 | 0.725 |
| Frogs | Original (log-linear) | 4.28 | 0.1338 | 0.658 (log-space) |
| Frogs | `sars` NLS | 4.98 | 0.1094 | 0.655 |
| Salamanders | Original (log-linear) | 8.95 | 0.1043 | 0.564 (log-space) |
| Salamanders | `sars` NLS | 7.85 | 0.1383 | 0.581 |

The direction of the c and z shifts is not uniform across groups, which differs from the pattern seen in the etn-sar analysis. For amphibians and salamanders, c drops and z rises under NLS — the expected direction, because log-linear regression overweights small-area sites and pulls c upward and z downward. For frogs, however, c rises slightly (4.28 → 4.98) and z drops (0.1338 → 0.1094). This reversal is attributable to the influence of UT Arboretum (1.01 km², 9 frog species), an outlier site that is highly frog-rich for its size — the NLS objective, which weighs all residuals equally in arithmetic space, reacts more strongly to this point than the log-linear objective does.

The original R² values are computed in log-space and are not directly comparable to the arithmetic-space R² from `sars`. Both, however, tell a consistent story: all three groups show moderate fits (R² ≈ 0.58–0.76), substantially lower than the near-perfect amphibian fit seen in the etn-sar project. This difference is expected: the appalachian-sar dataset spans a much wider geographic range, mixes nested and island sites, and includes localities with independent extinction/colonization dynamics — conditions that increase biological variance around any simple area-based prediction.

## z-Value Interpretation

All six z-values (original and NLS, across all three groups) fall within the canonical mainland/nested range of 0.10–0.20 (Rosenzweig 1995, Drakare et al. 2006).

| Group | z (Original) | z (NLS) | Range Context |
|---|---|---|---|
| Total Amphibians | 0.1079 | 0.1272 | Mainland/nested (0.10–0.20) |
| Frogs | 0.1338 | 0.1094 | Mainland/nested |
| Salamanders | 0.1043 | 0.1383 | Mainland/nested |

This is notable because the dataset includes true island sites (isolated bogs, small preserves, and nature reserves separated from larger continuous forest), which might be expected to push z toward the island range (0.20–0.35). The mainland-range z-values suggest that even the island sites in this dataset behave more like samples of a continuous regional biota than like truly isolated faunas — consistent with the Southern Appalachians being a species-rich source pool where dispersal limitation is rarely severe for amphibians.

## The Mixed Nested + Island Design

Unlike the etn-sar dataset, which comprised four strictly nested areas (Steele Creek Park ⊂ Sullivan County ⊂ NE Tennessee ⊂ Eastern Tennessee), the appalachian-sar dataset mixes:

- **Nested sites:** counties and ecoregional units embedded within the broader Southern Appalachian region
- **Island sites:** isolated bogs, small nature reserves, and management areas that are discrete habitat patches rather than subsets of a larger political or geographic unit

This mixed design has two important consequences. First, it increases the biological variance around the power-law prediction because island and nested sites have different underlying dynamics (Preston 1960, MacArthur & Wilson 1967). Second, it increases the effective range of areas available for regression (0.0081 km² to 2,114 km²), which generally improves parameter precision despite the added variance. The net result is that the appalachian-sar models are less precise per data point but well-constrained across the full area range.

## Outliers and Residual Patterns

Three localities stand out in residual diagnostics across multiple taxonomic groups:

**Whitfield County, GA (753 km²)** is the single largest negative outlier for both total amphibians (observed 20, predicted 30, std. residual = −2.03) and salamanders (observed 8, predicted 20, std. residual = −2.54). Its amphibian count is consistent with sites 10–50× smaller. Whitfield County sits in the Ridge and Valley physiographic province at the southern edge of the study area and is one of the most heavily developed and agriculturally modified counties in the dataset, which likely suppresses salamander richness well below area-based expectations.

**Fannin County, GA (1,015 km²)** shows large negative residuals for total amphibians (std. = −2.47) and salamanders (std. = −2.07). Like Whitfield County, it is located at the southern margin of the region. Its species counts are similar to sites 20–100× smaller in area, suggesting that the local biota is undersurveyed, or that the county's landscape configuration (with significant agricultural lowlands intermixed with mountain forest) reduces effective habitat relative to gross area.

**Great Smoky Mountains NP (2,114 km²)** is the largest positive outlier for both total amphibians (+2.24 std. residuals) and salamanders (+2.26 std. residuals). GSMNP is globally recognized as a salamander diversity hotspot, and its observed counts (45 amphibian species, 33 salamanders) exceed model predictions (34 and 23, respectively) by a substantial margin. This overprediction is ecologically meaningful: GSMNP benefits from exceptional habitat quality, long-term protection, and an unusually high proportion of endemic plethodontid salamanders that reflect deep evolutionary history rather than simple area dynamics.

These patterns suggest that the power model captures the broad area-richness relationship but is sensitive to local deviations driven by land use, habitat quality, and biogeographic distinctiveness.

## Model Selection with n = 26

A key advantage of this dataset over etn-sar is sample size. With n = 26 observations and k ≤ 3 estimated parameters (c, z, σ), the AICc correction term 2k(k+1)/(n−k−1) remains finite for all standard SAR models, enabling full AICc-based model selection, Akaike weights, model averaging, and bootstrap confidence intervals — none of which were applicable in the etn-sar analysis.

The `sars` multi-model comparison tests 20 model forms. Given the moderate R² values for all three groups, it is worth examining whether 3-parameter models (epm1, epm2, powerR) provide meaningful improvement over the standard 2-parameter power model, and whether any non-power-law models (linear, logarithmic, asymptotic) are competitive. The full results and AICc rankings are in [`analysis_updates.ipynb`](analysis_updates.ipynb).

## Threshold Analysis

`sars.sar_threshold()` tests for breakpoints (small-island effect) by comparing piecewise and simple linear models. Unlike etn-sar (where n = 4 made threshold detection meaningless), n = 26 provides enough power for a meaningful test. The wide range of small-area sites (0.0081–50 km²) in the dataset gives the threshold test a plausible opportunity to detect a left-horizontal or steep-slope segment at small areas. Results are reported in the notebook.

## Recommendations for Further Work

**Separate nested and island analyses.** The mixed design inflates residual variance relative to either a purely nested or purely island dataset. Fitting the power model separately to nested sites and island sites, then comparing z-values, would clarify whether the two site types behave differently — a fundamental question in SAR theory (Preston 1960).

**Land-use covariates.** The large residuals for Whitfield and Fannin counties point to anthropogenic modification as a driver of under-richness. Adding impervious surface percentage, forest cover, or a human footprint index as a covariate in a modified power model could substantially reduce residual variance.

**Elevational decomposition.** The Southern Appalachians span a steep elevational gradient from valley floors (~200 m) to high ridges (>2,000 m), and salamander communities in particular are strongly elevation-stratified. Stratifying the dataset by elevation band or physiographic province might reveal tighter SARs within strata.

**Integration with iNaturalist.** The `sars` library supports `sars.from_pyinaturalist()` for direct ingestion of citizen-science observation data, which could update species counts at poorly-sampled localities and potentially add new sites.

## Summary of Key Findings

1. The NLS fitting method in `sars` produces parameter shifts consistent with the different error assumptions of the two approaches, though the direction of change is not uniform: amphibians and salamanders show the expected c decrease and z increase, while frogs show the opposite, driven by the UT Arboretum outlier.

2. All z-values fall within the mainland/nested range (0.10–0.20), despite the inclusion of island sites, suggesting that Southern Appalachian amphibians experience the region as a continuous species pool rather than a true archipelago.

3. The moderate R² values (0.58–0.73) reflect genuine biological variance arising from the mixed nested/island design, wide geographic scope, and local drivers (land use, biogeographic history) that the power model cannot capture. They do not indicate a failed analysis, but rather the limits of a simple univariate area model applied across a complex landscape.

4. With n = 26, AICc-based model selection and bootstrap confidence intervals are fully available, enabling rigorous model comparison that was impossible in the etn-sar dataset.

5. Three localities are notable outliers: Whitfield County and Fannin County (both in GA) fall far below predictions, likely due to land-use modification; Great Smoky Mountains NP far exceeds predictions, reflecting its status as a global salamander diversity hotspot.

## References

- Drakare, S., Lennon, J.J. & Hillebrand, H. (2006). The imprint of the geographical, evolutionary and ecological context on species–area relationships. *Ecology Letters* 9:215–227.
- MacArthur, R.H. & Wilson, E.O. (1967). *The Theory of Island Biogeography*. Princeton University Press.
- Preston, F.W. (1960). Time and space and the variation of species. *Ecology* 41:611–627.
- Rosenzweig, M.L. (1995). *Species Diversity in Space and Time*. Cambridge University Press.
- Stout, J.B., Jessee, L.D. & McMeen, J.N. (2025). Nested and island models for determining the species-area relationship of southern Appalachian amphibians. *Journal of North American Herpetology* 2025(1):1–7.
- Tjørve, E. & Tjørve, K.M.C. (2021). Mathematical expressions for the species–area relationship and the assumptions behind the models. In: Matthews, T.J., Triantis, K.A. & Whittaker, R.J. (Eds.), *The Species–Area Relationship: Theory and Application*. Cambridge University Press, pp. 157–184.
