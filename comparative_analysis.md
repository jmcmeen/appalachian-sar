# Southern Appalachian Amphibians: Updated Comparative Analysis using NLS vs. Log-Linear SAR Fitting

A companion to [`analysis_updates.ipynb`](https://github.com/jmcmeen/appalachian-sar/blob/main/analysis_updates.ipynb), comparing the nonlinear least squares (NLS) results from the [`sars`](https://pypi.org/project/sars/) Python library against the original log-linear OLS calculations published in Stout, Jessee & McMeen (2025). Both the **nested model** (n = 7) and the **island model** (n = 26) are analyzed.

---

## Fitting Method Differences

The original calculations used log₁₀-transformed OLS regression:

> log₁₀(S) = z · log₁₀(A) + log₁₀(C)

This minimizes squared residuals in log-space and implicitly assumes multiplicative error — i.e., error proportional to the predicted value. The C and z parameters are recovered as C = 10^(intercept) and z = slope.

The `sars` library fits the power model directly in arithmetic space via nonlinear least squares (`scipy.optimize.least_squares`):

> S = cA^z

This minimizes squared residuals in the original species-count scale and assumes additive error. The NLS approach is generally preferred in modern SAR literature because it does not impose the variance structure that log-transformation assumes (Tjørve & Tjørve 2021).

## Parameter Estimates

### Nested Model (n = 7)

| Group | Method | c / C | z | R² |
|---|---|---|---|---|
| Amphibians | Original (log-linear) | 18.87 | 0.1222 | 0.9667 |
| Amphibians | `sars` NLS | 11.11 | 0.1779 | 0.9826 |
| Frogs | Original (log-linear) | 9.34 | 0.0598 | 0.9279 |
| Frogs | `sars` NLS | 8.22 | 0.0753 | 0.9281 |
| Salamanders | Original (log-linear) | 8.36 | 0.1734 | 0.9924 |
| Salamanders | `sars` NLS | 5.33 | 0.2192 | 0.9951 |

### Island Model (n = 26)

| Group | Method | c / C | z | R² |
|---|---|---|---|---|
| Amphibians | Original (log-linear) | 13.92 | 0.1079 | 0.7572 |
| Amphibians | `sars` NLS | 12.88 | 0.1272 | 0.7253 |
| Frogs | Original (log-linear) | 4.28 | 0.1338 | 0.6582 |
| Frogs | `sars` NLS | 4.98 | 0.1094 | 0.6549 |
| Salamanders | Original (log-linear) | 8.95 | 0.1043 | 0.5641 |
| Salamanders | `sars` NLS | 7.85 | 0.1383 | 0.5812 |

### Direction of Parameter Shifts

For the nested model, the expected pattern is clear: NLS produces lower c and higher z in all three groups. The effect is dramatic for amphibians (c drops from 18.87 to 11.11, a 41% decrease; z rises from 0.1222 to 0.1779, a 46% increase) and for salamanders (c drops 36%; z rises 26%). For frogs the shift is more modest, consistent with the already-strong log-linear fit.

For the island model, the shift direction is **not uniform**. Amphibians and salamanders follow the expected pattern (lower c, higher z), but frogs go the opposite way: NLS produces higher c (4.98 vs. 4.28) and lower z (0.1094 vs. 0.1338). This reversal is driven by the high scatter in the frog data and the influence of outlier sites — particularly the UT Arboretum (1.01 km², 9 frogs), which has unusually high frog richness for its size and exerts strong leverage in log-space.

The original R² values are unadjusted R² in log-space (the paper reports R² from `linregress`, not adjusted R²), while the `sars` values are unadjusted R² in arithmetic space. These are not directly comparable. In several cases, the NLS R² is slightly lower than the log-linear R² despite NLS producing a visually better fit — this is expected because NLS minimizes absolute residuals rather than proportional ones, and R² in arithmetic space penalizes large-area misfit more heavily.

## z-Value Interpretation

Canonical z-value ranges from Rosenzweig (1995) and Drakare et al. (2006) are 0.10–0.20 for mainland/nested designs and 0.25–0.35 for true island biogeography.

### Nested Model

| Group | z (Original) | z (NLS) | Range Context |
|---|---|---|---|
| Amphibians | 0.1222 | 0.1779 | Mainland/nested (0.10–0.20) |
| Frogs | 0.0598 | 0.0753 | Below canonical range |
| Salamanders | 0.1734 | 0.2192 | NLS pushes above 0.20 |

The nested model's z-values are broadly consistent with the mainland/nested range, as expected from the concentric sampling design (wetland pond ⊂ Steele Creek Park ⊂ Sullivan Co. ⊂ NE Tennessee ⊂ Eastern Tennessee ⊂ broader southern Appalachia). The frog z is notably low under both methods, reflecting the relatively flat frog richness across scales — frog diversity saturates quickly and adds few species even as area increases by orders of magnitude. Salamander z under NLS (0.2192) exceeds the canonical mainland ceiling of 0.20, approaching true-island territory. This likely reflects the strong elevational partitioning and high endemism among plethodontid salamanders, which makes even nested areas behave somewhat like ecological islands for this group.

### Island Model

| Group | z (Original) | z (NLS) | Range Context |
|---|---|---|---|
| Amphibians | 0.1079 | 0.1272 | Mainland/nested (0.10–0.20) |
| Frogs | 0.1338 | 0.1094 | Mainland/nested |
| Salamanders | 0.1043 | 0.1383 | Mainland/nested |

All island z-values fall within the mainland/nested range under both methods. This is somewhat surprising given that the island sites are independent localities (not concentric nested areas), which typically produce higher z-values. However, the island sites span a broad geographic and elevational range within a single ecoregion, and many share species pools. The result is consistent with the "mixed nested + island" design noted by the original authors.

## The Salamander Anomaly in the Island Model

The most analytically significant finding for the island dataset is the **weak power-model fit for salamanders** (NLS R² = 0.58). Examining the data reveals two patterns driving this:

1. **Outlier counties with suppressed richness.** Whitfield Co., GA (753.69 km², 8 salamanders) and Fannin Co., GA (1,015 km², 11 salamanders) have dramatically fewer salamanders than other areas of comparable size. Both are Georgia counties at the southern edge of the study area, where salamander diversity drops off due to physiographic and climatic transitions. These pull the power curve downward at intermediate scales.

2. **Exceptionally rich small sites.** Elk Knob SP (17.9 km², 13 salamanders) and Gorges SP (31.2 km², 14 salamanders) are montane parks with high salamander diversity relative to their size, reflecting the concentration of plethodontids in moist, high-elevation habitats.

Together, these produce a dataset that is poorly fit by any single monotonic curve. The linear model (R² = 0.6519) actually outperforms the power model (R² = 0.5812) for island salamanders by AIC and BIC, suggesting that the area–richness relationship for salamanders at this spatial scale is better approximated by a simple linear trend than by any concave power function.

This contrasts sharply with the nested model, where salamanders show the strongest power fit of any group (R² = 0.9951). The difference illustrates a key methodological point: the nested design, by using concentric areas that share species pools, eliminates between-site habitat variability and produces cleaner area–richness curves. The island model, using independent sites, exposes the ecological noise that habitat heterogeneity, edge effects, and geographic position introduce.

## Multi-Model Rankings

### AICc Considerations

With n = 7 (nested), the AICc correction term 2k(k+1)/(n−k−1) is large for 3-parameter models (k = 4 including σ) but still finite since n > k + 1. Results for 3-parameter models on the nested data should be interpreted cautiously — the additional parameter may produce artificially high R² with only 3 residual degrees of freedom.

With n = 26 (island), AICc is well-defined for all models, and model averaging and bootstrap confidence intervals are fully available.

### Nested Rankings

For all three nested groups, the power model provides an excellent fit. The 3-parameter powerR model slightly outperforms by BIC (amphibians: powerR R² = 0.9965 vs. power R² = 0.9826), but the improvement is marginal and comes at the cost of an additional parameter on only 7 data points. The 2-parameter power model is the most defensible choice.

### Island Rankings

For island **amphibians**, the linear model (R² = 0.7360, BIC = 164.86) and power model (R² = 0.7253, BIC = 165.89) perform nearly identically. The powerR model (R² = 0.7492) provides marginal improvement with a third parameter.

For island **frogs**, the power model dominates (R² = 0.6549, BIC = 115.43). The logarithmic model is a close second (R² = 0.6418, BIC = 116.40). The linear model performs poorly (R² = 0.5198).

For island **salamanders**, the linear model (R² = 0.6519, BIC = 157.83) outperforms the power model (R² = 0.5812, BIC = 162.65). This is the most striking departure from the standard SAR framework and is discussed in the salamander anomaly section above.

## Threshold Analysis

`sars.sar_threshold()` tests for breakpoints (small-island effect) by comparing continuous two-slope (ContOne), left-horizontal + right slope (ZslopeOne), and simple linear (no breakpoint) models.

### Nested

No breakpoints are detected for any group — the analysis selects the simple linear model in all cases. With n = 7, there is insufficient data to support piecewise structure.

### Island

For amphibians and salamanders, the threshold analysis selects the **ContOne** (continuous two-slope) model with a breakpoint at approximately log₁₀(area) ≈ 6.9 (corresponding to ~8,000,000 km², well outside the observed range). This is an extrapolation artifact and should not be interpreted as evidence of a genuine small-island effect. For frogs, the simple linear model is selected.

At n = 26, the island dataset is approaching the sample size where meaningful breakpoint detection becomes possible, but the current data do not support a small-island effect for any group.

## Model Averaging & Bootstrap CI

### Nested (n = 7)

With n = 7, model averaging via `sars.sar_average()` may produce Akaike weights that are dominated by a single model due to the large AICc corrections. Results should be treated as preliminary.

### Island (n = 26)

With n = 26, AICc-based model averaging and bootstrap confidence intervals are fully applicable. The notebook attempts both for each group. The Akaike weights indicate how much each of the 20 candidate models contributes to the information-theoretic ensemble prediction. Bootstrap CIs (500 replicates, 95% confidence) provide uncertainty bands around the model-averaged curve.

These tools represent a substantial analytical advance over the original log-linear approach, which produced only a single best-fit line with no formal uncertainty quantification.

## Nested vs. Island: Structural Comparison

The paper's central comparison — nested vs. island models — reveals consistent patterns under both the original and NLS methods:

1. **The nested model always produces higher c (baseline richness) and higher R² than the island model.** This is expected: the nested design eliminates between-site variability, and the largest nested areas encompass the species pools of all smaller areas.

2. **The z-values differ between models in group-specific ways.** For amphibians, the nested z is higher (0.1779 vs. 0.1272), indicating steeper scaling — the concentric areas capture species turnover more efficiently. For salamanders, nested z is also higher (0.2192 vs. 0.1383). For frogs, the pattern reverses: nested z (0.0753) is lower than island z (0.1094), driven by the high frog richness at the smallest nested site (a wetland pond with 7 frog species at 0.0018 km²), which anchors the low end of the curve and flattens the slope.

3. **The paper's recommendation to average the two models for practical prediction remains sound under NLS.** The nested model over-predicts (reflecting maximum regional diversity) while the island model under-predicts (reflecting observational gaps at individual sites). Their average provides a reasonable estimate at most spatial scales.

## Recommendations for Further Work

**Expand the island dataset.** Additional well-surveyed sites, particularly at intermediate spatial scales (50–500 km²), would strengthen the island model and potentially resolve the salamander anomaly. Even 5–10 additional sites would substantially improve multi-model discrimination.

**Expand the nested model.** The paper notes that the nested model is weakest at intermediate scales. Adding 2–3 additional nesting levels (e.g., individual watersheds, multi-county regions) would improve the model and reduce the AICc correction issue.

**Taxonomic decomposition.** Disaggregating frogs into anuran families (e.g., Hylidae, Ranidae) and salamanders into plethodontids vs. non-plethodontids could reveal whether the weak island fits are driven by specific clades.

**Integration with iNaturalist.** The `sars` library supports `sars.from_pyinaturalist()` for direct ingestion of citizen-science observation data, which could update species counts at all sites since the original data collection.

**GIS integration.** If site boundary shapefiles are available, `sars.from_geodataframe()` can ingest GeoDataFrames directly with automatic area calculation, streamlining the addition of new spatial scales.

## Summary of Key Findings

1. The NLS fitting method in `sars` produces systematically lower c and higher z values for the nested model, as expected from the different error assumptions. The direction of shift is not uniform for the island model, where the frog group reverses due to outlier leverage.

2. The qualitative conclusions of the original publication are fully supported: both nested and island models show that land area is a significant predictor of amphibian species richness in the southern Appalachians.

3. The nested model provides excellent power-law fits (R² = 0.93–0.99 across groups under NLS), confirming the strength of the concentric sampling design. The island model produces weaker but still significant fits (R² = 0.58–0.73).

4. Island salamander richness is better described by a linear model (R² = 0.65) than by the power model (R² = 0.58), driven by geographic outliers at the southern edge of the study area and exceptionally rich montane sites.

5. All z-values fall within or near the expected mainland/nested range (0.10–0.20), except for nested salamanders (NLS z = 0.22) which approach true-island values — consistent with the high endemism and elevational partitioning of plethodontid salamanders.

6. The island dataset (n = 26) supports full AICc-based model averaging and bootstrap confidence intervals, representing a significant analytical advance over the original log-linear approach.

7. No small-island effect (breakpoint) is statistically supported for either dataset at current sample sizes.

8. Averaging the nested and island model predictions remains a sound practical approach for estimating amphibian diversity at unstudied localities within the southern Appalachian ecoregion.

## References

- Drakare, S., Lennon, J.J. & Hillebrand, H. (2006). The imprint of the geographical, evolutionary and ecological context on species–area relationships. *Ecology Letters* 9:215–227.
- Rosenzweig, M.L. (1995). *Species Diversity in Space and Time*. Cambridge University Press.
- Stout, J.B., Jessee, L.D. & McMeen, J.N. (2025). Nested and island models for determining the species-area relationship of southern Appalachian amphibians. *Journal of North American Herpetology* 2025(1):1–7.
- Tjørve, E. & Tjørve, K.M.C. (2021). Mathematical expressions for the species–area relationship and the assumptions behind the models. In: Matthews, T.J., Triantis, K.A. & Whittaker, R.J. (Eds.), *The Species–Area Relationship: Theory and Application*. Cambridge University Press, pp. 157–184.
