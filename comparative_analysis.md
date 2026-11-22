# Southern Appalachian Amphibians: Updated Comparative Analysis Using NLS vs. Log-Linear SAR Fitting

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

| Group       | Method                | c / C | z      | R²     |
| ----------- | --------------------- | ----- | ------ | ------ |
| Amphibians  | Original (log-linear) | 18.87 | 0.1222 | 0.9667 |
| Amphibians  | `sars` NLS            | 11.11 | 0.1779 | 0.9826 |
| Frogs       | Original (log-linear) | 9.34  | 0.0598 | 0.9279 |
| Frogs       | `sars` NLS            | 8.22  | 0.0753 | 0.9281 |
| Salamanders | Original (log-linear) | 8.36  | 0.1734 | 0.9924 |
| Salamanders | `sars` NLS            | 5.33  | 0.2192 | 0.9951 |

### Island Model (n = 26)

| Group       | Method                | c / C | z      | R²     |
| ----------- | --------------------- | ----- | ------ | ------ |
| Amphibians  | Original (log-linear) | 13.92 | 0.1079 | 0.7572 |
| Amphibians  | `sars` NLS            | 12.88 | 0.1272 | 0.7253 |
| Frogs       | Original (log-linear) | 4.28  | 0.1338 | 0.6582 |
| Frogs       | `sars` NLS            | 4.98  | 0.1094 | 0.6549 |
| Salamanders | Original (log-linear) | 8.95  | 0.1043 | 0.5641 |
| Salamanders | `sars` NLS            | 7.85  | 0.1383 | 0.5812 |

### Direction of Parameter Shifts

For the nested model, the expected pattern is clear: NLS produces lower c and higher z in all three groups. The effect is dramatic for amphibians (c drops from 18.87 to 11.11, a 41% decrease; z rises from 0.1222 to 0.1779, a 46% increase) and for salamanders (c drops 36%; z rises 26%). For frogs the shift is more modest, consistent with the already-strong log-linear fit.

For the island model, the shift direction is **not uniform**. Amphibians and salamanders follow the expected pattern (lower c, higher z), but frogs go the opposite way: NLS produces higher c (4.98 vs. 4.28) and lower z (0.1094 vs. 0.1338). This reversal is driven by the high scatter in the frog data and the influence of outlier sites — particularly the UT Arboretum (1.01 km², 9 frogs), which has unusually high frog richness for its size and exerts strong leverage in log-space.

The original R² values are unadjusted R² in log-space (the paper reports R² from `linregress`, not adjusted R²), while the `sars` values are unadjusted R² in arithmetic space. These are not directly comparable. In several cases, the NLS R² is slightly lower than the log-linear R² despite NLS producing a visually better fit — this is expected because NLS minimizes absolute residuals rather than proportional ones, and R² in arithmetic space penalizes large-area misfit more heavily.

## z-Value Interpretation

Canonical z-value ranges from Rosenzweig (1995) and Drakare et al. (2006) are 0.10–0.20 for mainland/nested designs and 0.25–0.35 for true island biogeography.

### Nested Model

| Group       | z (Original) | z (NLS) | Range Context               |
| ----------- | ------------ | ------- | --------------------------- |
| Amphibians  | 0.1222       | 0.1779  | Mainland/nested (0.10–0.20) |
| Frogs       | 0.0598       | 0.0753  | Below canonical range       |
| Salamanders | 0.1734       | 0.2192  | NLS pushes above 0.20       |

The nested model's z-values are broadly consistent with the mainland/nested range, as expected from the concentric sampling design (wetland pond ⊂ Steele Creek Park ⊂ Sullivan Co. ⊂ NE Tennessee ⊂ Eastern Tennessee ⊂ broader southern Appalachia). The frog z is notably low under both methods, reflecting the relatively flat frog richness across scales — frog diversity saturates quickly and adds few species even as area increases by orders of magnitude. Salamander z under NLS (0.2192) exceeds the canonical mainland ceiling of 0.20, approaching true-island territory. This likely reflects the strong elevational partitioning and high endemism among plethodontid salamanders, which makes even nested areas behave somewhat like ecological islands for this group.

### Island Model

| Group       | z (Original) | z (NLS) | Range Context               |
| ----------- | ------------ | ------- | --------------------------- |
| Amphibians  | 0.1079       | 0.1272  | Mainland/nested (0.10–0.20) |
| Frogs       | 0.1338       | 0.1094  | Mainland/nested             |
| Salamanders | 0.1043       | 0.1383  | Mainland/nested             |

All island z-values fall within the mainland/nested range under both methods. This is somewhat surprising given that the island sites are independent localities (not concentric nested areas), which typically produce higher z-values. However, the island sites span a broad geographic and elevational range within a single ecoregion, and many share species pools. The result is consistent with the "mixed nested + island" design noted by the original authors.

## The Salamander Anomaly in the Island Model

The most analytically significant finding for the island dataset is the **weak power-model fit for salamanders** (NLS R² = 0.58). Examining the data reveals two patterns driving this:

1. **Outlier counties with suppressed richness.** Whitfield Co., GA (753.69 km², 8 salamanders) and Fannin Co., GA (1,015 km², 11 salamanders) have dramatically fewer salamanders than other areas of comparable size. Both are Georgia counties at the southern edge of the study area, where salamander diversity drops off due to physiographic and climatic transitions. These pull the power curve downward at intermediate scales.

2. **Exceptionally rich small sites.** Elk Knob SP (17.9 km², 13 salamanders) and Gorges SP (31.2 km², 14 salamanders) are montane parks with high salamander diversity relative to their size, reflecting the concentration of plethodontids in moist, high-elevation habitats.

Together, these produce a dataset that is poorly fit by any single monotonic curve. The linear model (R² = 0.6519) actually outperforms the power model (R² = 0.5812) for island salamanders by AIC and BIC, suggesting that the area–richness relationship for salamanders at this spatial scale is better approximated by a simple linear trend than by any concave power function.

This contrasts sharply with the nested model, where salamanders show the strongest power fit of any group (R² = 0.9951). The difference illustrates a key methodological point: the nested design, by using concentric areas that share species pools, eliminates between-site habitat variability and produces cleaner area–richness curves. The island model, using independent sites, exposes the ecological noise that habitat heterogeneity, edge effects, and geographic position introduce.

## Multi-Model Rankings

### AICc Considerations

With n = 7 (nested), the AICc correction term 2k(k+1)/(n−k−1) is large for 3-parameter models (k = 4 including σ) but still finite since n > k + 1. Results for 3-parameter models on the nested data should be interpreted cautiously — the additional parameter may produce artificially high R² with only 4 residual degrees of freedom (7 points less 3 parameters; σ is estimated from the residuals).

With n = 26 (island), AICc is well-defined for all models, and model averaging is fully available; the bootstrap intervals come with a caveat about independence (see below).

### Nested Rankings

For all three nested groups, the power model provides an excellent fit (R² = 0.93–0.995). Which model ranks first depends on the criterion. By AICc, the power model leads for frogs and salamanders, and the three-parameter epm1 for amphibians (R² = 0.998 against the power model's 0.983). By BIC, whose penalty for an extra parameter is much smaller than AICc's at n = 7 (about 2 against about 14), epm1 leads for all three groups, by a wide margin. The choice between them is a choice of penalty: on 7 points, a third parameter raises R² only slightly (to 0.988–0.999), and AICc's small-sample correction, designed for exactly this case, judges that not worth it except for amphibians. The 2-parameter power model is the more defensible choice for frogs and salamanders, and the evidence for amphibians is split.

### Island Rankings

For island **amphibians**, the linear model (R² = 0.7360, BIC = 164.86) and power model (R² = 0.7253, BIC = 165.89) perform nearly identically. The powerR model (R² = 0.7492) provides marginal improvement with a third parameter.

For island **frogs**, the power model dominates (R² = 0.6549, BIC = 115.43). The logarithmic model is a close second (R² = 0.6418, BIC = 116.40). The linear model performs poorly (R² = 0.5198).

For island **salamanders**, the linear model (R² = 0.6519, BIC = 157.83) outperforms the power model (R² = 0.5812, BIC = 162.65). This is the most striking departure from the standard SAR framework and is discussed in the salamander anomaly section above.

## Threshold Analysis

`sars.sar_threshold()` tests for breakpoints (small-island effect) by comparing continuous two-slope (ContOne), left-horizontal + right slope (ZslopeOne), and simple linear (no breakpoint) models.

### Nested

The analysis selects the simple linear model for every group, by AICc. That choice comes from AICc's small-sample correction rather than from the fits: by AIC and BIC the one-breakpoint ContOne model is far ahead (for amphibians, AIC 42.2 against 63.4 for the line), but with n = 7 the correction for its extra parameters is large enough to reverse the order (AICc 102.2 against 71.4). Seven concentric areas are too few to test piecewise structure, and the nested design isn't suited to it: no breakpoint is detected, which is not evidence that none exists.

### Island

For amphibians and salamanders, the threshold analysis selects the **ContOne** (continuous two-slope) model, with the breakpoint at ln(area) ≈ 6.92 for both groups. `sar_threshold()` fits on the natural log of area, so that is about 1,000 km², inside the observed range (0.008–2,114 km²): six of the 26 sites lie above it (Fannin Co. GA, Polk Co. TN, Hamilton Co. TN, Sevier Co. TN, Buncombe Co. NC and Great Smokies NP). ContOne beats the linear model by about 11 AICc units for both groups (157.3 vs 168.9 for amphibians, 153.3 vs 164.3 for salamanders), so the breakpoint has real support.

It is not a small-island effect. The slope *steepens* above the breakpoint, from 1.7 to 24 amphibian species per unit of ln(area) and from 1.0 to 21 for salamanders, whereas a small-island effect flattens the curve below a threshold; that is what the ZslopeOne model tests for, and it is not selected for any group. The steep upper segment rests on six county- and park-scale sites. Five of them are the five richest sites for amphibians as a whole (for salamanders, four of the five richest; Hamilton Co. TN, with 19, ranks below Unicoi Co. TN, Watauga Co. NC and New River Gorge NP WV, all three with 21), and the sixth, Fannin Co. GA, is the suppressed southern-edge county noted under the salamander anomaly, so the segment is better read as the signature of those sites than as a general property of the relationship. For frogs, the simple linear model is selected.

With n = 26 the island dataset can support a one-breakpoint model, as it does here; what it does not show, for any group, is a small-island effect.

## Model Averaging & Bootstrap CI

### Nested (n = 7)

Not every group's weights are dominated by a single model, and the weights depend heavily on the criterion, because AICc's small-sample correction is large at n = 7. The three largest weights under each:

| Group | Akaike (AICc) weights | BIC (Schwarz) weights |
| --- | --- | --- |
| Amphibians | epm1 0.62, power 0.30, powerR 0.07 | epm1 0.89, powerR 0.11, p2 0.002 |
| Frogs | power 0.51, epm1 0.27, powerR 0.16 | epm1 0.61, powerR 0.36, p2 0.01 |
| Salamanders | power 0.91, epm1 0.07, powerR 0.01 | epm1 0.74, powerR 0.13, epm2 0.05 |

The power model's BIC weight is 0.001–0.02. Two cases come close to a single model, salamanders under AICc (power 0.91) and amphibians under BIC (epm1 0.89), and they are different models. Where AICc and BIC disagree this sharply, the averaged curve depends on how hard extra parameters are penalised, so the nested averages should be treated as preliminary.

### Island (n = 26)

With n = 26 the two criteria broadly agree, and the weights are spread across many models (the top model's AICc weight is 0.20–0.34: linear for amphibians and salamanders, power for frogs). Model averaging is fully applicable here, and the notebook shows the averaged curve for each group.

The bootstrap intervals (500 resamples, 95% confidence) need a caveat. `sars.bootstrap_ci()` resamples sites with replacement, treating them as independent draws, and they aren't fully: some sites lie within others (UT Arboretum within Anderson Co. TN, Elk Knob SP largely within Watauga Co. NC, part of Great Smokies NP within Sevier Co. TN), and the sites draw on overlapping regional species pools. The intervals are therefore likely too narrow and should be read as rough. (The nested areas, concentric by design, can't be bootstrapped in the usual sense at all, and the notebook doesn't try.) A design-aware method, resampling by region for instance, would be needed for intervals with their usual meaning. Some resamples were dropped because their model average couldn't be computed: 3 of 500 for amphibians, 49 for frogs and 4 for salamanders (the notebook reports the counts). Those are not a random tenth of the frog resamples but the ones whose fits failed, which may narrow the frog interval further.

Model averaging and these intervals go beyond the original log-linear approach, which produced a single best-fit line with no formal uncertainty quantification, but the intervals' assumptions limit how far they can be taken.

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

6. The island dataset (n = 26) supports AICc-based model averaging; its bootstrap intervals assume independent sites, which the data only partly satisfy, so they are likely too narrow. On the nested data (n = 7), the averaging weights change sharply between AICc and BIC.

7. No small-island effect is supported for either dataset. The island amphibian and salamander fits do support a breakpoint near 1,000 km², above which richness climbs steeply across the six largest, county- and park-scale sites: a large-area steepening, not a small-island flattening.

8. Averaging the nested and island model predictions remains a sound practical approach for estimating amphibian diversity at unstudied localities within the southern Appalachian ecoregion.

## References

- Drakare, S., Lennon, J.J. & Hillebrand, H. (2006). The imprint of the geographical, evolutionary and ecological context on species–area relationships. *Ecology Letters* 9:215–227.
- Rosenzweig, M.L. (1995). *Species Diversity in Space and Time*. Cambridge University Press.
- Stout, J.B., Jessee, L.D. & McMeen, J.N. (2025). Nested and island models for determining the species-area relationship of southern Appalachian amphibians. *Journal of North American Herpetology* 2025(1):1–7.
- Tjørve, E. & Tjørve, K.M.C. (2021). Mathematical expressions for the species–area relationship and the assumptions behind the models. In: Matthews, T.J., Triantis, K.A. & Whittaker, R.J. (Eds.), *The Species–Area Relationship: Theory and Application*. Cambridge University Press, pp. 157–184.
