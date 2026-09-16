# Final reviewer-requested composition ablation

Completed 2026-09-12. This analysis resolves a specific comparator mismatch,
not a search over architectures or favorable splits.

The final reviewer noted that composition candidates in the repeated familiar-family
test still included family identity, even though that feature harmed the baseline.
The added comparison removes family identity from the candidate as well. Its inputs
are exactly the metadata-only baseline inputs plus 20 amino-acid frequencies.
All six saved assignments, supported test rows, full covered training cohort,
tree settings and aggregation rules are unchanged. Thirty additional fits completed.

| Assignment | Composition without family | Metadata conditions | Difference, conditional 95% interval |
| --- | ---: | ---: | --- |
| Original | 0.2397 | 0.1936 | +0.0461, -0.0084 to +0.1188 |
| Seed 47 | 0.2390 | 0.2388 | +0.0001, -0.0645 to +0.0796 |
| Seed 17 | 0.2160 | 0.1746 | +0.0414, -0.0006 to +0.0932 |
| Seed 23 | 0.1961 | 0.1731 | +0.0230, -0.0587 to +0.0988 |
| Seed 71 | 0.2314 | 0.1895 | +0.0419, +0.0161 to +0.0832 |
| Seed 101 | 0.2332 | 0.2566 | -0.0234, -0.0860 to +0.0867 |

Per-sequence averaging across available assignments gives composition 0.2263
versus metadata conditions 0.2021, a gain of 0.0242 with interval -0.0089 to 0.0824.
There are 80 sequences, eight families and four publication-connected uncertainty
components. Seventy-two sequences have six supported assignments and eight have
five. Repeated fits are not independent observations.

## Secondary decision utility

| Metric | Composition without family | Metadata conditions | Difference, conditional 95% interval |
| --- | ---: | ---: | --- |
| Pairwise ordering accuracy | 59.68% | 58.39% | +1.29 percentage points, +0.04 to +3.48 |
| Top-choice normalized rank | 0.6236 | 0.6050 | +0.0186, -0.0013 to +0.0536 |
| Best-solvent hit rate | 26.60% | 25.02% | +1.57 percentage points, -3.50 to +6.23 |

This ablation improves the case for a weak within-family composition signal.
It must not be summarized as every metric being null: the secondary pairwise
interval is just above zero. However, the primary matched-solvent correlation
and best-choice utility remain uncertain, and the small component count,
exploratory sequence and multiplicity preclude treating that secondary result
as a confirmed practical predictor. It does not establish a surface mechanism
or transfer to unfamiliar families.

The operational recommendation remains to withhold predictor-validation bench work
on this evidence. Preserve the narrower signal as a lead, not a demonstrated
absence of protein information. `analysis_plan.json`, `results.json`, per-family
scores, exact prediction tables and all signed fit caches are saved here.
