# Decision-facing secondary metrics

Completed 2026-09-11 from saved out-of-fold predictions, without model fitting or
model selection on these metrics. Tests cover duplicate solvents, observed ties,
prediction ties, unscoreable panels and averaging order.

Pairwise accuracy counts correctly ordered solvent pairs with unequal observed
responses, giving prediction ties half credit. Top-choice rank is the normalized
observed midrank of the predicted best solvent. Best-choice hit rate counts
selection of an observed optimum; exactly tied top predictions receive their
uniform expected score. Panels are averaged within sequence, then sequences.

The first familiar-family composition result looked useful against its original
plain conditions comparator: pairwise accuracy 60.51% versus 56.79%, and best-choice
hit rate 29.03% versus 24.17%. These selected-assignment results do not survive the
stronger metadata-aware comparator and fixed repeated-split analysis.

## Fixed six-assignment familiar-family summary

Each sequence is averaged across its five or six supported assignments and counted
once. There are 80 sequences and four family/publication-connected uncertainty
components. Runs are not treated as independent experiments.

| Metric | Composition plus metadata/family | Metadata conditions | Difference, conditional 95% interval |
| --- | ---: | ---: | --- |
| Pairwise accuracy | 58.45% | 58.39% | +0.06 percentage points, -2.46 to +2.40 |
| Top-choice normalized rank | 0.6156 | 0.6050 | +0.0106, -0.0231 to +0.0383 |
| Best-choice hit rate | 25.73% | 25.02% | +0.71 percentage points, -5.96 to +4.50 |

No stable decision advantage is established. The small number of independent
components limits inference; these conditional intervals are not proof of
equivalence or of the absence of a useful biological signal.

For the metadata-qualified unfamiliar-family-plus-publication test, ESM650 also
has no reliable advantage over metadata conditions on any of these three metrics:
pairwise difference +0.42 percentage points, top-rank difference +0.0195, and
best-choice-hit difference +0.40 percentage points, all intervals crossing zero.

The isolated positive secondary comparisons remain saved, including the original
familiar split and alternate split against weaker comparators. They are not used
to replace the primary matched-solvent outcome after seeing results.

See `repeated_familiar/results.json` for the aggregate, and named subfolders for
individual source-hashed prediction analyses. The reusable commands are
`python3 -m enzsolv.run_decision_utility` and
`python3 -m enzsolv.run_repeated_utility` with `PYTHONPATH=src` from the project root.
