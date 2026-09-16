# Direct ranking-target sensitivity

Completed 2026-09-11. Specified after the activity-regression results, so exploratory.

Training targets are within-experiment solvent midranks, normalized to [0,1].
Duplicate solvent rows are averaged first. Series with fewer than three solvents
or no response variation are excluded. Matched series cannot cross the strict
sequence-family-plus-publication folds. The target normalization uses no other
series. Fixed tree settings, cached representations and metadata flags are retained.
Evaluation still uses the original measured activity rankings.

The generic prediction helper applies a monotone `expm1` transform to its output.
That does not change the rank metric. These outputs are ranking scores, never
predicted activity percentages; no percentage-scale error is reported.

## Results and the essential comparator

| Cohort | Ranking-trained conditions | Composition | Surface | ESM650 | Earlier activity-trained metadata baseline |
| --- | ---: | ---: | ---: | ---: | ---: |
| Full covered | 0.0967 | 0.1854 | 0.1926 | 0.1887 | 0.1910 |
| Metadata-qualified | 0.1788 | 0.1966 | 0.2186 | 0.2358 | 0.2587 |

The full-cohort protein models beat their own ranking-trained comparator:
composition gain 0.0886, surface 0.0959, ESM 0.0920, all with conditional intervals
above zero. This is not a useful predictive advance over the already available
activity-trained metadata baseline. Against that stronger comparator, the gains
are -0.0056, 0.0016 and -0.0023, all with intervals crossing zero. The surface
comparison is 0.0016, interval -0.0538 to 0.0573.

In the qualified cohort, none beats either comparator reliably. Against the
activity-trained metadata baseline, composition is -0.0621, surface -0.0401 and
ESM -0.0229; all intervals cross zero. Surface also has no reliable gain over
composition in either cohort.

Conclusion: this bounded ranking-aware alternative does not rescue unfamiliar-family
solvent choice. The positive comparison against its weakened internal baseline is
retained transparently, not promoted into a claim of practical utility. The
alternative is a normalized-rank regression, not an exhaustive test of ranking
losses or neural architectures. Full-precision comparisons, including the cross-objective
paired scores on identical rows and folds, are saved in `results.json` and subfolders.
