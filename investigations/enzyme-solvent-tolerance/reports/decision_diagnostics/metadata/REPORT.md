# Source-motivated metadata sensitivity

Completed 2026-09-11. Rules were fixed from the source audit before these scores.
Original labels and earlier benchmark outputs were not changed.

## Cohorts and predictors

Literal exclusion of final-comment assay/control uncertainty leaves 2,328 rows
before organic-reference exclusions. Also excluding explicit organic-reference
comments leaves 2,316 rows from 110 sequences. Requiring an explicit aqueous-phase
procedure leaves 2,307 rows from 109 sequences and 98 publications. Exact structural
coverage reduces this strict subset to 2,293 rows, 107 sequences and 97 publications.
Its matched metric scores 320 series from 97 sequences.

The full covered cohort remains 5,675 rows and 289 sequences, with 763 scoreable
series from 270 sequences. All comparisons use identical rows and original folds.

Seven literal metadata flags were added: non-incubated control, water-incubated
control, unclassified control, aqueous-phase procedure, organic-solvent procedure,
mode uncertainty and explicit organic reference. Flags do not recover quantitative
salt, buffer, solvent carry-over or enzyme preparation. No uncertainty flag is not
proof of an unambiguous original assay. `inclusion_ledger.csv` records every choice.

## Primary results

| Cohort and holdout | Conditions | Conditions plus metadata | Composition plus metadata | Surface plus metadata | ESM650 plus metadata |
| --- | ---: | ---: | ---: | ---: | ---: |
| Full, family plus publication | 0.1805 | 0.1910 | 0.1870 | 0.2067 | 0.2199 |
| Qualified, family plus publication | 0.2403 | 0.2587 | 0.2572 | 0.2690 | 0.2719 |
| Qualified, family | 0.2935 | 0.2864 | 0.2559 | 0.2601 | 0.2705 |

Values are mean matched-solvent Spearman correlations, not accuracies. Protein
models do not reliably outperform the metadata comparator in any of these jobs.
In the qualified family-plus-publication split, surface gain is 0.0103 with
conditional 95% interval -0.0598 to 0.0753; ESM gain is 0.0132 with interval
-0.0583 to 0.0791. Full-cohort ESM gain is 0.0289, interval -0.0150 to 0.0756.

Separate scoring of the same predictions for water-incubated and non-incubated
controls also yields no protein-feature improvement interval entirely above zero.
The qualified water-control stratum scores only 18 sequences. These are descriptive
strata, not separately trained models or independent replications.

The qualified subset has higher absolute scores, but changing the test population
prevents interpreting that difference as a causal estimate of annotation quality.
The justified conclusion is narrower: these explicit metadata qualifications did
not produce a reliable protein-feature advantage under the tested strict transfers.
Full-precision results, strata, score tables and signed fit caches are saved here.
