# Surface features and transferable solvent ranking

Completed 2026-09-11. All planned fits and supplementary checks have finished.

## Main surface result

Global surface descriptors do not reliably improve solvent ranking over the
conditions-only baseline in the current structure-covered cohort. All 18
surface-versus-conditions intervals cross zero across the main analysis,
exposure/confidence sensitivities and additional random seeds.

The primary comparison uses the same 5,675 measurements from 289 proteins for
every model. The ranking metric is evaluated on 763 matched-solvent series from
270 proteins. Values below are mean within-protein Spearman correlations, not
fractions of correct predictions.

| Holdout | Conditions only | Whole-sequence composition | Surface | Surface improvement, 95% interval |
| --- | ---: | ---: | ---: | --- |
| 30% identity plus publication | 0.181 | 0.214 | 0.195 | 0.014, -0.030 to 0.063 |
| 30% identity | 0.212 | 0.229 | 0.205 | -0.007, -0.054 to 0.040 |
| 40% identity | 0.215 | 0.226 | 0.206 | -0.009, -0.054 to 0.036 |
| Exact sequence | 0.284 | 0.267 | 0.244 | -0.041, -0.083 to 0.001 |

Surface also fails to reliably outperform whole-sequence composition in any of
these four holdouts. Comparisons are based on paired predictions on identical
rows, not on previously reported scores from the larger full cohort.

## Sensitivities

Surface minus conditions; none of these intervals is entirely above zero.

| Variant | Publication-linked 30% | 30% identity | 40% identity | Exact sequence |
| --- | ---: | ---: | ---: | ---: |
| Primary, RSA at least 0.20 | 0.014 | -0.007 | -0.009 | -0.041 |
| RSA at least 0.10 | 0.027 | 0.006 | -0.005 | -0.041 |
| RSA at least 0.30 | 0.021 | 0.005 | -0.017 | -0.038 |
| Confidence-restricted | 0.013 | -0.021 | -0.005 | -0.013 |

Confidence restriction keeps 284 proteins with mean pLDDT at least 70, covering
5,597 measurements, and computes summaries using only residues with pLDDT at
least 70. Physical accessibility still uses the full provided geometry. Its
baseline and composition comparator are refitted on that same subset.

Additional seeds 23 and 47 in the primary publication-linked split give surface
improvements of 0.029 (-0.012 to 0.074) and 0.023 (-0.027 to 0.071), respectively.
The main seed 11 gives 0.014 (-0.030 to 0.063). These are seed checks on the same
folds, not independent replication or repeated outer split selection.

## Sequence and surface combination

All four ESM650 and combined-model comparisons completed on the same cohort.

| Holdout | Conditions only | ESM650 | ESM650 plus surface | Surface addition, paired 95% interval |
| --- | ---: | ---: | ---: | --- |
| 30% identity plus publication | 0.181 | 0.199 | 0.204 | 0.006, -0.012 to 0.022 |
| 30% identity | 0.212 | 0.221 | 0.208 | -0.013, -0.033 to 0.005 |
| 40% identity | 0.215 | 0.210 | 0.193 | -0.018, -0.041 to 0.004 |
| Exact sequence | 0.284 | 0.251 | 0.240 | -0.011, -0.034 to 0.011 |

Neither ESM650 nor the combination reliably beats conditions alone in these
four primary comparisons. The final column tests the combination against ESM650,
not against conditions. It provides no reliable evidence of added surface value.

## Secondary outcomes

Pooled ranking across all conditions, including concentration variation, does
improve in some comparisons. In the publication-linked split, improvements over
conditions are 0.056 for composition (0.020 to 0.094), 0.067 for ESM650 (0.032 to
0.105), and 0.069 for ESM650 plus surface (0.035 to 0.106). These secondary results
do not establish improved choice among solvents at matching conditions. Surface
alone does not reliably improve even this pooled metric in the four main splits.

Conditions-only also has lower percentage-point MAE than every protein-aware
model in all four main holdouts. Full-precision absolute-error and R2 results are
retained in `results.json`, not replaced by the more favorable secondary ranking.

## Second predictor family

A supplementary fixed-alpha ridge model explicitly represents protein-by-solvent
and protein-by-condition interactions. This was specified after seeing the main
surface-tree results and is reported as exploratory. It also fails to establish
a surface gain in all four holdouts. The baseline includes solvent-by-condition
interactions, so candidates are not rewarded merely for adding that capability.

ESM features are reduced to 16 PCs using only unique training proteins. No
held-out alpha or dimension selection is performed. Full results and targeted
numerical audits are in `../surface_interactions/REPORT.md`.

## Methods and verification

The encoder, target definition and cohort are unchanged. Every predictor uses
the seven existing condition features plus solvent identity. Extra Trees uses
250 trees, minimum leaf size five, all features considered at each split, and
seed 11 unless stated. Models fit log(1 + measured percent); activity above
100 percent is retained. Imputation and solvent categories are fitted on training
folds only. Protein features are 20 whole-sequence frequencies, 27 surface
descriptors or 1,280 frozen ESM650 features, depending on the model.

The same five-fold assignments are inherited from the original 295-sequence
cohort after selecting structure-covered rows. Components and complete sequences
remain separated across folds. In the publication-linked split, publications
also remain separated. Components are operational identity/coverage groups,
not guaranteed biological enzyme families.

The primary score ranks different solvents at matched recorded conditions,
averages series within each protein, then averages proteins. At least three
solvents and variable outcomes are required. A constant predicted ranking scores
zero. Paired 95% intervals use 5,000 component-bootstrap draws, conditional on
the existing fits. They are not multiplicity-adjusted and do not capture
uncertainty from retraining or choosing different splits. Do not select the most
favorable cutoff after inspecting these results.

The complete run contains 62 forest model configurations and 16 supplementary
ridge configurations, each with five held-out folds. Repeated controls across
cutoffs are included in these counts; these are not 78 independent scientific
hypotheses or replications.

Tests cover cohort alignment, duplicate-ID rejection, feature selection,
training-only imputation, group/publication separation, prediction caching,
pairwise scoring, interaction capability and training-only scaling. Thirty-three
tests pass. Checks of saved artifacts find finite predictions and correct group
separation. Repeated baseline/composition controls across exposure thresholds
agree within 4.55e-13 percentage points. Unseen-solvent rows are below 1% in
every holdout (39 to 53 of 5,675 rows).

## Interpretation and decision

There is currently no positive computational gate for prospective bench
validation of these models. This does not establish that surface chemistry or
sequence is biologically unrelated to solvent tolerance. Global descriptors,
mean-pooled sequence features, unverified solution assemblies and heterogeneous
literature measurements limit what this benchmark tests. Modest benefits remain
compatible with some intervals.

The defensible output is a reproducible feasibility benchmark with explicit
negative/inconclusive results, not a validated solvent-tolerance predictor. A
paper claim would need a justified narrower hypothesis, stronger data or another
genuinely supported result; this analysis alone does not guarantee publication.

## Reproduce

```bash
PYTHONPATH=src python3 -m enzsolv.run_surface_benchmark --phase all
PYTHONPATH=src python3 -m enzsolv.run_interactions
python3 -m pytest -q
```

`analysis_plan.json` records the cohort and choices. `results.json` stores full
precision. Each variant/split/seed directory contains out-of-fold predictions,
paired scores and input-signature-checked model prediction caches. Completed
folds resume without refitting. No new model or structure download is required.
