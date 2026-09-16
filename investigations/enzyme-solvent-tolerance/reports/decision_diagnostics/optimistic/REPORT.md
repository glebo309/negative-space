# Training fit and deliberately leaky row-split diagnostic

Completed 2026-09-11. Descriptive only, with no transfer-valid uncertainty interval.
The same 5,675 covered rows, fixed 250-tree settings and seven condition features
were used. This explicitly tests the difference between fitting a familiar design
and predicting a new enzyme. No scientific conclusion relies on the leaky scores.

| Evaluation | Model | Log-activity R2 | Mean matched-solvent Spearman |
| --- | --- | ---: | ---: |
| Random rows | Conditions | 0.492 | 0.259 |
| Random rows | Surface | 0.581 | 0.101 |
| Random rows | ESM650 | 0.575 | 0.072 |
| Training fit | Conditions | 0.705 | 0.557 |
| Training fit | Surface | 0.768 | 0.667 |
| Training fit | ESM650 | 0.763 | 0.725 |

In each random-row fold, at least 99.65% of test rows have a sequence already in
training, and at least 99.65% have a training publication. Matched-series predictions
can come from different fold-trained models, unlike the valid grouped analyses.
These row-split rankings are therefore not a clean solvent-panel deployment test.

Protein representations improve the global log-target fit here, but that gain
does not translate into the intended solvent-choice ranking. Training-set ranking
can be much higher than held-out performance. Thus a pleasant training plot or
row-level R2 would misrepresent the present predictive utility. This diagnostic
does not isolate the biological mechanism or quantify assay noise.

All folds, predictions, overlap fractions and settings are saved alongside the report.
