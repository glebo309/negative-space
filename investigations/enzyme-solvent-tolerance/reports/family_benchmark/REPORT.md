# Sequence cluster transfer benchmark

Date: 2026-09-11

## Finding

Whole-sequence amino-acid composition has no reliable incremental benefit for
ranking solvents at matching recorded conditions in this benchmark. This weakens
the initial composition result; it does not test or rule out protein language
models or explicit surface structure.

The primary comparison below uses the same data, transformations, folds and model
settings for both feature sets. Mean Spearman correlations first average over
eligible experimental series within each sequence and then over sequences.

| Held-out unit | Groups | Conditions only | Plus composition | Difference | Paired 95% bootstrap interval |
| --- | ---: | ---: | ---: | ---: | --- |
| Exact sequence | 295 | 0.286 | 0.284 | -0.002 | -0.044 to 0.039 |
| 40% identity component | 193 | 0.231 | 0.234 | 0.004 | -0.042 to 0.046 |
| 30% identity component | 160 | 0.220 | 0.211 | -0.010 | -0.048 to 0.031 |
| 30% identity component plus publication | 144 | 0.197 | 0.215 | 0.019 | -0.029 to 0.069 |

All four intervals include zero. There are 772 eligible matched series across
275 sequences. Group counts in the table refer to the full training cohort;
the bootstrap uses the groups represented among those 275 ranked sequences.

## Cohort and cleaning

Start from the previous 6,433-row primary cohort. Retain only uppercase standard
amino-acid sequences without unknown characters (6,357 rows), then measurements
with a numeric aqueous control of exactly 100, nonnegative outcome, solvent
concentration greater than zero and no more than 100, a DOI and solvent name
(5,801 rows). Remove three identical records, comparing all original fields except
measurement identifier. The resulting cohort has 5,798 rows, 295 unique sequences,
72 solvents and 272 publications.

Original source sequences and the downloaded CSV remain unchanged. Excluding
nonstandard strings is conservative and may omit recoverable formatting cases.
The control-100 restriction avoids assuming how differently normalized rows
should be combined. It does not certify that all source assay protocols are
equivalent. No row is excluded because of model error or a large positive outcome.

Incubation and assay temperatures and pH are parsed separately. Ranges, room
temperature and values with unspecified stages remain missing. The old parser
copied unlabelled values to both stages; this benchmark uses its own conservative
parser. A numeric incubation duration is extracted from the procedure field for
99.48% of retained rows. Units are converted to hours; missing values are imputed
using training-fold medians and represented by missingness indicators. Eight
random distinct procedure/time pairs were manually checked against their strings;
this is not a full source-paper audit.

Outcomes above 100% are retained as potential activation. Training predicts
log(1 + measured percentage) to reduce domination by extreme outcomes. Predictions
are transformed back for percentage-scale errors. This is different from the
initial raw-target model, so changes from the initial report cannot be attributed
solely to removal of homologues.

## Alignment and split definition

MMseqs2 all-against-all search uses no prefilter, explicit identity calculation,
minimum 30% sequence identity over alignment length, at least 80% coverage of both
query and target, and E-value at most 0.001. All 295 sequences have a self hit;
the resulting table contains 1,078 directed hits including self hits.

Connected components of the qualifying similarity edges are assigned to five
folds using GroupKFold. Connected components prevent a detected qualifying pair
from crossing folds, including through transitive links. In the publication
sensitivity analysis, proteins reported in the same DOI are joined before folds
are assigned. Assertions verify no qualifying edge crosses folds and no DOI
crosses folds in that analysis.

These are operational sequence similarity groups, not curated biological
families. Remote homologues and shared domains failing the coverage threshold may
remain across folds. No fold or structural-family separation is claimed.

## Models and metrics

Both models use 250 Extra Trees, minimum leaf size 5, all features available at
each split, seed 11 and four CPU threads. Baseline predictors are solvent identity,
concentration, solvent logP, incubation and assay temperature, incubation and assay
pH, and incubation time. The composition model adds twenty residue frequencies.
Preprocessing is fitted within each training fold. There is no hyperparameter
search and no use of DOI, sequence identifier or group label as predictors.

For matched-solvent ranking, rows share sequence, DOI, solvent concentration,
procedure, recorded temperature, pH, assay solution, substrates, substrate
concentrations, cofactor, shaking and comments. Replicated values within one
solvent are averaged for scoring. A series requires at least three solvents and
variable measured outcomes. Constant predictions score zero ranking information.
Series scores are averaged within sequence. Matching concerns recorded fields;
unreported experimental differences can still remain.

Uncertainty uses 5,000 paired bootstrap samples of entire held-out grouping units,
with sequence-weighted means. Intervals are conditional on these fitted models
and this fold assignment. They do not include retraining or fold-selection
uncertainty and are exploratory, not confirmatory hypothesis tests.

## Secondary results

For ranking all conditions within each sequence, including concentration changes,
the 30% split gives 0.338 for conditions and 0.331 with composition. The split
joining publications gives 0.292 versus 0.338, difference 0.046 (interval 0.004
to 0.089). This isolated positive result does not persist in the matched-solvent
metric. Composition may encode some general response differences without
reliably improving solvent choice.

Percentage-scale R2 remains negative. On the 30% split, percentage MAE is 41.94
for conditions and 44.07 with composition; log-scale R2 is 0.024 versus 0.088.
The ranking and transformed regression metrics answer different questions and
should not be substituted for one another to claim success.

## Interpretation and next experiment

The original 0.325 to 0.378 result is a preliminary observation, not established
transferable protein signal. The present analysis provides no convincing reason
to start prospective bench validation yet. A bounded next computational test
would compare a frozen sequence embedding with these same saved groups and
metrics. Whole-sequence composition omits residue order and surface exposure,
so its failure alone is not a scientific reason to rule out those representations.

Dataset heterogeneity remains substantial, including water-miscible and
water-immiscible solvent systems. A solvent-class restriction and source checks
remain needed before a publication claim. These runs are exploratory development
results; a final confirmatory benchmark needs a frozen evaluation plan.

## Reproduction and artifacts

From the project root run `PYTHONPATH=src python3 -m enzsolv.run_benchmark`.

- `results.json`: exact results, data checksum, installed package versions and alignment command.
- `clusters.csv`: sequence identifiers, original sequences and group assignments.
- `all_pairs.tsv`: alignment evidence for group construction.
- `*_predictions.csv`: every out-of-fold prediction and its fold assignment.
- `*_scores.csv`: paired per-sequence scores for both ranking metrics.

Twelve tests pass, including transitive grouping, conservative condition parsing,
duration units, cohort exclusions, paired cluster bootstrap and condition-matched
ranking. Output checks verified complete prediction coverage and group separation.

MMseqs2 was downloaded as an approximately 21 MB archive using the official
[installation instructions](https://github.com/soedinglab/MMseqs2#installation).
No pretrained protein model, structure database or GPU download was required.
