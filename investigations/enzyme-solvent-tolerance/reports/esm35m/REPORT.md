# Frozen ESM2 35M benchmark

Date: 2026-09-11

## Result

Frozen mean-pooled ESM2 35M embeddings, supplied to the existing Extra Trees
predictor, do not show a reliable benefit over solvent and condition features
for choosing solvents at matching recorded conditions. All four paired 95 percent
intervals cross zero. This is one representation and predictor combination; it
does not establish that sequence information is intrinsically uninformative.

| Holdout | Conditions | Composition | ESM2 35M | ESM difference from conditions | 95 percent interval |
| --- | ---: | ---: | ---: | ---: | --- |
| Exact sequence | 0.286 | 0.284 | 0.260 | -0.026 | -0.073 to 0.021 |
| 40 percent identity components | 0.231 | 0.234 | 0.197 | -0.034 | -0.085 to 0.017 |
| 30 percent identity components | 0.220 | 0.211 | 0.222 | 0.001 | -0.041 to 0.043 |
| 30 percent identity plus publication components | 0.197 | 0.215 | 0.232 | 0.036 | -0.010 to 0.085 |

Values are sequence-averaged Spearman correlations for 772 matched-solvent series
across 275 eligible sequences. Composition and conditions columns come from the
saved benchmark. Identical row order, source data hash, sequence IDs, condition
features, outcomes, clusters and fold assignments were verified before training.
Bootstrap intervals resample entire grouping units in pairs, 5,000 times. They
are conditional on the fitted models and fold assignments and remain exploratory.

## Computation

Model: `facebook/esm2_t12_35M_UR50D`, revision
`6fbf070e65b0b7291e7bbcd451118c216cff79d8`. Only the approximately 136 MB safetensors
checkpoint and small configuration/tokenizer files were downloaded. Existing
PyTorch, Transformers and safetensors installations were reused.

Hardware: Apple M4 Max, 36 GB unified memory, PyTorch MPS backend. After model
loading and a CPU/GPU comparison, feature extraction for all 295 sequences took
11.17 seconds. This excludes download, model loading and downstream model fitting.
The CPU/GPU check on one 128-residue input differed by at most 9.54e-7 per pooled
feature. GPU memory and energy consumption were not profiled.

Each protein produces 480 final-layer features. The encoder is frozen, evaluated
without gradients, and processed one sequence chunk at a time. Special and padding
tokens are excluded from the mean. Five proteins exceed 1,022 residues; these use
nonoverlapping chunks whose vectors are averaged with residue-count weighting.
No residues are truncated, but context across chunk boundaries is absent.

The resulting compressed cache is approximately 530 KiB. Sequence hashes and model
revision are verified on cache reuse. All retained sequences were checked for
complete tokenization without unknown tokens, and all 295 vectors are finite and
distinct.

## Predictor and evaluation

The small predictor is Extra Trees with 250 trees, minimum leaf size 5, all
features considered per split and seed 11. It receives the seven existing numeric
condition features, solvent identity and all 480 ESM coordinates. Training predicts
log(1 + measured percentage). Median imputation and solvent encoding are fitted
on training rows only. No hyperparameter search, encoder fine-tuning or feature
selection was performed. The conditions and composition reference predictions
were reused without refitting or changing folds.

The unchanged cohort contains 5,798 rows, 295 sequences, 72 solvents and 272
publications. Similarity and publication groups and primary scoring follow the
[cluster benchmark](../family_benchmark/REPORT.md). The exact same source and
recorded-condition limitations apply, including heterogeneous solvent phases.

On the strict 30 percent split, percentage MAE is 45.65 with ESM versus 41.94 for
conditions only. Percentage-scale R2 is -0.148; log-scale R2 is 0.070. These do
not support useful absolute activity prediction.

There is a positive secondary result for ranking all conditions together in the
publication-linked split: 0.292 to 0.359, difference 0.067 with interval 0.031 to
0.108. Because this includes concentration changes and does not reliably persist
for matched-solvent choice, it is insufficient for the primary claim. It should
not be selected retrospectively as a replacement success metric.

## Artifact correction and validation

The earlier cluster export reused `sequence` as both an amino-acid column and an
exact-sequence grouping name. This overwrote the strings in that CSV with sequence
IDs. The original CSV, alignment FASTA, group assignments and prediction runs
were unaffected. For ESM extraction, sequences were rebuilt from the immutable
source and every string was checked against the alignment FASTA.

The export code now uses `sequence_group` for the exact-sequence grouping. The
cluster CSV was repaired with verified strings while preserving assignments;
`clusters_legacy_ids.csv` retains the old export. A regression test prevents this
collision. Sixteen tests pass across data parsing, grouping, matching, pooling,
chunk coverage and ID-based feature alignment.

## Reproduction

Run `PYTHONPATH=src python3 -m enzsolv.run_esm` from the project root. The command
requires the pinned checkpoint in the Hugging Face cache and uses local files
only. It reuses validated embeddings and recomputes the downstream predictions.

- `extraction.json`: model revision, runtime, device, pooling, chunked IDs and versions.
- `embeddings.npz`: cached vectors, sequence IDs and hashes.
- `source_sequences.csv`: verified original amino-acid strings used by ESM.
- `results.json`: exact metrics, intervals, settings and source checksum.
- `*_predictions.csv`: baseline, composition and ESM predictions with saved folds.
- `*_scores.csv`: per-sequence paired scores for conditions versus ESM.

The computational experiment is complete. It does not pass the gate for bench
validation. A larger encoder, surface features or a different predictor remain
untested; none is guaranteed to improve the result.
