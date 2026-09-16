# Independent evaluation review

Date: 2026-09-11

Scope: read-only review of the existing evaluation, synthetic controls and proposed easier tasks. Only this report was created by the reviewer. No raw data, model code, existing predictions or vault files were changed.

## Bottom line

No material implementation bug, target leakage or row/fold misalignment was found in the inspected core benchmark. Independent calculations reproduce the primary surface comparison across all four splits, including its paired cluster bootstrap intervals. The planted scalar signal is recovered strongly, whereas null and shuffled-feature checks show no clear benefit. These observations support interpreting the existing result as an actual limitation of the tested predictive methods on this dataset, not as a demonstrated absence of sequence-dependent solvent tolerance.

Remaining concerns primarily restrict the scientific interpretation: heterogeneous measurement context, sparse independent repeats, limited family support, a simple positive control rather than comprehensive model power analysis, and conditional rather than full training-procedure uncertainty.

## Code and artifacts inspected

- `src/enzsolv/benchmark.py` and `run_benchmark.py`.
- `surface_benchmark.py` and `run_surface_benchmark.py`.
- `interaction_model.py` and `run_interactions.py`.
- `decision_diagnostics.py` and `run_diagnostic_controls.py`.
- Corresponding benchmark, surface, interaction and diagnostic tests.
- Primary surface predictions, results and original filtered cohort.
- Synthetic-control results and coverage summaries.
- The canonical decision plan in the vault.

The review did not independently rerun every model fit or independently verify every cached structural descriptor. It verifies saved predictions against source IDs/targets and independently recalculates scoring. Earlier structural and numerical audits remain separate evidence.

## Independently verified calculations

For each primary surface split, the reviewer reconstructed the filtered cohort, joined source values by measurement ID, verified unique IDs and matching sequence IDs and targets, checked whole-sequence and cluster fold separation, and checked publication separation for the strict publication split. Matched-series Spearman correlations were recalculated directly with SciPy, then averaged within each sequence. The cluster bootstrap was independently implemented as 5,000 sampled cluster-index sets with seed 11.

All four files contain 5,675 rows. Each comparison scores 270 sequences and 763 matched series.

| Split | Surface minus conditions mean Spearman | Independently reproduced 95% interval |
| --- | ---: | --- |
| Identity 30% plus publication | 0.0143825937 | [-0.0301672039, 0.0626152754] |
| Identity 30% | -0.0070873386 | [-0.0535899018, 0.0402408677] |
| Identity 40% | -0.0093772068 | [-0.0539506782, 0.0355031628] |
| Sequence | -0.0407865747 | [-0.0832907215, 0.0006510287] |

Differences from stored metrics are below numerical tolerances of 1e-12. Source raw SHA-256 remains `400823131c996a3899d95cdb0c82e37a6aa5833d32f0daa23e15a696045b3a06`.

## Implementation findings

### No material defect found

- Imputation and solvent categories are fitted inside training folds. The interaction model also fits context scaling, unique-training-protein scaling and ESM PCA on training data only.
- Numeric and protein feature selection are explicit. Target values, fold IDs, sequence IDs and publication IDs are not accidentally exposed as predictors.
- Surface joins enforce one descriptor row per sequence, and saved embedding IDs and sequence hashes are checked before use.
- The strict family-plus-publication grouping takes connected components over both relations, preventing publication leakage across folds. The weaker family and sequence splits intentionally do not guarantee publication independence.
- The paired metric compares candidates on identical series. Its simultaneous reassignment of comparison columns avoids overwriting the composition candidate.
- Averaging solvents within a matched series, series within a sequence, and then sequences avoids giving a heavily measured protein proportionally greater weight in the primary score. Model fitting itself remains row-weighted.
- The cache signatures bind row identities, target values, selected features, fold assignments and relevant implementation/settings. Core cache use does not silently reuse different targets.

### Limitations, not discovered implementation bugs

1. **Conditional confidence intervals.** Cluster resampling accounts for the specified test-unit dependence, conditional on the saved fits. It does not refit the full learning procedure, vary fold construction or include source curation uncertainty. Overlapping training sets also mean it is not a general confidence guarantee for the entire model-development procedure. The existing reports appropriately call these conditional intervals.
2. **Exploratory multiplicity.** Many models, splits and sensitivities have been examined. A selected positive interval would require a fresh validation strategy. Conversely, consistently small gains remain useful diagnostic evidence without claiming formal equivalence or biological absence.
3. **Endpoint and metadata adequacy.** The seven condition predictors omit important experimentally meaningful distinctions found by the source auditor. Correct code cannot recover information that is absent or inconsistently normalized. Exact matching on recorded strings also cannot prove identical experiments.
4. **Training objective differs from primary evaluation.** Trees and ridge optimize row-level log-target errors, while the main scientific endpoint is equally weighted per-protein solvent ranking. This is a defensible fixed baseline, not proof that ranking-aware training could never do better.
5. **Representation limitations.** Global surface summaries and mean-pooled embeddings are not exhaustive descriptions of local active-site, oligomeric, dynamic or solvent-specific effects. The supplementary 16-component PCA representation is a distinct constrained representation, not an exhaustive test of all ESM information.
6. **Operational families.** Sequence-identity connected components are operational homology groups under the reported coverage/search settings. They are not experimentally verified biochemical families, and this review does not independently validate the completeness of all possible homology edges.
7. **Fixed learner settings.** Avoiding held-out tuning is good practice, but a negative fixed configuration does not establish that the learner class has no predictive capacity. The synthetic controls address only a bounded portion of that concern.

## Synthetic-control interpretation

The primary planted relation uses a real acidic-surface fraction as a single standardized scalar, interacting with a smooth solvent logP term. The target generator uses the descriptor distribution globally to define the synthetic experiment; it is not fitting real target information from test folds. This global generator normalization is not real-data target leakage.

| Split | Strong candidate Spearman | Strong gain | Null gain | Shuffled gain |
| --- | ---: | ---: | ---: | ---: |
| Identity 30% plus publication | 0.7945 | 0.7813 | 0.0026 | 0.0095 |
| Identity 30% | 0.7886 | 0.7242 | 0.0005 | -0.0069 |

Both strong controls exceed the recorded criteria of candidate Spearman above 0.7 and gain above 0.3. The moderate controls also improve ranking. Both null and both shuffled gains have intervals spanning zero.

The control establishes that the real folds, condition preprocessing, tree learner and scoring can detect a simple, strong, transferable protein-by-solvent interaction. It does not calibrate detection power for weak biological effects, source-specific noise, high-dimensional embedding interactions or arbitrary nonlinear mechanisms. The one globally permuted protein-feature assignment is an illustrative negative control, not a formal permutation p-value and not a family-preserving exchangeability argument.

The synthetic outcomes make seven originally constant matched series variable. Synthetic scoring therefore covers 770 series, compared with 763 for the real-data primary endpoint, though both score 270 sequences. This is not a scoring defect. Recommended bounded sensitivity: use the original real-data eligible-series IDs and recalculate synthetic scores on those same 763 series, without refitting or changing the target generator. Disclose both results. Main agent owns this follow-up.

## Exact design for easier real-data diagnostics

These designs and support counts were established before fitting their task-specific models. They are exploratory, not preregistered confirmatory analyses. Thresholds below are practical support safeguards, not a formal statistical power calculation.

### Known enzyme in a held-out publication

1. Hold out one complete publication at a time, including every enzyme in that publication. Train on all other publications.
2. Score only held-out rows whose exact sequence appears in at least one training publication. Training may use the full eligible training cohort; do not train only on the repeated-sequence subset unless separately labeled.
3. Retain the established matched-series definition and at least three distinct solvents with a nonconstant observed response.
4. Compare conditions, conditions plus sequence identity as a diagnostic baseline, and available composition/surface/ESM descriptors on identical scored rows for any paired comparison.
5. A sequence-identity baseline must permit sequence-by-solvent effects. An additive identity intercept alone cannot change a fixed-condition solvent ranking. Trees with identity categories can express interactions, or an explicit regularized identity-by-solvent term can be used.
6. Aggregate multiple series and test publications within each sequence. For uncertainty, use sequence/publication-connected components where shared test publications connect multiple sequences. Report the number of such independent components; do not pretend that 60 series are 60 independent enzymes.
7. Do not reuse `fit_oof` unchanged: its whole-sequence fold assertion is correct for the original transfer task but incompatible with intentionally evaluating known sequences in a new publication.

Verified support before structural coverage filtering:

- 493 eligible rows from 15 repeated sequences in 31 publications.
- 392 rows contribute to 60 scoreable matched series.
- Those series cover only 13 sequences in 25 publications.

Run this task as an exploratory diagnostic. A negative result on 13 scoreable sequences cannot establish general impossibility of known-enzyme cross-study prediction. Any positive result needs per-sequence and per-publication inspection, a meaningful identity baseline and independent confirmation.

### Unseen sequence in a familiar family, with publication holdout

1. Build split groups as connected components of exact sequence identity and shared publication, without adding homology-family edges. Thus no sequence or publication crosses folds, but related training sequences can exist.
2. Use five-fold GroupKFold on these groups, with deterministic assignments made before modeling.
3. Prequalify operational 30%-identity families containing at least five distinct sequences and five distinct publications in the original cohort, based on metadata only.
4. Within each fold, score a test row only if its family has at least three distinct training sequences and three distinct training publications. Assert exact test-sequence absence from training.
5. Use the same matched-series eligibility as the original task. Compare conditions, a family-aware conditions baseline, composition, surface and ESM on common rows. Family identity is appropriate because the intended test family is already known.
6. Keep fold-group IDs separate from uncertainty-group IDs. A family deliberately spans folds in this task, so reusing the original `cluster` field for both split validation and family-level bootstrap would be incorrect. Bootstrap family/publication-connected evaluation components, report their count, and present per-family effects and leave-one-family-out influence checks.
7. A same-publication variant, if run, is a weaker diagnostic and must remain clearly separate. It is not a substitute for this study-independent test.

Independent support calculation used the above deterministic five-fold design and thresholds on the full 295-sequence cohort:

- 1,578 eligible rows from 87 sequences in 83 publications and eight families.
- 1,512 rows contribute to 187 scoreable series.
- Those series cover 82 sequences in 77 publications, still only eight families.

| Family ID | Scoreable sequences | Scoreable publications | Series |
| --- | ---: | ---: | ---: |
| seq0001 | 13 | 12 | 24 |
| seq0005 | 14 | 13 | 27 |
| seq0006 | 5 | 5 | 12 |
| seq0008 | 9 | 7 | 33 |
| seq0009 | 6 | 6 | 9 |
| seq0024 | 11 | 14 | 26 |
| seq0036 | 15 | 15 | 39 |
| seq0081 | 9 | 9 | 17 |

Surface coverage or a justified metadata filter will reduce these counts; recalculate and report that support explicitly rather than silently retaining these denominators. Eight families permit a useful pooled diagnostic but weak broad inference across enzyme families. A claimed successful narrow use case should show which families support the gain and whether it survives omission of any one family.

## Decision implications

Proceed with the metadata-qualified sensitivity and these two easier tasks. There is no code-review reason to restart the existing benchmark or download a larger encoder. If the additional tasks remain negative, the defensible decision is to stop the present predictive approach on this database, rather than claim that protein sequence or surface has no relationship to organic-solvent tolerance. If a narrower task improves, require a relevant baseline and study-independent support before choosing bench validation targets.
