# Easier real-data tasks: final diagnostic report

Date: 2026-09-11. All runs complete. Original measurements, sequences and prior
benchmarks were not modified.

## Decision-relevant conclusion

The easier tasks do not establish a reliable use case for the current predictor.
An initially promising within-family composition result does not survive the
combination of an adequate metadata-only baseline and repeated grouped split
assignments. Across six assignments, with each sequence counted once, composition
scores 0.1988 versus 0.2021 for metadata-only conditions: gain -0.0032, conditional
95% interval [-0.0620, 0.0570]. The known-enzyme cross-publication task is too small
for a strong negative inference and provides no demonstrated gain.

This is evidence against advancing these fitted predictors to bench validation,
not evidence that sequence or surface is biologically irrelevant. The initial
positive, the comparator weakness and all follow-up outcomes are retained below.

## Design and support

All models use the same structure-covered cohort of 5,675 rows from 289 sequences
for training eligibility. Tests score only task-supported subsets. The underlying
cohort is the existing database reference-construct incubation endpoint, not a
verified natural-WT, purified-soluble, identically normalized assay collection.
The source audit's endpoint and formulation limitations still apply.

Unseen-sequence, familiar-family task:

- Assign five folds to connected components of exact sequences and shared
  publications. Neither a sequence nor a publication crosses train/test.
- Assign on the original 295-sequence cohort before structural filtering.
- Prequalify operational 30%-identity families with at least five distinct
  sequences and five distinct publications in the original cohort.
- Score a test sequence only when its family has at least three other training
  sequences from at least three training publications. Train models on all
  eligible training rows, not just the supported test families.
- Bootstrap family/publication-connected evaluation components, distinct from
  the sequence/publication split groups. A family intentionally spans folds.

The full-cohort calculation independently reproduces the reviewer's 82 scoreable
sequences and 187 series. Structural coverage reduces the primary task to:

| Support | Count |
| --- | ---: |
| Eligible measurement rows | 1,517 |
| Eligible sequences / publications | 85 / 79 |
| Scoreable measurement rows | 1,452 |
| Scoreable sequences / publications | 80 / 73 |
| Matched solvent series | 182 |
| Operational families | 8 |
| Family/publication uncertainty components | 4 |

Shared papers connect families seq0001/seq0006/seq0036, seq0005/seq0009 and
seq0024/seq0081. Treating these as eight independent families would ignore those
shared studies. Four components remain a weak basis for general inference.

Known-enzyme task holds out an entire publication and predicts only exact
sequences represented in another training publication. It requires separate
split checks because sequence overlap is intentional. Structural coverage leaves
432 eligible rows from 13 repeated sequences in 27 publications; 332 rows form
55 scoreable series covering only 11 sequences in 21 publications. The scoring
uncertainty calculation has 11 sequence/publication components.

## Fixed learner and scoring

ExtraTrees uses 250 trees, minimum leaf size 5, all features available to each
split, seed 11 and log1p(percent) targets. Imputation and categorical encoding fit
only training data. The original seven numeric conditions and solvent identity
are retained. Composition contributes 20 residue fractions; surface contributes
the existing 27 RSA-0.20 descriptors; ESM uses the cached 650M representation.

Within-family protein candidates also receive family identity. The family-aware
conditions baseline permits family-by-solvent interactions. The known-enzyme
identity tree similarly permits sequence-by-solvent interactions, not just an
additive identity intercept.

Matched-series solvent Spearman is averaged within each sequence, then across
sequences. Paired intervals use 5,000 component draws conditional on the saved
fits. They do not include curation uncertainty, full learning-procedure uncertainty
or multiplicity correction. Follow-ups were explicitly documented after earlier
outcomes, not presented as confirmatory preregistration.

## Primary familiar-family result

| Model | Mean matched Spearman |
| --- | ---: |
| Conditions | 0.1600 |
| Conditions plus family identity | 0.1431 |
| Composition plus conditions and family | 0.2441 |
| Surface plus conditions and family | 0.2380 |
| ESM650 plus conditions and family | 0.1818 |

Composition improves over the family-aware baseline by +0.1011 [0.0719, 0.1177]
and over plain conditions by +0.0842 [0.0206, 0.1351]. The latter gain is positive
in six of eight families and remains +0.0572 to +0.1183 after omitting any one
family from the score summary. This was a real positive exploratory result, not
discarded because later checks were less favorable.

Surface improves over the family-aware baseline by +0.0949 [0.0260, 0.1447], but
not clearly over the stronger plain conditions comparator: +0.0780
[-0.0253, 0.1635]. Surface does not improve over composition: -0.0062
[-0.0587, 0.0430]. ESM likewise does not clearly beat conditions and is below
composition by -0.0623 [-0.1198, -0.0266] in this task.

## Metadata and split robustness

The first amendment adds the seven literal control/assay features from
`metadata_diagnostics.META_FEATURES`, including normalization type, assay phase
and curator mode uncertainty. It does not filter away uncertainty or claim full
metadata correction. On the original folds, composition still improves over the
metadata-plus-family baseline by +0.0710 [0.0391, 0.1000]; surface improves by
+0.0594 [0.0263, 0.0782].

On shuffled GroupKFold assignment seed 47, training-support safeguards leave 72
scoreable sequences and 165 series in eight families, connected into five
uncertainty components. Composition gain over metadata-plus-family is +0.0526
[-0.0659, 0.1347], and surface gain is +0.0258 [-0.0576, 0.0540]. Some individual
families then have very few supported test sequences: one in seq0006, two in
seq0009. Training-support safeguards are not a guarantee of precise family tests.

On identical shared measurement IDs and the same 165 matched series, composition
gain is +0.0828 in the original assignment versus +0.0526 in seed 47. Both shared
row calculations use the same five uncertainty components. Thus changed cohort
membership does not fully explain the attenuation. This comparison is score-only,
not a refit or independent replication.

Crucially, adding family identity also weakens the metadata-aware comparator.
Without family identity, metadata conditions scores 0.1936 in the original
assignment and 0.2388 in seed 47. Composition then has gains +0.0271
[-0.0379, 0.0982] and -0.0350 [-0.0788, 0.0089]. The apparent benefit against the
family-aware comparator cannot be treated as benefit against the stronger
available control.

## Final bounded six-assignment check

A second amendment fixed four additional shuffled group assignments, seeds 17,
23, 71 and 101, before their scores were inspected. Only three models were run:
metadata conditions, metadata conditions plus family, and composition plus both.
No seed or family was selected after results.

| Assignment | Scored sequences | Composition gain over metadata conditions | Conditional 95% interval |
| --- | ---: | ---: | --- |
| Original deterministic | 80 | +0.0271 | [-0.0379, 0.0982] |
| Shuffle 47 | 72 | -0.0350 | [-0.0788, 0.0089] |
| Shuffle 17 | 80 | +0.0367 | [0.0082, 0.0589] |
| Shuffle 23 | 80 | -0.0099 | [-0.1580, 0.1194] |
| Shuffle 71 | 80 | -0.0269 | [-0.0725, 0.0458] |
| Shuffle 101 | 80 | -0.0314 | [-0.0928, 0.0470] |

For the aggregate, each sequence's matched score is averaged over supported
assignments, then each sequence is counted once. Seventy-two sequences have six
supported assignments and eight have five. The six runs are not independent
biological replicates and are never counted as 472 independent sequences.

| Aggregate comparator | Composition mean | Comparator mean | Gain | Conditional 95% interval |
| --- | ---: | ---: | ---: | --- |
| Metadata conditions | 0.1988 | 0.2021 | -0.0032 | [-0.0620, 0.0570] |
| Metadata conditions plus family | 0.1988 | 0.1435 | +0.0554 | [-0.0010, 0.0933] |

The aggregate has 80 sequences, eight families and four uncertainty components.
Against metadata-only conditions, per-family mean gains are positive in four
families and negative in four. Leaving one family out gives aggregate gains from
-0.0264 to +0.0143, changing sign. There is no stable pooled advantage over the
strong comparator. This does not prove that every within-family architecture or
composition model without explicit family encoding must fail; it evaluates the
fixed, documented candidates and provides no validated narrower predictor.

## Known-enzyme held-out-publication result

| Model | Mean matched Spearman | Gain over conditions | Conditional 95% interval |
| --- | ---: | ---: | --- |
| Conditions | 0.2267 | Reference | |
| Sequence-identity interaction tree | 0.1690 | -0.0578 | [-0.2182, 0.0661] |
| Surface plus conditions | 0.1033 | -0.1234 | [-0.4022, 0.1145] |

All three models were fit for each of 27 held-out publications, using the same
eligible test rows. This supplies no positive result. Eleven scoreable sequences
are insufficient to adjudicate general known-enzyme cross-study learnability.
No increasingly expensive known-enzyme ESM search was undertaken.

## Verification and artifacts

Ten new tests were written and observed failing before implementation, then
passed. The complete shared test suite passed 54 tests at the final check.
Checks cover distinct-unit support, split leakage, train-only categories,
target-bound prediction caching, deterministic alternative folds and sequence-once
repeated-score aggregation.

Final artifact audit checked 16 prediction tables against raw measurement IDs,
sequence IDs and values, and checked all 206 cached model fits for finite
predictions. All sequence/publication split and training-support assertions pass.
Primary prediction files remain unchanged during follow-ups. Raw SHA-256 remains
`400823131c996a3899d95cdb0c82e37a6aa5833d32f0daa23e15a696045b3a06`.

Reproduction from the project directory:

```bash
python3 -m pytest -q
PYTHONPATH=src python3 -m enzsolv.run_easier_tasks --phase primary
PYTHONPATH=src python3 -m enzsolv.run_easier_tasks --phase robustness
PYTHONPATH=src python3 -m enzsolv.run_easier_tasks --phase repeated
```

The primary plan, two amendments, full-cohort folds, exact assignment and support
tables, prediction caches, sequence scores, per-family scores, score-only
leave-one-family-out tables and machine-readable results are retained in this
directory. The final repeated-assignment summary is
`repeated_split_stability/results.json`; per-sequence support frequency is in
`composition_vs_conditions_aggregate_scores.csv` beside it.

No individual seed, successful family, isolated interval or alternate comparator
should be presented as independent confirmation. These easier-task diagnostics
support withholding predictor-driven bench validation on the present evidence.
