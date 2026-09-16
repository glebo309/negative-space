# Enzyme solvent tolerance: computational decision

Completed 2026-09-12. All fits, artifact checks and four sequential independent
audits are complete. The final reviewer independently reproduced the last
ablation and endorsed the qualified operational recommendation below.

## Recommendation

Stop the present general sequence/surface predictor development on this dataset,
and do not initiate a robotic bench panel as validation of these predictions.
The tested methods have not established useful protein-specific solvent choice
beyond experimental-condition baselines. An initially large familiar-family gain
shrunk under stronger controls and fixed repeated group splits. The final fairer
ablation retains a modest hint, insufficient to establish useful solvent choice.

This is a decision about the current project configuration. It is not a claim
that protein sequence is unrelated to solvent tolerance, that all models must
fail, or that OrganiZymeDB is generally erroneous. The work has produced a
reproducible benchmark and a reusable surface/representation pipeline, not a
validated solvent-selection model or a demonstrated surface mechanism. A final
composition ablation preserves a weak familiar-family signal, described below;
the recommendation is not based on claiming every comparison was negative.

## The decisive comparisons

Mean within-sequence matched-solvent Spearman correlations. All comparisons use
identical evaluation rows within their stated job.

| Test | Protein model | Relevant comparator | Difference and conditional 95% interval |
| --- | ---: | ---: | --- |
| Full covered, unfamiliar family plus publication, ESM650 with metadata | 0.2199 | 0.1910 | +0.0289, -0.0150 to +0.0756 |
| Metadata-qualified, unfamiliar family plus publication, ESM650 | 0.2719 | 0.2587 | +0.0132, -0.0583 to +0.0791 |
| Same qualified test, surface descriptors | 0.2690 | 0.2587 | +0.0103, -0.0598 to +0.0753 |
| Familiar family, six fixed assignments, composition with metadata, no family label | 0.2263 | 0.2021 | +0.0242, -0.0089 to +0.0824 |
| Known enzyme in a held-out paper, surface | 0.1033 | 0.2267 | -0.1234, -0.4022 to +0.1145 |

The familiar-family aggregate counts each of 80 sequences once after averaging
its available assignments. Eight families connect into four uncertainty groups;
the assignments are not independent replications. The known-enzyme test scores
only 11 sequences, so its negative result is particularly imprecise.

The final repeated familiar-family composition model selects an observed best
solvent 26.60% of the time versus 25.02% for metadata conditions. The difference is
1.57 percentage points, interval -3.50 to +6.23. Pairwise-ordering accuracy is
59.68% versus 58.39%, a small positive secondary difference of 1.29 percentage
points with interval +0.04 to +3.48. This is a weak within-family lead, not a
confirmed best-solvent-selection advantage. The earlier candidate that retained
family identity scored 0.1988 versus 0.2021 on the primary metric; the final
ablation improves on it and is not omitted from the decision.

## What was tested, and what each test resolves

### Data support and original-source audit

The complete downloaded resource contains 28,173 measurements. Its incubation
stability subset contains 9,091 measurements, not 9,091 independent proteins.
The conservative analysis cohort contains 5,798 measurements from 295 exact
sequences. Exact structural matches cover 289 sequences and 5,675 measurements.
Only 15 sequences appear in multiple publications before structural filtering.

A purposive independent audit examined 25 records from eight publications:
11 exact numerical agreements, one immaterial 0.3 percentage-point discrepancy,
one source non-detection encoded as zero, and transparently marked context-only
or inaccessible cases. Genuine large solvent activation was confirmed. This
sample cannot estimate a population error rate.

The important issue is measurement context. Database `WT` can mean an engineered
reference construct, and numeric control 100 does not guarantee an aqueous
reference. Final curator comments flag assay/control uncertainty in 3,470 rows.
Control type, salt, buffer, solvent carry-over and enzyme preparation are not
fully represented by the original seven numeric conditions. The resource itself
acknowledges such ambiguity; this is not evidence of concealed or widespread error.

There are no cross-paper repeats matching all required complete numeric conditions
under the stated rules. Partial matches exist, but their differences cannot
estimate experimental noise or a prediction ceiling. Bibliographic normalization
did not collapse the 272 publication identities.

### Pipeline capability and independent code review

Strong synthetic protein-by-solvent signals are recovered through both strict
fold designs, with candidate correlations approximately 0.79. Moderate signals
are also recovered. Null and shuffled-feature controls show no clear advantage.
Results persist when scoring exactly the original 763 eligible series.

An independent reviewer reproduced all four primary surface comparisons and
bootstrap intervals, checked source targets/IDs/folds and found no material core
leakage, alignment or scoring defect. The positive control establishes capability
for a simple planted relation, not power for arbitrary biology or high-dimensional
embeddings. The prior numerical and exact structure-sequence audits are preserved.

### Representations and justified sensitivities

The earlier package compared amino-acid composition, frozen ESM35/650 embeddings,
27 explicit surface descriptors and combinations. Surface exposure cutoffs,
confidence filtering, learner seeds and a supplementary interaction-ridge model
did not establish unfamiliar-family solvent-choice gains.

The source-motivated qualified subset retains 2,307 rows from 109 sequences before
structural coverage, or 2,293 rows and 107 sequences after it. Its metric scores
320 panels from 97 sequences. Refitting all comparators and adding explicit
control/assay flags does not rescue protein-specific prediction. Water-incubated
and non-incubated control strata are also inconclusive.

Direct within-experiment rank-target training was tested to remove reliance on
comparable absolute activity scales. Protein models beat that method's weak
internal baseline but do not beat the existing activity-trained metadata baseline.
This positive-looking internal comparison is retained, not used to claim progress.

Training-fit and deliberately permissive row-split tests show why fit statistics
are not enough. Protein features improve global log-activity fit, but do not
thereby improve the intended held-out solvent decision. These leaky diagnostics
are descriptive only.

### The apparent familiar-family success and its challenge

The easier task holds out exact sequences and publications while allowing at
least three other training sequences from three publications in a prequalified
family. Families were selected by support, not outcome.

Initially composition reached 0.2441 versus 0.1600 for plain conditions, with a
positive conditional interval. Gains remained positive after omitting any one
family. However, explicit metadata improved the relevant baseline, and adding
family identity was itself degrading it. Composition failed to beat metadata-only
conditions reliably in the original or first alternate assignment.

Four further seeds were fixed before their scores were inspected. Across all six
assignments, metadata-only conditions beat composition in four point estimates.
The per-sequence aggregate is essentially zero improvement. One favorable seed
does not override the aggregate. The full trajectory, every seed and both stronger
and weaker comparators are saved rather than selecting a favorable result.

The final independent reviewer identified one remaining mismatch: the composition
candidate still contained the harmful family-label feature. Removing that feature
from the candidate, while retaining every assignment and comparator, raises its
aggregate score to 0.2263 versus 0.2021. Five of six point estimates are positive,
but the aggregate primary interval still crosses zero. A secondary pairwise gain
has a narrowly positive interval; best-choice utility remains uncertain. The
appropriate conclusion is limited evidence for a small familiar-family signal,
not a validated narrower predictor or a blanket claim of no signal. No additional
architecture search was undertaken to erase that uncertainty.

The final reviewer independently reproduced all final primary and utility
calculations. Six of eight families have positive point gains, but only two of
four publication-connected groups do. Omitting one six-sequence family reduces
the pooled gain to 0.0036. This influence analysis does not justify deleting that
family; it limits confidence that the modest hint is broadly useful. The reviewer
found no remaining implementation must-fix and requested no further model run.

## Scientific interpretation and publication boundary

The evidence supports an operational no-go, not a biological null. It cannot
cleanly apportion failure among missing context, limited independent examples,
representation limitations and genuinely weak cross-family transfer. Both data
heterogeneity and model limitations remain plausible. Larger encoder weights alone
are not the next justified step.

A positive predictor paper or a new surface-mechanism claim is not supported.
A modest benchmark/limitations paper is conceivable, but not automatically earned
by negative models. It would need a substantive reusable benchmark, careful
positioning against the database manuscript, and a venue interested in the result.
The narrower composition/surface idea already has prior art. See
[the 2018 esterase composition study](https://pubmed.ncbi.nlm.nih.gov/30176710/)
and the separate `LITERATURE_SCOPE.md` check.

The [resource record](https://zenodo.org/records/21874133) names Romain Debruyne
and Fabrizio Pucci. Its 10 August 2026 creation date is not the first study of
solvent-tolerance prediction. The associated resource manuscript remained pending
at this check. No authors were contacted or private data uploaded.

Reopening this direction would require materially new evidence: substantially
more comparable independent proteins/studies, source-supported condition repair
that improves a frozen held-out test, or a representation/target justified by a
specific mechanism and tested against the same strong controls. Building a new
standardized dataset could be a different project; it is not validation of the
current unproven predictor and is not recommended here as an automatic next step.

## Reproducibility and artifact map

Project root: `investigations/enzyme-solvent-tolerance/`.
Raw SHA-256: `400823131c996a3899d95cdb0c82e37a6aa5833d32f0daa23e15a696045b3a06`.
Original sequences, labels, structures and prior results were not overwritten.
All comparisons are exploratory. Bootstrap intervals are conditional on fitted
models and the stated grouping, not full training-procedure uncertainty, and
are not multiplicity adjusted. A zero-crossing interval is not an equivalence test.

| Folder | Contents |
| --- | --- |
| `source_audit/` | Original-source ledger, exact metadata predicates and citations |
| `evaluation_review/` | Independent code, metric and fold review |
| `coverage/` | Unique-sequence/family/publication support and partial repeat tables |
| `controls/` | Planted/null/shuffled signals and original-eligible-series check |
| `metadata/` | Inclusion ledger, strict subset and metadata comparators |
| `direct_ranking/` | Rank-target fits and paired cross-objective comparisons |
| `optimistic/` | Training-fit and deliberately leaky row-split diagnostics |
| `easier_tasks/` | Known-enzyme and familiar-family fits, amendments and six-split challenge |
| `utility/` | Pairwise and top-choice metrics from saved predictions |
| `composition_ablation/` | Final comparison omitting family identity from candidate and comparator |
| `decision_review/` | Independent final challenge, reproduced results and qualified decision |

Run from the project root with `PYTHONPATH=src`. Relevant modules:

```sh
python3 -m enzsolv.run_diagnostic_controls
python3 -m enzsolv.run_metadata_diagnostics
python3 -m enzsolv.run_rank_diagnostics
python3 -m enzsolv.run_optimistic_diagnostics
python3 -m enzsolv.run_easier_tasks --phase primary
python3 -m enzsolv.run_easier_tasks --phase robustness
python3 -m enzsolv.run_easier_tasks --phase repeated
python3 -m enzsolv.run_repeated_utility
python3 -m enzsolv.run_composition_ablation
python3 -m pytest -q
```

Use the preserved plans/amendments when reproducing the exploratory sequence.
Caches validate their bound inputs before reuse. New numerical comparisons used
test-first development for grouping, alignment, targets, metadata flags and metrics.
The final check passes 57 tests and verifies 37 prediction tables and 449 cached
diagnostic fits, including synthetic controls and repeated comparators. These
counts are computational inventory, not independent scientific replications.
`completion_checks.json` records file hashes, source identity/target checks,
finite outputs and appropriate fold separation. The earlier 390-fold benchmark
package remains separately preserved.
