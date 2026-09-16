# Independent decision challenge

Review dates: 2026-09-11 to 2026-09-12. Scope: decision logic, new diagnostic code and independent
recalculation from saved predictions. Raw data, code and vault notes were not
modified by this reviewer. Only this report is reviewer-owned.

## Status and recommendation

Stop or pause development of the current general predictor on this dataset and
withhold a predictor-validation bench panel. Retain a narrower exploratory
composition lead, not a validated narrow predictor. The fairest final ablation
reveals a small positive ranking hint, including a weak positive secondary
pairwise interval, but does not establish useful improvement in choosing the
best solvent. The recommendation is resource-bounded, not a biological null.

One bounded omitted ablation was identified before endorsing this decision:
remove explicit family identity from the familiar-family composition candidate,
because that same feature substantially degrades the metadata-only comparator.
The main agent accepted the challenge and completed it on all six existing
assignments without changing their rows, folds, support rules or comparator.
The reviewer independently verified the resulting scores below. This omission
is resolved; no implementation must-fix remains from this review.

This is a model-specification limitation, not an implementation bug. A worse
composition-plus-family result against metadata without family does not isolate
the incremental effect of composition. It was more actionable than downloading
another encoder or searching many new model classes. It is explicitly a
review-triggered exploratory check, not independent confirmation.

## Independently recomputed findings

The reviewer recalculated per-panel Spearman correlations and decision metrics
directly from all six saved familiar-family prediction tables using SciPy and
elementary pair comparisons, without calling the project's score functions.
Duplicate solvents were averaged within each matched panel; panels were averaged
within sequence; each sequence was averaged over its available assignments, then
counted once. Stored per-sequence scores agree within 2e-16.

| Metric | Composition with family and metadata | Metadata without family | Difference | Conditional 95% interval |
| --- | ---: | ---: | ---: | --- |
| Matched Spearman | 0.1988295120 | 0.2020669241 | -0.0032374121 | [-0.0619686961, 0.0569833452] |
| Best observed solvent hit | 0.2573247354 | 0.2502314815 | 0.0070932540 | [-0.0595511768, 0.0449891068] |
| Unequal-pair ordering accuracy | 0.5845456830 | 0.5839112971 | 0.0006343859 | [-0.0245859110, 0.0240432563] |

The bootstrap was separately implemented with 5,000 component draws and seed 11.
The saved grouping order was used to reproduce its numerical draws after
independently checking the grouping itself.

An independently implemented union-find over shared-publication family edges
recovers four evaluation components containing 31, 20, 9 and 20 scored sequences.
There are 80 scored sequences, of which 72 have six assignments and eight have
five. Neither 472 sequence-assignment records nor six CV runs are independent
biological sample sizes. Component labels can differ without changing membership.

All individual assignment deltas are reproduced: original +0.0270564; seed 47
-0.0349806; seed 17 +0.0367351; seed 23 -0.0099459; seed 71 -0.0269165;
seed 101 -0.0314149. The original favorable assignment has not been numerically
misreported or silently excluded.

The cross-objective surface comparison was independently recalculated too.
On 763 panels from 270 sequences, direct-rank surface prediction scores
0.1926005281 against 0.1910061979 for activity-trained metadata, a gain of
0.0015943302. On the qualified 320 panels from 97 sequences it scores
0.2185729480 against 0.2586799445, a difference of -0.0401069965. The joined
measurement IDs are unique and their original observed targets agree exactly.
Thus the rank-target method's advantage over its weaker internal baseline does
not supply a practical rescue against the already available comparator.

## Code review and interpretation challenges

Inspected code includes the familiar/known task support and split construction,
training-only preprocessing, prediction-cache signatures, repeated aggregation,
metadata predicates, normalized-rank targets and secondary utility metrics.
No material numerical, alignment, split or scoring defect was found in these
inspected paths. This does not certify every possible code path or model refit.

1. **Comparator selection:** using a stronger available metadata comparator is
   correct for deciding whether a predictor adds practical value. Selecting it
   after earlier outcomes does not make this a confirmatory study. Conversely,
   adding a known harmful family feature only to the candidate leaves an obvious
   ablation unanswered. The main agent resolved precisely that issue.
2. **Initial positives:** the original familiar-family composition gain, positive
   leave-one-family-out score summaries, and rank-trained protein gains are real
   exploratory observations. Their appropriate status is a hypothesis pursued
   by further checks, not independent validations and not fabricated successes.
   Later results must be reported with that trajectory rather than replacing it.
3. **Outcome-dependent scope:** the easier-task family support thresholds use
   sequence/publication counts rather than response values. Metadata filtering
   follows source-defined literal rules rather than selecting high-performing
   proteins. Repeated seeds were fixed before their scores. Those safeguards
   reduce outcome cherry-picking, but the overall development sequence remains
   exploratory because later questions were motivated by earlier results.
4. **Uncertainty:** four family/publication components are a weak inferential
   base. The exact reported bootstrap can be reproduced but cannot be treated
   as a precise population-level confidence guarantee. It conditions on fitted
   models, shares training data among predictions, excludes annotation
   uncertainty and has no multiplicity adjustment. A zero-crossing interval is
   not an equivalence result. The negative known-enzyme diagnostic covers only
   11 scoreable sequences and cannot adjudicate general cross-study learnability.
5. **Target mismatch:** log-activity training does not directly optimize the
   primary per-protein ranking objective. The normalized-rank sensitivity is a
   useful bounded challenge, not an exhaustive ranking-loss search. Its rankings
   are evaluated on original observed activity, and monotone transformed output
   is correctly not described as predicted percentages.
6. **Practical utility:** the secondary hit and pairwise scores correctly handle
   prediction ties and unequal observed pairs. They concern the solvent panels
   present in historical studies, not an arbitrary candidate solvent library or
   an externally calibrated probability of success. They ignore uncertainty in
   the observed optimum and do not establish a useful gain for the current fits.
7. **Source quality:** the audit is purposive, not a population error-rate study.
   It identifies genuine context heterogeneity and reference-construct naming,
   not wholesale transcription failure. The qualified subset remains only
   metadata-qualified; it is not certified homogeneous aqueous stability.
8. **Synthetic controls:** recovery of a strong planted scalar interaction is
   evidence of bounded pipeline capability. It does not prove power to extract
   weak mechanisms from frozen high-dimensional embeddings or rule out model
   limitations. Null and shuffled results are diagnostic checks, not a formal
   biological null or a calibrated permutation significance test.

## Claims the evidence can and cannot support

Supported now: the fitted predictors have not earned predictor-driven bench
validation; the current database/model combination has not demonstrated a stable
protein-specific solvent-choice advantage; richer experimental context and many
more independently comparable proteins could change the result. The resource
and exact structure/sequence pipeline remain reusable.

Unsupported: no sequence effect exists; solvent tolerance is fundamentally
enzyme-specific; the database is unlearnable or generally wrong; a larger model
must fail; the measured intervals bound all possible model improvement; a
negative benchmark automatically earns a paper. No finite model search can
establish those universal claims, and no such search is needed for a resource-
bounded project decision.

A modest limitations/benchmark manuscript is a separate publication decision,
requiring a substantive reusable contribution and accurate positioning against
prior work and the resource manuscript. It cannot currently be framed as a
validated predictor or discovered surface mechanism.

## Completed review-triggered ablation

The candidate now differs from metadata-only conditions solely by adding the 20
amino-acid frequencies. Neither candidate nor comparator receives family identity.
The inspected new code uses complete saved training-fold assignments, checks
sequence/publication separation and retains all covered training rows. Prediction
tables contain exactly the same ordered evaluation IDs and unchanged comparator
predictions. The reviewer independently recalculated all six assignments from
their raw prediction panels, without calling the project's score functions.

| Metric | Composition without family | Metadata without family | Difference | Conditional 95% interval |
| --- | ---: | ---: | ---: | --- |
| Matched Spearman | 0.2262640598 | 0.2020669241 | 0.0241971357 | [-0.0088900838, 0.0823850045] |
| Unequal-pair ordering accuracy | 0.5967955922 | 0.5839112971 | 0.0128842951 | [0.0004172181, 0.0348091546] |
| Predicted top solvent's normalized observed rank | 0.6235786944 | 0.6050048563 | 0.0185738381 | [-0.0013489404, 0.0536104486] |
| Best observed solvent hit | 0.2659771825 | 0.2502314815 | 0.0157457011 | [-0.0350437876, 0.0622507123] |

The sequence-once aggregate remains 80 sequences in four components. The
individual matched-rank gains are original +0.0461005; seed 47 +0.0001352;
seed 17 +0.0413633; seed 23 +0.0230219; seed 71 +0.0418631; and seed 101
-0.0233895. Five of six assignments and six of eight families have positive
point estimates. This is a directional exploratory signal and must be preserved.

Independent score-only influence calculations show positive mean rank gains in
two of the four uncertainty components: +0.1005881 and +0.0419336, versus
-0.0124142 and -0.0034278 in the others. Omitting any one family leaves a positive
pooled point estimate, ranging from +0.0035643 to +0.0463777. The smallest gain
follows omission of family seq0009, which contains six scored sequences. This
describes concentrated influence; it is not a refit, family selection rule or
independent robustness replication.

The pairwise interval narrowly excludes zero. It is wrong to summarize every
decision metric as null, or say the fair composition signal disappeared. It is
also wrong to promote this exploratory secondary result into independent
confirmation after many development choices. Best-choice hit improvement is
1.57 percentage points, with conditional interval -3.50 to +6.23 points. No useful
best-solvent improvement is established; no equivalence or maximum possible gain
is established either. No minimum practically useful gain was specified in
advance, so the conclusion is not a formal utility-threshold test.

## Final decision boundary

An operational stop on the present broad predictor and no predictor-validation
bench panel remain justified. The evidence now includes a modest familiar-family
composition hint, so the report must use that qualified conclusion rather than
the earlier stronger assertion that no narrower signal survived. Preserve this
lead for a materially better supported or independently validated study; it
does not automatically warrant more model search or immediate bench work.

This review does not demand infinitely many models to establish absence. The
bounded ablation directly addressed a documented harmful feature and materially
improved the candidate. It is now complete. Further attempts would be new
exploratory development, not a remaining correctness requirement. The scientific
question about transferable sequence information remains open even though the
current project configuration has not earned advancement to predictive validation.
