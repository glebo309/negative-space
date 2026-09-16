# Independent source and label audit

Audit date: 2026-09-11. Raw measurements were not modified.

## Conclusion

This audit does not establish a high transcription-error rate or an intrinsically
unlearnable database. It does establish that the current cohort definition is
less homogeneous than its earlier description suggests: a numeric control of
100 does not specify the reference experiment, and `mutation == WT` can designate
an engineered reference construct. Curators already flag assay/control ambiguity
in many records. These are concrete reasons for a metadata-qualified sensitivity
analysis, not grounds to discard all high values or declare that sequence has no
biological relevance.

## Sampling and limits

Purposive, non-random sample of 25 measurement IDs from eight publications, chosen
to cover very large activation, ordinary values, an explicit non-detection,
engineered reference constructs, inclusion bodies, salt-rich buffer, the largest
publication by cohort rows, and a highly represented buffer-comparison study.
This selection cannot estimate population error prevalence. Three papers were
available as full primary text and tables through Europe PMC XML; three additional
primary abstracts identified formulation/construct or buffer context; two were
limited to abstracts without the relevant original numbers. No inaccessible record
is counted as verified. The ledger distinguishes value checks from context checks.
Repeated independent-paper agreement is outside this audit and requires the main
coverage analysis; these source checks are not experimental replication.

## Direct source findings

### LipC12: activation is real, with one minor transcription discrepancy

DOI 10.1186/1475-2859-10-54, Table 2 and Methods. Seven selected records were
compared to numerical table entries. Six values agree exactly. ID 7979 has 1690.6%
in the database versus 1690.9% in the table. This 0.3 percentage-point discrepancy
is immaterial to solvent ranking here. The extreme methanol/THF values are not
decimal-place errors. Incubation is 48 h at 4°C, 15/30% v/v; enzyme is diluted
100-fold and solvent concentrations equalized before assay. The construct is
His-tagged, with no cloning mutation reported. Table control is described as the
reference condition set to 100%; the database's finer non-incubated interpretation
is not independently resolved by that footnote alone. Source:
[primary full text](https://link.springer.com/article/10.1186/1475-2859-10-54),
[primary XML and Table 2](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3161859/fullTextXML).

### EstA1: engineered reference constructs and censored non-detection

DOI 10.4014/jmb.1907.07057, Table 1 and Methods. Six selected records were
compared: five numerical values agree, and ID 15594 is encoded as zero although
the source says ND, explicitly defined as not detected. Thus it is a censored
measurement, not a proven exact zero. Both Δ20 and Δ54 enzymes are deletion
mutants; full-length WT failed expression. All 48 cohort rows are nevertheless
`mutation=WT`, with the truncations correctly disclosed in their names. This is
a reference-construct convention, not evidence that the database silently changed
the sequences. Incubation concentrations (15/30/50%), 30 min at 30°C, subsequent
assay pH 8 and 40°C agree. No-solvent activity is 100%; solvent carry-over is not
quantified in that methods description. Source:
[primary full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC9728341/),
[primary XML and Table 1](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9728341/fullTextXML).

### ADH/A1a: salt and normalization are essential context

DOI 10.1002/2211-5463.12557, Figure 2B caption and Methods. Three selected time
points were context-audited, not independently redigitized. The source explicitly
uses enzyme incubated without solvent as control and incubates the enzyme in
10 mM HEPES/NaOH pH 7.5 with 2 M NaCl at room temperature. These conditions are
recorded in the database. The existing seven numeric predictor features do not
include salt concentration, although `assay_solution` retains it and matched-series
scoring separates different strings. Room temperature correctly remains missing
in the numeric temperature parser. Source:
[primary full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC6356862/),
[primary XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6356862/fullTextXML).

### Other targeted context checks and inaccessible numeric results

- DOI 10.1007/s12010-012-0028-7, IDs 9245 and 9264: primary abstract confirms
  active inclusion bodies and solvent activation of rSML. The 21 cohort rows
  therefore should not silently be described as soluble purified enzymes. The
  exact two numerical values, fractions and normalization were not verified from
  full methods. [Primary abstract](https://link.springer.com/article/10.1007/s12010-012-0028-7).
- DOI 10.1007/s12010-012-0014-0, IDs 9266 and 9268: primary abstract explicitly
  identifies a chimeric BTL2 bearing A112G/H113E/Q115A. The five cohort rows are
  `mutation=WT` but their names disclose the chimera. Numerical values were not
  independently verified. [Primary abstract](https://link.springer.com/article/10.1007/s12010-012-0014-0).
- DOI 10.1016/j.chembiol.2007.08.010, ID 5512: the database explicitly says the
  reference and assay both contain solvent despite `aqueous_control=100`.
  Primary abstract confirms an engineering study, but does not verify this
  specific MtLT2 datum or its normalization. This is a database-internal cohort
  interpretation flag, not a source-confirmed numerical error.
  [Primary abstract](https://pubmed.ncbi.nlm.nih.gov/17884637/).
- DOI 10.1016/0141-0229(93)90102-8, IDs 2510 and 2511: primary abstract confirms
  that buffer identity changes stability at the same pH, and discusses protective
  additives. Full numerical tables were inaccessible, so these values are not
  verified. Buffer type is retained as free text but absent from predictors.
  [Primary abstract](https://pubmed.ncbi.nlm.nih.gov/7764108/).
- DOI 10.1016/j.chemosphere.2019.06.113, IDs 15217 and 15218: selected because
  this publication contributes 180 cohort measurements. Primary metadata match;
  numerical results and relevant methods were inaccessible.
  [Primary record](https://pubmed.ncbi.nlm.nih.gov/31234090/).

## Full-cohort metadata counts and exact suggested rules

Counts below are direct computations on `prepare_cohort` from the unchanged CSV:
5,798 rows, 295 exact sequences, 272 publication strings. They are metadata counts,
not source-verified counts. Case-insensitive matching was used.

| Predicate | Rows | Sequences | Publications |
| --- | ---: | ---: | ---: |
| Final comments segment contains `uncertain` | 3480 | 188 | 172 |
| Above, and same segment contains `assay`, `activity`, or `control` | 3470 | 187 | 171 |
| Other uncertainty only | 10 | 1 | 1 |
| No assay/control uncertainty under that rule | 2328 | 112 | 101 |
| No assay/control uncertainty plus procedure says `in aqueous phase` | 2307 | 109 | 98 |
| Comments match `control[^|]*presence of organic solvent` | 204 | 12 | 11 |

The final segment means `comments.fillna('').str.split('|').str[-1]`.
The ten other-only rows say `Uncertainty about substrate, here taken from the
reference`; do not misclassify them as uncertain assay mode. Some mode-uncertain
rows also mention sequence identity, buffer, substrate or temperature, so flags
are overlapping and not a complete ontology.

Control vocabulary: 4,837 rows mention `Non-incubated control`, 943 mention
`Water-incubated control`, and 18 mention neither. All 18 instead say
`Unclear control (in %), activity measured in aqueous phase | High uncertainty
about whether the assay was measured with or without organic solvent, uncertainty
about control`. Procedures contain `in aqueous phase` in 5,470 rows and
`in the presence of organic solvent` in 152. These literal matches are deliberately
conservative and do not exhaust synonyms or resolve carry-over.

Representative exact final comments segments and counts:

- `Uncertainty about whether assay was in the presence of organic solvent or in aqueous phase`: 386.
- `Uncertainty about whether assay really was in aqueous phase`: 307.
- `High uncertainty about whether assay in aqueous phase or in the presence of organic solvent`: 279.
- `Uncertainty about control`: 101.
- `Uncertainty about substrate, here taken from the reference`: 10.

Proposed sensitivity, without changing the main cohort or pretending it is
fully curated:

1. Exclude the 3,470 explicitly mode-uncertain rows; separately show a stricter
   subset whose procedure explicitly says aqueous phase (2,307 rows before any
   structure requirement). Preserve original folds and fit all comparators on
   identical remaining rows. Absence of a flag is not proof of correctness.
2. Exclude explicit organic-reference rows for an aqueous-reference question;
   retain them only as a separately named operational-stability endpoint.
3. Add explicit control type and assay mode as metadata, then inspect matched
   ranking separately by non-incubated versus water-incubated controls. They
   define different denominators even if both are numerically 100.
4. Flag formulation/construct examples above by DOI, not ad hoc value thresholds.
   If seeking natural soluble WT, a real construct/preparation audit is needed;
   the current mutation column cannot enforce it alone. Do not blanket-exclude
   truncated mature proteins without examining what was actually assayed.
5. Do not cap activation at 100 or delete all large values. The checked extreme
   values are genuine. Treat ND-as-zero as censoring in interpretation; this
   small audit cannot establish its prevalence elsewhere.

The audit provides one tiny numerical correction and several measurement-context
qualifications. It does not yet prove those qualifications explain the failed
benchmarks. Only the matched sensitivity runs can test that proposition.
