# Independent support and repeated-condition coverage

Completed 2026-09-11 on the original 5,798-row, 295-sequence cohort.

There are 160 operational 30%-identity families, of which eight contain at least
five sequences from at least five publications. The 40%-identity grouping has
193 families and also eight meeting that threshold. Thresholds describe practical
support, not a formal power calculation or verified biochemical family taxonomy.

Only 15 sequences occur in more than one publication, accounting for 493 rows.
The easier-task reports further count which rows form scoreable solvent panels.
Independent-study sample size is therefore much smaller than the measurement count.

## Cross-publication numerical matches

| Required common recorded fields, in addition to sequence and solvent | Matched groups | Sequences |
| --- | ---: | ---: |
| Solvent fraction | 70 | 9 |
| Plus incubation duration | 40 | 6 |
| Plus incubation temperature | 40 | 6 |
| Plus incubation pH | 0 | 0 |
| Duration, both temperatures and both pH values | 0 | 0 |

The last complete-condition definition is available in 2,579 rows, but none forms
a matching cross-publication group. The incubation-pH definition is available in
2,626 rows. Missing values and genuinely differing conditions both limit matching.
The forty partially matched groups have median cross-paper range 11.59 percentage
points, but they are not certified replicates: pH, buffer, preparation and other
context may differ or be missing. This range cannot estimate measurement noise,
reproducibility or an achievable prediction ceiling.

Basic publication-string normalization (case, whitespace, DOI URL prefixes and
trailing slash) leaves the 272 publication identities unchanged. Two publication
fields are a PDF URL and a PMID reference rather than DOI strings. No identity
repair was required from this simple check; it is not a full bibliographic audit.

All row-level tables and exact match definitions are saved alongside this report.
No absent repeat was interpreted as evidence that the labels are intrinsically noisy.
