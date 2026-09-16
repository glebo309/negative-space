# Planted-signal controls

Completed 2026-09-11. The raw measurements remain unchanged.

The same real condition design, structural coverage, strict folds, train-only
preprocessing, fixed tree settings and matched-solvent metric were used on a
synthetic target. Its protein term is a standardized acidic surface-area fraction
interacting with solvent logP. It deliberately represents a simple transferable
relationship, not a realistic generative model of every assay.

| Split | Signal | Candidate mean rho | Improvement over conditions | Conditional 95% interval |
| --- | --- | ---: | ---: | --- |
| Family plus publication | Strong | 0.7945 | 0.7813 | 0.6755 to 0.8780 |
| Family plus publication | Moderate | 0.6402 | 0.3775 | 0.2934 to 0.4504 |
| Family plus publication | None | 0.6197 | 0.0026 | -0.0089 to 0.0121 |
| Family plus publication | Shuffled feature | 0.0227 | 0.0095 | -0.0182 to 0.0381 |
| Family | Strong | 0.7886 | 0.7242 | 0.6298 to 0.8117 |
| Family | Moderate | 0.6386 | 0.4114 | 0.3147 to 0.5094 |
| Family | None | 0.6103 | 0.0005 | -0.0108 to 0.0110 |
| Family | Shuffled feature | 0.0575 | -0.0069 | -0.0503 to 0.0355 |

Both strong controls pass the recorded criterion: candidate rho above 0.7 and
improvement above 0.3. Both null and both shuffled-feature intervals contain zero.
The one feature permutation is an illustrative negative control, not a formal
permutation p-value or a family-preserving exchangeability test.

Synthetic noise makes seven originally constant series variable, so the initial
metric scores 770 series. The independent review requested scoring only the
original 763 eligible real-data series without refitting. That check changes no
conclusion: strong gains are 0.7819 and 0.7243; moderate gains are 0.3778 and 0.4117;
null and shuffled intervals still include zero. Both analyses score 270 sequences.

Conclusion: the pipeline can recover this simple global protein-by-solvent
signal. This is not a statistical-power claim for subtle biological effects,
high-dimensional embeddings, local structural mechanisms or missing assay context.
See `analysis_plan.json`, `results.json` and
`original_eligible_series_results.json` for exact settings and full precision.
