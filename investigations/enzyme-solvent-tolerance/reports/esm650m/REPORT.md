# ESM2 650M solvent-ranking benchmark

Completed 2026-09-11. The checkpoint SHA256 was verified before loading. All 295
sequences were encoded on MPS in 34.77 seconds after model loading. Features have
1,280 dimensions. CPU versus GPU maximum absolute difference was 1.31e-6.

## Primary outcome: solvent ranking at matched recorded conditions

| Holdout | Conditions only | ESM2 650M | Difference | Paired 95% interval |
| --- | ---: | ---: | ---: | --- |
| Exact sequence | 0.286 | 0.266 | -0.020 | -0.067 to 0.027 |
| 40% identity components | 0.231 | 0.218 | -0.012 | -0.067 to 0.041 |
| 30% identity components | 0.220 | 0.225 | 0.005 | -0.039 to 0.052 |
| 30% identity plus publication | 0.197 | 0.215 | 0.019 | -0.022 to 0.065 |

All four intervals include zero. Increasing encoder size from 35M to 650M has not
produced a reliable improvement over the conditions-only baseline on this primary
outcome. This is not a formal paired comparison of 650M versus 35M and does not
establish that sequence information is biologically irrelevant.

The cohort, 5,798 measurements from 295 sequences, and all four saved five-fold
assignments were unchanged. The primary metric averages matched-series Spearman
correlations within sequence, then across 275 evaluable sequences (772 series).
Intervals use 5,000 paired component-bootstrap draws, conditional on the existing
fits and folds, rather than repeated-training uncertainty.

The encoder is frozen. Final-layer residue means exclude padding and special
tokens. Five long proteins use nonoverlapping 1,022-residue chunks with
length-weighted pooling; no residues are discarded, but cross-chunk context is
unavailable. All 1,280 features are supplied with solvent/condition features to
the same 250-tree Extra Trees predictor used in the smaller-model benchmark.

## Secondary outcome

Ranking across all conditions, including concentration variation, improves over
baseline in the 40% and publication-linked splits with intervals above zero.
These are secondary findings and do not demonstrate better selection among
solvents at matching conditions. Raw-scale R2 remains negative in all splits.

## Reproduction

From the project directory:

```bash
PYTHONPATH=src python3 -m enzsolv.run_esm --size 650m
```

Full-precision metrics: `results.json`. Per-measurement predictions, sequence-level
scores, immutable sequence hashes, embedding cache and extraction metadata are
saved beside this report. Model: `facebook/esm2_t33_650M_UR50D`, pinned revision
`08e4846e537177426273712802403f7ba8261b6c`.
