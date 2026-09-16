# Supplementary solvent-interaction model

Completed 2026-09-11. This check was specified after seeing the primary surface-tree
results. It tests a second predictor family, not an independently preregistered
confirmatory hypothesis. No parameter was selected using held-out scores.

## Results

Mean matched-solvent Spearman, evaluated on the same 5,675 measurements from
289 proteins as the primary surface benchmark. The scoring subset has 763 series
across 270 proteins.

| Holdout | Conditions | Composition | Surface | ESM650, 16 PCs |
| --- | ---: | ---: | ---: | ---: |
| 30% identity plus publication | 0.158 | 0.196 | 0.137 | 0.172 |
| 30% identity | 0.162 | 0.168 | 0.152 | 0.161 |
| 40% identity | 0.127 | 0.152 | 0.153 | 0.121 |
| Exact sequence | 0.178 | 0.188 | 0.140 | 0.184 |

Surface minus conditions, paired component-bootstrap 95% intervals:

| Holdout | Difference | Interval |
| --- | ---: | --- |
| 30% identity plus publication | -0.021 | -0.088 to 0.052 |
| 30% identity | -0.010 | -0.064 to 0.046 |
| 40% identity | 0.026 | -0.029 to 0.081 |
| Exact sequence | -0.037 | -0.093 to 0.016 |

All surface-versus-conditions intervals include zero. Composition and ESM also
have intervals including zero against their corresponding conditions baseline.
Surface underperforms composition in the strictest split, with an exploratory,
unadjusted interval of -0.108 to -0.007. These checks do not rescue a transferable
surface-prediction claim. They also do not exclude other representations or
predictors.

## Model specification

Ridge regression uses fixed alpha 100, LSQR, tolerance 1e-6 and maximum 2,000
iterations. All completed fits converged in at most 109 iterations. The target is
log(1 + measured percent). Predicted log values below zero are floored at zero
before converting back to nonnegative percent values; counts are recorded.

The conditions baseline includes the seven existing numeric condition features,
missingness indicators, solvent one-hot indicators, and solvent-by-condition
products. Candidates add standardized protein features, protein-by-solvent
products and protein-by-condition products. Thus the baseline already allows
different condition responses for different solvents.

Whole-sequence composition uses 20 features; surface uses 27. ESM650 uses 16
principal components, fitted only on distinct training proteins and standardized
using those proteins. Surface and composition scaling likewise uses distinct
training proteins, not duplicated measurement rows. Condition imputation/scaling
and solvent categories are learned on training rows only. Unknown solvent
categories receive zero one-hot indicators. No three-way interaction is included.

Folds are the unchanged sequence/component/publication assignments, restricted
to the same structure-covered rows. All paired comparisons use identical held-out
measurements. Intervals resample components 5,000 times, conditional on these fits.
They are not multiplicity-adjusted or repeated-training uncertainty intervals.

## Numerical validation

The local NumPy 2.2.5 build uses Apple Accelerate. It emitted divide-by-zero,
overflow and invalid-operation warnings from matrix multiplication, including on
a small random finite-valued example. That example agreed with non-BLAS summation
within 1.07e-14. Similar warnings are documented in
[NumPy issue 28687](https://github.com/numpy/numpy/issues/28687).

Warnings were not taken as proof either of numerical failure or of safety. The
first fold of each model in every split was refitted and checked independently
(16 fits): predictions computed by non-BLAS summation agreed with the normal
prediction within 5.33e-15; the largest independently computed relative ridge
gradient residual was 2.25e-5; and regenerated predictions matched the cached
values exactly. All exported predictions are finite. This is a targeted audit,
not a numerical proof of every internal operation in every fit.

See `numerical_audit.json`, `analysis_plan.json`, `results.json`, fold caches and
per-measurement predictions in this directory.

## Reproduce

```bash
PYTHONPATH=src python3 -m enzsolv.run_interactions
```
