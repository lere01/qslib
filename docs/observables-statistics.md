# Observables and statistics

Exact expectation and variance functions normalize by the supplied state norm,
so callers may provide an unnormalized exact vector. Zero-norm and non-finite
states are errors. Matrix rows and columns retain the canonical exact-basis
order, and variance is evaluated from `H^2` rather than from a sampled square.

`WeightedMoments` uses normalized non-negative weights and a parallel-moments
merge formula. This makes independently accumulated batches equivalent to one
serial accumulator without discarding weight information. Complex estimators
retain independent real and imaginary moments through
`ComplexWeightedMoments`.

`autocorrelation` uses the named initial-positive-sequence estimator and
reports both integrated autocorrelation time and effective sample size.
`r_hat` is the classic equal-length, unsplit Gelman-Rubin diagnostic. The
algorithm is part of the result contract and should be recorded with reports.

`disorder_average` keeps each realization identifier and computes the weighted
fixed-realization mean separately from between-realization variance. This
prevents sampling uncertainty within one realization from being silently
combined with finite-disorder ensemble uncertainty.
