# Observables and statistics

Exact expectation and variance functions normalize by the supplied state norm,
so callers may provide an unnormalized exact vector. Zero-norm and non-finite
states are errors. Matrix rows and columns retain the canonical exact-basis
order, and Hermitian variance is evaluated from the centered vector norm rather
than from a subtraction of two potentially ill-conditioned moments.

`WeightedMoments` uses normalized non-negative weights and a parallel-moments
merge formula. This makes independently accumulated batches equivalent to one
serial accumulator without discarding weight information. Complex estimators
retain independent real and imaginary moments through
`ComplexWeightedMoments`.

`autocorrelation` uses Geyer's initial-positive-sequence estimator with a
common `1/N` autocovariance normalization, summing adjacent pairs until the
first non-positive pair, and reports
both integrated autocorrelation time and effective sample size. The reported
time is clamped to at least one, so alternating or very short chains never
produce a negative effective sample size. `r_hat` is the
classic equal-length, unsplit Gelman-Rubin diagnostic and returns its within-
chain `W`, between-chain `B`, value, and algorithm identity. The algorithm is
part of the result contract and should be recorded with reports.

`disorder_average` keeps each realization identifier and computes the weighted
fixed-realization mean separately from between-realization variance. One
realization reports ensemble variance as unavailable, rather than zero. This
prevents sampling uncertainty within one realization from being silently
combined with finite-disorder ensemble uncertainty. When per-realization
sampling variances are supplied, the summary reports the uncertainty of the
weighted ensemble mean, `sum_r w_r^2 Var_r / (sum_r w_r)^2`, considering only
positive-weight records.

The exact backend also provides direct Pauli-string expectations, two-site
correlations, axis-labelled Pauli magnetization, spin magnetization using
`S^a = sigma^a/2`, and the exact `<S_tot^2>` observable from an explicit basis
and state.
Structure-factor
results retain their physical axis and raw versus connected convention. A
weighted sublattice helper returns per-configuration signed, absolute, and
squared values; callers average those components over configurations.
