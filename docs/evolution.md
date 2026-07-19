# Transactional variational evolution

`qslib-quantum-variational` provides model-independent Euler and Heun
integration over a checked flat parameter vector. The integrator does not own a
neural-network or model object. A caller supplies the current flat state and a
velocity callback, then applies the accepted state to its model after each
accepted boundary.

## State and callbacks

`FlatState` contains finite parameter values and a deterministic layout
fingerprint. `EvolutionMetadata` records the schema version, physical time,
proposed next step, accepted and rejected counts, seed, method, metric, and
layout fingerprint. Metadata is JSON round-trippable and cannot be restored
against a different layout fingerprint.

The callback receives `(parameters, time, stage_seed)` and returns a
`Velocity`. A velocity contains a real parameter direction and may include a
validated positive-semidefinite dense QGT. Stage seeds are derived from the
master seed, accepted step, and stage index. Retry attempts expose their
attempt count for diagnostics, but the derived accepted-node seed is invariant
under rejection, so a resumed accepted-boundary trajectory reproduces the same
callback stream.

## Euler, Heun, and adaptive acceptance

Euler uses one velocity evaluation. Heun evaluates a predictor and then a
corrector. When adaptive mode is enabled, the predictor-corrector difference
is measured in either the Euclidean norm or the QGT norm

\[
\|\delta\|_S = \sqrt{\delta^T S\delta}.
\]

An accepted step commits parameters, advances physical time, increments the
accepted count, and records a bounded controller-proposed next step. The
controller applies the safety factor to both rejection shrinkage and accepted
step growth; an exact zero error permits growth up to the configured maximum.
Adaptive Euler is rejected because it has no error estimator. A rejected attempt changes
none of those accepted quantities. Only the rejected count and proposed step
are updated before the retry. Observation callbacks run only after acceptance,
so rejected trial states cannot appear in physical trajectories.

The integrator owns no sampling or automatic differentiation. Python or other
frontends can compute local energies, QGT directions, and observables and pass
them through the callback without per-sample foreign-function calls.
