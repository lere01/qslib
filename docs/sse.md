# Canonical SSE interface

`qslib-quantum-sse` is the sign-safe finite-temperature stochastic-series
expansion backend. It uses qslib's canonical `BasisBit` directly. Bit zero has
Pauli-Z eigenvalue `+1`; bit one has eigenvalue `-1` and is the Rydberg
occupation. No model-independent `Spin` alias is introduced.

The decomposition contract is

\[
H = E_{\rm shift} - \sum_a B_a,
\]

where each sampled matrix element of `B_a` is non-negative. `LocalSseModel`
provides explicit TFIM and Rydberg constructors. TFIM accepts explicit bonds
and non-negative `J` and `h`. Rydberg accepts per-site detunings and explicit
pair interactions, preserving pair-dependent coefficients and applying the
canonical occupation convention.

`BasisSseState` stores a padded operator string and validates propagation and
trace closure. The sampler performs three explicit Metropolis update families:
diagonal insertion/removal, paired off-diagonal vertices, and boundary basis
state flips. The latter is required for ergodic trace sampling and is model
agnostic because it operates on canonical `BasisBit` values. Operator strings
grow geometrically when their identity headroom is low, preserving all existing
vertices. Thermodynamic estimators use expansion-order moments.
The result also reports a naive independent-sample energy standard error;
correlated chains must be combined through independent logical chains for a
confidence interval.

Cutoff growth is permitted during thermalization only. If identity headroom is
insufficient after thermalization, the run fails explicitly instead of mixing
finite-cutoff ensembles in one estimate.

`qslib_sse::legacy` contains opt-in adapters for old `Up`/`Down` labels. The
adapter requires a model family because TFIM and Rydberg occupation semantics
differ. `derive_chain_seed` and `logical_chain_seeds` provide deterministic
logical chain seeds suitable for sequential, threaded, or process-level
execution without coupling results to a particular thread count.
