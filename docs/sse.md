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
trace closure. `SseSampler` currently provides sign-safe diagonal
insertion/removal updates and thermodynamic expansion-order estimators. A
caller must validate unsupported update families before requesting them; no
legacy `Spin` semantics are silently inferred.
