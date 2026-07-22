# qslib-quantum-exact

Exact-basis and exact-solver capabilities for
[qslib](https://github.com/lere01/qslib), the convention-first quantum
simulation library.

Most users should depend on the facade with the `exact` feature rather than
on this crate directly:

```toml
[dependencies]
qslib-quantum = { version = "0.2.0", features = ["exact"] }
```

The Cargo package is `qslib-quantum-exact`; the Rust library target is
`qslib_exact`, re-exported by the facade as `qslib::exact`.

## What this crate owns

- Full and fixed-weight exact bases that preserve canonical qslib site and
  state order.
- Dense and CSR matrices assembled from canonical Hamiltonians, with
  Hermiticity validation.
- Deterministic reference eigensolvers: a real-Jacobi dense solver and a
  fully reorthogonalized Lanczos ground-state solver, both reporting
  residuals so callers can reject under-resolved results.
- Thermal sums with overflow reporting, and real- and imaginary-time
  evolution through full spectral decomposition.
- An observable suite: expectations, correlations, magnetizations,
  structure factors, entropies, and quantum Fisher information.

The backend is deliberately small-system: bases scale exponentially with
site count, and the solvers are pure-Rust references rather than BLAS
bindings. See the
[limitations](https://lere01.github.io/qslib/limitations.html) page.

## Documentation

- [User guide](https://lere01.github.io/qslib/), especially the
  [exact backend](https://lere01.github.io/qslib/exact.html) and
  [evolution](https://lere01.github.io/qslib/evolution.html) pages.

Licensed under Apache-2.0.
