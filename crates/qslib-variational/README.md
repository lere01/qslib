# qslib-quantum-variational

Variational statistics, TDVP, and integration kernels for
[qslib](https://github.com/lere01/qslib), the convention-first quantum
simulation library.

Most users should depend on the facade with the `variational` feature rather
than on this crate directly:

```toml
[dependencies]
qslib-quantum = { version = "0.2.0", features = ["variational"] }
```

The Cargo package is `qslib-quantum-variational`; the Rust library target is
`qslib_variational`, re-exported by the facade as `qslib::variational`.

## What this crate owns

- Online, mergeable weighted moments with overflow-atomic updates.
- Chain diagnostics: integrated autocorrelation time, effective sample
  size, and split-chain Gelman-Rubin R-hat.
- Disorder averaging with between-realization uncertainty.
- TDVP estimation from sampled ratios and local energies, fingerprinted
  parameter layouts, dense and streaming quantum geometric tensors,
  conjugate-gradient solves, and Tikhonov regularization with GCV scoring.
- Transactional Euler/Heun evolution drivers with adaptive stepping,
  deterministic per-stage seeds, and checkpoint restore.

The crate supplies statistics and solver infrastructure; it does not bundle
a wavefunction ansatz or optimizer. See the
[limitations](https://lere01.github.io/qslib/limitations.html) page.

## Documentation

- [User guide](https://lere01.github.io/qslib/), especially the
  [TDVP](https://lere01.github.io/qslib/tdvp.html),
  [evolution](https://lere01.github.io/qslib/evolution.html), and
  [observables and statistics](https://lere01.github.io/qslib/observables-statistics.html)
  pages.

Licensed under Apache-2.0.
