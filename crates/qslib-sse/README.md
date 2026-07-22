# qslib-quantum-sse

Sign-safe finite-temperature stochastic series expansion (SSE) quantum Monte
Carlo for [qslib](https://github.com/lere01/qslib), the convention-first
quantum simulation library.

Most users should depend on the facade with the `sse` feature rather than on
this crate directly:

```toml
[dependencies]
qslib-quantum = { version = "0.2.0", features = ["sse"] }
```

The Cargo package is `qslib-quantum-sse`; the Rust library target is
`qslib_sse`, re-exported by the facade as `qslib::sse`.

## What this crate owns

- Sign-safe local decompositions of transverse-field Ising and Rydberg
  Hamiltonians with explicit non-negative shifts, built from tuples or from
  canonical resolved interactions.
- Trace-validated padded operator strings over canonical `BasisBit` states,
  with adaptive cutoff growth.
- Update families selected per run through `UpdateScheme`: local
  (diagonal insertion/removal, paired off-diagonal vertices, boundary
  flips), the deterministic TFIM linked-cluster update, production
  site-local Rydberg world-line moves, and a Metropolis-corrected global
  cluster correctness reference.
- Expansion-order thermodynamics (energy and heat capacity), optional
  recorded per-measurement series, and update statistics.
- Deterministic, worker-count-independent logical chain seeds under the
  versioned `qslib-seed-v1` scheme, plus explicit legacy adapters for the
  standalone SSE program's spin labels and seeds.

Statistical estimates require user-selected thermalization, chain counts,
and autocorrelation analysis; see the
[limitations](https://lere01.github.io/qslib/limitations.html) page.

## Documentation

- [User guide](https://lere01.github.io/qslib/), especially the
  [SSE](https://lere01.github.io/qslib/sse.html) and
  [migration from standalone SSE](https://lere01.github.io/qslib/migration-sse.html)
  pages.

Licensed under Apache-2.0.
