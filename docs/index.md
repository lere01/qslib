# qslib documentation

qslib is a quantum simulation library for finite lattice models. The project
defines physical conventions first, then exposes the same checked objects
through Rust, the command line, and Python.

## Start here

- [Scientific conventions](conventions.md) - canonical site order, bits,
  physical axes, bases, couplings, boundaries, and normalization.
- [CLI guide](cli.md) - run the four-site exact and tiny SSE examples.
- [Python guide](python.md) - use the `qslib_quantum` NumPy binding.
- [Geometry and interactions](geometry-interactions.md) - pair-dependent and
  disordered couplings.
- [Operators and models](operators-models.md) - Hamiltonian signs and model
  construction.
- [Exact methods](exact.md) - basis enumeration, spectra, evolution, and
  residual diagnostics.
- [Observables and statistics](observables-statistics.md) - totals, densities,
  axes, correlations, and uncertainty.
- [TDVP](tdvp.md) - real- and imaginary-time statistics and solves.
- [SSE](sse.md) - sign-safe finite-temperature sampling.
- [IO](io.md) - versioned configurations, checkpoints, and trajectories.
- [Symmetry](symmetry.md) - site permutations, characters, and projections.
- [Reproducibility](reproducibility.md) - provenance, seeds, checksums, and
  uncertainty.
- [Limitations](limitations.md) - finite-system scope and numerical caveats.
- [Migration from ncli](migration-ncli.md)
- [Migration from standalone SSE](migration-sse.md)

## API reference

Build the Rust API reference locally with:

```text
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --all-features
```

The generated entry point is `target/doc/qslib/index.html`. The API reference
is generated output and is intentionally not committed to the source tree.

## Architecture and governance

- [Architecture overview](architecture/README.md)
- [Architecture decisions](decisions/README.md)
- [qslib 1.0 execution plan](plans/qslib-v1.md)
- [Toolchain policy](toolchain-policy.md)
- [Security policy](../SECURITY.md)
- [Contribution policy](../CONTRIBUTING.md)

All worked examples in the CLI guide are executed by integration tests. The
Python wheel and its contract tests are built separately because Python is an
optional foreign-function boundary rather than a dependency of the Rust core.
