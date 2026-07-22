# qslib-quantum-core

Foundational scientific vocabulary for
[qslib](https://github.com/lere01/qslib), the convention-first quantum
simulation library. Every other qslib crate builds on the checked types
defined here.

Most users should depend on the facade rather than on this crate directly:

```toml
[dependencies]
qslib-quantum = "0.2.0"
```

The Cargo package is `qslib-quantum-core`; the Rust library target is
`qslib_core`, re-exported by the facade as `qslib::core` and at the facade
root.

## What this crate owns

- Canonical row-major geometry: rectangular and triangular lattices with
  independent per-axis boundary conditions, custom coordinate sets,
  minimum-image displacements, and shell-based pair selection.
- Little-endian basis states: `BasisBit`, packed words, full and
  fixed-weight sector bases.
- Pair-dependent weighted interactions with typed channels and pinned
  disorder provenance.
- Model resolution for TFIM, isotropic Heisenberg, J1-J2, disordered
  exchange, and Rydberg Hamiltonians with provenance-labelled coefficients.
- Pauli strings, Hermitian Hamiltonians, and local-energy evaluation.
- Symmetry utilities: permutations, translation and point groups, spin
  inversion, and diagonal gauges.
- The versioned `qslib-seed-v1` deterministic seeding scheme.

## Documentation

- [User guide](https://lere01.github.io/qslib/), especially the
  [geometry and interactions](https://lere01.github.io/qslib/geometry-interactions.html),
  [operators and models](https://lere01.github.io/qslib/operators-models.html), and
  [symmetry](https://lere01.github.io/qslib/symmetry.html) pages.
- [Scientific conventions](https://github.com/lere01/qslib/blob/main/docs/conventions.md)
  are the normative contract implementations are tested against.

Licensed under Apache-2.0.
