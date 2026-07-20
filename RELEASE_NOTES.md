# qslib 0.1.0 candidate

This is a local release-candidate snapshot for the convention-first qslib
quantum simulation library. It is not a published release and does not carry
1.0 stability or publication authorization.

## Included capabilities

- Canonical row-major geometry, little-endian basis states, explicit physical
  axes, and pair-dependent weighted interactions.
- TFIM, isotropic Heisenberg, J1-J2, disordered exchange, and Rydberg model
  resolution with provenance-labelled coefficients.
- Exact small-system matrices, spectra, ground states, thermal sums, evolution,
  observables, symmetry utilities, TDVP kernels, and sign-safe SSE sampling.
- Versioned YAML/JSON configuration, checkpoint, trajectory, Parquet, and
  provenance contracts.
- A physicist-facing CLI and an ABI3 Python module imported as
  `qslib_quantum`.

## Validation in this snapshot

The local evidence includes Rust 1.85 and stable workspace tests, strict
Clippy and rustdoc checks, cargo-deny and cargo-audit policy checks, bounded
CLI and state-conversion fuzzing, core Miri tests on nightly Rust, branch
coverage, semver comparison against the previous local candidate, clean Python
wheel/sdist smoke tests, and a checksum-verified local artifact bundle.

## Remaining 1.0 gates

Remote Linux/macOS/Windows CI execution remains outside this local-only run.
The separately owned ncli repository now contains an opt-in qslib exact backend
and focused parity tests. Do not tag, publish, sign, deploy, push, or pull this
candidate without explicit owner authorization. Consult
[`docs/plans/qslib-v1.md`](docs/plans/qslib-v1.md) and
[`docs/release-evidence.md`](docs/release-evidence.md) for exact evidence and
known limitations.
