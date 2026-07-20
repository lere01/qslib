# API stability review

This is the pre-1.0 public API freeze record for qslib. It complements the
Rust documentation: rustdoc explains what an item calculates, while this page
records whether a physicist or downstream integrator may rely on that item in
the 1.0 contract.

## Stability dispositions

| Surface | Disposition | Compatibility rule |
| --- | --- | --- |
| `qslib` facade and `qslib::core` | Candidate 1.0 stable | The default facade is core-only. Re-exports preserve canonical basis, geometry, interaction, model, operator, randomness, scalar, and symmetry meanings. |
| `qslib::exact` | Candidate 1.0 stable behind `exact` | Exact matrices, spectra, residuals, observables, thermodynamics, and evolution retain their documented conventions and quantity-specific errors. |
| `qslib::variational` | Candidate 1.0 stable behind `variational` | TDVP statistics, regularization, solver provenance, and transactional evolution retain their documented signs and metadata. |
| `qslib::sse` | Candidate 1.0 stable behind `sse` | SSE update semantics, chain seeds, measurements, and confidence metadata retain the canonical qslib conventions. |
| `qslib::io` | Candidate 1.0 stable behind `io` | Versioned schemas are separate compatibility contracts. Unknown fields and unsupported schema versions remain errors. |
| `qslib-quantum-cli` | Candidate 1.0 user interface | Commands, physically labelled JSON fields, schema identifiers, and error categories are tested semantic interfaces. Human-readable prose may improve without changing machine fields. |
| `qslib-quantum-python` | Candidate 1.0 Python interface | The ABI3 module owns returned arrays and exposes only documented coarse-grained kernels. NumPy shape, dtype, order, and exception contracts are tested. |
| `qslib-sse::legacy` | Stable migration adapter, not a new model API | Legacy spin labels and seed reproduction are explicit compatibility tools. They do not redefine canonical qslib behavior. |
| `qslib-test-support` | Internal test-support surface | Fixtures and oracle helpers are for conformance tests and are not part of the published scientific API. |

## Enum and extension policy

The binary-state, site-order, basis-axis, boundary, lattice, interaction,
Pauli, and model-family enums describe the currently supported scientific
contract. They are intentionally documented as closed for qslib 1.0: adding a
new physical meaning requires a convention review, independent reference
fixtures, and a compatibility decision.

Error enums and numerical-control enums are also candidate-stable, but their
variants are not a license for callers to infer undocumented implementation
details. A new error or backend algorithm requires an API review, a rustdoc
entry, conformance tests, and an explicit semver decision before it is added.
No experimental algorithm is re-exported from a stable module, and the
repository currently has no `experimental` module or feature.

## Freeze evidence

- Every production crate denies missing public documentation.
- The facade has an empty default feature set and additive capability features;
  feature-boundary tests verify the dependency policy.
- The current facade was compared with baseline commit `2584261` by
  `cargo-semver-checks 0.42.0` as an assumed patch release: 165 checks passed
  and 12 were skipped.
- The SSE sampler no longer suppresses deprecated RNG traits; it uses the
  current `rand_core::Rng` contract explicitly.
- Public JSON and serialized artifacts carry independent version identifiers;
  Rust storage layout is not treated as a schema guarantee.

Before a 1.0 tag, a maintainer must rerun the semver check against the selected
release baseline, review any newly exported item, and update this page if an
API moves from candidate-stable to an explicitly experimental or migration
status. This review does not authorize publication, tagging, or remote changes.
