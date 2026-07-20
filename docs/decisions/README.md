# Architectural decision records

This directory records consequential qslib design decisions and their
rationale. Use an ADR when a choice is difficult to reverse, changes public
semantics, introduces a major dependency, establishes a cross-cutting
abstraction, or affects interoperability with another project.

Copy [`0000-template.md`](0000-template.md) and assign the next four-digit
number. Keep one decision per file. Use one of these statuses:

- `proposed`
- `accepted`
- `superseded by ADR-NNNN`
- `rejected`

An ADR does not override `docs/conventions.md` silently. A decision that changes
a scientific convention must update the convention schema and the normative
document explicitly.

## Decision index

- [ADR-0001](0001-workspace-boundaries.md): layered Cargo workspace boundaries
  - accepted
- [ADR-0002](0002-linear-algebra-sparse-storage.md): linear algebra and sparse
  storage - accepted
- [ADR-0003](0003-randomness-reproducibility.md): deterministic randomness and
  seed derivation - accepted
- [ADR-0004](0004-artifacts-serialization.md): versioned serialization and
  columnar artifacts - accepted
- [ADR-0005](0005-python-ffi-ownership.md): Python FFI and array ownership -
  accepted
- [ADR-0006](0006-api-stability-and-features.md): public stability and feature
  policy - accepted
- [ADR-0007](0007-standalone-sse-migration.md): additive standalone SSE
  migration - accepted
- [ADR-0008](0008-repository-ownership.md): dedicated repository ownership
  boundary - accepted
- [ADR-0009](0009-registry-package-names.md): registry package identifiers -
  accepted
