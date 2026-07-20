# ADR-0001: Layered Cargo workspace boundaries

- Status: accepted
- Date: 2026-07-19
- Owners: qslib maintainers

## Context

qslib must serve lightweight model-building users as well as exact, variational,
thermal, command-line, and Python workflows. Putting every capability in one
crate would make foundational types depend on heavy numerical and interface
stacks. Parallel definitions of sites, basis states, or interactions would be
scientifically unsafe.

## Decision

Use one Cargo workspace with a small root facade and these responsibility
crates:

- `qslib-core` owns identifiers, basis meaning, geometry, weighted
  interactions, operators, canonical models, symmetry primitives, and
  observable definitions;
- `qslib-exact` owns exact bases, matrices, eigensolvers, thermodynamics, and
  exact evolution;
- `qslib-variational` owns local-energy, statistics, TDVP, and integration
  kernels that consume caller-supplied wavefunction data;
- `qslib-sse` owns sign-safe stochastic series expansion algorithms;
- `qslib-io` owns versioned configuration, artifacts, and checkpoints;
- `qslib-python` owns Python bindings and foreign-array conversion;
- `qslib-cli` owns the physicist-first executable; and
- unpublished `qslib-test-support` owns independent fixtures and test oracles.

Dependencies point from interfaces and algorithms toward `qslib-core`. Core
must not depend on another qslib crate. The facade re-exports stable APIs under
features and contains no second implementation. Algorithm crates do not depend
on `qslib-io`; IO converts public domain values to versioned data-transfer
objects. CLI and Python crates depend only on public library APIs. Production
crates never depend on test support. Workspace package versions are synchronized
for 1.0.

If the facade is published to crates.io, every non-development path dependency
must also be published under its approved prefixed package identifier. A facade
release may not depend on unpublished path-only capability crates.

## Alternatives considered

- One crate was rejected because it couples compile time, dependencies, and
  API stability across unrelated capabilities.
- Many model-specific crates were rejected because they encourage duplicate
  physical vocabulary and make cross-model algorithms harder to reuse.
- Keeping SSE outside the workspace was rejected because it would preserve an
  ambiguous spin convention and a second geometry vocabulary.

## Consequences

The core-only path remains small and suitable for future simulation programs.
Heavy capabilities have explicit dependency and feature costs. Cross-crate API
changes require coordination, but the ownership boundary is clear.

## Validation

- `cargo tree -p qslib-core` contains no solver, CLI, Python, Arrow, or SSE
  dependency.
- Every workspace dependency edge follows the documented direction.
- Core-only and all-feature builds run in CI.
- Public examples use the facade unless they explicitly teach a backend.
