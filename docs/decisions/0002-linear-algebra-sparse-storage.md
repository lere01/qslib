# ADR-0002: Linear algebra and sparse storage

- Status: accepted
- Date: 2026-07-19
- Owners: qslib maintainers

## Context

Exact diagonalization and TDVP need real and complex dense algebra, sparse
Hermitian matrices, and matrix-free solves. Scientific output must not expose a
backend layout as its meaning. Native BLAS dependencies complicate portability,
while a general sparse package does not by itself provide the residual and
degeneracy guarantees required by qslib.

## Decision

Own a small checked compressed-sparse-row representation at the qslib backend
boundary. Rows and columns use basis indices, entries in each row are sorted by
column, duplicates are combined deterministically, and exact zeros are removed.
Serialization uses `u64` row offsets and column indices and never serializes
`usize`. Deserialization checks every conversion to the host index width.
The public physical API exposes Hamiltonian action rather than backend matrix
storage, with explicit conversion functions for users who request matrices.

Use concrete public reference scalars `f64` and `Complex64`; do not make the 1.0
API generically scalar without a demonstrated use case. Use `faer` behind exact
and variational crates
for dense decompositions and suitable matrix-free kernels. Disable unneeded
default parallel behavior and let qslib own thread policy. A qslib Hermitian
Lanczos implementation is the sparse extremal reference path unless a
separately accepted dependency proves all required ordering, complex, residual,
and degeneracy behavior. Dense diagonalization remains the guarded small-system
fallback and oracle.

Dependency versions are pinned by `Cargo.lock` and must compile at the workspace
MSRV. Backend types do not cross the stable facade without an explicitly named
adapter feature.

## Alternatives considered

- `nalgebra` is strong for small statically shaped algebra but is not the best
  primary fit for large dynamic Hermitian problems.
- `sprs` supplies useful CSR and CSC data structures, but its pre-1.0 API and
  lack of the required Hermitian eigensolver make an owned narrow CSR boundary
  simpler.
- `ndarray-linalg` and direct LAPACK bindings were rejected as the default
  because they introduce native library and platform configuration into the
  normal build.
- Writing all dense decompositions in qslib was rejected as unnecessary risk.

## Consequences

qslib gets portable dense reference algebra and stable sparse semantics. It
must maintain a small CSR implementation and a carefully tested Lanczos path.
Backend upgrades are possible without changing scientific types.

The Milestone 0 dependency probe validated `faer 0.24.4` with default features
disabled and `std`, `linalg`, and `sparse` enabled. Its complete selected graph
compiled on Rust 1.85 and current stable on macOS and passed the configured
license and source policy. Cross-platform behavior remains an implementation CI
gate rather than an assumption in this record.

The graph contains unmaintained `paste 1.0.15` through `gemm`. RustSec
RUSTSEC-2024-0436 reports maintenance status, not a vulnerability, and no safe
upgrade exists. qslib records a monitored exception for this compile-time macro
and must remove it when upstream migrates.

## Validation

- Independently assembled dense and CSR matrices agree entry by entry and under
  matrix-vector action.
- Hermitian solvers report residuals and compare degenerate invariant
  subspaces, not arbitrary eigenvectors.
- Real and complex reference cases pass on Linux, macOS, Windows, and MSRV.
- A core-only dependency audit excludes `faer`.
- The selected `faer` release and complete locked graph pass license, source,
  Rust 1.85, and stable checks. Linux and Windows remain mandatory CI gates
  before the exact or variational capability is released.
