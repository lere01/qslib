# ADR-0005: Python FFI and array ownership

- Status: accepted
- Date: 2026-07-19
- Owners: qslib maintainers

## Context

ncli and physicist workflows need NumPy interoperability without moving model
semantics back into Python. Foreign buffers can be non-contiguous, mutable,
wrongly typed, or shorter-lived than a Rust view. Unsafe FFI must not enter
scientific crates.

## Decision

Use PyO3 and the Rust `numpy` crate in the dedicated `qslib-python` crate, built
with Maturin. Target Python 3.12 or later to match ncli and test `abi3-py312`
wheel compatibility before adopting that wheel policy. The binding crate is the
only owner of Python object conversion and any unavoidable FFI `unsafe` code.
Core scientific crates forbid unsafe code.

Borrow Python arrays only for the dynamic extent of one Python call while the
GIL token and owning object remain valid. Validate dtype, rank, shape, binary
values, strides, finite values, mutability, and canonical site order before
calling a scientific kernel. Copy into owned Rust memory when contiguity,
alignment, lifetime, or concurrent execution cannot be proven. Never retain a
borrowed NumPy pointer in a Rust model, checkpoint, iterator, or background
task.

Release the GIL around pure Rust work only after all inputs are safely owned or
immutably borrowed according to PyO3's contract. Convert structured Rust errors
to a documented Python exception hierarchy. Rust panics must not cross the FFI
boundary.

## Alternatives considered

- A C ABI was rejected for the primary interface because PyO3 provides safer
  ownership and exception integration.
- Embedding Python in scientific crates was rejected because it reverses the
  dependency direction.
- Always copying was rejected as the only policy because safe contiguous
  read-only inputs can support efficient call-scoped borrowing.
- Exposing Torch internals directly was rejected for 1.0. ncli converts or
  supplies validated CPU arrays at the binding boundary.

## Consequences

Python packaging is isolated and independently testable. Some calls copy data
for safety, and zero-copy behavior is a documented optimization rather than a
semantic guarantee.

The Milestone 0 dependency probe validated PyO3 `0.29.0` with `abi3-py312` and
the Rust NumPy crate `0.29.0`. Their complete selected graph compiled on Rust
1.85 and current stable and passed configured license and source checks. This
accepts the ownership boundary and candidate ABI; actual wheels remain subject
to Milestone 12 platform and Python tests.

## Validation

- Binding tests cover C and Fortran order, sliced and negative strides, wrong
  dtype, wrong shape, read-only and mutable arrays, and non-finite values.
- Repeated calls under Python garbage collection and multiple threads pass
  sanitizing and stress tests.
- Rust and Python exact results agree on neutral fixtures.
- An unsafe-code audit finds no unsafe block outside the FFI crate or reviewed
  dependencies.
- Python 3.12+ import, ABI, ownership, and supported-wheel validation remain
  mandatory before binding release.
