# Migration from standalone SSE

The standalone SSE implementation is represented by the `qslib-sse` crate and
the qslib facade's `sse` feature. Canonical `BasisBit`, weighted interactions,
logical chain seeds, operator-string trace closure, and thermodynamic moments
are owned by qslib types.

Migration should proceed in this order:

1. Resolve the model using qslib's canonical row-major sites and explicit
   per-pair coefficients.
2. Select the sign-safe `LocalSseModel` decomposition and verify the supported
   interaction channels.
3. Preserve logical chain identifiers and the versioned seed derivation when
   changing worker counts.
4. Compare tiny thermal energies against exact enumeration before scaling the
   system or changing sweep controls.
5. Store resolved inputs, controls, checksums, and uncertainty metadata in the
   versioned IO artifacts.

Legacy spin labels and x-major site order require named adapters. They must not
be inferred from an array shape or silently accepted as canonical qslib input.
