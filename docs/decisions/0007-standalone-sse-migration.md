# ADR-0007: Additive migration of standalone SSE

- Status: accepted
- Date: 2026-07-19
- Owners: qslib maintainers

## Context

The standalone `quantum-sse` crate contains valuable TFIM and Rydberg SSE
behavior, but its `Spin::Up` means Pauli-Z plus for TFIM and occupied for
Rydberg. Those meanings cannot share one canonical conversion. A destructive
copy or replacement would risk both physics and existing users.

## Decision

Build `qslib-sse` on canonical `qslib-core` basis, geometry, model, and weighted
interaction types. Port algorithms incrementally, beginning with neutral tests
of physical decomposition, shift restoration, update balance, measurements,
and deterministic chains. Do not import the ambiguous legacy spin enum.

Keep standalone SSE operational until qslib parity gates pass. Put every legacy
conversion in a model-aware adapter named for TFIM or Rydberg. Preserve current
sign-safe parameter restrictions unless a new decomposition is independently
proved. After parity, a separate owner-authorized release may turn the old crate
into a compatibility shim or deprecate it. This ADR does not authorize deleting
or publishing anything in the standalone repository.

## Alternatives considered

- Copying source wholesale was rejected because it would import conflicting
  conventions and duplicate shared types.
- Making legacy `Spin` canonical was rejected because no model-independent
  mapping exists.
- Rewriting SSE without parity fixtures was rejected because validated update
  behavior and deterministic execution would be lost.
- Immediately retiring the old crate was rejected as destructive and outside
  autonomous authority.

## Consequences

Migration takes place behind explicit adapters and duplicate paths remain for a
time. The result has one canonical physical vocabulary and a recoverable
transition for existing users.

## Validation

- Model-aware state vectors cover both legacy mappings in the convention
  specification.
- Resolved decompositions and restored physical energies match independently.
- Tiny systems agree with exact thermal traces within stated statistical
  criteria.
- Per-chain streams and results remain invariant under thread-count changes.
- Unsupported sign structures fail before sampling begins.
