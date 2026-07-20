# ADR-0003: Deterministic randomness and seed derivation

- Status: accepted
- Date: 2026-07-19
- Owners: qslib maintainers

## Context

SSE chains, disorder realizations, stochastic estimators, and adaptive stages
need independent random streams. Worker scheduling must not change scientific
streams. A seed alone is insufficient unless the algorithm, derivation, and
logical stream identity are fixed and recorded.

## Decision

Use `ChaCha20Rng` as the qslib 1.x portable pseudorandom generator. Persist its
algorithm identifier and qslib seed-scheme version. Represent a master seed as
32 bytes. Convenience constructors may expand a `u64` master value through the
same versioned derivation.

Derive each 32-byte child seed with BLAKE3 keyed by the 32-byte master seed. The
hashed message has this exact canonical byte encoding:

1. ASCII bytes `qslib-seed-v1` followed by one zero byte;
2. the domain byte length as `u32` little-endian, followed by the UTF-8 domain
   such as `disorder`, `sse_chain`, or `integrator_stage`;
3. the number of logical indices as `u32` little-endian; and
4. every stable logical index or accepted-state identifier as `u64`
   little-endian in documented order.

Expand a convenience `u64` master seed with `blake3::derive_key` using context
string `qslib master seed v1` and the seed's little-endian bytes. Worker number,
thread count, process identifier, hash-map order, and rejected-attempt ordinal
are not stream identities. A retry at the same physical stage may use a
documented retry index, while rejection must not consume streams belonging to
later accepted states. Persist realized disorder coefficients in addition to
RNG provenance.

## Alternatives considered

- `StdRng` was rejected because its algorithm is intentionally not a stable
  reproducibility contract.
- Sequentially splitting one mutable master stream was rejected because task
  scheduling would affect child assignment.
- A custom SplitMix-style derivation was retained only for explicit legacy SSE
  parity because it lacks domain separation and collision-resistant framing.
- Operating-system randomness remains available only for creating a new master
  seed, not for replayable execution.

## Consequences

Logical jobs reproduce independently of parallel scheduling and can be replayed
from complete metadata. Exact sample sequences may change only with a named
seed-scheme or algorithm change.

The Milestone 0 dependency probe validated `rand_chacha 0.10.0`,
`rand_core 0.10.1`, and `blake3 1.8.5`. The selected graph compiled on Rust
1.85 and current stable and passed the configured license, advisory, and source
checks.

## Validation

- Fixed derivation vectors cover domains, indices, and `u64` expansion.
- Chain and realization streams match across thread counts and job order.
- Resume produces the same accepted trajectory as uninterrupted execution.
- Serialized provenance names the RNG and seed scheme, and disorder artifacts
  contain every realized coefficient.
- Fixed vectors and supported-platform determinism remain implementation gates
  before the RNG API is released.
