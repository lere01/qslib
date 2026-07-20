# ADR-0004: Versioned serialization and columnar artifacts

- Status: accepted
- Date: 2026-07-19
- Owners: qslib maintainers

## Context

Physicists need editable configurations, inspectable summaries, scalable
trajectories, and restartable checkpoints. Native Rust layouts and unlabelled
arrays are not durable scientific formats. One format is not optimal for every
role.

## Decision

Use a versioned, role-specific artifact set:

- YAML and JSON are accepted human-readable inputs and resolve to one canonical
  schema with unknown fields rejected by default;
- JSON stores manifests, resolved configurations, concise summaries, checksum
  indexes, and checkpoint metadata;
- an append-only dataset of immutable Apache Parquet part files stores
  trajectory and measurement tables, with an atomically replaced manifest; and
- a JSON checkpoint envelope plus semantically named NPY files stores typed
  array payloads.

Complex columns are represented by named real and imaginary children. Array
metadata records shape, dtype, and logical order. Every durable scientific
document records `qslib-conventions-v1` and its own schema identifier. BLAKE3
checksums cover exact file bytes and are listed in a manifest. A checkpoint is
written to a temporary sibling directory, flushed and synced where supported,
then atomically renamed. Each Parquet part is completed atomically before its
manifest entry appears. A completed-run marker is written last.

Schemas evolve additively within a version. A loader rejects an unknown major
schema and never guesses a legacy convention. Compatibility loaders resolve to
canonical documents before algorithms see the data.

## Alternatives considered

- JSON for all data was rejected because large numeric trajectories are bulky
  and weakly typed.
- CSV was rejected because it cannot reliably carry nested, complex, dtype, or
  schema metadata.
- A Rust-specific binary serializer was rejected as the durable public format
  because cross-language longevity is a core requirement.
- HDF5 was rejected as the default because native linking and global-library
  behavior complicate portable builds.

## Consequences

The IO crate carries substantial Arrow and Parquet dependencies, so it remains
feature-gated away from core. Artifacts are portable to Python and other data
tools, but schema code and atomic recovery tests are mandatory.

The Milestone 0 dependency probe validated Serde `1.0.229`, serde_json
`1.0.150`, serde_yaml_ng `0.10.0`, Parquet `59.1.0` with Arrow and Zstandard,
ndarray `0.17.2`, ndarray-npy `0.10.0`, and BLAKE3 `1.8.5`. The complete graph
compiled on Rust 1.85 and current stable and passed configured license and
source checks.

Parquet's graph also contains unmaintained compile-time `paste 1.0.15`.
RUSTSEC-2024-0436 has no safe upgrade and reports no vulnerability. The same
monitored exception as ADR-0002 applies until Arrow migrates.

## Validation

- Golden schema round trips preserve resolved couplings, conventions, dtype,
  shape, tolerances, algorithms, and fingerprints.
- Unknown fields, versions, truncated files, and checksum mismatches fail with
  structured context.
- Parquet and NPY payloads round trip through independent Python readers.
- Failure injection proves that partial writes are never reported complete.
- Independent Python round trips and Linux and Windows builds remain mandatory
  implementation gates before IO release.
