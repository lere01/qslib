# Versioned scientific IO

`qslib-quantum-io` stores scientific meaning alongside numerical payloads. A
configuration records the convention schema, site and byte order, simulation
basis, resolved coefficients, scalar dtype, backend, tolerances, and the
versioned RNG scheme. JSON and YAML are strict: unknown fields and unsupported
schema versions are rejected.

Artifact manifests bind every immutable file to a BLAKE3 checksum and byte
length and record the convention schema. Checkpoints bind an accepted step,
configuration checksum, parameter-layout fingerprint, RNG algorithm and state
schema, complete evolution controls, accepted-state schema, and checksums for
the typed payload and every named little-endian C-order `f64` NPY array.
ChaCha20 positions are recorded in complete 64-byte blocks and expose their
equivalent 16-word stream position. `atomic_write` completes a sibling
temporary file before replacing the target, so an interrupted write cannot
masquerade as a complete artifact.

Accepted trajectory rows are represented as equal-length typed columns and can
be written as immutable Apache Parquet parts, matching ADR-0004. The dataset
manifest is bound to the resolved configuration and is the transaction boundary
for a set of parts. Completion is published only after every part is readable,
the complete manifest is atomically written, and the exact `COMPLETE` marker is
written last. Readers must validate checksums before interpreting arrays. A
seed is provenance, not a substitute for storing realized disorder or resolved
coefficients.

`ParquetDatasetManifest::load` may explicitly recover a durable complete
manifest whose marker publication was interrupted. Read-only audit tools use
`ParquetDatasetManifest::inspect`, which validates the exact marker and parts
without rewriting user data.

The independent standard-library verifier in
[`tools/verify_io_artifacts.py`](../tools/verify_io_artifacts.py) checks
checkpoint JSON/NPY structure and completed Parquet framing without importing
qslib. Install `pyarrow` separately when full Parquet column decoding is
desired.

To produce a fixture for that check locally:

```text
cargo run -p qslib-quantum-io --example io_artifacts -- /tmp/qslib-io-fixture
python tools/verify_io_artifacts.py \
  --checkpoint /tmp/qslib-io-fixture/checkpoint \
  --dataset /tmp/qslib-io-fixture/dataset
```
