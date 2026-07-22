# qslib-quantum-io

Versioned scientific configuration and artifacts for
[qslib](https://github.com/lere01/qslib), the convention-first quantum
simulation library.

Most users should depend on the facade with the `io` feature rather than on
this crate directly:

```toml
[dependencies]
qslib-quantum = { version = "0.2.0", features = ["io"] }
```

The Cargo package is `qslib-quantum-io`; the Rust library target is
`qslib_io`, re-exported by the facade as `qslib::io`.

## What this crate owns

- Versioned YAML and JSON schemas for scientific configuration, run
  specifications, manifests, and summaries, with strict round-trips and
  semantic validation on load.
- blake3-checksummed artifact manifests bound to their configuration
  checksum.
- Checkpoint bundles with named NPY arrays and captured RNG state for exact
  resume.
- JSON and Apache Parquet trajectory storage, including append-only Parquet
  datasets with immutable parts, a completion marker, and crash recovery.
- Atomic write primitives (temp file, rename, directory fsync) used by
  every artifact path.

Every schema carries an explicit `schema_version`; unknown fields and
unsupported versions are errors, never silent reinterpretation.

## Documentation

- [User guide](https://lere01.github.io/qslib/), especially the
  [IO](https://lere01.github.io/qslib/io.html) and
  [reproducibility](https://lere01.github.io/qslib/reproducibility.html)
  pages.

Licensed under Apache-2.0.
