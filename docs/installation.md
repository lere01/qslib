# Installation

qslib is currently distributed as source, a Rust workspace, and an optional
Python wheel. The release bundle contains the `qslib` command, API docs, and
wheel artifacts when those platform builds are available.

## Run the command locally

From a source checkout with Rust 1.85 or newer:

```text
cargo run --locked -p qslib-quantum-cli -- inspect conventions
```

For a local optimized binary:

```text
cargo build --locked --release -p qslib-quantum-cli
target/release/qslib inspect environment
```

## Rust library

Add the facade to a Rust application using a local checkout while the v1 crate
remains unpublished:

```toml
[dependencies]
qslib-quantum = { path = "../qslib" }
```

Enable only the capabilities needed by the application, such as `exact`,
`variational`, `sse`, or `io`. The public Rust API is linked from the combined
[Rust API reference](api/qslib/index.html).

## Python binding

Install a platform wheel from a release bundle with:

```text
python -m pip install qslib_quantum-*.whl
```

The import name is `qslib_quantum`; the Python package exposes NumPy-backed
exact, observable, and TDVP contracts documented in [Python bindings](python.md).

## Verify an installation

Every release bundle includes `SHA256SUMS`. Verify the files before use, then
run:

```text
qslib inspect conventions --json
qslib conformance self-test --json
```

The self-test is a one-fixture smoke test. Scientific production runs must
also execute the full conformance and model-specific validation suites.
