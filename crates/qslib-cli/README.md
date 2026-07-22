# qslib-quantum-cli

Physicist-first command line interface for
[qslib](https://github.com/lere01/qslib), the convention-first quantum
simulation library. The installed command is `qslib`.

Most users should download the platform archive from the
[releases page](https://github.com/lere01/qslib/releases) or install with
Cargo:

```bash
cargo install qslib-quantum-cli
qslib inspect conventions
```

Rust programs should depend on the `qslib-quantum` facade instead of this
crate; the CLI is a thin interface over the validated library kernels.

## Commands

```text
qslib inspect conventions [--json]
qslib inspect environment [--json]
qslib model validate CONFIG [--json]
qslib exact ground-state CONFIG [--json]
qslib exact evolve CONFIG --t-max TIME [--imaginary] [--json]
qslib sse run CONFIG [--json]
qslib artifacts inspect PATH [--json]
qslib conformance self-test [--json]
```

Configuration files are YAML or JSON against the versioned
`qslib-model-input-v1` schema; unknown fields are rejected, and JSON results
carry a provenance object naming the resolved conventions. Errors identify
the physical field that failed validation and exit with status 2.

## Documentation

- [CLI guide](https://lere01.github.io/qslib/cli.html) for configuration
  conventions and per-model field contracts.
- [User guide](https://lere01.github.io/qslib/) for the underlying physics
  and library documentation.

Licensed under Apache-2.0.
