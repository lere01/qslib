# qslib

`qslib` is a quantum simulation library written in Rust. It is intended to
provide shared scientific definitions and numerical building blocks for exact
diagonalization, variational Monte Carlo, TDVP, stochastic series expansion,
symmetry analysis, homogeneous and disordered lattice models, observables, and
future simulation programs.

The project is at the convention-first stage. Its first public contract is the
[scientific convention specification](docs/conventions.md). Implementations
must be tested against that specification rather than inheriting implicit
conventions from an earlier program.

The current dependency direction and subsystem boundaries are documented in
the [architecture overview](docs/architecture/README.md). Consequential design
choices are recorded as [architectural decision records](docs/decisions/README.md).
The qslib 1.0 release criteria are recorded in the project's internal
execution plan; the public documentation describes the supported behavior and
quality guarantees.
The [documentation index](docs/index.md) collects the physicist-first guides
and the locally generated Rust API reference entry point.

The current Cargo workspace exposes the `qslib` facade, while capability crates
isolate core, exact, variational, SSE, IO, CLI, Python, and test-support
responsibilities. The Rust backends implement checked basis and geometry,
heterogeneous model assembly, exact solvers and observables, TDVP, SSE, and
versioned scientific artifacts. The [CLI guide](docs/cli.md) and
[Python guide](docs/python.md) describe the current user-facing boundaries.

## Workspace and features

The default facade contains only `qslib-core`. Optional additive features are
`exact`, `variational`, `sse`, and `io`; `full` enables all four Rust library
capabilities. CLI and Python are separate packages, so core users do not compile
interface or heavy-backend dependencies. Cargo package names use the
`qslib-quantum-` prefix, while Rust library targets use concise names such as
`qslib_core` and `qslib_exact`.

The language-neutral fixtures under
[`fixtures/conformance/v1/`](fixtures/conformance/v1/README.md) record the
small analytic systems that later implementations must reproduce. They include
explicit basis order, resolved coefficients, matrix layout, oracle provenance,
comparison policy, and a checksummed manifest. Passing the current harness
means the evidence is internally valid. It does not claim that scientific
algorithms have already been implemented.

The local CI contract is defined in
[`.github/workflows/ci.yml`](.github/workflows/ci.yml). Remote workflow
execution remains pending owner-authorized push activity.

For a physicist-first introduction to site order, boundaries, pair-dependent
couplings, and disorder provenance, see the
[geometry and interactions guide](docs/geometry-interactions.md).
For operator action, Hamiltonian signs, basis rotations, and model constructors,
see the [operators and models guide](docs/operators-models.md).
For gather-direction site permutations, lattice groups, spin inversion, and
projection, see the [symmetry guide](docs/symmetry.md).
For command-line configuration and physically labelled output, see the
[CLI guide](docs/cli.md).

## Project policies

- [Scientific conventions](docs/conventions.md)
- [Architecture](docs/architecture/README.md)
- [Architectural decisions](docs/decisions/README.md)
- [Rust toolchain and MSRV](docs/toolchain-policy.md)
- [Contribution policy](docs/contribution-policy.md)
- [Security policy](docs/security-policy.md)
- [Changelog](CHANGELOG.md)
- [Local release notes](RELEASE_NOTES.md)
- [Apache-2.0 license](LICENSE)

External publication is disabled. The `qslib` crates.io and PyPI distribution
names belong to an unrelated project. The approved Rust and Python distribution
name is `qslib-quantum`; Rust imports use `qslib`, while Python imports use the
collision-safe `qslib_quantum`. Naming does not authorize publication. See
[ADR-0009](docs/decisions/0009-registry-package-names.md).

The project does not currently accept unsolicited external contributions. See
the [contribution policy](docs/contribution-policy.md) for the current status.
