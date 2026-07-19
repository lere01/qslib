# qslib architecture

This directory describes qslib's current architecture. Scientific meanings are
defined by [`../conventions.md`](../conventions.md), while consequential design
choices and their rationale are recorded in [`../decisions/`](../decisions/).

## Architectural direction

qslib should develop as a layered workspace whose dependency direction follows
physical foundations toward algorithms and applications:

```text
foundations
  geometry, basis, scalar policy, identifiers
      |
operators and interactions
      |
physical models and observables
      |
state representations and numerical kernels
      |
exact, variational, TDVP, thermal, and SSE algorithms
      |
adapters, command-line programs, and external interfaces
```

The dependency direction will be implemented in Milestone 1 through the
workspace topology approved in
[ADR-0001](../decisions/0001-workspace-boundaries.md):

```text
qslib-core
  ^       ^              ^
  |       |              |
exact  variational      sse
  ^       ^              ^
  +-------+-------+------+
                  |
               facade

qslib-io converts public domain and result values to versioned DTOs.
qslib-cli and qslib-python consume coarse public APIs.
qslib-test-support is never a production dependency.
```

Arrows point from a consumer toward its foundational dependency. Algorithm
crates do not depend on IO, CLI, Python, or test support. The facade owns no
second implementation. ADR-0001, ADR-0008, and ADR-0009 approve this topology,
the dedicated repository boundary, and package names. No competing production
structure may be introduced without a superseding decision record.

Cargo uses collision-safe distribution names while Rust code retains concise
imports:

| Responsibility | Cargo package | Rust target |
| --- | --- | --- |
| Facade | `qslib-quantum` | `qslib` |
| Foundations | `qslib-quantum-core` | `qslib_core` |
| Exact methods | `qslib-quantum-exact` | `qslib_exact` |
| Variational methods | `qslib-quantum-variational` | `qslib_variational` |
| Thermal SSE | `qslib-quantum-sse` | `qslib_sse` |
| Scientific IO | `qslib-quantum-io` | `qslib_io` |
| Python boundary | `qslib-quantum-python` | `qslib_quantum` |
| Command line | `qslib-quantum-cli` | `qslib` binary |
| Test support | `qslib-test-support` | `qslib_test_support` |

The facade has no default optional features. Core is always present. The
additive `exact`, `variational`, `sse`, and `io` features expose capability
crates, and `full` enables those four Rust capabilities. CLI and Python remain
separate packages.

## Stable boundaries

- Geometry defines sites, coordinates, boundaries, bonds, and distances. It
  does not know about Hamiltonians or solvers.
- Basis types define state meaning and layout. They do not attach model-specific
  meanings such as Rydberg occupation to a generic spin label.
- Operators and weighted interactions define mathematical action with resolved
  per-term coefficients.
- Models assemble physical operators without depending on how a state is
  sampled, optimized, diagonalized, or thermally expanded.
- Observables state normalization and physical axes independently of the
  estimator used to evaluate them.
- Algorithms consume validated abstractions and expose their approximations and
  numerical diagnostics.
- Adapters own every conversion to legacy ncli, legacy SSE, Python arrays,
  serialized formats, and future accelerator backends.
- qslib-owned core scientific crates forbid unsafe code. Reviewed dependency
  internals and the isolated Python FFI boundary follow their own audited
  policies.

## Development architecture

qslib uses specification-driven test-first development for supported behavior:

```text
scientific convention or accepted decision
                  |
independent oracle and failing acceptance test
                  |
          production implementation
                  |
       refactoring under a green suite
```

The test oracle must be independent of the production path being tested. An
oracle may be an analytic result, explicit small matrix, exact enumeration,
unitary equivalence, manufactured numerical problem, or separately implemented
reference backend. Tests must demonstrate the initial failure for the intended
reason before implementation proceeds.

The test architecture has complementary layers:

- conformance tests encode normative conventions and required reference vectors;
- model tests verify matrix elements, Hermiticity, conserved quantities,
  heterogeneous interactions, and analytic limits;
- property tests verify general invariants across generated valid inputs;
- backend tests compare equivalent representations after explicit conversion;
- numerical tests verify residuals, convergence, and precision-specific error;
- stochastic tests use reproducible streams and robust statistical criteria;
- serialization tests preserve physical meaning and provenance through a round
  trip;
- doctests keep public physicist-first examples executable.

Production code must not depend on test fixtures or a test-only oracle. Shared
reference functionality that is scientifically useful outside testing belongs
in an explicit reference or exact backend with its own public contract.
Regression fixes begin with a failing reproduction. Tests should bind to public
behavior and scientific invariants, leaving private implementation structure
free to evolve.

Neutral conformance data lives under `fixtures/conformance/v1/`. Every fixture
states one physical or representational claim, its resolved conventions, an
independent derivation, authorship and review provenance, and a
quantity-specific comparison policy. Matrices record basis masks, complex
`f64` entries, shape, row-major layout, and the matrix-element definition. A
sorted manifest binds all eight required fixture kinds to their BLAKE3 digests.
The test-support harness validates this evidence without depending on a
production qslib crate.

## Validation architecture

Every optimized path should have a small-system reference route that is simple
enough to inspect independently. Exact enumeration and explicit matrices are
scientific oracles at small sizes, not production scalability strategies.
Cross-backend comparisons must resolve both sides to the same convention,
basis, site order, coefficients, and observable normalization before numerical
comparison.

## Evolution

Update this document when layer responsibilities change. Create an
architectural decision record when a decision is difficult to reverse, affects
multiple layers, or constrains future projects.
