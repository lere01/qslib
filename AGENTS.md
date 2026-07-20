# qslib agent instructions

## Project purpose

qslib is a long-lived Rust library for quantum simulation. Optimize technical
decisions for scientific correctness, clear physical meaning, robustness,
scalability, and maintainability. The primary audience includes physicists who
may not write Rust regularly.

## Sources of authority

Apply project decisions in this order:

1. `docs/conventions.md` for normative scientific and data conventions.
2. Accepted records in `docs/decisions/` for architectural decisions.
3. `docs/architecture/` for the current design and dependency boundaries.
4. Public API documentation and implementation.

Do not silently resolve contradictions between these sources. Identify the
conflict and update the authoritative source as part of an approved change.
Compatibility behavior must be explicit and must never redefine a canonical
qslib convention.

## Scientific requirements

- Separate physical definitions from representations and algorithms.
- Represent every resolved Hamiltonian term with its own coefficient. Treat a
  global or shell coupling as constructor shorthand, not an evaluator
  assumption.
- Keep physical axes, simulation bases, bit states, occupation values, and
  legacy representations distinct in types and names.
- Use canonical row-major site order and little-endian packed states unless an
  explicitly named adapter performs a conversion.
- Distinguish total quantities, intensive quantities, raw correlations, and
  connected correlations in public names.
- Preserve algorithmic shifts, regularizers, gauges, cutoffs, and tolerances as
  numerical metadata. Do not fold them silently into the physical model.
- Store realized per-interaction couplings for disordered models. Do not rely
  on a seed alone to reproduce a disorder realization.
- Reject invalid or ambiguous scientific input with contextual errors.

## Architecture and API rules

- Keep foundational types independent of solver and model implementations.
- Prefer checked domain types such as `SiteId`, `BasisState`, `Boundary`, and
  `WeightedInteraction` over primitive values with implicit meaning.
- Keep model construction, operator action, observables, and numerical solvers
  in separate modules with one-way dependencies.
- Put ncli, legacy SSE, file-format, and foreign-array conversions in explicit
  adapter modules.
- Avoid public APIs that expose storage layout unless layout is the purpose of
  the type.
- Return structured errors. Do not panic for invalid user input in library
  code.
- Introduce dependencies only when their maintenance, numerical behavior, and
  license are suitable for a foundational library.

## Documentation rules

- Lead with the physical meaning before Rust implementation details.
- Define every symbol, unit, sign, basis, normalization, and boundary policy
  needed to interpret a result.
- Include a small worked physical example for each public model or algorithm.
- Give public items useful Rust documentation and keep examples compatible with
  doctests where practical.
- Explain failure modes and numerical limitations without assuming familiarity
  with Rust ownership or trait terminology.

Use `$qslib-write-physics-docs` for substantial public documentation work.

## Specification-driven test-first development

Use test-driven development for supported qslib behavior. Derive tests from the
scientific specification, an analytic calculation, or an independent reference
method before writing the production implementation. Do not derive expected
results from the code under test.

Follow this cycle for new behavior:

1. State the physical and public API contract, including invariants, units,
   conventions, error behavior, and numerical guarantees.
2. Write comprehensive acceptance and conformance tests for that contract.
3. Run the new tests and confirm that they fail for the intended reason. A
   missing API that fails to compile is a valid initial failure.
4. Implement the smallest complete behavior that makes the tests pass.
5. Run the focused tests and then the complete relevant suite.
6. Refactor only while the suite remains green.

Write a failing reproduction test before fixing a defect. For existing untested
code that must be changed, add characterization tests and independently derived
correctness tests before modifying behavior.

Tests should protect public behavior and scientific invariants rather than
incidental private structure. Include analytic limits, explicit small matrices,
heterogeneous couplings, invalid inputs, serialization, and representation
parity as applicable. Numerical tests must use quantity-specific tolerances.
Stochastic tests must use reproducible streams, report their statistical
criterion, and avoid fragile exact-sample assertions. Separate stochastic,
floating-point, approximation, and finite-disorder-ensemble uncertainty.

Exploratory prototypes may temporarily precede tests only when clearly marked
unsupported and isolated from the public API. They must acquire a test-first
contract before being merged into supported production code.

## Verification

Changes are complete only after verification proportional to their scientific
risk. When the corresponding code exists, the standard checks are:

```text
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-targets --all-features
RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --all-features
```

Model and convention changes must also include small independently derived
reference tests. Prefer full matrices, exact enumeration, analytic limits, and
cross-representation parity at sizes where those checks are inexpensive.
Never update an expected value merely to match the implementation.

Use `$qslib-audit-model` for model changes and `$qslib-check-conformance` for
convention, backend, serialization, or compatibility changes.

## Decisions and scope

Record durable, consequential architectural choices in `docs/decisions/`.
Keep personal preferences and temporary work notes out of normative documents.
Do not use generated agent memory as the sole source for a project requirement.

For changes spanning multiple subsystems, scientific conventions, migrations,
or more than one focused working session, use an ExecPlan governed by
`PLANS.md`. The active qslib 1.0 plan is `docs/plans/qslib-v1.md`. Read the
entire governing plan before implementation, keep its progress and decision
sections current, and continue to the next milestone without requesting routine
next steps. Stop only at the owner gates defined by the plan or when completion
criteria pass.

The `qslib-architect` custom agent is the preferred reviewer for changes that
cross scientific domains or module boundaries. Use parallel subagents only for
independent, bounded work such as read-only exploration, test execution, or
separate audits. Avoid concurrent edits to overlapping files.
