# Contributing to qslib

qslib does not currently accept unsolicited external contributions. Pull
requests and contribution proposals may be closed without review. Opening the
project requires a later owner decision, a private conduct-reporting channel,
and an adopted code of conduct.

The current owner decision is recorded in
[`docs/governance/code-of-conduct-decision.md`](docs/governance/code-of-conduct-decision.md).

The remaining sections define the scientific and review standard used by the
owner and authorized agents.

## Start from the scientific contract

Read `AGENTS.md`, `docs/conventions.md`, the accepted records in
`docs/decisions/`, and the relevant architecture document. A change that alters
a sign, basis, bit meaning, site order, normalization, unit, boundary policy,
or serialized meaning must update the normative specification deliberately. It
must never hide the change in implementation code.

Use an ExecPlan governed by `PLANS.md` for work spanning multiple subsystems,
public models, migrations, or more than one focused session.

## Test-first workflow

For supported behavior:

1. state the physical and public API contract;
2. write an independent failing acceptance or conformance test;
3. run it and confirm the failure is for the intended missing behavior;
4. implement the smallest complete behavior;
5. run focused and affected suites; and
6. refactor only while the suite remains green.

Expected values must come from analytic calculations, explicit small matrices,
exact enumeration, manufactured numerical problems, or neutral fixtures. Do
not calculate an expected value by calling the production path under test.

## Documentation

Write for physicists first. Define symbols, signs, axes, basis, units,
normalization, boundaries, limitations, and failure modes before explaining
Rust types. Every public model or algorithm needs a small worked example.

## Required local checks

Once the workspace exists, run:

```text
cargo fmt --check
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo test --workspace --all-features
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --all-features
```

Run the applicable model audit and convention-conformance workflows described
in `AGENTS.md`. Never weaken a reference value, tolerance, lint, or test merely
to make a change pass.

## Changes and reviews

Keep changes focused and preserve unrelated work. Explain the physical meaning,
test oracle, numerical tolerances, compatibility impact, and user-visible
documentation. Consequential or difficult-to-reverse choices require an ADR.

All contributions are licensed under Apache-2.0. By submitting a contribution,
you represent that you have the right to do so under that license.
