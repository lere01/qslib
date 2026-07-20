---
name: qslib-check-conformance
description: Verify qslib code, backends, serialized artifacts, and compatibility adapters against the normative scientific conventions through specification-driven test-first development and required reference vectors. Use when planning tests for new qslib behavior, performing parity testing, changing basis or site-order conversions, comparing backends, changing schemas, implementing legacy ncli or SSE interoperability, fixing regressions, or changing any convention-sensitive behavior.
---

# Check qslib conformance

Establish that two implementations or representations describe the same
physical object before comparing their numerical output.

## Test-first sequence

1. State the claim and applicable scientific or API contract before changing
   production code.
2. Select an expected result derived from the specification, an analytic
   calculation, or an independent reference implementation.
3. Encode comprehensive acceptance and conformance tests, including invalid and
   boundary cases appropriate to the claim.
4. Run the tests and verify that they fail for the intended reason. Preserve the
   failing reproduction for a defect.
5. Implement the smallest complete change that satisfies the contract.
6. Run focused tests, the relevant conformance vectors, and the complete
   affected suite before refactoring.

## Workflow

1. Read the nearest `AGENTS.md`, `docs/conventions.md` section 24, and every
   convention section governing the feature under test.
2. Record each side's site order, bit convention, physical and simulation axes,
   boundary and multiplicity policy, resolved coefficients, dtype, and output
   normalization.
3. Apply explicit adapters before comparison. Never reinterpret a flat state,
   spin enum, or coupling table under a different convention.
4. Select an independent oracle appropriate to the claim:
   - analytic one-site or one-bond calculation;
   - explicit small Hamiltonian matrix;
   - exact basis enumeration;
   - unitary basis transformation;
   - a separately implemented reference backend.
5. Test structure before floating-point values: dimensions, canonical ordering,
   sparsity pattern, Hermiticity, conserved sectors, and interaction identity.
6. Test matrix elements, spectra or invariant subspaces, local energies,
   observables, and serialized round trips as applicable.
7. Use tolerances tied to the tested quantity and dtype. Report absolute and
   relative error where scale varies. Do not use one unexplained global epsilon.
8. Run deterministic tests with complete provenance. For disorder, persist and
   compare the realized coupling table rather than regenerating it from a seed.
9. Run the repository's format, lint, test, and documentation checks when code
   changes are in scope.

## Integrity rules

- Do not derive expected values from the implementation being tested.
- Do not update a golden result solely because the implementation changed.
- Do not write production behavior before its contract and failing test unless
  the work is an isolated, explicitly unsupported exploration.
- Do not compare individual eigenvectors inside a degenerate eigenspace. Compare
  eigenvalues and invariant-subspace projectors.
- Separate stochastic error, numerical approximation error, and finite disorder
  ensemble uncertainty.
- Treat any unexplained convention conversion as a failed check even if selected
  scalar outputs happen to agree.

Report the claim tested, oracle, resolved conventions, commands, tolerances,
results, and any untested surface.
