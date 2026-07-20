---
name: qslib-audit-model
description: Audit a qslib quantum model or Hamiltonian for scientific semantics, weighted couplings, basis action, Hermiticity, conserved quantities, observables, serialization, and API correctness. Use when adding or reviewing a model family, interaction graph, local-energy implementation, disordered coupling scheme, basis transform, or model compatibility adapter.
---

# Audit a qslib model

Perform an evidence-based physics and API audit. Treat the repository's
`docs/conventions.md` as normative and follow the nearest `AGENTS.md`.

## Workflow

1. Locate the model specification, constructors, operator action, observables,
   serialization, tests, and relevant adapters.
2. Write the implemented physical Hamiltonian explicitly. Define every sign,
   unit, physical axis, simulation basis, boundary condition, and normalization.
3. Compare that definition with the relevant convention sections. Report a
   convention conflict rather than silently choosing one interpretation.
4. Resolve every uniform, shell, matrix, or generated coupling into conceptual
   per-interaction coefficients. Verify endpoint ordering, multiplicity,
   duplicate behavior, signed values, zero terms, and disorder provenance.
5. Derive diagonal and connected-state matrix elements independently. Check
   Hermiticity, constant energy terms, basis rotations, and double counting.
6. Identify exact and conditional invariants: conserved sectors, lattice
   symmetries, spin inversion, realness, positivity assumptions, and algorithmic
   sign restrictions.
7. Verify that observables and outputs name physical axes and distinguish totals,
   densities, raw correlations, connected correlations, and disorder averages.
8. Build or require the smallest falsifying tests. Prefer analytic one-bond
   cases, explicit small matrices, basis-rotation spectra, and heterogeneous
   coupling vectors.
9. Run applicable format, lint, test, and documentation checks when changes are
   requested.

## Required audit output

Lead with concrete findings ordered by scientific severity. For each finding,
give the affected location, physical consequence, derivation or test evidence,
and the smallest robust correction. Then list verified properties, missing
evidence, and residual risks.

Do not approve a model merely because tests pass. Check that tests encode an
independently derived physical result. Do not call bond disorder frustration
without showing incompatible preferred constraints.
