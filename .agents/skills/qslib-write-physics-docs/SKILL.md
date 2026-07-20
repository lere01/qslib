---
name: qslib-write-physics-docs
description: Create or revise qslib crate docs, Rust API docs, guides, examples, and CLI help for a physicist-first audience that may not know Rust. Use when documenting models, observables, solvers, configuration, errors, conventions, examples, tutorials, or public releases.
---

# Write qslib physics documentation

Explain what physical calculation the user is requesting before explaining Rust
types or implementation machinery. Follow the nearest `AGENTS.md` and treat
`docs/conventions.md` as normative.

## Workflow

1. Identify the audience and task: choosing a model, constructing a system,
   running an algorithm, interpreting output, extending Rust code, or diagnosing
   a failure.
2. Inspect the implementation, tests, and applicable convention sections. Do not
   document intended behavior that the code does not implement without labeling
   it as planned.
3. Start with the physical object or result. Give the defining equation when it
   removes ambiguity.
4. Define every symbol and state:
   - units and natural-unit assumptions;
   - sign conventions and physical axes;
   - simulation basis and state encoding;
   - geometry, boundaries, and interaction multiplicity;
   - coupling scope, including per-pair values;
   - observable normalization and uncertainty meaning.
5. Translate those concepts into the smallest public API or CLI example that
   answers the user's task. Prefer descriptive domain names over Rust jargon.
6. Explain returned quantities, errors, numerical approximations, convergence
   diagnostics, and scientifically invalid parameter combinations.
7. Add links from overview material to detailed API documentation rather than
   duplicating long reference text.
8. Compile doctests or run documentation checks when the relevant crate exists.

## Public item checklist

For a public model, algorithm, observable, configuration type, or command,
include as applicable:

- one-sentence physical purpose;
- mathematical definition;
- parameter meanings, units, signs, defaults, and allowed ranges;
- basis, indexing, boundary, coupling, and normalization conventions;
- a minimal realistic example;
- output interpretation;
- errors, panics, safety constraints, and numerical limitations;
- links to related higher-level guides and lower-level Rust APIs.

Write for a physicist first, but do not hide information needed by a Rust
integrator. Use precise language and short examples. Avoid unexplained trait,
ownership, lifetime, and generic-type terminology in introductory material.
