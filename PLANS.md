# qslib execution plans

An execution plan, or ExecPlan, is a living implementation specification that
an agent or a new contributor can follow without relying on chat history or
generated memory. Use an ExecPlan for every qslib change that spans more than
one subsystem, changes a scientific convention, introduces a public model or
algorithm, performs a migration, or is expected to take more than one focused
working session.

The active qslib 1.0 plan is `docs/plans/qslib-v1.md`.

## Executor contract

Read the applicable `AGENTS.md`, this file, the entire active ExecPlan, and the
scientific convention sections named by the plan before editing production
code. Treat the ExecPlan as the durable state of the work. Do not assume access
to an earlier conversation.

Proceed from one milestone to the next without asking the user for routine next
steps. Resolve reversible technical ambiguity using the priorities and defaults
in the plan, then record the decision and evidence. Stop for user input only
when work requires new external authority, a destructive or irreversible
operation, a change to the accepted scientific contract, unavailable private
credentials, or a product choice that the plan explicitly reserves for the
owner.

Keep one primary agent responsible for overlapping edits. Read-only exploration,
independent physics audits, test execution, and documentation review may be
delegated. Parallel workers must not modify the same files or share an
uncommitted working directory unless the plan explicitly coordinates them.

Do not publish crates, wheels, binaries, documentation, tags, releases, or
remote branches unless the user explicitly authorizes that external action.
Producing verified local release artifacts is allowed when required by the
plan.

## Test-first execution

For each supported behavior, state the contract and write tests before the
production implementation. Run the new tests and confirm that they fail for the
intended reason. Then implement the smallest complete behavior, run focused
tests, run the affected suite, and refactor only while tests remain green.

Expected values must come from conventions, analytic calculations, independent
reference code, or neutral fixtures. A defect fix starts with a failing
reproduction. Numerical tolerances belong to individual quantities. Stochastic
tests must state their statistical criterion and must not depend on fragile
exact sample sequences unless the sequence itself is the public contract.

An exploratory spike may precede tests only when the plan labels it as a spike,
keeps it outside the supported public API, and states promotion or deletion
criteria.

## Living-document requirements

Update the active ExecPlan at every meaningful stopping point. It must always be
possible for a fresh agent to determine what is complete, what failed, why a
decision was made, which commands were run, and what should happen next.

Every ExecPlan must contain these sections:

- `Purpose and user-visible outcome`
- `Progress`
- `Surprises and discoveries`
- `Decision log`
- `Outcomes and retrospective`
- `Context and orientation`
- `Scope and non-goals`
- `Architecture and interfaces`
- `Milestones`
- `Concrete commands`
- `Validation and acceptance`
- `Idempotence and recovery`
- `External authority and owner gates`

Only `Progress` uses checkboxes. Give completed entries a UTC timestamp. Record
short evidence in `Surprises and discoveries`, not just conclusions. Record the
rationale and date for every decision that changes the implementation route.
At each phase gate, update `Outcomes and retrospective` with what users can now
do and what remains.

## Milestone quality gates

A milestone is complete only when its promised behavior is observable and all
of its acceptance tests pass. Compilation alone is insufficient. Each milestone
must leave the workspace buildable, documented, and safer to resume than it was
at the start.

Before marking a milestone complete:

1. Run formatting, lint, unit, integration, doctest, and relevant parity tests.
2. Run `$qslib-audit-model` for changed model semantics.
3. Run `$qslib-check-conformance` for convention-sensitive behavior.
4. Run `$qslib-write-physics-docs` for new public surfaces.
5. Review the diff for accidental API expansion, hidden convention changes,
   unchecked numerical assumptions, and unrelated user files.
6. Update the ExecPlan with commands and concise results.

Commit only files within the milestone scope. Never stage unrelated changes in
a dirty parent repository. A milestone commit is recommended when the project
has a clean dedicated repository, but a commit is not a substitute for the
tests and evidence recorded in the plan.

## Recovery and blockers

Prefer additive, idempotent changes. Preserve a working legacy path while a new
adapter or implementation is being validated. When a command fails, diagnose
the cause, record material evidence, and try a safe alternative. Do not weaken a
test, tolerance, lint, or scientific requirement merely to advance progress.

Mark work blocked only after safe in-scope alternatives have been exhausted and
the remaining condition needs user authority or an external state change. The
blocker report must name the exact decision or permission required and the
verified state from which work can resume.

## Completion

An ExecPlan is complete only when every promised user-visible behavior and
quality gate passes, local release artifacts can be built from a clean checkout,
the public documentation describes the implemented behavior, and no required
work remains. Move optional future ideas to a separately named roadmap rather
than leaving the active plan indefinitely incomplete.
