# ADR-0008: Dedicated repository ownership boundary

- Status: accepted
- Date: 2026-07-19
- Owners: qslib owner

## Context

Before this decision, `qslib/` resolved to the parent ncli Git root. The parent
worktree contained extensive unrelated modifications and untracked simulation
outputs.
qslib has its own version, CI, release artifacts, issue history, Python wheels,
Cargo packages, agent instructions, and consumers. Milestone commits cannot be
isolated safely while qslib remains an untracked subtree.

## Decision

Make `qslib/` a dedicated Git repository before Milestone 1. The owner approved
this boundary on 2026-07-19, and the repository was initialized with default
branch `main`. Project agent configuration is present inside the repository as
`.codex/agents/qslib-architect.toml` and `.agents/skills/qslib-*`. Configure its
`origin` remote as `https://github.com/lere01/qslib.git`. The owner approved
local commits but prohibited push and pull, so configuring the remote does not
authorize network contact.

The parent ncli repository should consume released qslib packages or an
explicitly configured local path during development. Removing the current
untracked parent entries or rewriting parent history is not part of this ADR.

## Alternatives considered

- Keeping qslib as an ncli subtree through 1.0 avoids a repository setup step,
  but entangles versioning and makes clean milestone commits risky in the
  current dirty worktree.
- A Cargo workspace rooted at ncli was rejected because qslib has independent
  consumers and release policy.
- Git submodules or subtree history may become a consumer integration choice,
  but neither is required to establish qslib ownership.

## Consequences

The repository root, agent assets, CI, tags, and release history become
unambiguous. The nested repository now provides independent ownership metadata
and a configured `origin`. It has no commit at the time this ADR is accepted
because the initial commit is the final Milestone 0 action after clean-tree
validation. The later ncli integration mechanism remains a separate decision.

## Validation

- `git rev-parse --show-toplevel` from any qslib directory returns the qslib
  root.
- A clean checkout contains all governing docs and agent assets.
- `git status --short` contains only qslib-scoped changes.
- No ncli file is staged, moved, or removed by the transition.
