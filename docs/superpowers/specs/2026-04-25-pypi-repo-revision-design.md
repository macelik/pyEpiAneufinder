# pyEpiAneufinder PyPI Repo Revision Design

**Date:** 2026-04-25

## Goal

Prepare the `pyEpiAneufinder` repository for a credible first PyPI release by revising packaging, documentation, tests, and hosted docs workflow without making algorithmic changes.

## Constraints

- No algorithmic changes are allowed in this revision.
- Any code change required for packaging, documentation, or testability must preserve current computational behavior.
- Every intentional repository change in this effort must be documented in user-facing or maintainer-facing docs.
- The documented public API for the first release is limited to:
  - `epiAneufinder`
  - `split_subclones`
  - `karyo_gainloss`
  - `plot_single_cell_profile`
  - `compute_aneuploidy_across_sample`
  - `compute_aneuploidy_by_chr`
  - `compute_heterogeneity_across_sample`
  - `compute_heterogeneity_by_chr`
  - `compute_cnv_burden_cell`

## Scope

This revision is a balanced package-foundation pass. It should make the repository installable, testable, documentable, and ready for hosted documentation, while avoiding release-engineering work that is unnecessary for the first publishable version.

In scope:

- PyPI-facing package metadata cleanup
- separation of runtime, test, and docs dependencies
- Sphinx documentation
- Read the Docs configuration
- CI for install, tests, and docs build
- public API documentation for the agreed functions
- targeted tests around documented public behavior
- explicit documentation of all meaningful repo changes made during this revision

Out of scope:

- algorithm redesign or output changes
- large refactors of internal implementation
- CLI design unless a minimal entry point is already required by the current package
- automated release publishing beyond what is needed to validate packaging locally and in CI

## Current-State Summary

The repository already contains:

- a package directory with `pyproject.toml`
- a top-level `README.md`
- a small existing `tests/` suite
- sample data assets
- exported functions in `pyEpiAneufinder/__init__.py`

The repository does not yet appear to contain:

- a `docs/` tree for Sphinx
- Read the Docs configuration
- CI workflow files for tests and docs
- a formal change log for this packaging revision

There is also a mismatch between the README and the actual current plotting API surface: the README mentions `plot_karyo_annotated`, while the release boundary for this revision is to document `karyo_gainloss` instead.

## Design

### 1. Packaging

The `pyproject.toml` file should be revised to reflect a real distributable Python package.

Expected packaging changes:

- keep build backend simple unless the current layout requires more
- enrich project metadata:
  - project description
  - readme declaration
  - homepage/repository/documentation URLs
  - classifiers
  - license expression compatible with modern packaging metadata
- move non-runtime tooling out of `dependencies`
  - `pytest` must not remain a runtime dependency
  - docs tools must live in a docs extra or equivalent optional dependency group
- verify package discovery/package data rules so the intended Python package installs cleanly
- define supported Python versions explicitly

Packaging changes must not alter algorithm behavior. If any import-path cleanup is needed, it should preserve public behavior and be documented.

### 2. Public API Boundary

The first PyPI release should document and support only the agreed public functions already exposed to users through the README workflow.

Rules for the API boundary:

- only the agreed public functions are documented in the API reference
- everything else is treated as internal implementation unless later promoted deliberately
- the README and docs must stop advertising functions outside this agreed boundary
- the plotting documentation should describe `karyo_gainloss` rather than `plot_karyo_annotated`

This keeps the maintenance burden low and reduces accidental compatibility promises.

### 3. Documentation

Documentation should be built with Sphinx and hosted on Read the Docs.

Recommended docs structure:

- `docs/index.*`
- `docs/installation.*`
- `docs/quickstart.*`
- `docs/input-output.*`
- `docs/api.*`
- `docs/revision-log.*` or `docs/changelog.*`

Documentation responsibilities:

- `index`: package overview, links to main sections, citation/reference links
- `installation`: pip installation, supported Python versions, dependency notes
- `quickstart`: minimal end-to-end example using the documented main workflow
- `input-output`: required inputs, important parameters, output files, caveats
- `api`: reference material for the agreed public functions only
- `revision-log` or `changelog`: record packaging/repo-revision changes made during this effort

The top-level `README.md` should remain concise and act as the front door. Deeper instructions should live in Sphinx pages and be linked from the README.

### 4. Tests

Tests for this revision should protect documented public behavior rather than internal implementation details.

Testing priorities:

- import smoke tests for the public API
- focused tests for documented helper functions where existing behavior is stable
- plotting tests that validate expected outputs or validation behavior without relying on manual visual inspection
- packaging-oriented tests where appropriate, such as ensuring documented imports resolve

Testing constraints:

- do not rewrite algorithms to satisfy new tests
- do not invent new behavior beyond what the current package already intends
- if a function’s current behavior is awkward but acceptable for this release, document it instead of silently changing it

### 5. Automation

Automation should validate the package as a package, not just as local source code.

Required automation:

- GitHub Actions workflow for test execution
- GitHub Actions workflow step or job for Sphinx docs build
- Read the Docs configuration checked into the repository
- install steps that use the repository’s declared dependencies rather than ad hoc environment setup

CI should answer three questions on every change:

1. Does the package install?
2. Do the tests pass?
3. Do the docs build?

### 6. Change Documentation

This revision must document every intentional repo change.

Minimum documentation standard:

- maintain a dedicated revision log or changelog page for this packaging effort
- update relevant user-facing docs whenever setup, imports, or documented workflows change
- keep commit messages specific to the repo revision work

If a code change is necessary for packaging hygiene, docs generation, or test stability, the reason and expected impact should be captured in the revision log.

## Risks and Guardrails

### Risks

- package metadata cleanup may expose unclear dependency boundaries
- current docstrings may be incomplete for API reference generation
- existing tests may not cover the public API evenly
- README examples may drift from actual callable behavior

### Guardrails

- prefer documentation alignment over behavior changes
- treat any behavior-changing code edit as out of scope unless explicitly approved later
- keep the public API small
- document mismatches instead of “fixing” scientific logic during packaging work

## Acceptance Criteria

This design is successful when:

- `pyEpiAneufinder` has PyPI-ready metadata and dependency separation
- the agreed public API is clearly documented
- Sphinx docs build locally and in CI
- Read the Docs can build from the repository configuration
- the test suite covers documented public behavior at a reasonable first-release level
- the README is aligned with the actual public API
- all meaningful repo changes in this revision are documented
- no algorithmic behavior has been intentionally changed

## Implementation Notes for the Follow-Up Plan

The implementation plan should decompose work into small tasks covering:

- packaging metadata and dependency cleanup
- README alignment
- Sphinx scaffolding and page authoring
- API reference generation strategy
- test additions/adjustments focused on public behavior
- GitHub Actions configuration
- Read the Docs configuration
- revision-log/changelog updates

The follow-up plan must preserve the no-algorithm-change constraint in every task.
