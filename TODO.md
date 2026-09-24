# Running TODO — Generalized Semi-Clifford Tester

Last updated: 2026-09-24

This is the operational queue for incremental development. Keep
`docs/status-and-roadmap.md` as the higher-level scientific roadmap. Each
completed change should be small enough to test, document, commit, and push as
one coherent unit.

## Current focus

- [ ] **Add a machine-readable command-line interface.**
  Emit dense and circuit-test results as versioned JSON while preserving
  finite-shot evidence, numerical tolerances, and `UNKNOWN` states.

## Next

- [ ] Investigate an exact algebraic backend and write down the supported gate
  domain before implementation.
- [ ] Investigate the prime-qudit generalization only after the qubit API and
  correctness contract stabilize.

## Maintenance

- [ ] Keep README examples executable and synchronized with the public API.
- [ ] Keep package version and `CITATION.cff` version synchronized for releases.
- [ ] Require the full test suite, Ruff, example execution, and package build
  before pushing an automated change.
- [ ] Record measured performance claims rather than estimated ones.
- [ ] Preserve `UNKNOWN` or candidate-style outcomes whenever a search or
  statistical conclusion is incomplete.

## Completed

- [x] 2026-09-24 — Compare sparse-transfer leakage intervals with ideal and
  Aer depolarizing Bell sampling, keeping empirical leakage and simultaneous
  finite-shot confidence bounds distinct from exact claims.

- [x] 2026-09-23 — Add exhaustive small-qubit Pauli-transfer extraction with a
  thresholded sparse representation, explicit arithmetic/tolerance and
  discarded-mass metadata, and dense one- and two-qubit cross-checks.
- [x] 2026-09-22 — Replace positive noisy all-pairs output-Lagrangian scoring
  with a threshold-complete observed-support index, cross-checked against
  exhaustive scoring while preserving exact negative-result diagnostics.
- [x] 2026-09-21 — Add a minimal three-qubit GSC-but-not-semi-Clifford basis
  permutation with an exact bit-level obstruction, a minimality argument, and
  a deliberately rejected semi-Clifford classification.
- [x] 2026-09-20 — Add parameterized one- and two-qubit semi-Clifford/GSC
  fixture families, a provable one-qubit non-GSC rotation, and deliberately
  wrong Pauli-MASA witnesses.
- [x] 2026-09-19 — Separate exact benchmark workload metrics from host
  measurements, record source-circuit depth and operation counts, and check in
  reproducible one- and two-qubit identity baselines.
- [x] 2026-09-18 — Add a reproducible one-qubit benchmark runner and versioned
  JSON schema for Lagrangians, circuits, shots, candidate pairs, runtime, and
  Python peak memory.
- [x] 2026-09-17 — Cover every validation branch in `statistics.py`, including
  all-leakage observations and invalid confidence and shot-planning inputs.
- [x] 2026-09-16 — Add discovery followed by fresh fixed-witness verification,
  including separate shot costs and explicit acceptance/rejection results.
- [x] 2026-09-15 — Add fixed-witness Clopper–Pearson leakage bounds with a
  familywise correction across generators.
- [x] 2026-09-15 — Add zero-event shot planning for a target leakage and
  confidence level.
- [x] 2026-09-15 — Add Qiskit GSC witness testing and small-qubit discovery.
- [x] 2026-09-15 — Extend canonical Lagrangian discovery through four qubits.

## Rules for automated maintenance

1. Work from the first unchecked, currently feasible item unless a failing
   build or correctness bug has higher priority.
2. Split work that cannot be completed safely in one run into a concrete next
   milestone and update this file.
3. Move finished work to `Completed` with the date and a concise result.
4. Do not mark theoretical claims complete merely because an implementation
   exists; record the proof or precise empirical scope.
5. Do not push unless all relevant validation passes.
