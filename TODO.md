# Running TODO — Generalized Semi-Clifford Tester

Last updated: 2026-09-16

This is the operational queue for incremental development. Keep
`docs/status-and-roadmap.md` as the higher-level scientific roadmap. Each
completed change should be small enough to test, document, commit, and push as
one coherent unit.

## Current focus

- [ ] **Add unit tests for every validation branch in `statistics.py`.** Cover
  all-leakage observations and invalid confidence/shot-planning arguments with
  focused deterministic cases.

## Next

- [ ] Add benchmark tooling that records qubits, Lagrangian count, circuits,
  shots, candidate pairs, classical runtime, and peak memory.
- [ ] Add analytically justified fixtures for GSC, semi-Clifford,
  GSC-but-not-semi-Clifford, and deliberately wrong witnesses.
- [ ] Replace the noisy all-pairs output-Lagrangian fallback with a
  support-aware candidate search and compare it against exhaustive scoring.
- [ ] Add a sparse Pauli-transfer representation with dense cross-checks for
  small qubit counts.
- [ ] Demonstrate transpilation and execution with a documented Qiskit noise
  model, including the effect on leakage confidence bounds.
- [ ] Add a CLI that emits machine-readable JSON results without overstating
  finite-shot conclusions.
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
