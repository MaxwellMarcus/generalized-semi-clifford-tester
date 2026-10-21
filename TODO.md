# Running TODO — Generalized Semi-Clifford Tester

Last updated: 2026-10-21

This is the operational queue for incremental development. Keep
`docs/status-and-roadmap.md` as the higher-level scientific roadmap. Each
completed change should be small enough to test, document, commit, and push as
one coherent unit.

## Current focus

- [ ] **Define a bounded fourth two-qutrit rejection prefix.**
  Measure the same exact non-Clifford fixture through four canonical input
  Lagrangians with explicit coefficient and candidate-pair budgets, preserving
  `UNKNOWN` unless the complete search establishes a classification.

## Next

- [ ] Add a checked fourth-prefix record only if the measured host cost remains
  suitable for the example suite; otherwise document it without turning
  machine-dependent timing into a pass/fail invariant.

## Maintenance

- [ ] Keep README examples executable and synchronized with the public API.
- [ ] Keep package version and `CITATION.cff` version synchronized for releases.
- [ ] Require the full test suite, Ruff, example execution, and package build
  before pushing an automated change.
- [ ] Record measured performance claims rather than estimated ones.
- [ ] Preserve `UNKNOWN` or candidate-style outcomes whenever a search or
  statistical conclusion is incomplete.

## Completed

- [x] 2026-10-21 — Measure and check in the three-input exact two-qutrit
  non-Clifford rejection prefix: 486 coefficients and 120 candidate pairs
  reject all three canonical inputs before the cap returns global `UNKNOWN`,
  with host timing retained as observational metadata only.

- [x] 2026-10-20 — Measure and check in the two-input exact two-qutrit
  non-Clifford rejection prefix: 324 coefficients and 80 candidate pairs
  reject both canonical inputs before the cap returns global `UNKNOWN`, with
  host timing retained as observational metadata only.

- [x] 2026-10-19 — Add an exact two-qutrit non-Clifford rejection-prefix
  record that computes one input Lagrangian's two expansions, rejects all 40
  output candidates, and returns `UNKNOWN` at the next coefficient boundary
  rather than promoting the bounded prefix to a global negative result.

- [x] 2026-10-18 — Add a bounded exact two-qutrit SUM positive record with a
  324-coefficient budget for two search expansions plus independent witness
  verification and a one-pair prefix that finds the deterministic first
  candidate; retain `UNKNOWN` semantics for any prefix that does not find a
  witness and keep host timing observational.

- [x] 2026-10-16 — Cache each exact field's immutable primitive root after the
  81-coefficient profile identified repeated SymPy reconstruction as a material
  hotspot; retain the analytic SUM support and `UNKNOWN` result while reducing
  the same `max_coefficients=161` measurement to 0.22 seconds and 203,686 bytes
  of Python peak memory on the development host.

- [x] 2026-10-15 — Replace dense output-Weyl adjoints and matrix products in
  odd-prime coefficient extraction with the same exact trace inner product on
  the Weyl matrix's monomial support, independently cross-check it against the
  dense formula, and re-measure the first 81-coefficient two-qutrit expansion.

- [x] 2026-10-14 — Add the exact two-qutrit SUM permutation to the checked
  qutrit benchmark suite with an 80-coefficient pre-expansion cap, preserving
  `UNKNOWN` and asserting only portable exact-work counters and semantics.

- [x] 2026-10-13 — Specify the exact two-qutrit SUM permutation fixture and an
  80-coefficient checked-run budget that finishes before the first indivisible
  81-coefficient expansion, preserving a measured, reproducible `UNKNOWN`
  scaling baseline without implying a completed classification.

- [x] 2026-10-12 — Add an exact odd-prime CLI accepting versioned rational
  cyclotomic JSON and emitting the existing exact result schema unchanged,
  preserving resource-capped `UNKNOWN` and exact-only work counters.

- [x] 2026-10-09 — Benchmark qutrit positive, completed-negative, and
  candidate-capped `UNKNOWN` exact searches in a versioned suite, keeping
  portable exact work counters separate from host runtime and Python memory.

- [x] 2026-10-08 — Add bounded odd-prime input/output Lagrangian search with
  a versioned result schema, exact completed negatives, independently
  reverified positive witnesses, and qudit/coefficient/candidate-pair caps
  that preserve `UNKNOWN`.

- [x] 2026-10-07 — Add exact odd-prime Weyl matrices and fixed-witness
  verification with independently reconstructed finite-field spans, a separate
  verified/rejected/`UNKNOWN` result, explicit qudit and coefficient caps, and
  an exact qutrit Fourier regression and example.

- [x] 2026-10-06 — Enumerate canonical odd-prime Lagrangians through two
  qudits for `p <= 5`, with the four one-qutrit lines and all 40 two-qutrit
  planes independently cross-checked for count, uniqueness, span size, and
  isotropy.

- [x] 2026-10-05 — Fix the odd-prime `(x | z)` Weyl phase and commutation
  convention, add exact finite-field symplectic primitives, and exhaustively
  cross-check the multiplication laws against all one-qutrit matrix pairs.

- [x] 2026-10-04 — Benchmark positive bounded exact searches through three
  qubits plus a completed one-qubit negative and capped three-qubit `UNKNOWN`,
  keeping exact work counts separate from host runtime and Python peak memory.

- [x] 2026-10-03 — Fix exact witness revalidation to canonicalize valid
  enumerated multi-qubit bases while still requiring their stored elements to
  equal the independently reconstructed span; add a two-qubit CNOT regression.

- [x] 2026-10-02 — Add bounded exhaustive exact Pauli-MASA search with qubit,
  coefficient, and candidate-pair caps; preserve `UNKNOWN` on interruption;
  independently revalidate positive witnesses; and export a versioned JSON-
  serializable result schema.

- [x] 2026-10-01 — Compute complete exact Pauli-conjugation expansions over
  the declared cyclotomic field and independently revalidate supplied
  input/output Lagrangians, with a separate no-tolerance fixed-witness result
  that makes no exhaustive-search claim.

- [x] 2026-09-30 — Implement canonical cyclotomic scalars and matrices over a
  declared `Q(zeta_m)`, common-field embedding, exact adjoints, products,
  tensors and unitarity, plus Clifford+T, permutation, and controlled-phase
  constructors behind an optional SymPy dependency.

- [x] 2026-09-29 — Specify a cyclotomic `Q(zeta_m)` exact-backend domain,
  supported and rejected gate inputs, representation tradeoffs, a separate
  exact result contract, exact/`UNKNOWN` claim boundaries, and independently
  testable implementation milestones.

- [x] 2026-09-28 — Formalize the fixed-witness tolerant-testing contract with
  explicit sampling hypotheses, a familywise soundness proof, a
  margin-dependent finite-sample completeness bound, and the exact zero-leakage
  endpoint kept separate from finite-shot certification.

- [x] 2026-09-27 — Add a machine-readable CLI for dense and circuit tests,
  with versioned deterministic JSON preserving numerical tolerances,
  `UNKNOWN` outcomes, raw finite-shot evidence, and fresh-verification bounds.

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
