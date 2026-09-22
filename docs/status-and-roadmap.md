# Status and roadmap

Last audited: 2026-09-22.

## What works now

| Layer | Capability | Practical boundary |
| --- | --- | --- |
| Binary algebra | Exact symplectic arithmetic and canonical Lagrangian enumeration | 1--4 qubits enumerated: 3, 15, 135, 2,295 Lagrangians |
| Dense GSC tester | Exhaustive input/output Pauli-MASA search with independently checkable witnesses | Defaults to `UNKNOWN` above 3 qubits; complex128 with a stated tolerance |
| Semi-Clifford circuits | Bell-sampling conjugation circuits, exhaustive discovery, and `n`-circuit candidate verification | Finite-shot candidate evidence; exhaustive discovery uses `4**n - 1` circuits |
| GSC circuits | Full sampled Pauli-support decoding, exhaustive discovery, and `n`-circuit fixed-witness verification | Discovery through 4 qubits; support-aware fast path with an exhaustive noisy fallback |
| Statistics | Exact one-sided binomial leakage bounds with familywise confidence for a fixed GSC witness | Valid for a witness fixed independently of the verification samples |
| Benchmark records | Versioned JSON separates exact source-circuit depth/operation counts from wall time and Python peak memory | Checked-in identity workloads cover one and two qubits; depth is pre-transpilation, while runtime and memory remain machine-dependent |
| Analytic fixtures | Parameterized phase and controlled-phase families, a one-qubit non-GSC rotation, and rejected witnesses | GSC-but-not-semi-Clifford still needs a separately justified higher-qubit fixture |
| Conjugation sampling | Short words at arbitrary bounded depth and recursive finite-subgroup escape without previous-group enumeration | Exact-oracle containment theorem counts group operations, not expanded U queries; dense frontend is numerical and limited to 3 qubits |
| SC-distance evidence | Replayable hierarchy paths and an exhaustive unrestricted SC-distance fallback | Restricted and unrestricted targets are distinguished; floating-point margins are not interval certificates |

The key conceptual distinction is now represented in the API:

- `has_candidate_witness` reports that empirical leakage passed the chosen
  threshold;
- `has_confidence_certified_witness` reports that the one-sided simultaneous
  upper bound also passed it;
- exhaustive discovery never claims post-selection confidence from the same
  samples. A discovered witness should be checked with fresh calls to
  `run_gsc_witness_test`.

## Current scaling

- Dense numerical testing grows over every ordered pair of Lagrangians and is
  intended for at most three qubits.
- Circuit GSC discovery constructs 3, 13, 47, and 165 distinct conjugation
  experiments for one through four qubits, respectively.
- A supplied input/output Pauli-MASA witness always needs `n` circuits.
- Exact-support data usually constructs the forced output Lagrangian directly.
  The finite-shot fallback can still score every input/output pair and is the
  main classical bottleneck at four qubits.

## Prioritized to-do list

### P0: make experimental claims defensible

- [x] Add fixed-witness, familywise finite-shot leakage bounds.
- [x] Add a zero-observed-leakage shot-planning helper that inverts the
  confidence bound for a target leakage and confidence level.
- [ ] Formalize the tolerant testing promise and prove completeness/soundness
  for fixed-witness verification.
- [x] Provide discovery followed by held-out fixed-witness verification, with
  separate results and shot accounting for both phases.

### P1: improve the algorithm rather than only brute force

- [x] Implement recursive conjugation-word sampling without full preceding-group access.
- [x] Document a finite-subgroup escape argument for the Clifford-target route,
  with exact finite-group regression tests and explicit work caps.
- [x] Derive the current group-operation and expanded-query recurrences, and
  distinguish these from runtime and finite-precision membership costs.
- [ ] Independently review the new probability argument and verify the paper's
  precise sufficient-criterion hypotheses against the implemented indexing.
- [ ] Replace exponentially expanded conjugation words by an efficient compact
  representation, or prove a short-word detection bound; group-operation
  complexity alone does not establish a polynomial-query tester.
- [ ] Supply a general exact or rigorously error-bounded membership backend;
  floating-point results do not inherit the exact-oracle theorem automatically.
- [ ] Replace noisy all-pairs output-MASA scoring with a support-aware search
  or optimization routine.
- [ ] Add a sparse Pauli-transfer representation and compare it with Bell
  sampling on simulated circuits.
- [ ] Add a stabilizer-tableau path for Clifford-heavy circuits.
- [x] Add source-circuit depth recording and checked-in one- and two-qubit
  identity workload baselines, with host runtime and memory kept separately.

### P2: broaden evidence and usability

- [ ] Add an analytically derived GSC-but-not-semi-Clifford fixture; one- and
  two-qubit positive families, a one-qubit negative, and wrong witnesses are
  already covered.
- [ ] Demonstrate transpilation and execution against an IBM-compatible noisy
  backend or noise model.
- [ ] Add an exact or symbolic backend for algebraic gate sets.
- [ ] Extend the binary implementation to prime-dimensional qudits.
- [ ] Add a command-line interface and machine-readable result export.

## Honest interpretation

The dense tester decides its bounded numerical search problem, not exact
symbolic GSC membership. The circuit tester measures Pauli coefficient mass
and can rigorously bound approximate leakage for a pre-specified witness, but
it does not yet constitute a complete property tester for exact GSC membership.
Negative finite-shot results are therefore evidence against the tested witness,
not a general proof that the unitary is not generalized semi-Clifford.
