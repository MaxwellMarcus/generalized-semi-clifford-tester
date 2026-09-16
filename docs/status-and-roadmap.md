# Status and roadmap

Last audited: 2026-09-16.

## What works now

| Layer | Capability | Practical boundary |
| --- | --- | --- |
| Binary algebra | Exact symplectic arithmetic and canonical Lagrangian enumeration | 1--4 qubits enumerated: 3, 15, 135, 2,295 Lagrangians |
| Dense GSC tester | Exhaustive input/output Pauli-MASA search with independently checkable witnesses | Defaults to `UNKNOWN` above 3 qubits; complex128 with a stated tolerance |
| Semi-Clifford circuits | Bell-sampling conjugation circuits, exhaustive discovery, and `n`-circuit candidate verification | Finite-shot candidate evidence; exhaustive discovery uses `4**n - 1` circuits |
| GSC circuits | Full sampled Pauli-support decoding, exhaustive discovery, and `n`-circuit fixed-witness verification | Discovery through 4 qubits; support-aware fast path with an exhaustive noisy fallback |
| Statistics | Exact one-sided binomial leakage bounds with familywise confidence for a fixed GSC witness | Valid for a witness fixed independently of the verification samples |

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

- [ ] Replace noisy all-pairs output-MASA scoring with a support-aware search
  or optimization routine.
- [ ] Add a sparse Pauli-transfer representation and compare it with Bell
  sampling on simulated circuits.
- [ ] Add a stabilizer-tableau path for Clifford-heavy circuits.
- [ ] Separate and benchmark circuit count, circuit depth, shots, classical
  runtime, and memory.

### P2: broaden evidence and usability

- [ ] Add analytically derived GSC, semi-Clifford, non-semi-Clifford, and
  non-GSC fixture families rather than isolated examples.
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
