# Generalized Semi-Clifford Tester

[![CI](https://github.com/MaxwellMarcus/generalized-semi-clifford-tester/actions/workflows/ci.yml/badge.svg)](https://github.com/MaxwellMarcus/generalized-semi-clifford-tester/actions/workflows/ci.yml)
![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)

An early-stage research software project for deciding whether a unitary has
generalized semi-Clifford structure and producing a checkable Pauli-MASA
certificate when possible.

> **Status:** a definition-first exhaustive tester is implemented for small qubit
> systems. It uses complex floating-point arithmetic and defaults to `UNKNOWN`
> above three qubits. It must not be presented as an exact symbolic proof.

The tester checks the defining condition directly: it searches for binary
Lagrangians (L,S) such that (U\mathcal A_LU^\dagger=\mathcal A_S). Every
result records its tolerance and arithmetic, and positive results contain a
witness that can be independently recomputed.

## Implemented now

- Exact matrix arithmetic over \(\mathbb F_2\).
- Construction of the standard binary symplectic form.
- Validation of candidate symplectic matrices.
- Canonical RREF enumeration of the 3, 15, 135, and 2,295 Lagrangians for one
  through four qubits; the dense numerical tester retains its three-qubit
  default limit.
- Dense Pauli matrices in a documented binary \((x\mid z)\) convention.
- A bounded naïve GSC tester with `GSC`, `NOT_GSC`, and `UNKNOWN` outcomes.
- Independently verifiable numerical witnesses.
- Replayable conjugation-word sampling without enumerating earlier groups,
  an exact-oracle randomized Clifford-containment algorithm, and a bounded
  numerical GSC/SC-distance frontend (see below).
- Analytically justified one- and two-qubit fixture families, including
  explicit valid and deliberately invalid Pauli-MASA witnesses.
- An experimental Qiskit semi-Clifford tester based on Pauli-conjugation
  circuits and Bell-basis sampling.
- Qiskit GSC discovery and \(n\)-circuit candidate-witness verification based
  on the full observed Pauli support of each conjugate.
- Packaging, tests, continuous integration, and contribution guidance.

## Quick start

Python 3.10 or newer is required.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
python examples/check_symplectic.py
python examples/naive_tester.py
python -m pytest
```

Install the optional Qiskit implementation with:

```bash
python -m pip install -e ".[qiskit]"
python examples/qiskit_semi_clifford.py
python examples/qiskit_gsc.py
```

## Python API

```python
import numpy as np
from generalized_semi_clifford import check_gsc_naive, verify_gsc_witness

t_gate = np.diag([1, np.exp(1j * np.pi / 4)])
result = check_gsc_naive(t_gate)

print(result.status, result.best_leakage)
assert result.witness is not None
assert verify_gsc_witness(t_gate, result.witness)
```

The reusable fixture catalog keeps its mathematical claim separate from the
floating-point check. It includes arbitrary one-qubit phase gates and
two-qubit controlled-phase gates with fixed-Z witnesses, plus a one-qubit
rotation whose Bloch action has three nonzero Pauli coefficients in every
column and is therefore not GSC:

```python
from generalized_semi_clifford import analytic_gsc_fixtures

for fixture in analytic_gsc_fixtures():
    result = check_gsc_naive(fixture.unitary)
    assert result.status is fixture.expected_status
```

These fixtures support regression testing; their concise algebraic
justifications do not turn the package's floating-point result into an exact
symbolic proof.

The default three-qubit limit is a complexity guard, not a mathematical
restriction on the definition. See `docs/algorithm-design.md` for the exact
criterion and scaling discussion.

## Experimental Qiskit semi-Clifford tester

```python
from qiskit import QuantumCircuit

from generalized_semi_clifford import (
    run_semi_clifford_sampling_test,
    run_semi_clifford_witness_test,
)

unitary = QuantumCircuit(2)
unitary.h(0)
unitary.cx(0, 1)
unitary.t(1)

result = run_semi_clifford_sampling_test(unitary, shots=1024, seed=7)
print(result.has_candidate_witness)
print(result.witness)

# If theory supplies a candidate input Lagrangian, verify it with only n circuits.
candidate = ((1, 0, 0, 0), (0, 1, 0, 0))  # X on qubits 0 and 1
verification = run_semi_clifford_witness_test(unitary, candidate, shots=1024)
```

For every nonidentity Pauli \(P\), this prototype prepares the Choi state of
\(UPU^\dagger\) and measures it in the Pauli Bell basis. It then searches the
high-probability Pauli conjugates for \(n\) independent commuting input/output
pairs. Such a basis is a candidate Lagrangian witness for semi-Cliffordness.

This is a finite-shot experiment, not an exact proof or a completed property-
testing theorem. Exhaustive discovery uses \(4^n-1\) circuits, while checking a
proposed Lagrangian basis uses only \(n\). Sampler jobs can be batched for
backend limits. Negative results from a nonexhaustive input set are explicitly
marked as incomplete.

## Experimental Qiskit GSC tester

The GSC tester reuses the same circuits but retains every Bell outcome rather
than only the dominant Pauli. Given input and output Lagrangian bases, it checks
whether the empirical Pauli support of each conjugated input generator stays
inside the output Lagrangian:

```python
from generalized_semi_clifford import run_gsc_witness_test, zero_event_shots_required

shots = zero_event_shots_required(
    0.01,
    confidence_level=0.95,
    simultaneous_tests=len(input_basis),
)

result = run_gsc_witness_test(
    unitary,
    input_basis,
    output_basis,
    shots=shots,
    leakage_threshold=0.01,
)

print(result.has_candidate_witness)
print(result.maximum_leakage_upper_bound)
print(result.has_confidence_certified_witness)
```

Candidate verification uses \(n\) circuits. Small-qubit discovery enumerates
all input/output Lagrangians, but reuses each Pauli-conjugation experiment and
only constructs circuits for Paulis appearing in the canonical input bases: 3,
13, 47, and 165 circuits for one through four qubits, respectively. Exact-
support cases construct the forced output Lagrangian directly; noisy cases fall
back to exhaustive output-Lagrangian scoring.

Fixed-witness verification also reports an exact one-sided binomial upper
bound on maximum leakage, Bonferroni-corrected across the `n` generators. This
confidence statement is valid when the witness was chosen independently of
the verification samples. After exploratory discovery, verification must use
fresh samples; reusing the discovery counts would introduce selection bias.
The convenience workflow enforces this separation and reports each phase's
total shot cost:

```python
from generalized_semi_clifford import run_gsc_discovery_then_verification

result = run_gsc_discovery_then_verification(
    unitary,
    discovery_shots=256,
    verification_shots=2048,
    leakage_threshold=0.01,
    confidence_level=0.95,
)
print(result.discovery_shot_cost, result.verification_shot_cost)
print(result.has_confidence_certified_witness)
```

The discovery phase may select a candidate but never receives a confidence
claim. Only the separate verification run can set
`has_confidence_certified_witness`; a failed verification is evidence against
that candidate rather than a general non-GSC conclusion.

For measured small-qubit scaling records, `run_gsc_sampling_benchmark` wraps a
discovery run and emits a versioned JSON record containing the qubit and
Lagrangian counts, circuits, shots, candidate pairs, exact untranspiled source-
circuit depth and operation counts, elapsed host time, and Python-managed peak
memory. Exact workload fields and machine-dependent measurements are separate
objects in the version-two schema. The checked-in one- and two-qubit identity
workloads in `benchmarks/identity-small-qubits.json` are regression baselines;
they intentionally omit runtime and memory. The fixed one-qubit workload in
`examples/benchmark_gsc.py` can be rerun with `python examples/benchmark_gsc.py`.
Runtime includes sampling and post-processing, and peak memory uses
`tracemalloc`; neither field is a hardware-independent complexity claim. The
source depth is pre-transpilation circuit structure, not a hardware depth.

## Conjugation sampling and SC-distance testing

The new experimental frontend samples either `Gamma_2(U)` against `C_(k-2)`
or `Gamma_(k-2)(U)` against `C_2`, using replayable words rather than full
preceding groups:

```python
from generalized_semi_clifford import check_gsc_sandwich

result = check_gsc_sandwich(t_gate, hierarchy_level=3, epsilon=0.01, seed=7)
print(result.status, result.gsc_witness)
```

This three-qubit-limited numerical frontend returns a direct GSC witness,
a distance bound from semi-Clifford gates, or `INCONCLUSIVE`. It distinguishes
distance from **all** semi-Clifford gates from distance only within `C_k`.
Passing short sampled words never supplies its positive GSC certificate.

Separately, `test_clifford_conjugation_containment` implements a recursive
finite-subgroup escape algorithm with a false-acceptance bound under **exact
group operations and an exact Clifford membership oracle**. It needs no
enumeration of any previous conjugation group. Its polynomial group-operation
count at fixed depth does not yet imply polynomial physical queries to `U`:
shared words can expand exponentially, and exact arithmetic costs remain an
additional issue. Numerical backends do not receive an unconditional
probability guarantee.

See [the algorithm and proof notes](docs/conjugation-sampling.md) for the
precise claims, explicit complexity recurrences, comparison with the
property-testing literature, remaining research gaps, and API. At fixed
error probability, the current recursive construction uses `O(n^4)`,
`O(n^8)`, and `O(n^13)` abstract group operations for hierarchy levels
4, 5, and 6, respectively; literal expanded seed-call bounds are exponential.
These are not polynomial total-runtime or physical-query guarantees. Run
`python examples/conjugation_sampling.py` for a small numerical demonstration.

## Roadmap

1. Add an exact arithmetic backend for algebraic gate sets.
2. Replace the noisy-data exhaustive output-Lagrangian fallback with a more
   scalable support-aware optimization.
3. Prove completeness/soundness bounds for the circuit-based semi-Clifford
   tester and replace exhaustive Pauli enumeration with sampling where possible.
4. Add stabilizer-tableau and sparse Pauli-transfer representations.
5. Validate against a larger corpus of known examples and counterexamples.
6. Benchmark scaling and document practical complexity limits.

See `docs/algorithm-design.md` for the acceptance criteria and proposed module
boundaries, and `docs/status-and-roadmap.md` for the current audit and
prioritized work list.

## License and citation

The code is available under the MIT License. Citation metadata is provided in
`CITATION.cff`.
