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
- Exhaustive enumeration of the 3, 15, and 135 Lagrangians for one, two, and
  three qubits.
- Dense Pauli matrices in a documented binary \((x\mid z)\) convention.
- A bounded naïve GSC tester with `GSC`, `NOT_GSC`, and `UNKNOWN` outcomes.
- Independently verifiable numerical witnesses.
- An experimental Qiskit semi-Clifford tester based on Pauli-conjugation
  circuits and Bell-basis sampling.
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

## Roadmap

1. Add an exact arithmetic backend for algebraic gate sets.
2. Replace input/output-pair enumeration with a support-derived candidate search.
3. Prove completeness/soundness bounds for the circuit-based semi-Clifford
   tester and replace exhaustive Pauli enumeration with sampling where possible.
4. Add stabilizer-tableau and sparse Pauli-transfer representations.
5. Validate against a larger corpus of known examples and counterexamples.
6. Benchmark scaling and document practical complexity limits.

See `docs/algorithm-design.md` for the acceptance criteria and proposed module
boundaries.

## License and citation

The code is available under the MIT License. Citation metadata is provided in
`CITATION.cff`.
