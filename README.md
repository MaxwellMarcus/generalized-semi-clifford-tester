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

## Roadmap

1. Add an exact arithmetic backend for algebraic gate sets.
2. Replace input/output-pair enumeration with a support-derived candidate search.
3. Add stabilizer-tableau and sparse Pauli-transfer representations.
4. Validate against a larger corpus of known examples and counterexamples.
5. Benchmark scaling and document practical complexity limits.

See `docs/algorithm-design.md` for the acceptance criteria and proposed module
boundaries.

## License and citation

The code is available under the MIT License. Citation metadata is provided in
`CITATION.cff`.
