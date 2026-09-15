# Generalized Semi-Clifford Tester

[![CI](https://github.com/MaxwellMarcus/generalized-semi-clifford-tester/actions/workflows/ci.yml/badge.svg)](https://github.com/MaxwellMarcus/generalized-semi-clifford-tester/actions/workflows/ci.yml)
![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)

An early-stage research software project aimed at deciding whether a unitary has
generalized semi-Clifford structure and producing a checkable certificate when
possible.

> **Status:** foundational implementation, not yet a GSC decision procedure.
> The current release provides exact binary symplectic primitives, tests, and an
> implementation plan. It must not be cited as deciding GSC membership.

That status is deliberate: the mathematical criterion and certificate format
will be fixed before a public `is_gsc` API is exposed. This avoids turning a
numerical heuristic into an apparently authoritative yes/no answer.

## Implemented now

- Exact matrix arithmetic over \(\mathbb F_2\).
- Construction of the standard binary symplectic form.
- Validation of candidate symplectic matrices.
- Packaging, tests, continuous integration, and contribution guidance.

## Quick start

Python 3.10 or newer is required.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
python examples/check_symplectic.py
python -m pytest
```

## Roadmap

1. Freeze the mathematical membership criterion and supported input class.
2. Define machine-checkable positive and negative certificates.
3. Implement Pauli/symplectic normalization and invariant-subspace search.
4. Validate against analytically known examples and counterexamples.
5. Benchmark scaling and publish explicit complexity limits.

See `docs/algorithm-design.md` for the acceptance criteria and proposed module
boundaries.

## License and citation

The code is available under the MIT License. Citation metadata is provided in
`CITATION.cff`.
