"""Benchmark representative bounded exact searches through three qubits."""

from __future__ import annotations

import json
import platform

from generalized_semi_clifford.exact_benchmarking import run_exact_search_benchmark
from generalized_semi_clifford.exact_cyclotomic import (
    CyclotomicField,
    CyclotomicMatrix,
    cnot,
    hadamard,
    t_gate,
)

SUITE_SCHEMA = "generalized-semi-clifford/exact-search-benchmark-suite-v1"


def build_benchmark_suite() -> dict[str, object]:
    """Run the documented fixtures and return their versioned records."""

    field = CyclotomicField(8)
    identity = CyclotomicMatrix.identity(field, 2)
    fixtures = (
        ("one_qubit_h_positive", hadamard(field), {}),
        (
            "one_qubit_t_h_t_complete_negative",
            t_gate(field) @ hadamard(field) @ t_gate(field),
            {},
        ),
        ("two_qubit_cnot_positive", cnot(field), {}),
        ("three_qubit_cnot_identity_positive", cnot(field).tensor(identity), {}),
        (
            "three_qubit_h_identity_identity_capped_unknown",
            hadamard(field).tensor(identity).tensor(identity),
            {"max_candidate_pairs": 1},
        ),
    )
    records = [
        run_exact_search_benchmark(name, unitary, **limits).to_dict()
        for name, unitary, limits in fixtures
    ]
    return {
        "environment": {
            "implementation": platform.python_implementation(),
            "platform": platform.platform(),
            "python_version": platform.python_version(),
        },
        "records": records,
        "schema": SUITE_SCHEMA,
    }


def main() -> None:
    print(json.dumps(build_benchmark_suite(), allow_nan=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
