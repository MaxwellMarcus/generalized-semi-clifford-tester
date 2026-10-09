"""Benchmark positive, negative, and capped exact qutrit searches."""

from __future__ import annotations

import json
import platform

from sympy import sqrt

from generalized_semi_clifford.exact_cyclotomic import (
    CyclotomicField,
    CyclotomicMatrix,
)
from generalized_semi_clifford.exact_odd_prime_benchmarking import (
    run_exact_prime_search_benchmark,
)

SUITE_SCHEMA = "generalized-semi-clifford/odd-prime-exact-search-benchmark-suite-v1"


def build_benchmark_suite() -> dict[str, object]:
    """Run the documented qutrit fixtures and return versioned records."""

    field = CyclotomicField(24)
    omega = field.root_power(8)
    scale = field.scalar(1 / sqrt(3))
    fourier = CyclotomicMatrix.from_rows(
        field,
        (
            (scale * omega ** (row * column) for column in range(3))
            for row in range(3)
        ),
    )
    phase = field.root_power(1)
    diagonal = CyclotomicMatrix.from_rows(
        field,
        ((1, 0, 0), (0, 1, 0), (0, 0, phase)),
    )
    non_gsc = diagonal @ fourier @ diagonal
    fixtures = (
        ("qutrit_fourier_positive", fourier, {}),
        ("qutrit_phase_fourier_phase_complete_negative", non_gsc, {}),
        (
            "qutrit_phase_fourier_phase_candidate_capped_unknown",
            non_gsc,
            {"max_candidate_pairs": 2},
        ),
    )
    records = [
        run_exact_prime_search_benchmark(name, unitary, 3, **limits).to_dict()
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
