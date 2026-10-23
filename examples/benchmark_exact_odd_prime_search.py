"""Benchmark positive, negative, and capped exact qutrit searches."""

from __future__ import annotations

import json
import platform

from sympy import sqrt

from generalized_semi_clifford.exact_cyclotomic import (
    CyclotomicField,
    CyclotomicMatrix,
    permutation_gate,
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
    two_qutrit_non_clifford = non_gsc.tensor(CyclotomicMatrix.identity(field, 3))
    two_qutrit_sum = permutation_gate(field, (0, 1, 2, 4, 5, 3, 8, 6, 7))
    fixtures = (
        ("qutrit_fourier_positive", fourier, {}),
        ("qutrit_phase_fourier_phase_complete_negative", non_gsc, {}),
        (
            "qutrit_phase_fourier_phase_candidate_capped_unknown",
            non_gsc,
            {"max_candidate_pairs": 2},
        ),
        (
            "two_qutrit_sum_coefficient_capped_unknown",
            two_qutrit_sum,
            {"max_coefficients": 80, "max_candidate_pairs": 1},
        ),
        (
            "two_qutrit_sum_verified_positive",
            two_qutrit_sum,
            {"max_coefficients": 324, "max_candidate_pairs": 1},
        ),
        (
            "two_qutrit_non_clifford_first_input_rejected_unknown",
            two_qutrit_non_clifford,
            {"max_coefficients": 162, "max_candidate_pairs": 40},
        ),
        (
            "two_qutrit_non_clifford_first_two_inputs_rejected_unknown",
            two_qutrit_non_clifford,
            {"max_coefficients": 324, "max_candidate_pairs": 80},
        ),
        (
            "two_qutrit_non_clifford_first_three_inputs_rejected_unknown",
            two_qutrit_non_clifford,
            {"max_coefficients": 486, "max_candidate_pairs": 120},
        ),
        (
            "two_qutrit_non_clifford_first_four_inputs_rejected_unknown",
            two_qutrit_non_clifford,
            {"max_coefficients": 648, "max_candidate_pairs": 160},
        ),
        (
            "two_qutrit_non_clifford_first_five_inputs_rejected_unknown",
            two_qutrit_non_clifford,
            {"max_coefficients": 810, "max_candidate_pairs": 200},
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
