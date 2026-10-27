from __future__ import annotations

import json
from pathlib import Path

import pytest

sympy = pytest.importorskip("sympy")

from generalized_semi_clifford.exact_cyclotomic import (  # noqa: E402
    CyclotomicField,
    CyclotomicMatrix,
)
from generalized_semi_clifford.exact_odd_prime import ExactPrimeSearchStatus  # noqa: E402
from generalized_semi_clifford.exact_odd_prime_benchmarking import (  # noqa: E402
    EXACT_PRIME_SEARCH_BENCHMARK_SCHEMA,
    run_exact_prime_search_benchmark,
)


def _qutrit_fourier(field: CyclotomicField) -> CyclotomicMatrix:
    omega = field.root_power(field.order // 3)
    scale = field.scalar(1 / sympy.sqrt(3))
    return CyclotomicMatrix.from_rows(
        field,
        (
            (scale * omega ** (row * column) for column in range(3))
            for row in range(3)
        ),
    )


def test_prime_benchmark_separates_exact_work_from_host_measurements() -> None:
    record = run_exact_prime_search_benchmark(
        "qutrit_fourier_positive",
        _qutrit_fourier(CyclotomicField(24)),
        3,
    )

    payload = json.loads(record.to_json())
    assert payload["schema"] == EXACT_PRIME_SEARCH_BENCHMARK_SCHEMA
    assert payload["outcome"] == {
        "search_complete": False,
        "status": "gsc",
        "stop_reason": "witness_found",
        "witness_verified": True,
    }
    assert payload["exact_workload"]["prime"] == 3
    assert payload["exact_workload"]["num_qudits"] == 1
    assert payload["exact_workload"]["candidate_pairs_checked"] > 0
    assert payload["exact_workload"]["coefficients_computed"] >= 18
    assert payload["measurements"]["host_runtime_seconds"] >= 0
    assert payload["measurements"]["python_peak_memory_bytes"] > 0
    assert record.result.status is ExactPrimeSearchStatus.GSC


def test_prime_benchmark_rejects_ambiguous_measurement_sessions() -> None:
    import tracemalloc

    tracemalloc.start()
    try:
        with pytest.raises(RuntimeError, match="tracemalloc"):
            run_exact_prime_search_benchmark(
                "identity",
                CyclotomicMatrix.identity(CyclotomicField(24), 3),
                3,
            )
    finally:
        tracemalloc.stop()


def test_checked_in_qutrit_suite_covers_outcome_boundaries() -> None:
    path = Path(__file__).parents[1] / "benchmarks" / "exact-odd-prime-search-qutrit.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    records = payload["records"]

    assert payload["schema"] == (
        "generalized-semi-clifford/odd-prime-exact-search-benchmark-suite-v1"
    )
    assert {record["exact_workload"]["prime"] for record in records} == {3}
    assert {record["exact_workload"]["num_qudits"] for record in records} == {1, 2}
    assert {record["outcome"]["status"] for record in records} == {
        "gsc",
        "not_gsc",
        "unknown",
    }
    negative, = (record for record in records if record["outcome"]["status"] == "not_gsc")
    candidate_capped, = (
        record
        for record in records
        if record["fixture"] == "qutrit_phase_fourier_phase_candidate_capped_unknown"
    )
    coefficient_capped, = (
        record
        for record in records
        if record["fixture"] == "two_qutrit_sum_coefficient_capped_unknown"
    )
    two_qutrit_positive, = (
        record
        for record in records
        if record["fixture"] == "two_qutrit_sum_verified_positive"
    )
    two_qutrit_rejection_prefix, = (
        record
        for record in records
        if record["fixture"]
        == "two_qutrit_non_clifford_first_input_rejected_unknown"
    )
    two_qutrit_second_rejection_prefix, = (
        record
        for record in records
        if record["fixture"]
        == "two_qutrit_non_clifford_first_two_inputs_rejected_unknown"
    )
    two_qutrit_third_rejection_prefix, = (
        record
        for record in records
        if record["fixture"]
        == "two_qutrit_non_clifford_first_three_inputs_rejected_unknown"
    )
    two_qutrit_fourth_rejection_prefix, = (
        record
        for record in records
        if record["fixture"]
        == "two_qutrit_non_clifford_first_four_inputs_rejected_unknown"
    )
    two_qutrit_fifth_rejection_prefix, = (
        record
        for record in records
        if record["fixture"]
        == "two_qutrit_non_clifford_first_five_inputs_rejected_unknown"
    )
    two_qutrit_sixth_rejection_prefix, = (
        record
        for record in records
        if record["fixture"]
        == "two_qutrit_non_clifford_first_six_inputs_rejected_unknown"
    )
    two_qutrit_seventh_rejection_prefix, = (
        record
        for record in records
        if record["fixture"]
        == "two_qutrit_non_clifford_first_seven_inputs_rejected_unknown"
    )
    two_qutrit_eighth_rejection_prefix, = (
        record
        for record in records
        if record["fixture"]
        == "two_qutrit_non_clifford_first_eight_inputs_rejected_unknown"
    )
    two_qutrit_ninth_rejection_prefix, = (
        record
        for record in records
        if record["fixture"]
        == "two_qutrit_non_clifford_first_nine_inputs_rejected_unknown"
    )
    assert negative["outcome"]["search_complete"] is True
    assert negative["outcome"]["stop_reason"] == "exhausted"
    assert candidate_capped["outcome"] == {
        "search_complete": False,
        "status": "unknown",
        "stop_reason": "candidate_pair_cap",
        "witness_verified": None,
    }
    assert coefficient_capped["outcome"] == {
        "search_complete": False,
        "status": "unknown",
        "stop_reason": "coefficient_cap",
        "witness_verified": None,
    }
    assert coefficient_capped["exact_workload"] == {
        "arithmetic": "cyclotomic_exact_odd_prime",
        "candidate_pairs_checked": 0,
        "coefficients_computed": 0,
        "field_order": 24,
        "input_lagrangians_checked": 0,
        "lagrangians_total": 40,
        "limits": {
            "max_candidate_pairs": 1,
            "max_coefficients": 80,
            "max_qudits": 2,
        },
        "num_qudits": 2,
        "prime": 3,
        "witness_verification_coefficients": 0,
    }
    assert two_qutrit_positive["outcome"] == {
        "search_complete": False,
        "status": "gsc",
        "stop_reason": "witness_found",
        "witness_verified": True,
    }
    assert two_qutrit_positive["exact_workload"] == {
        "arithmetic": "cyclotomic_exact_odd_prime",
        "candidate_pairs_checked": 1,
        "coefficients_computed": 324,
        "field_order": 24,
        "input_lagrangians_checked": 1,
        "lagrangians_total": 40,
        "limits": {
            "max_candidate_pairs": 1,
            "max_coefficients": 324,
            "max_qudits": 2,
        },
        "num_qudits": 2,
        "prime": 3,
        "witness_verification_coefficients": 162,
    }
    assert two_qutrit_rejection_prefix["outcome"] == {
        "search_complete": False,
        "status": "unknown",
        "stop_reason": "coefficient_cap",
        "witness_verified": None,
    }
    assert two_qutrit_rejection_prefix["exact_workload"] == {
        "arithmetic": "cyclotomic_exact_odd_prime",
        "candidate_pairs_checked": 40,
        "coefficients_computed": 162,
        "field_order": 24,
        "input_lagrangians_checked": 1,
        "lagrangians_total": 40,
        "limits": {
            "max_candidate_pairs": 40,
            "max_coefficients": 162,
            "max_qudits": 2,
        },
        "num_qudits": 2,
        "prime": 3,
        "witness_verification_coefficients": 0,
    }
    assert two_qutrit_second_rejection_prefix["outcome"] == {
        "search_complete": False,
        "status": "unknown",
        "stop_reason": "coefficient_cap",
        "witness_verified": None,
    }
    assert two_qutrit_second_rejection_prefix["exact_workload"] == {
        "arithmetic": "cyclotomic_exact_odd_prime",
        "candidate_pairs_checked": 80,
        "coefficients_computed": 324,
        "field_order": 24,
        "input_lagrangians_checked": 2,
        "lagrangians_total": 40,
        "limits": {
            "max_candidate_pairs": 80,
            "max_coefficients": 324,
            "max_qudits": 2,
        },
        "num_qudits": 2,
        "prime": 3,
        "witness_verification_coefficients": 0,
    }
    assert two_qutrit_third_rejection_prefix["outcome"] == {
        "search_complete": False,
        "status": "unknown",
        "stop_reason": "coefficient_cap",
        "witness_verified": None,
    }
    assert two_qutrit_third_rejection_prefix["exact_workload"] == {
        "arithmetic": "cyclotomic_exact_odd_prime",
        "candidate_pairs_checked": 120,
        "coefficients_computed": 486,
        "field_order": 24,
        "input_lagrangians_checked": 3,
        "lagrangians_total": 40,
        "limits": {
            "max_candidate_pairs": 120,
            "max_coefficients": 486,
            "max_qudits": 2,
        },
        "num_qudits": 2,
        "prime": 3,
        "witness_verification_coefficients": 0,
    }
    assert two_qutrit_fourth_rejection_prefix["outcome"] == {
        "search_complete": False,
        "status": "unknown",
        "stop_reason": "coefficient_cap",
        "witness_verified": None,
    }
    assert two_qutrit_fourth_rejection_prefix["exact_workload"] == {
        "arithmetic": "cyclotomic_exact_odd_prime",
        "candidate_pairs_checked": 160,
        "coefficients_computed": 648,
        "field_order": 24,
        "input_lagrangians_checked": 4,
        "lagrangians_total": 40,
        "limits": {
            "max_candidate_pairs": 160,
            "max_coefficients": 648,
            "max_qudits": 2,
        },
        "num_qudits": 2,
        "prime": 3,
        "witness_verification_coefficients": 0,
    }
    assert two_qutrit_fifth_rejection_prefix["outcome"] == {
        "search_complete": False,
        "status": "unknown",
        "stop_reason": "coefficient_cap",
        "witness_verified": None,
    }
    assert two_qutrit_fifth_rejection_prefix["exact_workload"] == {
        "arithmetic": "cyclotomic_exact_odd_prime",
        "candidate_pairs_checked": 200,
        "coefficients_computed": 810,
        "field_order": 24,
        "input_lagrangians_checked": 5,
        "lagrangians_total": 40,
        "limits": {
            "max_candidate_pairs": 200,
            "max_coefficients": 810,
            "max_qudits": 2,
        },
        "num_qudits": 2,
        "prime": 3,
        "witness_verification_coefficients": 0,
    }
    assert two_qutrit_sixth_rejection_prefix["outcome"] == {
        "search_complete": False,
        "status": "unknown",
        "stop_reason": "coefficient_cap",
        "witness_verified": None,
    }
    assert two_qutrit_sixth_rejection_prefix["exact_workload"] == {
        "arithmetic": "cyclotomic_exact_odd_prime",
        "candidate_pairs_checked": 240,
        "coefficients_computed": 972,
        "field_order": 24,
        "input_lagrangians_checked": 6,
        "lagrangians_total": 40,
        "limits": {
            "max_candidate_pairs": 240,
            "max_coefficients": 972,
            "max_qudits": 2,
        },
        "num_qudits": 2,
        "prime": 3,
        "witness_verification_coefficients": 0,
    }
    assert two_qutrit_seventh_rejection_prefix["outcome"] == {
        "search_complete": False,
        "status": "unknown",
        "stop_reason": "coefficient_cap",
        "witness_verified": None,
    }
    assert two_qutrit_seventh_rejection_prefix["exact_workload"] == {
        "arithmetic": "cyclotomic_exact_odd_prime",
        "candidate_pairs_checked": 280,
        "coefficients_computed": 1134,
        "field_order": 24,
        "input_lagrangians_checked": 7,
        "lagrangians_total": 40,
        "limits": {
            "max_candidate_pairs": 280,
            "max_coefficients": 1134,
            "max_qudits": 2,
        },
        "num_qudits": 2,
        "prime": 3,
        "witness_verification_coefficients": 0,
    }
    assert two_qutrit_eighth_rejection_prefix["outcome"] == {
        "search_complete": False,
        "status": "unknown",
        "stop_reason": "coefficient_cap",
        "witness_verified": None,
    }
    assert two_qutrit_eighth_rejection_prefix["exact_workload"] == {
        "arithmetic": "cyclotomic_exact_odd_prime",
        "candidate_pairs_checked": 320,
        "coefficients_computed": 1296,
        "field_order": 24,
        "input_lagrangians_checked": 8,
        "lagrangians_total": 40,
        "limits": {
            "max_candidate_pairs": 320,
            "max_coefficients": 1296,
            "max_qudits": 2,
        },
        "num_qudits": 2,
        "prime": 3,
        "witness_verification_coefficients": 0,
    }
    assert two_qutrit_ninth_rejection_prefix["outcome"] == {
        "search_complete": False,
        "status": "unknown",
        "stop_reason": "coefficient_cap",
        "witness_verified": None,
    }
    assert two_qutrit_ninth_rejection_prefix["exact_workload"] == {
        "arithmetic": "cyclotomic_exact_odd_prime",
        "candidate_pairs_checked": 360,
        "coefficients_computed": 1458,
        "field_order": 24,
        "input_lagrangians_checked": 9,
        "lagrangians_total": 40,
        "limits": {
            "max_candidate_pairs": 360,
            "max_coefficients": 1458,
            "max_qudits": 2,
        },
        "num_qudits": 2,
        "prime": 3,
        "witness_verification_coefficients": 0,
    }
    assert all(record["measurements"]["host_runtime_seconds"] > 0 for record in records)
    assert all(record["measurements"]["python_peak_memory_bytes"] > 0 for record in records)
