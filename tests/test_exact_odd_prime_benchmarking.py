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
    assert {record["exact_workload"]["num_qudits"] for record in records} == {1}
    assert {record["outcome"]["status"] for record in records} == {
        "gsc",
        "not_gsc",
        "unknown",
    }
    negative, = (record for record in records if record["outcome"]["status"] == "not_gsc")
    unknown, = (record for record in records if record["outcome"]["status"] == "unknown")
    assert negative["outcome"]["search_complete"] is True
    assert negative["outcome"]["stop_reason"] == "exhausted"
    assert unknown["outcome"]["search_complete"] is False
    assert unknown["outcome"]["stop_reason"] == "candidate_pair_cap"
    assert all(record["measurements"]["host_runtime_seconds"] > 0 for record in records)
    assert all(record["measurements"]["python_peak_memory_bytes"] > 0 for record in records)
