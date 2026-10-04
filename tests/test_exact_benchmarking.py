from __future__ import annotations

import json
from pathlib import Path

import pytest

from generalized_semi_clifford.exact_benchmarking import (
    EXACT_SEARCH_BENCHMARK_SCHEMA,
    run_exact_search_benchmark,
)
from generalized_semi_clifford.exact_cyclotomic import (
    CyclotomicField,
    ExactSearchStatus,
    hadamard,
)


def test_exact_benchmark_separates_workload_from_host_measurements() -> None:
    record = run_exact_search_benchmark(
        "one_qubit_h_positive",
        hadamard(CyclotomicField(8)),
    )

    payload = json.loads(record.to_json())
    assert payload["schema"] == EXACT_SEARCH_BENCHMARK_SCHEMA
    assert payload["outcome"] == {
        "search_complete": False,
        "status": "gsc",
        "stop_reason": "witness_found",
        "witness_verified": True,
    }
    assert payload["exact_workload"]["candidate_pairs_checked"] == 3
    assert payload["exact_workload"]["coefficients_computed"] == 8
    assert payload["measurements"]["host_runtime_seconds"] >= 0
    assert payload["measurements"]["python_peak_memory_bytes"] > 0
    assert record.result.status is ExactSearchStatus.GSC


def test_exact_benchmark_rejects_ambiguous_measurement_sessions() -> None:
    import tracemalloc

    tracemalloc.start()
    try:
        with pytest.raises(RuntimeError, match="tracemalloc"):
            run_exact_search_benchmark("h", hadamard(CyclotomicField(8)))
    finally:
        tracemalloc.stop()


def test_checked_in_suite_covers_sizes_and_outcome_boundaries() -> None:
    path = Path(__file__).parents[1] / "benchmarks" / "exact-search-small-qubits.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    records = payload["records"]

    assert payload["schema"] == "generalized-semi-clifford/exact-search-benchmark-suite-v1"
    assert {record["exact_workload"]["num_qubits"] for record in records} == {1, 2, 3}
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
