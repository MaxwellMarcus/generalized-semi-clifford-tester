from __future__ import annotations

import json

import pytest

qiskit = pytest.importorskip("qiskit")

from generalized_semi_clifford import (  # noqa: E402
    GSC_BENCHMARK_SCHEMA,
    build_gsc_benchmark_record,
    run_gsc_sampling_benchmark,
    run_gsc_sampling_test,
)


def test_small_qubit_benchmark_json_schema_is_stable() -> None:
    from qiskit import QuantumCircuit

    result = run_gsc_sampling_test(QuantumCircuit(1), shots=16, seed=7)
    record = build_gsc_benchmark_record(
        result,
        classical_runtime_seconds=0.125,
        peak_memory_bytes=4096,
        seed=7,
    )

    payload = json.loads(record.to_json())
    assert payload == {
        "candidate_pairs_checked": 1,
        "circuits_executed": 3,
        "classical_runtime_seconds": 0.125,
        "has_candidate_witness": True,
        "lagrangian_count": 3,
        "num_qubits": 1,
        "peak_memory_bytes": 4096,
        "schema": GSC_BENCHMARK_SCHEMA,
        "search_mode": "exhaustive-lagrangian-discovery",
        "seed": 7,
        "shots_per_circuit": 16,
        "total_shots": 48,
    }


def test_benchmark_measures_runtime_and_peak_memory() -> None:
    from qiskit import QuantumCircuit

    benchmark = run_gsc_sampling_benchmark(QuantumCircuit(1), shots=8, seed=3)

    assert benchmark.result.has_candidate_witness
    assert benchmark.record.classical_runtime_seconds > 0
    assert benchmark.record.peak_memory_bytes > 0
    assert benchmark.record.total_shots == 24


@pytest.mark.parametrize(
    ("runtime", "memory", "message"),
    [(-0.1, 0, "classical_runtime_seconds"), (0.1, -1, "peak_memory_bytes")],
)
def test_benchmark_record_rejects_negative_measurements(
    runtime: float,
    memory: int,
    message: str,
) -> None:
    from qiskit import QuantumCircuit

    result = run_gsc_sampling_test(QuantumCircuit(1), shots=8, seed=1)
    with pytest.raises(ValueError, match=message):
        build_gsc_benchmark_record(
            result,
            classical_runtime_seconds=runtime,
            peak_memory_bytes=memory,
            seed=1,
        )
