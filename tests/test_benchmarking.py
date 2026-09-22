from __future__ import annotations

import json
from pathlib import Path

import pytest

qiskit = pytest.importorskip("qiskit")

from generalized_semi_clifford import (  # noqa: E402
    GSC_BENCHMARK_SCHEMA,
    build_gsc_benchmark_record,
    build_semi_clifford_test_circuits,
    run_gsc_sampling_benchmark,
    run_gsc_sampling_test,
)


def test_small_qubit_benchmark_json_schema_is_stable() -> None:
    from qiskit import QuantumCircuit

    result = run_gsc_sampling_test(QuantumCircuit(1), shots=16, seed=7)
    source_circuits = tuple(
        circuit
        for _, circuit in build_semi_clifford_test_circuits(QuantumCircuit(1))
    )
    record = build_gsc_benchmark_record(
        result,
        source_circuits=source_circuits,
        classical_runtime_seconds=0.125,
        peak_memory_bytes=4096,
        seed=7,
    )

    payload = json.loads(record.to_json())
    assert payload == {
        "has_candidate_witness": True,
        "measurements": {
            "classical_runtime_seconds": 0.125,
            "peak_memory_bytes": 4096,
        },
        "schema": GSC_BENCHMARK_SCHEMA,
        "search_mode": "support-aware-lagrangian-discovery",
        "seed": 7,
        "workload": {
            "candidate_pairs_checked": 1,
            "circuits_executed": 3,
            "lagrangian_count": 3,
            "num_qubits": 1,
            "shots_per_circuit": 16,
            "source_circuit_depth_histogram": {"8": 3},
            "source_circuit_depth_max": 8,
            "source_circuit_depth_total": 24,
            "source_operation_counts": {
                "cx": 6,
                "h": 6,
                "measure": 6,
                "queried_unitary": 3,
                "queried_unitary_dagger": 3,
                "x": 1,
                "y": 1,
                "z": 1,
            },
            "total_shots": 48,
        },
    }


def test_benchmark_measures_runtime_and_peak_memory() -> None:
    from qiskit import QuantumCircuit

    benchmark = run_gsc_sampling_benchmark(QuantumCircuit(1), shots=8, seed=3)

    assert benchmark.result.has_candidate_witness
    assert benchmark.record.measurements.classical_runtime_seconds > 0
    assert benchmark.record.measurements.peak_memory_bytes > 0
    assert benchmark.record.workload.total_shots == 24
    assert benchmark.record.workload.source_circuit_depth_max == 8


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
            source_circuits=(
                circuit
                for _, circuit in build_semi_clifford_test_circuits(QuantumCircuit(1))
            ),
            classical_runtime_seconds=runtime,
            peak_memory_bytes=memory,
            seed=1,
        )


def test_benchmark_record_requires_one_source_circuit_per_observation() -> None:
    from qiskit import QuantumCircuit

    result = run_gsc_sampling_test(QuantumCircuit(1), shots=8, seed=1)
    with pytest.raises(ValueError, match="source_circuits"):
        build_gsc_benchmark_record(
            result,
            source_circuits=(),
            classical_runtime_seconds=0.1,
            peak_memory_bytes=1,
            seed=1,
        )


def test_checked_in_identity_workload_baselines_are_reproducible() -> None:
    from qiskit import QuantumCircuit

    baseline_path = Path(__file__).parents[1] / "benchmarks" / "identity-small-qubits.json"
    payload = json.loads(baseline_path.read_text())

    assert payload["schema"] == "generalized-semi-clifford/gsc-workload-baselines-v1"
    assert payload["unitary"] == "identity"
    reproduced = []
    for expected in payload["baselines"]:
        circuit = QuantumCircuit(
            expected["num_qubits"],
            name=f"identity_{expected['num_qubits']}q",
        )
        benchmark = run_gsc_sampling_benchmark(
            circuit,
            shots=expected["shots_per_circuit"],
            seed=payload["seed"],
        )
        reproduced.append(benchmark.record.workload.to_dict())

    assert reproduced == payload["baselines"]
