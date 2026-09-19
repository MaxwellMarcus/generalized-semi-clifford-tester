"""Reproducible records for small-qubit GSC sampling benchmarks."""

from __future__ import annotations

import json
import math
import time
import tracemalloc
from collections import Counter
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from typing import Any

from .gsc_qiskit import GSCSamplingResult, run_gsc_sampling_test
from .lagrangian import enumerate_lagrangians, lagrangian_count
from .semi_clifford_qiskit import build_semi_clifford_test_circuits

GSC_BENCHMARK_SCHEMA = "generalized-semi-clifford/gsc-sampling-benchmark-v2"


@dataclass(frozen=True)
class GSCBenchmarkWorkload:
    """Exact workload dimensions and source-circuit operation counts.

    Circuit metrics describe the untranspiled Qiskit circuits submitted by the
    benchmark. They are exact counts for those source circuits, not hardware
    depths or estimates of backend cost.
    """

    num_qubits: int
    lagrangian_count: int
    circuits_executed: int
    shots_per_circuit: int
    total_shots: int
    candidate_pairs_checked: int
    source_circuit_depth_max: int
    source_circuit_depth_total: int
    source_circuit_depth_histogram: tuple[tuple[int, int], ...]
    source_operation_counts: tuple[tuple[str, int], ...]

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-compatible exact-workload payload."""

        payload = asdict(self)
        payload["source_circuit_depth_histogram"] = {
            str(depth): count for depth, count in self.source_circuit_depth_histogram
        }
        payload["source_operation_counts"] = dict(self.source_operation_counts)
        return payload


@dataclass(frozen=True)
class GSCBenchmarkMeasurements:
    """Machine-dependent measurements for one benchmark run."""

    classical_runtime_seconds: float
    peak_memory_bytes: int


@dataclass(frozen=True)
class GSCBenchmarkRecord:
    """Machine-readable resource record for one GSC discovery run.

    Exact source-circuit counts are kept in ``workload``. Runtime is elapsed
    host wall time for sampling plus classical post-processing. Peak memory is
    measured with :mod:`tracemalloc`, so it covers Python-managed allocations
    rather than total process memory. These host measurements live separately
    in ``measurements`` and must not be treated as exact complexity claims.
    """

    schema: str
    workload: GSCBenchmarkWorkload
    measurements: GSCBenchmarkMeasurements
    seed: int | None
    search_mode: str
    has_candidate_witness: bool

    def to_dict(self) -> dict[str, object]:
        """Return the versioned JSON-compatible record."""

        return {
            "has_candidate_witness": self.has_candidate_witness,
            "measurements": asdict(self.measurements),
            "schema": self.schema,
            "search_mode": self.search_mode,
            "seed": self.seed,
            "workload": self.workload.to_dict(),
        }

    def to_json(self, *, indent: int | None = 2) -> str:
        """Serialize the record deterministically."""

        return json.dumps(self.to_dict(), allow_nan=False, indent=indent, sort_keys=True)


@dataclass(frozen=True)
class GSCBenchmarkRun:
    """Pair the sampling result with its resource record."""

    result: GSCSamplingResult
    record: GSCBenchmarkRecord


def build_gsc_benchmark_record(
    result: GSCSamplingResult,
    *,
    source_circuits: Iterable[Any],
    classical_runtime_seconds: float,
    peak_memory_bytes: int,
    seed: int | None,
) -> GSCBenchmarkRecord:
    """Build a validated benchmark record from an existing sampling result."""

    if not math.isfinite(classical_runtime_seconds) or classical_runtime_seconds < 0:
        raise ValueError("classical_runtime_seconds must be finite and nonnegative")
    if peak_memory_bytes < 0:
        raise ValueError("peak_memory_bytes must be nonnegative")
    circuits = tuple(source_circuits)
    if len(circuits) != len(result.observations):
        raise ValueError("source_circuits must match the executed observations")
    depths: list[int] = []
    operation_counts: Counter[str] = Counter()
    for circuit, observation in zip(circuits, result.observations, strict=True):
        depth = circuit.depth()
        if not isinstance(depth, int) or depth < 0:
            raise ValueError("source circuit depths must be nonnegative integers")
        depths.append(depth)
        num_qubits = result.num_qubits
        pauli_weight = sum(
            x_bit or z_bit
            for x_bit, z_bit in zip(
                observation.input_pauli[:num_qubits],
                observation.input_pauli[num_qubits:],
                strict=True,
            )
        )
        dagger_index = 2 * num_qubits
        unitary_index = dagger_index + pauli_weight + 1
        for index, instruction in enumerate(circuit.data):
            if index == dagger_index:
                operation_counts["queried_unitary_dagger"] += 1
            elif index == unitary_index:
                operation_counts["queried_unitary"] += 1
            else:
                operation_counts[instruction.operation.name] += 1
    depth_histogram = Counter(depths)
    workload = GSCBenchmarkWorkload(
        num_qubits=result.num_qubits,
        lagrangian_count=lagrangian_count(result.num_qubits),
        circuits_executed=len(circuits),
        shots_per_circuit=result.shots_per_pauli,
        total_shots=sum(observation.shots for observation in result.observations),
        candidate_pairs_checked=result.candidate_pairs_checked,
        source_circuit_depth_max=max(depths),
        source_circuit_depth_total=sum(depths),
        source_circuit_depth_histogram=tuple(sorted(depth_histogram.items())),
        source_operation_counts=tuple(sorted(operation_counts.items())),
    )
    return GSCBenchmarkRecord(
        schema=GSC_BENCHMARK_SCHEMA,
        workload=workload,
        measurements=GSCBenchmarkMeasurements(
            classical_runtime_seconds=float(classical_runtime_seconds),
            peak_memory_bytes=peak_memory_bytes,
        ),
        seed=seed,
        search_mode=result.search_mode,
        has_candidate_witness=result.has_candidate_witness,
    )


def run_gsc_sampling_benchmark(
    unitary: Any,
    *,
    unitary_dagger: Any | None = None,
    sampler: Any | None = None,
    shots: int = 1024,
    leakage_threshold: float = 0.01,
    seed: int | None = None,
    batch_size: int | None = None,
) -> GSCBenchmarkRun:
    """Run GSC discovery and measure host runtime and Python peak memory.

    A benchmark refuses to run while another ``tracemalloc`` session is active
    because resetting or stopping that shared process-global session would make
    either measurement ambiguous.
    """

    if tracemalloc.is_tracing():
        raise RuntimeError("cannot benchmark while tracemalloc is already active")
    num_qubits = getattr(unitary, "num_qubits", None)
    if not isinstance(num_qubits, int) or num_qubits < 1:
        raise TypeError("unitary must be a Qiskit circuit or instruction acting on qubits")
    lagrangians = enumerate_lagrangians(num_qubits)
    required_inputs = tuple(
        dict.fromkeys(label for lagrangian in lagrangians for label in lagrangian.basis)
    )
    experiments = build_semi_clifford_test_circuits(
        unitary,
        unitary_dagger=unitary_dagger,
        input_paulis=required_inputs,
    )
    source_circuits = tuple(circuit for _, circuit in experiments)
    tracemalloc.start()
    started = time.perf_counter()
    try:
        result = run_gsc_sampling_test(
            unitary,
            unitary_dagger=unitary_dagger,
            sampler=sampler,
            shots=shots,
            leakage_threshold=leakage_threshold,
            seed=seed,
            batch_size=batch_size,
        )
        elapsed = time.perf_counter() - started
        _, peak_memory_bytes = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    record = build_gsc_benchmark_record(
        result,
        source_circuits=source_circuits,
        classical_runtime_seconds=elapsed,
        peak_memory_bytes=peak_memory_bytes,
        seed=seed,
    )
    return GSCBenchmarkRun(result=result, record=record)
