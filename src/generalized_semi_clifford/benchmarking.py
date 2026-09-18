"""Reproducible records for small-qubit GSC sampling benchmarks."""

from __future__ import annotations

import json
import math
import time
import tracemalloc
from dataclasses import asdict, dataclass
from typing import Any

from .gsc_qiskit import GSCSamplingResult, run_gsc_sampling_test
from .lagrangian import lagrangian_count

GSC_BENCHMARK_SCHEMA = "generalized-semi-clifford/gsc-sampling-benchmark-v1"


@dataclass(frozen=True)
class GSCBenchmarkRecord:
    """Machine-readable resource record for one GSC discovery run.

    Runtime is elapsed host wall time for sampling plus classical
    post-processing. Peak memory is measured with :mod:`tracemalloc`, so it
    covers Python-managed allocations rather than total process memory.
    """

    schema: str
    num_qubits: int
    lagrangian_count: int
    circuits_executed: int
    shots_per_circuit: int
    total_shots: int
    candidate_pairs_checked: int
    classical_runtime_seconds: float
    peak_memory_bytes: int
    seed: int | None
    search_mode: str
    has_candidate_witness: bool

    def to_dict(self) -> dict[str, object]:
        """Return the versioned JSON-compatible record."""

        return asdict(self)

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
    classical_runtime_seconds: float,
    peak_memory_bytes: int,
    seed: int | None,
) -> GSCBenchmarkRecord:
    """Build a validated benchmark record from an existing sampling result."""

    if not math.isfinite(classical_runtime_seconds) or classical_runtime_seconds < 0:
        raise ValueError("classical_runtime_seconds must be finite and nonnegative")
    if peak_memory_bytes < 0:
        raise ValueError("peak_memory_bytes must be nonnegative")
    circuits = len(result.observations)
    return GSCBenchmarkRecord(
        schema=GSC_BENCHMARK_SCHEMA,
        num_qubits=result.num_qubits,
        lagrangian_count=lagrangian_count(result.num_qubits),
        circuits_executed=circuits,
        shots_per_circuit=result.shots_per_pauli,
        total_shots=sum(observation.shots for observation in result.observations),
        candidate_pairs_checked=result.candidate_pairs_checked,
        classical_runtime_seconds=float(classical_runtime_seconds),
        peak_memory_bytes=peak_memory_bytes,
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
        classical_runtime_seconds=elapsed,
        peak_memory_bytes=peak_memory_bytes,
        seed=seed,
    )
    return GSCBenchmarkRun(result=result, record=record)
