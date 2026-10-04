"""Reproducible host measurements for bounded exact GSC searches."""

from __future__ import annotations

import json
import math
import time
import tracemalloc
from dataclasses import asdict, dataclass

from .exact_cyclotomic import CyclotomicMatrix, ExactSearchResult, search_exact_lagrangian_witness

EXACT_SEARCH_BENCHMARK_SCHEMA = "generalized-semi-clifford/exact-search-benchmark-v1"


@dataclass(frozen=True)
class ExactSearchBenchmarkMeasurements:
    """Machine-dependent measurements for one exact-search run."""

    host_runtime_seconds: float
    python_peak_memory_bytes: int


@dataclass(frozen=True)
class ExactSearchBenchmarkRecord:
    """One benchmark record with exact work and host measurements separated."""

    fixture: str
    result: ExactSearchResult
    measurements: ExactSearchBenchmarkMeasurements
    schema: str = EXACT_SEARCH_BENCHMARK_SCHEMA

    def to_dict(self) -> dict[str, object]:
        """Return a deterministic, JSON-compatible benchmark payload."""

        result = self.result
        witness_coefficients = (
            result.witness.coefficients_computed if result.witness is not None else 0
        )
        return {
            "exact_workload": {
                "arithmetic": result.arithmetic,
                "candidate_pairs_checked": result.candidate_pairs_checked,
                "coefficients_computed": result.coefficients_computed,
                "field_order": result.field_order,
                "input_lagrangians_checked": result.input_lagrangians_checked,
                "lagrangians_total": result.lagrangians_total,
                "limits": {
                    "max_candidate_pairs": result.max_candidate_pairs,
                    "max_coefficients": result.max_coefficients,
                    "max_qubits": result.max_qubits,
                },
                "num_qubits": result.num_qubits,
                "witness_verification_coefficients": witness_coefficients,
            },
            "fixture": self.fixture,
            "measurements": asdict(self.measurements),
            "outcome": {
                "search_complete": result.search_complete,
                "status": result.status.value,
                "stop_reason": result.stop_reason,
                "witness_verified": result.witness.verified if result.witness else None,
            },
            "schema": self.schema,
        }

    def to_json(self, *, indent: int | None = 2) -> str:
        """Serialize the benchmark record deterministically."""

        return json.dumps(self.to_dict(), allow_nan=False, indent=indent, sort_keys=True)


def run_exact_search_benchmark(
    fixture: str,
    unitary: CyclotomicMatrix,
    *,
    max_qubits: int = 3,
    max_coefficients: int = 100_000,
    max_candidate_pairs: int = 100_000,
) -> ExactSearchBenchmarkRecord:
    """Run one bounded exact search and measure host runtime and Python memory.

    Work counters come from the exact search itself. Wall time and
    :mod:`tracemalloc` peak memory are machine-dependent observations and are
    intentionally serialized under a separate key.
    """

    if not isinstance(fixture, str) or not fixture.strip():
        raise ValueError("fixture must be a nonempty string")
    if tracemalloc.is_tracing():
        raise RuntimeError("cannot benchmark while tracemalloc is already active")

    tracemalloc.start()
    started = time.perf_counter()
    try:
        result = search_exact_lagrangian_witness(
            unitary,
            max_qubits=max_qubits,
            max_coefficients=max_coefficients,
            max_candidate_pairs=max_candidate_pairs,
        )
        elapsed = time.perf_counter() - started
        _, peak_memory_bytes = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()

    if not math.isfinite(elapsed) or elapsed < 0:  # pragma: no cover - defensive clock check
        raise RuntimeError("host runtime measurement must be finite and nonnegative")
    return ExactSearchBenchmarkRecord(
        fixture=fixture.strip(),
        result=result,
        measurements=ExactSearchBenchmarkMeasurements(
            host_runtime_seconds=elapsed,
            python_peak_memory_bytes=peak_memory_bytes,
        ),
    )
