"""Reproducible host measurements for bounded odd-prime exact searches."""

from __future__ import annotations

import json
import math
import time
import tracemalloc
from dataclasses import asdict, dataclass

from .exact_cyclotomic import CyclotomicMatrix
from .exact_odd_prime import ExactPrimeSearchResult, search_exact_prime_lagrangian_witness

EXACT_PRIME_SEARCH_BENCHMARK_SCHEMA = (
    "generalized-semi-clifford/odd-prime-exact-search-benchmark-v1"
)


@dataclass(frozen=True)
class ExactPrimeSearchBenchmarkMeasurements:
    """Machine-dependent measurements for one odd-prime exact-search run."""

    host_runtime_seconds: float
    python_peak_memory_bytes: int


@dataclass(frozen=True)
class ExactPrimeSearchBenchmarkRecord:
    """One odd-prime benchmark with exact work and host data separated."""

    fixture: str
    result: ExactPrimeSearchResult
    measurements: ExactPrimeSearchBenchmarkMeasurements
    schema: str = EXACT_PRIME_SEARCH_BENCHMARK_SCHEMA

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
                    "max_qudits": result.max_qudits,
                },
                "num_qudits": result.num_qudits,
                "prime": result.prime,
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


def run_exact_prime_search_benchmark(
    fixture: str,
    unitary: CyclotomicMatrix,
    prime: int,
    *,
    max_qudits: int = 2,
    max_coefficients: int = 100_000,
    max_candidate_pairs: int = 100_000,
) -> ExactPrimeSearchBenchmarkRecord:
    """Run one bounded odd-prime exact search with host measurements.

    Portable work counters come from the exact search. Wall time and
    :mod:`tracemalloc` peak memory remain explicitly machine-dependent.
    """

    if not isinstance(fixture, str) or not fixture.strip():
        raise ValueError("fixture must be a nonempty string")
    if tracemalloc.is_tracing():
        raise RuntimeError("cannot benchmark while tracemalloc is already active")

    tracemalloc.start()
    started = time.perf_counter()
    try:
        result = search_exact_prime_lagrangian_witness(
            unitary,
            prime,
            max_qudits=max_qudits,
            max_coefficients=max_coefficients,
            max_candidate_pairs=max_candidate_pairs,
        )
        elapsed = time.perf_counter() - started
        _, peak_memory_bytes = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()

    if not math.isfinite(elapsed) or elapsed < 0:  # pragma: no cover - defensive clock check
        raise RuntimeError("host runtime measurement must be finite and nonnegative")
    return ExactPrimeSearchBenchmarkRecord(
        fixture=fixture.strip(),
        result=result,
        measurements=ExactPrimeSearchBenchmarkMeasurements(
            host_runtime_seconds=elapsed,
            python_peak_memory_bytes=peak_memory_bytes,
        ),
    )
