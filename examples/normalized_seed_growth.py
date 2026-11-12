"""Record deterministic normalized seed-program growth in exact arithmetic."""

from __future__ import annotations

import json
from dataclasses import dataclass

from generalized_semi_clifford import (
    evaluate_conjugation_seed_program,
    resynthesize_conjugation_word,
    sample_conjugation_word,
)


@dataclass(frozen=True)
class DihedralBackend:
    """Exact projective coordinates for ``<Rz(2 pi / modulus), X>``."""

    modulus: int = 256
    num_qubits: int = 1
    exact: bool = True
    arithmetic: str = "exact dihedral group coordinates"
    seed: tuple[int, int] = (1, 0)

    def identity(self) -> tuple[int, int]:
        return (0, 0)

    def pauli(self, label: tuple[int, int]) -> tuple[int, int]:
        x, z = label
        return (self.modulus // 2 * z, x)

    def multiply(
        self, left: tuple[int, int], right: tuple[int, int]
    ) -> tuple[int, int]:
        a, b = left
        c, d = right
        return ((a + (-1) ** b * c) % self.modulus, b ^ d)

    def adjoint(self, value: tuple[int, int]) -> tuple[int, int]:
        a, b = value
        return ((-((-1) ** b) * a) % self.modulus, b)


def benchmark_records() -> list[dict[str, int | str]]:
    """Return maxima over a fixed parameter grid without extrapolating a bound."""

    backend = DihedralBackend()
    records: list[dict[str, int | str]] = []
    for depth in range(1, 8):
        maximum = -1
        first_maximizer: tuple[int, int, int] | None = None
        maximizer_count = 0
        for word_length in range(1, 5):
            for seed in range(32):
                sample = sample_conjugation_word(
                    backend,
                    depth=depth,
                    word_length=word_length,
                    seed=seed,
                    max_operations=200_000,
                )
                program = resynthesize_conjugation_word(
                    sample.word, max_materialized_factors=100_000_000
                )
                assert (
                    evaluate_conjugation_seed_program(
                        program, backend, max_operations=1_000_000
                    )
                    == sample.value
                )
                program_length = len(program.factors)
                if program_length > maximum:
                    maximum = program_length
                    first_maximizer = (word_length, seed, program.seed_uses)
                    maximizer_count = 1
                elif program_length == maximum:
                    maximizer_count += 1
        assert first_maximizer is not None
        word_length, seed, seed_calls = first_maximizer
        records.append(
            {
                "depth": depth,
                "first_maximizing_seed": seed,
                "first_maximizing_word_length": word_length,
                "max_normalized_program_length": maximum,
                "maximizer_count": maximizer_count,
                "normalized_seed_calls_at_first_maximizer": seed_calls,
            }
        )
    return records


def main() -> None:
    report = {
        "arithmetic": DihedralBackend.arithmetic,
        "cases_scanned": 7 * 4 * 32,
        "depths": [1, 7],
        "maxima": benchmark_records(),
        "resource_caps": {
            "group_operations_per_case": 200_000,
            "materialized_factors_per_case": 100_000_000,
            "replay_operations_per_case": 1_000_000,
        },
        "schema_version": 2,
        "seeds": [0, 31],
        "word_lengths": [1, 4],
    }
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
