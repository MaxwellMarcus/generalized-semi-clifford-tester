"""Record deterministic normalized seed-program growth in exact arithmetic."""

from __future__ import annotations

import json
from dataclasses import dataclass

from generalized_semi_clifford import (
    compile_conjugation_word,
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
    """Return a fixed depth prefix without extrapolating a growth bound."""

    backend = DihedralBackend()
    records: list[dict[str, int | str]] = []
    for depth in range(1, 8):
        sample = sample_conjugation_word(backend, depth=depth, word_length=2, seed=19)
        plan = compile_conjugation_word(sample.word)
        program = resynthesize_conjugation_word(
            sample.word, max_materialized_factors=1_000_000
        )
        assert evaluate_conjugation_seed_program(program, backend) == sample.value
        records.append(
            {
                "abstract_operations": plan.group_operations,
                "depth": depth,
                "expanded_seed_calls": plan.expanded_seed_uses,
                "normalized_program_length": len(program.factors),
                "normalized_seed_calls": program.seed_uses,
                "stored_dag_nodes": plan.stored_nodes,
            }
        )
    return records


def main() -> None:
    report = {
        "arithmetic": DihedralBackend.arithmetic,
        "cases": benchmark_records(),
        "schema_version": 1,
        "seed": 19,
        "word_length": 2,
    }
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
