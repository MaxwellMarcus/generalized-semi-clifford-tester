"""Reproduce exact Clifford-containment outcomes without numerical tolerances."""

from __future__ import annotations

import json

from generalized_semi_clifford.conjugation_sampling import (
    test_clifford_conjugation_containment,
)
from generalized_semi_clifford.exact_cyclotomic import (
    CyclotomicField,
    CyclotomicMatrix,
    ExactCyclotomicConjugationBackend,
    t_gate,
)


def _record(name: str, result: object) -> dict[str, object]:
    return {
        "name": name,
        "status": result.status.value,
        "arithmetic": result.arithmetic,
        "rounds_completed": result.rounds_completed,
        "rounds_required": result.rounds_required,
        "oracle_calls": result.oracle_calls,
        "rigorous_false_accept_bound": result.rigorous_false_accept_bound,
        "has_violation_witness": result.witness is not None,
        "message": result.message,
    }


def main() -> None:
    field8 = CyclotomicField(8)
    supported_backend = ExactCyclotomicConjugationBackend(t_gate(field8))
    supported = test_clifford_conjugation_containment(
        supported_backend,
        supported_backend.is_clifford,
        depth=1,
        oracle_is_exact=True,
        seed=7,
    )

    field16 = CyclotomicField(16)
    eighth_turn = CyclotomicMatrix.from_rows(field16, ((1, 0), (0, field16.zeta)))
    violating_backend = ExactCyclotomicConjugationBackend(eighth_turn)
    violating = test_clifford_conjugation_containment(
        violating_backend,
        violating_backend.is_clifford,
        depth=1,
        oracle_is_exact=True,
        seed=7,
    )

    # An interrupted or resource-capped exact oracle returns None.  The
    # containment test must preserve that as INCONCLUSIVE with no error bound.
    inconclusive = test_clifford_conjugation_containment(
        supported_backend,
        lambda _: None,
        depth=1,
        oracle_is_exact=True,
        seed=7,
    )

    records = (
        _record("t_gate_supported", supported),
        _record("eighth_turn_violation", violating),
        _record("indeterminate_oracle", inconclusive),
    )
    print(json.dumps(records, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
