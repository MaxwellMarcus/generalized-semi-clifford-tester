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


def _record(
    name: str,
    result: object,
    *,
    oracle_stop_reason: str | None = None,
) -> dict[str, object]:
    record = {
        "name": name,
        "status": result.status.value,
        "arithmetic": result.arithmetic,
        "rounds_completed": result.rounds_completed,
        "rounds_required": result.rounds_required,
        "oracle_calls": result.oracle_calls,
        "exact_oracle_coefficients_computed": result.exact_oracle_coefficients_computed,
        "rigorous_false_accept_bound": result.rigorous_false_accept_bound,
        "has_violation_witness": result.witness is not None,
        "message": result.message,
    }
    if result.exact_oracle_coefficients_computed is not None:
        record["oracle_stop_reason"] = oracle_stop_reason
    return record


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

    # Each one-qubit generator expansion costs four exact coefficients.  This
    # budget funds the first expansion but not the second, so the stateful
    # exact oracle returns None and containment must preserve INCONCLUSIVE.
    capped_oracle = supported_backend.clifford_oracle(max_coefficients=4)
    inconclusive = test_clifford_conjugation_containment(
        supported_backend,
        capped_oracle,
        depth=1,
        oracle_is_exact=True,
        seed=7,
    )

    # A two-qubit Pauli expansion costs 16 exact coefficients.  T tensor I is
    # still a supported exact seed, but this budget funds only the first of the
    # four generator expansions required for a complete Clifford decision.
    # The bounded prefix is resource evidence, not a completed classification.
    two_qubit_seed = t_gate(field8).tensor(CyclotomicMatrix.identity(field8, 2))
    two_qubit_backend = ExactCyclotomicConjugationBackend(two_qubit_seed)
    two_qubit_capped_oracle = two_qubit_backend.clifford_oracle(max_coefficients=16)
    two_qubit_inconclusive = test_clifford_conjugation_containment(
        two_qubit_backend,
        two_qubit_capped_oracle,
        depth=1,
        oracle_is_exact=True,
        seed=7,
    )

    records = (
        _record("t_gate_supported", supported),
        _record("eighth_turn_violation", violating),
        _record(
            "coefficient_capped_oracle",
            inconclusive,
            oracle_stop_reason=capped_oracle.stop_reason,
        ),
        _record(
            "two_qubit_coefficient_capped_oracle",
            two_qubit_inconclusive,
            oracle_stop_reason=two_qubit_capped_oracle.stop_reason,
        ),
    )
    print(json.dumps(records, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
