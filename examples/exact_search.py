"""Run bounded exact Pauli-MASA searches with explicit work caps."""

from __future__ import annotations

import json

from generalized_semi_clifford.exact_cyclotomic import (
    CyclotomicField,
    hadamard,
    search_exact_lagrangian_witness,
    t_gate,
)


def main() -> None:
    field = CyclotomicField(8)
    positive = search_exact_lagrangian_witness(hadamard(field))
    negative_gate = t_gate(field) @ hadamard(field) @ t_gate(field)
    negative = search_exact_lagrangian_witness(negative_gate)
    capped = search_exact_lagrangian_witness(negative_gate, max_candidate_pairs=2)

    for name, result in (
        ("hadamard", positive),
        ("t_h_t", negative),
        ("t_h_t_capped", capped),
    ):
        print(name)
        print(json.dumps(result.to_dict(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
