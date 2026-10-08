"""Verify and search exact qutrit Weyl-Lagrangian witnesses."""

import json

from sympy import sqrt

from generalized_semi_clifford import prime_lagrangian_from_basis
from generalized_semi_clifford.exact_cyclotomic import (
    CyclotomicField,
    CyclotomicMatrix,
)
from generalized_semi_clifford.exact_odd_prime import (
    search_exact_prime_lagrangian_witness,
    verify_exact_prime_lagrangian_witness,
)

field = CyclotomicField(24)
omega = field.root_power(8)
scale = field.scalar(1 / sqrt(3))
qutrit_fourier = CyclotomicMatrix.from_rows(
    field,
    (
        (scale * omega ** (row * column) for column in range(3))
        for row in range(3)
    ),
)
z_line = prime_lagrangian_from_basis(((0, 1),), 1, 3)
x_line = prime_lagrangian_from_basis(((1, 0),), 1, 3)

result = verify_exact_prime_lagrangian_witness(qutrit_fourier, z_line, x_line)
print("fixed qutrit Fourier witness")
print(
    json.dumps(
        {
            "status": result.status.value,
            "complete": result.complete,
            "prime": result.prime,
            "num_qudits": result.num_qudits,
            "coefficients_computed": result.coefficients_computed,
            "stop_reason": result.stop_reason,
        },
        indent=2,
        sort_keys=True,
    )
)

phase = field.root_power(1)
non_gsc = (
    CyclotomicMatrix.from_rows(
        field,
        ((1, 0, 0), (0, 1, 0), (0, 0, phase)),
    )
    @ qutrit_fourier
    @ CyclotomicMatrix.from_rows(
        field,
        ((1, 0, 0), (0, 1, 0), (0, 0, phase)),
    )
)

for name, search in (
    ("qutrit_fourier", search_exact_prime_lagrangian_witness(qutrit_fourier, 3)),
    ("qutrit_completed_negative", search_exact_prime_lagrangian_witness(non_gsc, 3)),
    (
        "qutrit_capped",
        search_exact_prime_lagrangian_witness(non_gsc, 3, max_candidate_pairs=2),
    ),
):
    print(name)
    print(json.dumps(search.to_dict(), indent=2, sort_keys=True))
