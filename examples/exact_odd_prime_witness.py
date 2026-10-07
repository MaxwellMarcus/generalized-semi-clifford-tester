"""Verify one exact qutrit input/output Weyl-Lagrangian witness."""

from sympy import sqrt

from generalized_semi_clifford import prime_lagrangian_from_basis
from generalized_semi_clifford.exact_cyclotomic import (
    CyclotomicField,
    CyclotomicMatrix,
)
from generalized_semi_clifford.exact_odd_prime import verify_exact_prime_lagrangian_witness

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
print(
    {
        "status": result.status.value,
        "complete": result.complete,
        "prime": result.prime,
        "num_qudits": result.num_qudits,
        "coefficients_computed": result.coefficients_computed,
        "stop_reason": result.stop_reason,
    }
)
