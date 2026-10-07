import pytest

sympy = pytest.importorskip("sympy")

from generalized_semi_clifford.exact_cyclotomic import (  # noqa: E402
    CyclotomicField,
    CyclotomicMatrix,
)
from generalized_semi_clifford.exact_odd_prime import (  # noqa: E402
    ExactPrimeWitnessStatus,
    exact_prime_weyl_matrix,
    verify_exact_prime_lagrangian_witness,
)
from generalized_semi_clifford.odd_prime import (  # noqa: E402
    PrimeLagrangian,
    prime_lagrangian_from_basis,
)


def _qutrit_fourier(field: CyclotomicField) -> CyclotomicMatrix:
    omega = field.root_power(field.order // 3)
    scale = field.scalar(1 / sympy.sqrt(3))
    return CyclotomicMatrix.from_rows(
        field,
        (
            (scale * omega ** (row * column) for column in range(3))
            for row in range(3)
        ),
    )


def test_exact_qutrit_weyl_product_uses_documented_phase() -> None:
    field = CyclotomicField(24)
    omega = field.root_power(8)

    left = exact_prime_weyl_matrix(field, (1, 0), 3)
    right = exact_prime_weyl_matrix(field, (0, 1), 3)
    total = exact_prime_weyl_matrix(field, (1, 1), 3)

    assert left @ right == CyclotomicMatrix(
        field,
        tuple(tuple(omega * entry for entry in row) for row in total.rows),
    )


def test_exact_qutrit_fourier_witness_is_verified_or_rejected_per_pair() -> None:
    field = CyclotomicField(24)
    z_line = prime_lagrangian_from_basis(((0, 1),), 1, 3)
    x_line = prime_lagrangian_from_basis(((1, 0),), 1, 3)

    accepted = verify_exact_prime_lagrangian_witness(
        _qutrit_fourier(field), z_line, x_line
    )
    rejected = verify_exact_prime_lagrangian_witness(
        _qutrit_fourier(field), z_line, z_line
    )

    assert accepted.status is ExactPrimeWitnessStatus.VERIFIED
    assert accepted.verified is True
    assert accepted.complete
    assert accepted.coefficients_computed == accepted.coefficients_required == 9
    assert accepted.arithmetic == "cyclotomic_exact_odd_prime"
    assert rejected.status is ExactPrimeWitnessStatus.REJECTED
    assert rejected.verified is False
    assert rejected.complete


def test_exact_qutrit_witness_caps_preserve_unknown() -> None:
    field = CyclotomicField(24)
    line = prime_lagrangian_from_basis(((0, 1),), 1, 3)

    coefficient_capped = verify_exact_prime_lagrangian_witness(
        CyclotomicMatrix.identity(field, 3),
        line,
        line,
        max_coefficients=8,
    )
    qudit_capped = verify_exact_prime_lagrangian_witness(
        CyclotomicMatrix.identity(field, 3),
        line,
        line,
        max_qudits=0,
    )

    assert coefficient_capped.status is ExactPrimeWitnessStatus.UNKNOWN
    assert coefficient_capped.verified is None
    assert not coefficient_capped.complete
    assert coefficient_capped.stop_reason == "coefficient_cap"
    assert qudit_capped.status is ExactPrimeWitnessStatus.UNKNOWN
    assert qudit_capped.stop_reason == "qudit_cap"


def test_exact_qutrit_witness_revalidates_spans_and_field() -> None:
    field = CyclotomicField(24)
    line = prime_lagrangian_from_basis(((0, 1),), 1, 3)
    malformed = PrimeLagrangian(3, 1, line.basis, ((0, 0),))

    with pytest.raises(ValueError, match="canonical span"):
        verify_exact_prime_lagrangian_witness(
            CyclotomicMatrix.identity(field, 3), malformed, line
        )
    with pytest.raises(ValueError, match="divisible"):
        exact_prime_weyl_matrix(CyclotomicField(8), (1, 0), 3)
