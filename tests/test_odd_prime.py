from itertools import product

import numpy as np
import pytest

from generalized_semi_clifford import (
    enumerate_prime_lagrangians,
    is_prime_symplectic,
    normalize_prime_label,
    prime_lagrangian_count,
    prime_standard_form,
    prime_symplectic_pairing,
    weyl_product_phase_exponent,
)


def _qutrit_weyl(label: tuple[int, int]) -> np.ndarray:
    prime = 3
    omega = np.exp(2j * np.pi / prime)
    shift = np.roll(np.eye(prime, dtype=np.complex128), 1, axis=0)
    phase = np.diag(omega ** np.arange(prime))
    x, z = label
    scalar = omega ** (pow(2, -1, prime) * x * z % prime)
    return scalar * np.linalg.matrix_power(shift, x) @ np.linalg.matrix_power(phase, z)


def test_qutrit_weyl_product_and_commutation_phases() -> None:
    prime = 3
    omega = np.exp(2j * np.pi / prime)
    labels = tuple(product(range(prime), repeat=2))

    for left, right in product(labels, repeat=2):
        total = tuple((a + b) % prime for a, b in zip(left, right, strict=True))
        product_phase = weyl_product_phase_exponent(left, right, prime)
        pairing = prime_symplectic_pairing(left, right, prime)
        np.testing.assert_allclose(
            _qutrit_weyl(left) @ _qutrit_weyl(right),
            omega**product_phase * _qutrit_weyl(total),
            atol=1e-12,
        )
        np.testing.assert_allclose(
            _qutrit_weyl(left) @ _qutrit_weyl(right),
            omega**pairing * _qutrit_weyl(right) @ _qutrit_weyl(left),
            atol=1e-12,
        )


def test_qutrit_shift_and_phase_have_documented_pairing() -> None:
    shift = (1, 0)
    phase = (0, 1)

    assert prime_symplectic_pairing(shift, phase, 3) == 2
    assert weyl_product_phase_exponent(shift, phase, 3) == 1
    assert prime_symplectic_pairing(phase, shift, 3) == 1


def test_prime_standard_form_and_symplectic_maps() -> None:
    assert prime_standard_form(1, 3) == ((0, 2), (1, 0))
    assert is_prime_symplectic(((1, 0), (0, 1)), 3)
    assert is_prime_symplectic(((0, 2), (1, 0)), 3)
    assert not is_prime_symplectic(((1, 0), (0, 2)), 3)


def test_labels_are_canonicalized_modulo_prime() -> None:
    assert normalize_prime_label((4, -1, 8, -6), 5) == (4, 4, 3, 4)
    assert prime_symplectic_pairing((1, 0), (0, 1), 5) == 4


@pytest.mark.parametrize("prime", [2, 9, 15])
def test_invalid_odd_prime_is_rejected(prime: int) -> None:
    with pytest.raises(ValueError, match="odd prime"):
        normalize_prime_label((1, 0), prime)


def test_invalid_label_and_matrix_shapes_are_rejected() -> None:
    with pytest.raises(ValueError, match="positive even width"):
        normalize_prime_label((1, 0, 2), 3)
    with pytest.raises(ValueError, match="same width"):
        prime_symplectic_pairing((1, 0), (1, 0, 0, 0), 3)
    with pytest.raises(ValueError, match="common nonzero width"):
        is_prime_symplectic(((1, 0), (0,)), 3)
    with pytest.raises(TypeError, match="integers"):
        is_prime_symplectic(((1.0, 0), (0, 1)), 3)


def test_one_qutrit_lagrangians_are_the_four_projective_lines() -> None:
    lagrangians = enumerate_prime_lagrangians(1, 3)

    assert tuple(item.basis for item in lagrangians) == (
        ((1, 0),),
        ((1, 1),),
        ((1, 2),),
        ((0, 1),),
    )
    assert lagrangians[1].elements == ((0, 0), (1, 1), (2, 2))
    assert all(item.prime == 3 and item.num_qudits == 1 for item in lagrangians)


def test_two_qutrit_enumeration_has_independent_count_and_isotropy_checks() -> None:
    lagrangians = enumerate_prime_lagrangians(2, 3)

    assert len(lagrangians) == 40
    assert len({item.basis for item in lagrangians}) == 40
    assert all(len(item.basis) == 2 and len(item.elements) == 9 for item in lagrangians)
    assert all(
        prime_symplectic_pairing(left, right, 3) == 0
        for item in lagrangians
        for left, right in product(item.elements, repeat=2)
    )


def test_prime_lagrangian_count_and_enumeration_caps() -> None:
    assert prime_lagrangian_count(1, 3) == 4
    assert prime_lagrangian_count(2, 5) == 156

    with pytest.raises(ValueError, match="primes at most 5"):
        enumerate_prime_lagrangians(1, 7)
    with pytest.raises(ValueError, match="at most 2 qudits"):
        enumerate_prime_lagrangians(3, 3)
