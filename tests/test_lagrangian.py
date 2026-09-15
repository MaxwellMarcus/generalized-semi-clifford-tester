import pytest

from generalized_semi_clifford import enumerate_lagrangians, is_symplectic, lagrangian_count
from generalized_semi_clifford.lagrangian import (
    is_lagrangian,
    lagrangian_containing,
    lagrangian_from_basis,
)


@pytest.mark.parametrize(
    ("num_qubits", "expected"),
    [(1, 3), (2, 15), (3, 135), (4, 2295)],
)
def test_lagrangian_counts(num_qubits: int, expected: int) -> None:
    lagrangians = enumerate_lagrangians(num_qubits)

    assert lagrangian_count(num_qubits) == expected
    assert len(lagrangians) == expected
    assert all(is_lagrangian(item.elements, num_qubits) for item in lagrangians)
    assert all(len(item.basis) == num_qubits for item in lagrangians)


def test_invalid_lagrangian_inputs() -> None:
    with pytest.raises(ValueError, match="positive"):
        enumerate_lagrangians(0)
    with pytest.raises(ValueError, match="limited to 4"):
        enumerate_lagrangians(5)
    assert not is_lagrangian(((0, 0),), 1)
    assert not is_lagrangian(((0, 2), (0, 0)), 1)


def test_existing_symplectic_api_remains_available() -> None:
    assert is_symplectic(((0, 1), (1, 0)))


def test_lagrangian_from_basis_validates_and_builds_span() -> None:
    lagrangian = lagrangian_from_basis(((1, 0, 0, 0), (0, 1, 0, 0)), 2)
    assert len(lagrangian.basis) == 2
    assert len(lagrangian.elements) == 4
    assert is_lagrangian(lagrangian.elements, 2)


def test_lagrangian_from_basis_rejects_anticommuting_basis() -> None:
    with pytest.raises(ValueError, match="commute"):
        lagrangian_from_basis(((1, 0, 0, 0), (0, 0, 1, 0)), 2)


def test_lagrangian_containing_extends_isotropic_span() -> None:
    x0 = (1, 0, 0, 0)
    lagrangian = lagrangian_containing((x0,), 2)
    assert lagrangian is not None
    assert x0 in lagrangian.elements
    assert is_lagrangian(lagrangian.elements, 2)


def test_lagrangian_containing_rejects_anticommuting_span() -> None:
    assert lagrangian_containing(((1, 0), (0, 1)), 1) is None
