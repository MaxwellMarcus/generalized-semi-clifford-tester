import pytest

from generalized_semi_clifford import enumerate_lagrangians, is_symplectic, lagrangian_count
from generalized_semi_clifford.lagrangian import is_lagrangian


@pytest.mark.parametrize(("num_qubits", "expected"), [(1, 3), (2, 15), (3, 135)])
def test_lagrangian_counts(num_qubits: int, expected: int) -> None:
    lagrangians = enumerate_lagrangians(num_qubits)

    assert lagrangian_count(num_qubits) == expected
    assert len(lagrangians) == expected
    assert all(is_lagrangian(item.elements, num_qubits) for item in lagrangians)
    assert all(len(item.basis) == num_qubits for item in lagrangians)


def test_invalid_lagrangian_inputs() -> None:
    with pytest.raises(ValueError, match="positive"):
        enumerate_lagrangians(0)
    with pytest.raises(ValueError, match="limited to 3"):
        enumerate_lagrangians(4)
    assert not is_lagrangian(((0, 0),), 1)
    assert not is_lagrangian(((0, 2), (0, 0)), 1)


def test_existing_symplectic_api_remains_available() -> None:
    assert is_symplectic(((0, 1), (1, 0)))
