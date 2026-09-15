import pytest

from generalized_semi_clifford import identity, is_symplectic, matmul, standard_form, transpose


def test_identity_and_hadamard_actions_are_symplectic() -> None:
    assert is_symplectic(identity(2))
    assert is_symplectic(((0, 1), (1, 0)))


def test_singular_and_odd_dimensional_maps_are_not_symplectic() -> None:
    assert not is_symplectic(((1, 0), (1, 0)))
    assert not is_symplectic(identity(3))


def test_standard_form_is_involutory() -> None:
    form = standard_form(3)
    assert matmul(form, form) == identity(6)
    assert transpose(form) == form


def test_matrix_validation_rejects_nonbinary_and_ragged_input() -> None:
    with pytest.raises(ValueError, match="0 or 1"):
        matmul(((2,),), ((1,),))
    with pytest.raises(ValueError, match="common nonzero width"):
        transpose(((1, 0), (1,)))


def test_matmul_rejects_incompatible_dimensions() -> None:
    with pytest.raises(ValueError, match="inner matrix dimensions"):
        matmul(((1, 0),), ((1, 0),))
