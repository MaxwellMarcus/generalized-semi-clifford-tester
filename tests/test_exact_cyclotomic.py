from fractions import Fraction

import pytest

sympy = pytest.importorskip("sympy")

from generalized_semi_clifford.exact_cyclotomic import (  # noqa: E402
    CyclotomicField,
    CyclotomicMatrix,
    ccz,
    cnot,
    common_cyclotomic_field,
    controlled_phase,
    cz,
    hadamard,
    pauli_x,
    pauli_y,
    pauli_z,
    permutation_gate,
    phase_s,
    swap,
    t_gate,
    toffoli,
)


def test_field_rejects_implicit_approximate_inputs() -> None:
    field = CyclotomicField(8)

    with pytest.raises(TypeError, match="floating-point"):
        field.scalar(0.5)
    with pytest.raises(TypeError, match="floating-point"):
        field.scalar(1j)
    with pytest.raises(TypeError, match="floating-point"):
        field.scalar(sympy.Float("0.5"))
    with pytest.raises(ValueError, match="multiple of eight"):
        CyclotomicField(12)

    assert field.scalar(Fraction(1, 2)) + field.scalar(Fraction(1, 2)) == field.one


def test_common_field_embedding_is_exact_and_declared() -> None:
    field8 = CyclotomicField(8)
    field24 = CyclotomicField(24)
    common = common_cyclotomic_field(field8, field24)

    assert common.order == 24
    assert field8.zeta.embed(common) == common.root_power(3)
    assert t_gate(field8).embed(common) == t_gate(common)


@pytest.mark.parametrize(
    "constructor",
    [
        pauli_x,
        pauli_y,
        pauli_z,
        hadamard,
        phase_s,
        t_gate,
        cnot,
        cz,
        swap,
        toffoli,
        ccz,
    ],
)
def test_fixed_gate_constructors_are_exactly_unitary(constructor) -> None:
    matrix = constructor(CyclotomicField(8))

    assert matrix.is_unitary()
    assert matrix.adjoint() @ matrix == CyclotomicMatrix.identity(matrix.field, matrix.shape[0])


def test_controlled_root_phase_and_tensor_products() -> None:
    field = CyclotomicField(16)
    gate = controlled_phase(field, 3, 1)
    product = hadamard(field).tensor(t_gate(field))

    assert gate.rows[-1][-1] == field.zeta
    assert gate.is_unitary()
    assert product.shape == (4, 4)
    assert product.is_unitary()


def test_permutation_validation_and_deliberately_nonunitary_matrix() -> None:
    field = CyclotomicField(8)
    nonunitary = CyclotomicMatrix.from_rows(field, ((1, 1), (0, 1)))

    assert not nonunitary.is_unitary()
    with pytest.raises(ValueError, match="every index"):
        permutation_gate(field, (0, 0))
    with pytest.raises(ValueError, match="equal length"):
        CyclotomicMatrix.from_rows(field, ((1, 0), (1,)))
