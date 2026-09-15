import numpy as np
import pytest

from generalized_semi_clifford import (
    GSCStatus,
    GSCWitness,
    Lagrangian,
    check_gsc_naive,
    verify_gsc_witness,
)

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def test_identity_has_a_verifiable_witness() -> None:
    result = check_gsc_naive(I2)

    assert result.status is GSCStatus.GSC
    assert result.is_gsc is True
    assert result.witness is not None
    assert result.best_leakage == pytest.approx(0.0)
    assert verify_gsc_witness(I2, result.witness)


def test_t_gate_is_gsc_because_it_preserves_the_z_masa() -> None:
    t_gate = np.diag([1, np.exp(1j * np.pi / 4)])
    result = check_gsc_naive(t_gate)

    assert result.status is GSCStatus.GSC
    assert result.witness is not None
    z_label = (0, 1)
    assert z_label in result.witness.input_lagrangian.elements
    assert z_label in result.witness.output_lagrangian.elements


def test_generic_one_qubit_rotation_is_not_gsc() -> None:
    angle = 0.37
    axis = (X + Y + Z) / np.sqrt(3)
    unitary = np.cos(angle) * I2 - 1j * np.sin(angle) * axis
    result = check_gsc_naive(unitary)

    assert result.status is GSCStatus.NOT_GSC
    assert result.is_gsc is False
    assert result.lagrangians_checked == 3
    assert result.candidate_pairs_checked == 9
    assert result.best_leakage is not None and result.best_leakage > result.tolerance


def test_large_input_returns_unknown_without_searching() -> None:
    result = check_gsc_naive(np.eye(16), max_qubits=3)

    assert result.status is GSCStatus.UNKNOWN
    assert result.is_gsc is None
    assert result.lagrangians_checked == 0
    assert result.best_leakage is None


def test_rejects_invalid_matrices_and_options() -> None:
    with pytest.raises(ValueError, match="not unitary"):
        check_gsc_naive(np.ones((2, 2)))
    with pytest.raises(ValueError, match="power of two"):
        check_gsc_naive(np.eye(3))
    with pytest.raises(ValueError, match="positive"):
        check_gsc_naive(I2, tolerance=0)


def test_witness_verifier_rejects_malformed_lagrangians() -> None:
    valid_result = check_gsc_naive(I2)
    assert valid_result.witness is not None
    malformed_output = Lagrangian(
        num_qubits=1,
        basis=((1, 0),),
        elements=((0, 0),),
    )
    malformed_witness = GSCWitness(
        input_lagrangian=valid_result.witness.input_lagrangian,
        output_lagrangian=malformed_output,
        max_leakage=0.0,
    )

    assert not verify_gsc_witness(I2, malformed_witness)
