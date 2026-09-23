import numpy as np
import pytest

from generalized_semi_clifford import (
    dense_pauli_transfer,
    extract_sparse_pauli_transfer,
)

IDENTITY = np.eye(2, dtype=complex)
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)


def test_sparse_clifford_transfer_matches_dense_reference() -> None:
    transfer = extract_sparse_pauli_transfer(H)
    dense = dense_pauli_transfer(H)

    assert transfer.num_qubits == 1
    assert transfer.retained_entries == 4
    assert transfer.total_entries == 16
    assert transfer.max_discarded_l2_mass == pytest.approx(0.0, abs=1e-28)
    assert np.allclose(transfer.to_dense(), dense, atol=1e-15)
    assert not transfer.to_dense().flags.writeable
    assert transfer.coefficient((1, 0), (0, 1)) == pytest.approx(1.0)
    assert transfer.coefficient((0, 1), (1, 0)) == pytest.approx(1.0)


def test_two_qubit_sparse_transfer_cross_checks_every_coefficient() -> None:
    cnot = np.array(
        [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]],
        dtype=complex,
    )
    transfer = extract_sparse_pauli_transfer(cnot)

    assert transfer.retained_entries == 16
    assert transfer.total_entries == 256
    assert np.array_equal(transfer.to_dense(), dense_pauli_transfer(cnot))


def test_threshold_metadata_quantifies_omitted_mass() -> None:
    angle = 1e-3
    pauli_x = np.array([[0, 1], [1, 0]], dtype=complex)
    rotation = np.cos(angle) * IDENTITY - 1j * np.sin(angle) * pauli_x
    transfer = extract_sparse_pauli_transfer(rotation, coefficient_tolerance=0.01)
    dense = dense_pauli_transfer(rotation)

    difference = dense - transfer.to_dense()
    per_column_mass = np.sum(np.square(difference), axis=0)
    assert np.allclose(transfer.discarded_l2_mass, per_column_mass)
    assert transfer.max_discarded_l2_mass > 0
    assert transfer.coefficient_tolerance == 0.01
    assert transfer.unitary_tolerance == 1e-9
    assert transfer.unitary_residual < transfer.unitary_tolerance
    assert "exhaustive" in transfer.arithmetic


def test_sparse_transfer_validates_limits_tolerances_and_labels() -> None:
    with pytest.raises(ValueError, match="coefficient_tolerance"):
        extract_sparse_pauli_transfer(IDENTITY, coefficient_tolerance=0)
    with pytest.raises(ValueError, match="dense Pauli-transfer limit"):
        extract_sparse_pauli_transfer(np.eye(4), max_qubits=1)
    with pytest.raises(ValueError, match="not unitary"):
        dense_pauli_transfer(2 * IDENTITY)

    transfer = extract_sparse_pauli_transfer(IDENTITY)
    with pytest.raises(KeyError, match="input"):
        transfer.coefficient((1, 0, 0, 0), (0, 0))
    with pytest.raises(KeyError, match="output"):
        transfer.coefficient((0, 0), (1, 0, 0, 0))
