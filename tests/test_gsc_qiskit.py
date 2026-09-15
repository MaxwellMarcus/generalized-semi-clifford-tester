from __future__ import annotations

import numpy as np
import pytest

qiskit = pytest.importorskip("qiskit")

from generalized_semi_clifford import (  # noqa: E402
    run_gsc_sampling_test,
    run_gsc_witness_test,
    run_semi_clifford_sampling_test,
)


def _permutation_circuit(permutation: list[int]):
    from qiskit import QuantumCircuit
    from qiskit.circuit.library import UnitaryGate

    dimension = len(permutation)
    matrix = np.zeros((dimension, dimension), dtype=np.complex128)
    matrix[permutation, np.arange(dimension)] = 1
    circuit = QuantumCircuit(dimension.bit_length() - 1)
    circuit.append(UnitaryGate(matrix), circuit.qubits)
    return circuit


def test_gsc_discovery_finds_identity_masas() -> None:
    from qiskit import QuantumCircuit

    identity = QuantumCircuit(2)
    result = run_gsc_sampling_test(identity, shots=32, seed=4, batch_size=4)
    assert result.has_candidate_witness
    assert result.witness is not None
    assert result.witness.maximum_empirical_leakage == 0.0
    assert result.search_mode == "exhaustive-lagrangian-discovery"


def test_permutation_gate_has_diagonal_gsc_witness() -> None:
    # A computational-basis permutation normalizes the diagonal algebra, even
    # when it does not map individual Z Paulis to individual Paulis.
    unitary = _permutation_circuit([6, 2, 3, 4, 0, 5, 1, 7])
    z_basis = (
        (0, 0, 0, 1, 0, 0),
        (0, 0, 0, 0, 1, 0),
        (0, 0, 0, 0, 0, 1),
    )
    result = run_gsc_witness_test(
        unitary,
        z_basis,
        z_basis,
        shots=256,
        leakage_threshold=0.0,
        seed=8,
    )
    assert result.has_candidate_witness
    assert result.best_empirical_leakage == 0.0
    assert len(result.observations) == 3


def test_permutation_fixture_is_not_detected_as_semi_clifford() -> None:
    unitary = _permutation_circuit([6, 2, 3, 4, 0, 5, 1, 7])
    result = run_semi_clifford_sampling_test(
        unitary,
        shots=512,
        pauli_probability_threshold=0.99,
        seed=11,
    )
    assert not result.has_candidate_witness


def test_gsc_witness_rejects_wrong_output_masa() -> None:
    unitary = _permutation_circuit([6, 2, 3, 4, 0, 5, 1, 7])
    z_basis = (
        (0, 0, 0, 1, 0, 0),
        (0, 0, 0, 0, 1, 0),
        (0, 0, 0, 0, 0, 1),
    )
    x_basis = (
        (1, 0, 0, 0, 0, 0),
        (0, 1, 0, 0, 0, 0),
        (0, 0, 1, 0, 0, 0),
    )
    result = run_gsc_witness_test(
        unitary,
        z_basis,
        x_basis,
        shots=256,
        leakage_threshold=0.01,
        seed=8,
    )
    assert not result.has_candidate_witness
    assert result.best_empirical_leakage > 0.01
