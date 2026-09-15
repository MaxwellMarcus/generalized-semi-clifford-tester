from __future__ import annotations

from math import pi

import pytest

qiskit = pytest.importorskip("qiskit")

from generalized_semi_clifford.semi_clifford_qiskit import (  # noqa: E402
    PauliConjugationObservation,
    bell_counts_to_observation,
    build_pauli_conjugation_test_circuit,
    find_lagrangian_witness,
    run_semi_clifford_sampling_test,
    symplectic_pairing,
)


def _observation(
    input_pauli: tuple[int, ...],
    output_pauli: tuple[int, ...],
) -> PauliConjugationObservation:
    return PauliConjugationObservation(
        input_pauli=input_pauli,
        dominant_output_pauli=output_pauli,
        dominant_probability=1.0,
        shots=100,
        counts=(("0000", 100),),
    )


def test_symplectic_pairing_detects_anticommutation() -> None:
    x0 = (1, 0, 0, 0)
    x1 = (0, 1, 0, 0)
    z0 = (0, 0, 1, 0)
    assert symplectic_pairing(x0, x1) == 0
    assert symplectic_pairing(x0, z0) == 1


def test_find_lagrangian_witness_requires_commuting_independent_pairs() -> None:
    observations = (
        _observation((1, 0, 0, 0), (0, 0, 1, 0)),
        _observation((0, 1, 0, 0), (0, 0, 0, 1)),
    )
    witness = find_lagrangian_witness(
        observations,
        num_qubits=2,
        pauli_probability_threshold=0.99,
    )
    assert witness is not None
    assert len(witness.input_basis) == 2


def test_bell_counts_decode_x_and_z_coordinates() -> None:
    # For one qubit, Qiskit displays classical bits as xz.
    x_observation = bell_counts_to_observation((0, 1), {"10": 32})
    z_observation = bell_counts_to_observation((1, 0), {"01": 32})
    assert x_observation.dominant_output_pauli == (1, 0)
    assert z_observation.dominant_output_pauli == (0, 1)


def test_hadamard_circuit_maps_x_to_z_under_conjugation() -> None:
    from qiskit import QuantumCircuit
    from qiskit.primitives import StatevectorSampler

    hadamard = QuantumCircuit(1)
    hadamard.h(0)
    circuit = build_pauli_conjugation_test_circuit(hadamard, (1, 0))
    result = StatevectorSampler(seed=9).run([circuit], shots=32).result()[0]
    observation = bell_counts_to_observation((1, 0), result.data.bell.get_counts())
    assert observation.dominant_output_pauli == (0, 1)
    assert observation.dominant_probability == 1.0


def test_full_sampling_test_finds_identity_witness() -> None:
    from qiskit import QuantumCircuit

    identity = QuantumCircuit(2)
    result = run_semi_clifford_sampling_test(identity, shots=16, seed=4)
    assert result.has_candidate_witness
    assert result.witness is not None
    assert len(result.witness.input_basis) == 2


def test_full_sampling_test_rejects_generic_single_qubit_rotation() -> None:
    from qiskit import QuantumCircuit

    generic = QuantumCircuit(1)
    generic.u(pi / 3, pi / 5, pi / 7, 0)
    result = run_semi_clifford_sampling_test(
        generic,
        shots=2048,
        pauli_probability_threshold=0.95,
        seed=12,
    )
    assert not result.has_candidate_witness
