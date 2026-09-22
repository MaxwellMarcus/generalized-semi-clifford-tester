from __future__ import annotations

import numpy as np
import pytest

qiskit = pytest.importorskip("qiskit")

from generalized_semi_clifford import (  # noqa: E402
    binomial_proportion_upper_bound,
    enumerate_lagrangians,
    run_gsc_discovery_then_verification,
    run_gsc_sampling_test,
    run_gsc_witness_test,
    run_semi_clifford_sampling_test,
    zero_event_shots_required,
)
from generalized_semi_clifford.gsc_qiskit import (  # noqa: E402
    _output_support_index,
    _support_aware_output_indices,
    empirical_lagrangian_leakage,
    find_gsc_sampling_witness,
)
from generalized_semi_clifford.semi_clifford_qiskit import (  # noqa: E402
    bell_counts_to_observation,
)


def test_exact_binomial_upper_bound_for_zero_events() -> None:
    bound = binomial_proportion_upper_bound(0, 100, confidence_level=0.95)
    assert bound == pytest.approx(1.0 - 0.05 ** (1.0 / 100))


def test_zero_event_shot_planner_meets_familywise_target() -> None:
    shots = zero_event_shots_required(
        0.01,
        confidence_level=0.95,
        simultaneous_tests=3,
    )
    per_test_confidence = 1.0 - 0.05 / 3
    assert binomial_proportion_upper_bound(
        0,
        shots,
        confidence_level=per_test_confidence,
    ) <= 0.01
    assert binomial_proportion_upper_bound(
        0,
        shots - 1,
        confidence_level=per_test_confidence,
    ) > 0.01


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
    assert result.search_mode == "support-aware-lagrangian-discovery"
    assert result.maximum_leakage_upper_bound is None
    assert not result.has_confidence_certified_witness


def test_fixed_identity_witness_gets_familywise_confidence_bound() -> None:
    from qiskit import QuantumCircuit

    identity = QuantumCircuit(2)
    z_basis = ((0, 0, 1, 0), (0, 0, 0, 1))
    result = run_gsc_witness_test(
        identity,
        z_basis,
        z_basis,
        shots=1024,
        leakage_threshold=0.01,
        confidence_level=0.95,
        seed=4,
    )
    assert result.maximum_leakage_upper_bound is not None
    assert result.maximum_leakage_upper_bound < 0.01
    assert result.has_confidence_certified_witness


def test_discovery_then_verification_uses_fresh_samples_and_reports_costs() -> None:
    from qiskit import QuantumCircuit

    result = run_gsc_discovery_then_verification(
        QuantumCircuit(2),
        discovery_shots=32,
        verification_shots=1024,
        leakage_threshold=0.01,
        confidence_level=0.95,
        discovery_seed=4,
        verification_seed=5,
        discovery_batch_size=4,
    )

    assert result.has_discovered_candidate
    assert result.verification is not None
    assert result.has_confidence_certified_witness
    assert result.discovery.search_mode == "support-aware-lagrangian-discovery"
    assert result.verification.search_mode == "candidate-witness-verification"
    assert result.discovery.observations is not result.verification.observations
    assert result.discovery_shot_cost == 13 * 32
    assert result.verification_shot_cost == 2 * 1024


def test_discovery_then_verification_keeps_rejected_candidate_without_claim(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from qiskit import QuantumCircuit

    import generalized_semi_clifford.gsc_qiskit as gsc_qiskit

    original_verification = gsc_qiskit.run_gsc_witness_test

    def reject_with_wrong_output(unitary, input_basis, output_basis, **kwargs):
        del output_basis
        z_basis = ((0, 0, 1, 0), (0, 0, 0, 1))
        return original_verification(unitary, input_basis, z_basis, **kwargs)

    monkeypatch.setattr(gsc_qiskit, "run_gsc_witness_test", reject_with_wrong_output)
    result = run_gsc_discovery_then_verification(
        QuantumCircuit(2),
        discovery_shots=16,
        verification_shots=64,
        leakage_threshold=0.01,
        discovery_seed=2,
        verification_seed=3,
    )

    assert result.has_discovered_candidate
    assert result.verification is not None
    assert not result.verification.has_candidate_witness
    assert not result.has_confidence_certified_witness
    assert result.verification_shot_cost == 2 * 64


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


def _bitstring_for_pauli(label: tuple[int, ...]) -> str:
    num_qubits = len(label) // 2
    classical_bits = label[num_qubits:] + label[:num_qubits]
    return "".join(str(bit) for bit in reversed(classical_bits))


def test_support_aware_candidates_are_complete_against_exhaustive_scoring() -> None:
    rng = np.random.default_rng(20260922)
    outputs = enumerate_lagrangians(2)
    input_lagrangian = outputs[7]
    support_index = _output_support_index(outputs)
    paulis = tuple(
        tuple((value >> bit) & 1 for bit in range(4)) for value in range(16)
    )
    saw_reduction = False

    for _ in range(24):
        observations = []
        for input_pauli in input_lagrangian.basis:
            raw_counts = rng.integers(0, 9, size=len(paulis))
            if not raw_counts.any():
                raw_counts[0] = 1
            counts = {
                _bitstring_for_pauli(label): int(count)
                for label, count in zip(paulis, raw_counts, strict=True)
                if count
            }
            observations.append(bell_counts_to_observation(input_pauli, counts))
        samples = tuple(observations)

        for threshold in (0.0, 0.1, 0.25, 0.5):
            candidates = set(
                _support_aware_output_indices(
                    samples,
                    outputs,
                    leakage_threshold=threshold,
                    support_index=support_index,
                )
            )
            feasible = {
                index
                for index, output in enumerate(outputs)
                if max(
                    empirical_lagrangian_leakage(observation, output)
                    for observation in samples
                )
                <= threshold
            }
            assert feasible <= candidates
            saw_reduction |= len(candidates) < len(outputs)

            witness, best, _, _ = find_gsc_sampling_witness(
                samples,
                input_lagrangians=(input_lagrangian,),
                output_lagrangians=outputs,
                leakage_threshold=threshold,
            )
            assert (witness is not None) == bool(feasible)
            if witness is not None:
                assert witness.output_lagrangian in tuple(
                    outputs[index] for index in feasible
                )
            else:
                exhaustive_best = min(
                    max(
                        empirical_lagrangian_leakage(observation, output)
                        for observation in samples
                    )
                    for output in outputs
                )
                assert best == pytest.approx(exhaustive_best)

    assert saw_reduction
