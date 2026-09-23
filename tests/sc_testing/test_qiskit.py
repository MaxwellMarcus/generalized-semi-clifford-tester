from types import SimpleNamespace

import numpy as np
import pytest

pytest.importorskip("qiskit")

from qiskit import QuantumCircuit  # noqa: E402
from qiskit.circuit import Parameter  # noqa: E402
from qiskit.circuit.library import UnitaryGate  # noqa: E402
from qiskit.primitives import StatevectorSampler  # noqa: E402
from qiskit.quantum_info import Statevector  # noqa: E402

from generalized_semi_clifford.sc_testing import (  # noqa: E402
    SCDecision,
    build_choi_bell_sampling_circuit,
    decode_choi_bell_outcome,
    run_choi_sc_test,
)


def _t_gate():
    circuit = QuantumCircuit(1)
    circuit.t(0)
    return circuit


def test_circuit_uses_two_forward_queries_and_no_inverse():
    gate = UnitaryGate(np.diag([1, np.exp(1j * np.pi / 4)]))
    circuit = build_choi_bell_sampling_circuit(gate)
    assert circuit.num_qubits == 4
    assert circuit.num_clbits == 4
    assert circuit.count_ops()["unitary"] == 2
    assert circuit.metadata["inverse_queries_per_shot"] == 0
    assert circuit.cregs[0].name == "bell"
    assert build_choi_bell_sampling_circuit(gate, measure=False).num_clbits == 0


def test_circuit_distribution_matches_direct_complex_choi_bell_amplitudes():
    gate = _t_gate()
    circuit = build_choi_bell_sampling_circuit(gate, measure=False)
    measured = Statevector.from_instruction(circuit).probabilities_dict()
    from generalized_semi_clifford.pauli import pauli_matrix

    # For a one-qubit diagonal gate there is no output/reference endian
    # ambiguity in its Choi vector. Complex phases test transpose signs.
    choi = np.array([1, 0, 0, np.exp(1j * np.pi / 4)]) / np.sqrt(2)
    for outcome in range(16):
        bitstring = f"{outcome:04b}"
        label = decode_choi_bell_outcome(bitstring, 1)
        expected = abs(choi @ pauli_matrix(label) @ choi) ** 2 / 4
        assert measured.get(bitstring, 0) == pytest.approx(expected)


def test_t_gate_end_to_end_with_theorem_budget_and_query_accounting():
    result = run_choi_sc_test(
        _t_gate(),
        hierarchy_level=3,
        seed=12,
        assume_hierarchy_membership=True,
        assume_ideal_sampling=True,
        batch_shots=17,
    )
    assert result.decision is SCDecision.SC_SUPPORTED
    assert result.analysis.projected_input_basis == ((0, 1),)
    assert result.analysis.samples_used == result.plan.bell_difference_samples
    assert result.unitary_queries == result.plan.unitary_queries
    assert result.circuit_shots == 2 * result.analysis.samples_used
    assert result.sampler_jobs == (result.circuit_shots + 16) // 17
    assert result.inverse_queries == 0
    assert result.conditional_false_accept_bound == 0.01
    assert result.conditional_sc_distance_lower_bound is None


def test_asymmetric_two_qubit_gate_preserves_qiskit_reference_indexing():
    circuit = QuantumCircuit(2)
    circuit.x(0)
    circuit.t(1)
    result = run_choi_sc_test(
        circuit,
        hierarchy_level=3,
        seed=23,
        assume_hierarchy_membership=True,
        assume_ideal_sampling=True,
    )
    assert result.decision is SCDecision.SC_SUPPORTED
    assert set(result.analysis.projected_input_basis) == {
        (1, 0, 0, 0),
        (0, 0, 1, 0),
        (0, 0, 0, 1),
    }


def test_gsc_but_not_sc_is_never_reported_as_non_gsc():
    # Basis permutation: fixes 0..4 and cycles 5,6,7. GSC, not SC.
    # No claim is made that this fixture belongs to a specified C_k.
    matrix = np.zeros((8, 8))
    matrix[[0, 1, 2, 3, 4, 6, 7, 5], np.arange(8)] = 1
    result = run_choi_sc_test(
        UnitaryGate(matrix),
        hierarchy_level=3,
        difference_samples=128,
        assume_ideal_sampling=True,
        seed=31,
    )
    assert result.decision is SCDecision.NON_SC
    assert result.analysis.maximum_input_isotropic_dimension == 1
    assert result.conditional_sc_distance_lower_bound is None
    assert result.conditional_false_accept_bound is None
    assert "not a non-GSC" in result.message


@pytest.mark.parametrize("hierarchy,ideal", [(False, False), (True, False), (False, True)])
def test_unacknowledged_premises_do_not_get_positive_guarantees(hierarchy, ideal):
    result = run_choi_sc_test(
        _t_gate(),
        hierarchy_level=3,
        seed=1,
        assume_hierarchy_membership=hierarchy,
        assume_ideal_sampling=ideal,
    )
    assert result.decision is SCDecision.INCONCLUSIVE
    assert result.conditional_false_accept_bound is None


def test_short_run_and_caps_do_not_claim_success_or_rejection():
    short = run_choi_sc_test(
        _t_gate(),
        hierarchy_level=3,
        difference_samples=1,
        seed=1,
        assume_hierarchy_membership=True,
        assume_ideal_sampling=True,
    )
    assert short.decision is SCDecision.INCONCLUSIVE
    assert short.unitary_queries == 4
    assert short.conditional_false_accept_bound is None
    capped = run_choi_sc_test(_t_gate(), hierarchy_level=10, max_difference_samples=2)
    assert capped.analysis is None
    assert capped.unitary_queries == capped.sampler_jobs == 0
    large = run_choi_sc_test(QuantumCircuit(4), hierarchy_level=2)
    assert large.decision is SCDecision.INCONCLUSIVE
    assert large.analysis is None


def test_persistent_rng_makes_results_independent_of_batch_partition(monkeypatch):
    import generalized_semi_clifford.sc_testing.qiskit as implementation

    original = implementation.analyze_choi_bell_differences
    captured = []

    def capture(n, samples):
        values = tuple(samples)
        captured.append(values)
        return original(n, values)

    monkeypatch.setattr(implementation, "analyze_choi_bell_differences", capture)
    for batch_shots in (1, 17, 10000):
        run_choi_sc_test(_t_gate(), hierarchy_level=3, seed=55, batch_shots=batch_shots)
    assert captured[0] == captured[1] == captured[2]
    assert len(set(captured[0])) > 1  # Repeated seed=55 per one-shot job would give only zero.


class _RecordingSampler:
    def __init__(self):
        self.inner = StatevectorSampler(seed=np.random.default_rng(9))
        self.shots = []

    def run(self, circuits, *, shots):
        self.shots.append(shots)
        return self.inner.run(circuits, shots=shots)


def test_custom_sampler_bypasses_only_local_simulation_cap():
    sampler = _RecordingSampler()
    result = run_choi_sc_test(
        _t_gate(),
        hierarchy_level=3,
        sampler=sampler,
        max_simulation_qubits=1,
        difference_samples=5,
        batch_shots=3,
    )
    assert sampler.shots == [3, 3, 3, 1]
    assert result.analysis.samples_used == 5
    assert result.unitary_queries == 20


@pytest.mark.parametrize("failure", ["counts_only", "short", "malformed"])
def test_invalid_sampler_results_are_not_silently_accepted(failure):
    class BadSampler:
        def run(self, circuits, *, shots):
            if failure == "counts_only":
                data = SimpleNamespace(bell=SimpleNamespace(get_counts=lambda: {"0000": shots}))
            else:
                strings = ["0000"] * (shots - 1) if failure == "short" else ["001x"] * shots
                data = SimpleNamespace(bell=SimpleNamespace(get_bitstrings=lambda: strings))
            return SimpleNamespace(result=lambda: [SimpleNamespace(data=data)])

    with pytest.raises(ValueError):
        run_choi_sc_test(_t_gate(), hierarchy_level=3, sampler=BadSampler(), difference_samples=2)


def test_nonideal_negative_has_no_decision_or_distance():
    gate = QuantumCircuit(1)
    gate.rx(0.71, 0)
    gate.rz(0.43, 0)
    result = run_choi_sc_test(gate, hierarchy_level=3, difference_samples=100, seed=4)
    assert not result.analysis.semi_clifford_compatible
    assert result.decision is SCDecision.INCONCLUSIVE
    assert result.conditional_sc_distance_lower_bound is None


def test_invalid_circuits_parameters_and_labels():
    measured = QuantumCircuit(1, 1)
    measured.measure(0, 0)
    parameterized = QuantumCircuit(1)
    parameterized.rx(Parameter("theta"), 0)
    reset = QuantumCircuit(1)
    reset.reset(0)
    for circuit in (measured, parameterized, parameterized.to_gate(), reset):
        with pytest.raises(ValueError):
            build_choi_bell_sampling_circuit(circuit)
    with pytest.raises(TypeError):
        build_choi_bell_sampling_circuit(np.eye(2))
    for value in ("001", "001x", 1):
        with pytest.raises(ValueError):
            decode_choi_bell_outcome(value, 1)
    for options in ({"batch_shots": 0}, {"assume_ideal_sampling": "yes"}):
        with pytest.raises(ValueError):
            run_choi_sc_test(_t_gate(), hierarchy_level=3, **options)


def test_global_phase_does_not_change_bell_probabilities():
    first = _t_gate()
    second = _t_gate()
    second.global_phase = 0.37
    p = Statevector.from_instruction(build_choi_bell_sampling_circuit(first, measure=False))
    q = Statevector.from_instruction(build_choi_bell_sampling_circuit(second, measure=False))
    assert p.probabilities() == pytest.approx(q.probabilities())
