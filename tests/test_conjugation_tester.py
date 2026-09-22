import numpy as np
import pytest

from generalized_semi_clifford import (
    DenseConjugationBackend,
    SandwichStatus,
    check_gsc_sandwich,
    evaluate_conjugation_word,
    probe_gsc_conjugation_condition,
    probe_hierarchy_distance,
    semi_clifford_distance_lower_bound,
    verify_gsc_witness,
)

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]])
Z = np.diag([1, -1])
H = (X + Z) / np.sqrt(2)
T = np.diag([1, np.exp(1j * np.pi / 4)])


def rotation(angle):
    return np.cos(angle) * I2 - 1j * np.sin(angle) * (X + Y + Z) / np.sqrt(3)


def test_dense_clifford_oracle_and_drift_detection():
    backend = DenseConjugationBackend(H)
    assert backend.clifford_membership(H)
    assert backend.clifford_membership(1j * H)
    assert not backend.clifford_membership(T)
    assert backend.clifford_membership(1.01 * H) is None
    assert backend.clifford_membership(np.full((2, 2), np.nan)) is None


def test_hierarchy_probe_distance_and_phase_invariance():
    evidence = probe_hierarchy_distance(T, level=2)
    assert len(evidence.pauli_path) == 1
    assert evidence.distance_lower_bound > 0.1
    # Compare against a known member I of C_2 as an upper bound on the true distance.
    upper = np.linalg.norm(T - np.exp(1j * np.pi / 8) * I2, "fro") / np.sqrt(2)
    assert evidence.distance_lower_bound <= upper
    phased = probe_hierarchy_distance(np.exp(0.37j) * T, level=2)
    assert phased.distance_lower_bound == pytest.approx(evidence.distance_lower_bound)
    assert probe_hierarchy_distance(T, level=3, trials=100, seed=4).distance_lower_bound == 0
    assert probe_hierarchy_distance(1j * X, level=1).distance_lower_bound == 0


@pytest.mark.parametrize("route", ["shallow", "deep"])
def test_both_routes_pass_known_sc_hierarchy_members_without_claiming_containment(route):
    unitary = np.diag([1, np.exp(1j * np.pi / 16)])  # C_5
    result = probe_gsc_conjugation_condition(
        unitary, hierarchy_level=5, route=route, samples=20, seed=7
    )
    assert result.words_checked == 20
    assert result.evidence is None
    assert not result.budget_exhausted
    assert "do not establish containment" in result.message


@pytest.mark.parametrize("route", ["shallow", "deep"])
def test_violation_is_replayable_and_only_bounds_hierarchy_restricted_sc(route):
    # This gate is exactly SC, but outside the requested low hierarchy level.
    unitary = np.diag([1, np.exp(0.37j)])
    result = probe_gsc_conjugation_condition(
        unitary, hierarchy_level=4, route=route, samples=16, seed=7
    )
    assert result.evidence is not None
    evidence = result.evidence
    assert evidence.distance_from_sc_in_ck_lower_bound > 0
    replay = evaluate_conjugation_word(evidence.sample.word, DenseConjugationBackend(unitary))
    assert np.allclose(replay, evidence.sample.value)
    assert semi_clifford_distance_lower_bound(unitary) == 0


def test_unrestricted_sc_distance_is_phase_invariant_and_bounded_by_known_sc_distance():
    unitary = rotation(0.37)
    bound = semi_clifford_distance_lower_bound(unitary)
    assert bound > 0.05
    assert bound <= np.linalg.norm(unitary - I2, "fro") / np.sqrt(2)
    assert semi_clifford_distance_lower_bound(1j * unitary) == pytest.approx(bound)
    assert semi_clifford_distance_lower_bound(H @ T @ H) == 0
    # Near-SC inputs must not be declared farther away than their actual distance.
    assert semi_clifford_distance_lower_bound(rotation(1e-5)) < 1e-5


def test_sandwich_positive_uses_direct_witness_even_when_containment_fails():
    unitary = np.diag([1, np.exp(0.37j)])
    result = check_gsc_sandwich(unitary, hierarchy_level=4, samples=8, seed=7)
    assert result.status is SandwichStatus.GSC_WITNESS
    assert result.conjugation.evidence is not None
    assert result.gsc_witness is not None
    assert verify_gsc_witness(unitary, result.gsc_witness)
    assert "not interval-certified" in result.arithmetic


def test_sandwich_far_and_inconclusive_are_distinct():
    result = check_gsc_sandwich(rotation(0.37), hierarchy_level=4, samples=4, seed=7)
    assert result.status is SandwichStatus.FAR_FROM_SC
    assert result.distance_from_sc_lower_bound > result.epsilon
    inconclusive = check_gsc_sandwich(
        rotation(0.37), hierarchy_level=4, epsilon=1.0, samples=4, seed=7
    )
    assert inconclusive.status is SandwichStatus.INCONCLUSIVE
    large = check_gsc_sandwich(np.eye(16), hierarchy_level=4)
    assert large.status is SandwichStatus.INCONCLUSIVE
    assert large.conjugation is None


def test_budget_exhaustion_is_not_passing_the_test():
    result = probe_gsc_conjugation_condition(I2, hierarchy_level=5, route="deep", max_operations=4)
    assert result.budget_exhausted
    assert result.words_checked == 0
    assert result.evidence is None


def test_numerical_failure_is_not_mislabelled_as_budget_exhaustion(monkeypatch):
    from generalized_semi_clifford import conjugation_tester

    def fail_probe(*args, **kwargs):
        raise ArithmeticError("simulated numerical drift")

    monkeypatch.setattr(conjugation_tester, "probe_hierarchy_distance", fail_probe)
    result = probe_gsc_conjugation_condition(I2, hierarchy_level=4, samples=1)
    assert result.numerical_failure
    assert not result.budget_exhausted
    assert result.words_checked == 0
    assert result.evidence is None


def test_budget_boundary_after_a_completed_sample():
    # A length-two Gamma_2 sample uses 14 group operations in this implementation.
    result = probe_gsc_conjugation_condition(
        I2, hierarchy_level=4, samples=2, word_length=2, max_operations=14, seed=7
    )
    assert result.words_checked == 1
    assert result.budget_exhausted
    assert not result.numerical_failure


@pytest.mark.parametrize("num_qubits", [2, 3])
def test_multiqubit_diagonal_hierarchy_members(num_qubits):
    # A phase of pi/2 on |11...1> has diagonal hierarchy level n+1.
    phases = np.ones(1 << num_qubits, dtype=complex)
    phases[-1] = 1j
    result = probe_gsc_conjugation_condition(
        np.diag(phases), hierarchy_level=num_qubits + 1, samples=8, seed=13
    )
    assert result.evidence is None
    assert result.words_checked == 8


@pytest.mark.parametrize("bad", [np.ones((2, 2)), np.eye(3), np.full((2, 2), np.nan)])
def test_invalid_matrices_are_rejected(bad):
    with pytest.raises(ValueError):
        check_gsc_sandwich(bad, hierarchy_level=4)


def test_invalid_options_are_rejected():
    for kwargs in (
        {"route": "wrong"},
        {"samples": 0},
        {"word_length": 0},
        {"numerical_margin": float("nan")},
    ):
        with pytest.raises(ValueError):
            probe_gsc_conjugation_condition(I2, hierarchy_level=4, **kwargs)
    for level in (0, 9):
        with pytest.raises(ValueError):
            probe_hierarchy_distance(I2, level=level)
    with pytest.raises(ValueError):
        DenseConjugationBackend(np.eye(16))
    with pytest.raises(ValueError):
        DenseConjugationBackend(I2).pauli((1, 0, 0, 0))
    with pytest.raises(ValueError):
        check_gsc_sandwich(I2, hierarchy_level=4, epsilon=-1)
