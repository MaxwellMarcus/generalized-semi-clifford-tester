from itertools import product
from math import log

import numpy as np
import pytest

from generalized_semi_clifford.lagrangian import (
    _int_to_label,
    _label_to_int,
    _row_reduce,
    _symplectic_pairing,
)
from generalized_semi_clifford.pauli import all_pauli_labels, pauli_matrix
from generalized_semi_clifford.sc_testing import (
    analyze_choi_bell_differences,
    plan_choi_sc_samples,
)


def _bell_difference_distribution(gate):
    # Compute the actual identical-state Bell measurement probabilities,
    # then convolve them, rather than assuming the Fourier identity under test.
    n = len(gate).bit_length() - 1
    state = gate.reshape(-1) / np.sqrt(len(gate))
    paulis = [pauli_matrix(label) for label in all_pauli_labels(2 * n)]
    bell = np.array([abs(state @ p @ state) ** 2 / len(state) for p in paulis])
    indices = np.arange(len(bell))
    difference = np.array([np.sum(bell * bell[indices ^ i]) for i in indices])
    expectations = np.array([np.vdot(state, p @ state).real for p in paulis])
    return difference, expectations


@pytest.mark.parametrize("gate", [np.diag([1, np.exp(1j * np.pi / 4)]), np.diag([1, 1, 1, 1j])])
def test_actual_bell_difference_fourier_identity_and_recovered_stabilizers(gate):
    n = len(gate).bit_length() - 1
    distribution, expectations = _bell_difference_distribution(gate)
    assert distribution.sum() == pytest.approx(1)
    for a, expectation in enumerate(expectations):
        odd = [i for i in range(len(distribution)) if _symplectic_pairing(a, i, 2 * n)]
        assert distribution[odd].sum() == pytest.approx((1 - expectation**4) / 2)
    support = [_int_to_label(i, 4 * n) for i, weight in enumerate(distribution) if weight > 1e-12]
    result = analyze_choi_bell_differences(n, support)
    exact_stabilizers = [i for i, value in enumerate(expectations) if abs(value) > 1 - 1e-12]
    assert result.compatible_stabilizer_dimension == len(
        _row_reduce(tuple(exact_stabilizers), 4 * n)
    )
    assert result.semi_clifford_compatible
    assert result.maximum_input_isotropic_dimension == n


def test_finite_samples_recover_t_gate_with_planned_budget():
    gate = np.diag([1, np.exp(1j * np.pi / 4)])
    distribution, _ = _bell_difference_distribution(gate)
    distribution /= distribution.sum()
    plan = plan_choi_sc_samples(1, 3)
    rng = np.random.default_rng(52)
    draws = rng.choice(len(distribution), size=plan.bell_difference_samples, p=distribution)
    result = analyze_choi_bell_differences(1, (_int_to_label(int(i), 4) for i in draws))
    assert result.compatible_stabilizer_dimension == 1
    assert result.projected_input_basis == ((0, 1),)
    assert result.semi_clifford_compatible


def test_non_sc_without_hierarchy_promise_full_support_only():
    # This generic gate is not claimed to be in a hierarchy level. Exhaustive
    # support is used to test postprocessing, not a promised query bound.
    x = pauli_matrix((1, 0))
    y = pauli_matrix((1, 1))
    z = pauli_matrix((0, 1))
    gate = np.sqrt(3) / 2 * np.eye(2) - 0.5j * (x + y + z) / np.sqrt(3)
    distribution, _ = _bell_difference_distribution(gate)
    result = analyze_choi_bell_differences(
        1, (_int_to_label(i, 4) for i, p in enumerate(distribution) if p > 1e-12)
    )
    assert result.compatible_stabilizer_dimension == 0
    assert not result.semi_clifford_compatible


@pytest.mark.parametrize("n", [1, 3, 30])
def test_basis_only_processing_sc_and_non_sc_extremes(n):
    empty = analyze_choi_bell_differences(n, [])
    assert empty.semi_clifford_compatible
    full = [_int_to_label(1 << i, 4 * n) for i in range(4 * n)]
    result = analyze_choi_bell_differences(n, full)
    assert result.compatible_stabilizer_dimension == 0
    assert not result.semi_clifford_compatible


def test_projection_uses_reference_blocks_not_output_blocks():
    # Compatibility space is spanned by reference Z. Supply its entire
    # orthogonal hyperplane as samples to make the intended space explicit.
    reference_z = _label_to_int((0, 0, 0, 1))
    samples = [
        _int_to_label(i, 4) for i in range(16) if _symplectic_pairing(i, reference_z, 2) == 0
    ]
    result = analyze_choi_bell_differences(1, samples)
    assert result.projected_input_basis == ((0, 1),)


def test_isotropic_dimension_not_just_projected_rank():
    # A two-dimensional nondegenerate input subspace on two qubits has
    # isotropic dimension one, insufficient for a two-qubit SC witness.
    n = 2
    desired = (1 << n, 1 << (3 * n))  # X and Z on the first reference qubit.
    samples = [
        _int_to_label(i, 4 * n)
        for i in range(1 << (4 * n))
        if all(_symplectic_pairing(i, row, 2 * n) == 0 for row in desired)
    ]
    result = analyze_choi_bell_differences(n, samples)
    assert len(result.projected_input_basis) == 2
    assert result.maximum_input_isotropic_dimension == 1
    assert not result.semi_clifford_compatible


@pytest.mark.parametrize("n,k", [(1, 2), (3, 4), (100, 5)])
def test_budget_controls_union_bound(n, k):
    plan = plan_choi_sc_samples(n, k, failure_probability=0.01)
    beta = 1 / (2 * 4 ** (k - 2))
    assert 4 * n * log(2) - beta * plan.bell_difference_samples <= log(0.01)
    assert plan.unitary_queries == 4 * plan.bell_difference_samples
    assert plan.sc_distance_gap == pytest.approx(np.sqrt(2) / 2 ** (k - 1))


def test_bad_inputs():
    for n, k, delta in [(0, 3, 0.01), (1, 1, 0.01), (1, 3, 0), (1, 3, float("nan"))]:
        with pytest.raises(ValueError):
            plan_choi_sc_samples(n, k, failure_probability=delta)
    with pytest.raises(ValueError):
        analyze_choi_bell_differences(1, [(0, 1)])
    with pytest.raises(ValueError):
        analyze_choi_bell_differences(1, [(0, 1, 2, 0)])


def test_complement_and_projection_against_brute_force():
    rng = np.random.default_rng(19)
    for _ in range(10):
        vectors = tuple(int(i) for i in rng.integers(256, size=5))
        result = analyze_choi_bell_differences(2, (_int_to_label(v, 8) for v in vectors))
        complement = [
            i for i in range(256) if all(_symplectic_pairing(i, v, 4) == 0 for v in vectors)
        ]
        projected = tuple(((i >> 2) & 3) | (((i >> 6) & 3) << 2) for i in complement)
        expected = _row_reduce(projected, 4)
        actual = _row_reduce(tuple(map(_label_to_int, result.projected_input_basis)), 4)
        assert actual == expected
        assert all(_symplectic_pairing(v, a, 4) == 0 for v, a in product(vectors, complement))
