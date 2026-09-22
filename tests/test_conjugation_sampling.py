"""Exact finite-group tests, plus floating-point reference-backend checks."""

from dataclasses import dataclass

import numpy as np
import pytest

from generalized_semi_clifford import (
    ConjugationWord,
    ContainmentStatus,
    DenseConjugationBackend,
    SamplingLimitExceeded,
    evaluate_conjugation_word,
    sample_conjugation_word,
    sample_subgroup_escape_word,
)
from generalized_semi_clifford import (
    test_clifford_conjugation_containment as check_containment,
)


@dataclass
class DihedralBackend:
    """Exact projective one-qubit group <Rz(2pi/m), X>.

    (a,b) represents r**a s**b, r=diag(1, exp(2pi i/m)), s=X.
    This is an exact arithmetic backend, not a floating comparison oracle.
    """

    modulus: int = 16
    num_qubits: int = 1
    exact: bool = True
    arithmetic: str = "exact dihedral group coordinates"
    seed: tuple[int, int] = (1, 0)

    def identity(self):
        return (0, 0)

    def pauli(self, label):
        x, z = label
        return (self.modulus // 2 * z, x)

    def multiply(self, left, right):
        a, b = left
        c, d = right
        return ((a + (-1) ** b * c) % self.modulus, b ^ d)

    def adjoint(self, value):
        a, b = value
        return ((-((-1) ** b) * a) % self.modulus, b)

    def is_clifford(self, value):
        return value[0] % (self.modulus // 4) == 0


def closure(backend, generators):
    elements = {backend.identity()}
    frontier = list(elements)
    while frontier:
        element = frontier.pop()
        for generator in generators:
            candidate = backend.multiply(element, generator)
            if candidate not in elements:
                elements.add(candidate)
                frontier.append(candidate)
    return elements


def exact_gamma(backend, depth):
    previous = {backend.seed}
    probes = [backend.pauli((1, 0)), backend.pauli((0, 1))]
    for _ in range(depth):
        generators = [
            backend.multiply(backend.multiply(g, p), backend.adjoint(g))
            for g in previous
            for p in probes
        ]
        previous = closure(backend, generators)
    return previous


@pytest.mark.parametrize("depth", [1, 2, 3, 4])
def test_short_samples_are_actual_group_elements_and_replay(depth):
    backend = DihedralBackend(32)
    group = exact_gamma(backend, depth)
    for seed in range(10):
        sample = sample_conjugation_word(backend, depth=depth, word_length=2, seed=seed)
        assert sample.value in group
        assert evaluate_conjugation_word(sample.word, backend) == sample.value
        assert sample.expanded_seed_uses == 2 * 4 ** (depth - 1)


def test_shallow_sampler_reproducible_and_seed_dependent():
    backend = DihedralBackend(32)
    first = sample_conjugation_word(backend, seed=7)
    repeat = sample_conjugation_word(backend, seed=7)
    assert first.value == repeat.value
    assert first.group_operations == repeat.group_operations
    assert len({sample_conjugation_word(backend, seed=s).value for s in range(20)}) > 1


def test_escape_sample_does_not_require_enumeration_and_replays():
    backend = DihedralBackend(32)
    sample = sample_subgroup_escape_word(backend, depth=2, subgroup_log2_bound=5, seed=7)
    assert sample.value in exact_gamma(backend, 2)
    assert evaluate_conjugation_word(sample.word, backend) == sample.value
    # Compact group-operation count is not the number of expanded U queries.
    assert sample.expanded_seed_uses > sample.group_operations


def test_escape_sampler_small_target_uses_paulis_at_any_depth():
    backend = DihedralBackend()
    sample = sample_subgroup_escape_word(backend, depth=7, subgroup_log2_bound=1, seed=7)
    assert sample.value in exact_gamma(backend, 2)
    assert sample.expanded_seed_uses == 0


def test_exact_containment_accepts_c4_phase_and_rejects_c5_phase_at_depth_two():
    contained = DihedralBackend(16)
    result = check_containment(
        contained,
        contained.is_clifford,
        depth=2,
        seed=7,
        failure_probability=0.05,
        oracle_is_exact=True,
    )
    assert result.status is ContainmentStatus.SUPPORTED
    assert result.rigorous_false_accept_bound == 0.05
    assert result.rounds_completed == result.rounds_required
    assert all(contained.is_clifford(g) for g in exact_gamma(contained, 2))

    outside = DihedralBackend(32)
    result = check_containment(outside, outside.is_clifford, depth=2, seed=7, oracle_is_exact=True)
    assert result.status is ContainmentStatus.VIOLATION
    assert result.witness is not None
    assert result.witness.value in exact_gamma(outside, 2)
    assert not outside.is_clifford(result.witness.value)
    assert evaluate_conjugation_word(result.witness.word, outside) == result.witness.value
    assert result.rigorous_false_accept_bound is None


def test_depth_one_is_deterministic_generator_check():
    backend = DihedralBackend(8)
    result = check_containment(backend, backend.is_clifford, depth=1, oracle_is_exact=True)
    assert result.status is ContainmentStatus.SUPPORTED
    assert result.rigorous_false_accept_bound == 0
    assert result.oracle_calls == 2


def test_recursive_depth_three_detects_violation_without_previous_group_access():
    backend = DihedralBackend(64)
    result = check_containment(
        backend,
        backend.is_clifford,
        depth=3,
        seed=0,
        oracle_is_exact=True,
        max_operations=200_000,
    )
    assert result.status is ContainmentStatus.VIOLATION
    assert result.witness is not None
    # Closure is only the independent test oracle, never called by the sampler.
    assert result.witness.value in exact_gamma(backend, 3)
    assert not backend.is_clifford(result.witness.value)
    assert (
        evaluate_conjugation_word(result.witness.word, backend, max_operations=200_000)
        == result.witness.value
    )
    assert result.witness.expanded_seed_uses.bit_length() > 100


def test_probability_claim_requires_exact_backend_and_oracle():
    backend = DihedralBackend(8)
    result = check_containment(backend, backend.is_clifford, depth=1)
    assert result.ideal_false_accept_bound == 0
    assert result.rigorous_false_accept_bound is None
    backend.exact = False
    result = check_containment(backend, backend.is_clifford, depth=1, oracle_is_exact=True)
    assert result.rigorous_false_accept_bound is None


def test_budget_and_unknown_oracle_never_accept():
    backend = DihedralBackend()
    limited = check_containment(backend, backend.is_clifford, depth=3, max_operations=30)
    assert limited.status is ContainmentStatus.INCONCLUSIVE
    assert limited.group_operations == 30
    assert limited.ideal_false_accept_bound is None
    unknown = check_containment(backend, lambda _: None, depth=1)
    assert unknown.status is ContainmentStatus.INCONCLUSIVE
    assert unknown.oracle_calls == 1
    assert unknown.rounds_completed == 0
    with pytest.raises(TypeError, match="bool or None"):
        check_containment(backend, lambda _: "yes", depth=1)


def test_dense_replay_agrees_with_exact_backend():
    exact = DihedralBackend(32)
    unitary = np.diag([1, np.exp(2j * np.pi / 32)])
    dense = DenseConjugationBackend(unitary)
    sample = sample_conjugation_word(exact, depth=3, seed=19, word_length=2)
    matrix = evaluate_conjugation_word(sample.word, dense)
    a, b = sample.value
    expected = np.diag([1, np.exp(2j * np.pi * a / 32)])
    if b:
        expected = expected @ np.array([[0, 1], [1, 0]])
    overlap = np.vdot(expected, matrix)
    assert np.allclose(matrix, overlap / abs(overlap) * expected)


def test_sampling_and_replay_validation_and_caps():
    backend = DihedralBackend()
    for kwargs in ({"depth": 0}, {"depth": 9}, {"word_length": 0}, {"word_length": True}):
        with pytest.raises(ValueError):
            sample_conjugation_word(backend, **kwargs)
    with pytest.raises(ValueError):
        sample_subgroup_escape_word(backend, depth=2, subgroup_log2_bound=-1)
    with pytest.raises(SamplingLimitExceeded):
        sample_conjugation_word(backend, max_operations=1)
    sample = sample_conjugation_word(backend, seed=3)
    with pytest.raises(SamplingLimitExceeded):
        evaluate_conjugation_word(sample.word, backend, max_operations=1)
    with pytest.raises(ValueError, match="malformed"):
        evaluate_conjugation_word(ConjugationWord("invalid"), backend)
    for probability in (0, 1, float("nan"), float("inf")):
        with pytest.raises(ValueError):
            check_containment(
                backend, backend.is_clifford, depth=2, failure_probability=probability
            )
