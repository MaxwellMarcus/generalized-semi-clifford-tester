"""Small dense experiments with explicit, numerical distance evidence.

No result here is interval-certified. The sampler's exact-arithmetic theorem
must not be transferred to these floating-point calculations unconditionally.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from math import isfinite, sqrt

import numpy as np
from numpy.typing import ArrayLike

from .conjugation_sampling import (
    ConjugationSample,
    SamplingLimitExceeded,
    _integer,
    sample_conjugation_word,
)
from .lagrangian import PauliLabel, enumerate_lagrangians
from .pauli import ComplexMatrix, all_pauli_labels, infer_num_qubits, pauli_matrix
from .tester import GSCStatus, GSCWitness, check_gsc_naive


def _positive(value: float, name: str) -> None:
    if not isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")


class DenseConjugationBackend:
    """Numerical reference backend, deliberately capped at three qubits.

    No rounding to a Clifford or silent projection to a unitary is performed.
    Long shared words can amplify floating-point error as well as query cost.
    """

    exact = False
    arithmetic = "complex128 floating point; not interval-certified"

    def __init__(self, unitary: ArrayLike, *, tolerance: float = 1e-9) -> None:
        _positive(tolerance, "tolerance")
        matrix = np.asarray(unitary, dtype=np.complex128)
        self.num_qubits = infer_num_qubits(matrix)
        if self.num_qubits > 3:
            raise ValueError("the dense conjugation backend is limited to three qubits")
        if not np.all(np.isfinite(matrix)):
            raise ValueError("unitary entries must be finite")
        self.tolerance = tolerance
        if not self.is_numerically_unitary(matrix):
            raise ValueError("matrix is not unitary within tolerance")
        self.seed = np.array(matrix, copy=True)
        self.seed.flags.writeable = False

    def identity(self) -> ComplexMatrix:
        return np.eye(1 << self.num_qubits, dtype=np.complex128)

    def pauli(self, label: PauliLabel) -> ComplexMatrix:
        if len(label) != 2 * self.num_qubits:
            raise ValueError("Pauli label has the wrong number of qubits")
        return pauli_matrix(label)

    def multiply(self, left: ComplexMatrix, right: ComplexMatrix) -> ComplexMatrix:
        return left @ right

    def adjoint(self, value: ComplexMatrix) -> ComplexMatrix:
        return value.conj().T

    def is_numerically_unitary(self, matrix: ComplexMatrix) -> bool:
        if not np.all(np.isfinite(matrix)):
            return False
        residual = np.linalg.norm(matrix.conj().T @ matrix - self.identity(), ord=2)
        return bool(residual <= self.tolerance)

    def clifford_membership(self, matrix: ComplexMatrix) -> bool | None:
        """Check all 2n generators; return None on detectable numerical drift."""
        if not self.is_numerically_unitary(matrix):
            return None
        width = 2 * self.num_qubits
        for index in range(width):
            label = tuple(int(i == index) for i in range(width))
            image = matrix @ self.pauli(label) @ matrix.conj().T
            if _pauli_distance(image, self.num_qubits) > self.tolerance:
                return False
        return True


def _pauli_distance(matrix: ComplexMatrix, num_qubits: int) -> float:
    """min_phase,P ||matrix - phase P||_F / sqrt(d), evaluated stably."""
    best_overlap = -1.0
    best_pauli = None
    phase = 1.0 + 0.0j
    for label in all_pauli_labels(num_qubits):
        pauli = pauli_matrix(label)
        overlap = np.vdot(pauli, matrix)
        if abs(overlap) > best_overlap:
            best_overlap = float(abs(overlap))
            best_pauli = pauli
            phase = overlap / abs(overlap) if abs(overlap) else 1.0 + 0.0j
    assert best_pauli is not None
    # Unlike sqrt(2 - 2*max_overlap/d), this is stable near zero.
    return float(np.linalg.norm(matrix - phase * best_pauli, ord="fro") / sqrt(len(matrix)))


@dataclass(frozen=True)
class HierarchyDistanceEvidence:
    """One Pauli-conjugation path giving a lower bound on distance to C_level."""

    level: int
    pauli_path: tuple[PauliLabel, ...]
    leaf_pauli_distance: float
    numerical_margin: float
    distance_lower_bound: float


def probe_hierarchy_distance(
    unitary: ArrayLike,
    *,
    level: int,
    trials: int = 16,
    seed: int | None = None,
    numerical_margin: float = 1e-9,
) -> HierarchyDistanceEvidence:
    """Find a numerical lower bound, not an approximate-membership certificate.

    Descend level-1 times by Pauli conjugation. For a leaf at distance rho
    from the Paulis, distance from the original operator to C_level is at
    least rho / 2**(level-1). Zero evidence is inconclusive. For level two,
    all standard Pauli generators are checked instead of random paths.
    The margin is user-selected, not an automatic roundoff-error bound.
    """
    _integer(level, "level")
    if level > 8:
        raise ValueError("level exceeds the prototype's limit of 8")
    _integer(trials, "trials")
    _positive(numerical_margin, "numerical_margin")
    backend = DenseConjugationBackend(unitary, tolerance=numerical_margin)
    rng = np.random.default_rng(seed)
    width = 2 * backend.num_qubits
    best = HierarchyDistanceEvidence(level, (), 0.0, numerical_margin, 0.0)
    count = 1 if level == 1 else width if level == 2 else trials
    for trial in range(count):
        matrix = backend.seed
        path = []
        for _ in range(level - 1):
            label = (
                tuple(int(i == trial) for i in range(width))
                if level == 2
                else tuple(int(bit) for bit in rng.integers(2, size=width))
            )
            path.append(label)
            matrix = matrix @ backend.pauli(label) @ matrix.conj().T
        if not backend.is_numerically_unitary(matrix):
            raise ArithmeticError("numerical drift in a hierarchy probe; no distance claim")
        distance = _pauli_distance(matrix, backend.num_qubits)
        lower_bound = max(0.0, distance - numerical_margin) / 2 ** (level - 1)
        if trial == 0 or lower_bound > best.distance_lower_bound:
            best = HierarchyDistanceEvidence(
                level, tuple(path), distance, numerical_margin, lower_bound
            )
    return best


@dataclass(frozen=True)
class ConjugationDistanceEvidence:
    sample: ConjugationSample[ComplexMatrix]
    hierarchy: HierarchyDistanceEvidence
    distance_from_sc_in_ck_lower_bound: float


@dataclass(frozen=True)
class ConjugationProbeResult:
    hierarchy_level: int
    route: str
    words_checked: int
    group_operations: int
    evidence: ConjugationDistanceEvidence | None
    budget_exhausted: bool
    numerical_failure: bool
    message: str


def probe_gsc_conjugation_condition(
    unitary: ArrayLike,
    *,
    hierarchy_level: int,
    route: str = "shallow",
    samples: int = 32,
    word_length: int = 2,
    hierarchy_trials: int = 16,
    seed: int | None = None,
    numerical_margin: float = 1e-9,
    max_operations: int = 100_000,
) -> ConjugationProbeResult:
    """Probe Gamma_2 <= C_(k-2), or Gamma_(k-2) <= C_2 with short words.

    No detection probability is asserted for these short-word distributions.
    Any positive bound concerns SC intersect C_k, not all SC and not non-GSC.
    This interpretation uses SC intersect C_k => the containment conditions.
    """
    _integer(hierarchy_level, "hierarchy_level", 3)
    if hierarchy_level > 10:
        raise ValueError("hierarchy_level exceeds the prototype's limit of 10")
    _integer(samples, "samples")
    _integer(word_length, "word_length")
    _integer(hierarchy_trials, "hierarchy_trials")
    _integer(max_operations, "max_operations")
    if route not in ("shallow", "deep"):
        raise ValueError("route must be 'shallow' or 'deep'")
    backend = DenseConjugationBackend(unitary, tolerance=numerical_margin)
    depth, target = (2, hierarchy_level - 2) if route == "shallow" else (hierarchy_level - 2, 2)
    rng = np.random.default_rng(seed)
    best = None
    operations = 0
    checked = 0
    stopped = False
    numerical_failure = False
    message = "short-word exploration only; passing samples do not establish containment"
    for _ in range(samples):
        if operations >= max_operations:
            stopped = True
            message = "group-operation budget exhausted; no containment conclusion"
            break
        try:
            sample = sample_conjugation_word(
                backend,
                depth=depth,
                word_length=word_length,
                seed=int(rng.integers(2**63)),
                max_operations=max_operations - operations,
            )
        except SamplingLimitExceeded:
            operations = max_operations
            stopped = True
            message = "group-operation budget exhausted; no containment conclusion"
            break
        operations += sample.group_operations
        try:
            evidence = probe_hierarchy_distance(
                sample.value,
                level=target,
                trials=hierarchy_trials,
                seed=int(rng.integers(2**63)),
                numerical_margin=numerical_margin,
            )
        except (ArithmeticError, ValueError):
            message = "numerical drift prevented completion; no containment conclusion"
            numerical_failure = True
            break
        checked += 1
        uses = sample.expanded_seed_uses
        bound = evidence.distance_lower_bound / uses if uses else 0.0
        if bound > 0 and (best is None or bound > best.distance_from_sc_in_ck_lower_bound):
            best = ConjugationDistanceEvidence(sample, evidence, bound)
    return ConjugationProbeResult(
        hierarchy_level, route, checked, operations, best, stopped, numerical_failure, message
    )


def semi_clifford_distance_lower_bound(
    unitary: ArrayLike, *, numerical_margin: float = 1e-9
) -> float:
    """Small-system numerical bound on distance to ALL semi-Clifford gates.

    Return (min_L max_{P in L} d(U P U-dagger, Paulis) - margin)_+ / 2.
    This exhaustive fallback uses no hierarchy promise. The metric is
    phase-optimized Frobenius norm divided by sqrt(dimension).
    """
    backend = DenseConjugationBackend(unitary, tolerance=numerical_margin)
    residuals = {}
    for label in all_pauli_labels(backend.num_qubits):
        image = backend.seed @ backend.pauli(label) @ backend.seed.conj().T
        residuals[label] = _pauli_distance(image, backend.num_qubits)
    residual = min(
        max(residuals[label] for label in lagrangian.elements)
        for lagrangian in enumerate_lagrangians(backend.num_qubits)
    )
    return max(0.0, residual - numerical_margin) / 2


class SandwichStatus(str, Enum):
    GSC_WITNESS = "numerical_gsc_witness"
    FAR_FROM_SC = "numerically_far_from_sc"
    FAR_FROM_SC_IN_CK = "numerically_far_from_sc_in_ck"
    INCONCLUSIVE = "inconclusive"


@dataclass(frozen=True)
class GSCSandwichResult:
    status: SandwichStatus
    arithmetic: str
    epsilon: float
    gsc_witness: GSCWitness | None
    distance_from_sc_lower_bound: float | None
    conjugation: ConjugationProbeResult | None
    message: str


def check_gsc_sandwich(
    unitary: ArrayLike,
    *,
    hierarchy_level: int,
    epsilon: float = 0.01,
    route: str = "shallow",
    samples: int = 32,
    word_length: int = 2,
    seed: int | None = None,
    numerical_margin: float = 1e-9,
) -> GSCSandwichResult:
    """Bounded numerical first step toward the GSC/SC-distance tester.

    A positive GSC outcome always has a direct Pauli-MASA witness; it is never
    inferred from passing random words. A separate exhaustive small-system
    bound can justify the unrestricted FAR_FROM_SC outcome. Above three qubits
    the dense prototype returns INCONCLUSIVE rather than running the search.
    """
    _integer(hierarchy_level, "hierarchy_level", 3)
    if hierarchy_level > 10:
        raise ValueError("hierarchy_level exceeds the prototype's limit of 10")
    _positive(epsilon, "epsilon")
    _positive(numerical_margin, "numerical_margin")
    _integer(samples, "samples")
    _integer(word_length, "word_length")
    if route not in ("shallow", "deep"):
        raise ValueError("route must be 'shallow' or 'deep'")
    matrix = np.asarray(unitary, dtype=np.complex128)
    n = infer_num_qubits(matrix)
    if not np.all(np.isfinite(matrix)):
        raise ValueError("unitary entries must be finite")
    arithmetic = DenseConjugationBackend.arithmetic
    if n > 3:
        return GSCSandwichResult(
            SandwichStatus.INCONCLUSIVE,
            arithmetic,
            epsilon,
            None,
            None,
            None,
            "input exceeds the three-qubit dense-prototype limit",
        )
    DenseConjugationBackend(matrix, tolerance=numerical_margin)
    conjugation = probe_gsc_conjugation_condition(
        matrix,
        hierarchy_level=hierarchy_level,
        route=route,
        samples=samples,
        word_length=word_length,
        seed=seed,
        numerical_margin=numerical_margin,
    )
    direct = check_gsc_naive(matrix, tolerance=numerical_margin, unitary_tolerance=numerical_margin)
    if direct.status is GSCStatus.GSC:
        return GSCSandwichResult(
            SandwichStatus.GSC_WITNESS,
            arithmetic,
            epsilon,
            direct.witness,
            None,
            conjugation,
            "direct numerical Pauli-MASA witness; sampling is not the acceptance certificate",
        )
    sc_bound = semi_clifford_distance_lower_bound(matrix, numerical_margin=numerical_margin)
    if sc_bound > epsilon:
        status = SandwichStatus.FAR_FROM_SC
        message = "exhaustive numerical distance bound from all SC exceeds epsilon"
    elif conjugation.evidence is not None and (
        conjugation.evidence.distance_from_sc_in_ck_lower_bound > epsilon
    ):
        status = SandwichStatus.FAR_FROM_SC_IN_CK
        message = "only the distance bound from SC intersect C_k exceeds epsilon"
    else:
        status = SandwichStatus.INCONCLUSIVE
        message = "no GSC witness or distance bound exceeding epsilon was established"
    return GSCSandwichResult(status, arithmetic, epsilon, None, sc_bound, conjugation, message)
