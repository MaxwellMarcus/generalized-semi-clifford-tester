"""Polynomial postprocessing for a promised Choi/Bell-difference SC test.

The planner assumes exact U in C_k and ideal independent Bell-difference
samples from its Choi state. Positive finite-data compatibility alone is
not a certificate. This tests SC, not the full GSC class; see
docs/semi-clifford/choi-tester.md.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from math import ceil, isfinite, log

from ..lagrangian import _int_to_label, _label_to_int, _row_reduce, _symplectic_pairing


def _positive(value: int, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{name} must be a positive integer")


def _insert(basis: tuple[int, ...], vector: int) -> tuple[int, ...]:
    for row in basis:
        if vector >> (row.bit_length() - 1) & 1:
            vector ^= row
    if not vector:
        return basis
    pivot = vector.bit_length() - 1
    rows = [row ^ vector if row >> pivot & 1 else row for row in basis]
    rows.append(vector)
    return tuple(sorted(rows, key=int.bit_length, reverse=True))


def _symplectic_complement(basis: tuple[int, ...], qubits: int) -> tuple[int, ...]:
    mask = (1 << qubits) - 1
    constraints = tuple(((row & mask) << qubits) | (row >> qubits) for row in basis)
    rows = _row_reduce(constraints, 2 * qubits)
    pivots = {row.bit_length() - 1 for row in rows}
    answer = []
    for free in range(2 * qubits):
        if free in pivots:
            continue
        vector = 1 << free
        for row in rows:
            if (row & vector).bit_count() % 2:
                vector ^= 1 << (row.bit_length() - 1)
        answer.append(vector)
    return tuple(answer)


@dataclass(frozen=True)
class ChoiSCAnalysis:
    """Finite-data compatibility; positive claims need the sampling theorem."""

    num_qubits: int
    samples_used: int
    sampled_span_dimension: int
    compatible_stabilizer_dimension: int
    projected_input_basis: tuple[tuple[int, ...], ...]
    maximum_input_isotropic_dimension: int

    @property
    def semi_clifford_compatible(self) -> bool:
        return self.maximum_input_isotropic_dimension == self.num_qubits


def analyze_choi_bell_differences(
    num_qubits: int,
    samples: Iterable[Sequence[int]],
) -> ChoiSCAnalysis:
    """Process XORs of independent Bell outcomes on pairs of Choi states.

    Labels have 4n bits in (x_output, x_reference | z_output, z_reference)
    order. Each block has n bits. They are not the 2n-bit labels obtained
    by Bell-measuring the Choi state of U P U† in the other tester.

    The compatible stabilizer space contains the true unsigned stabilizer
    space under ideal sampling. Therefore a negative compatibility result
    excludes SC; a positive result needs a sufficient sample budget and
    the hierarchy promise for its false-acceptance guarantee.
    """
    _positive(num_qubits, "num_qubits")
    n = num_qubits
    sampled_basis: tuple[int, ...] = ()
    used = 0
    for label in samples:
        if len(label) != 4 * n or any(bit not in (0, 1) for bit in label):
            raise ValueError("Bell-difference labels must contain exactly 4n binary entries")
        sampled_basis = _insert(sampled_basis, _label_to_int(tuple(label)))
        used += 1
    compatible = _symplectic_complement(sampled_basis, 2 * n)
    mask = (1 << n) - 1
    projected: tuple[int, ...] = ()
    for vector in compatible:
        # Transposing a reference Pauli changes only its sign, not its label.
        source = ((vector >> n) & mask) | (((vector >> (3 * n)) & mask) << n)
        projected = _insert(projected, source)
    gram = tuple(
        sum(_symplectic_pairing(left, right, n) << j for j, right in enumerate(projected))
        for left in projected
    )
    gram_rank = len(_row_reduce(gram, len(projected)))
    isotropic_dimension = len(projected) - gram_rank // 2
    return ChoiSCAnalysis(
        n,
        used,
        len(sampled_basis),
        len(compatible),
        tuple(_int_to_label(row, 2 * n) for row in projected),
        isotropic_dimension,
    )


@dataclass(frozen=True)
class ChoiSCSamplePlan:
    num_qubits: int
    hierarchy_level: int
    failure_probability: float
    bell_difference_samples: int

    @property
    def unitary_queries(self) -> int:
        """Four fresh Choi states per sample; inverse queries are unnecessary."""
        return 4 * self.bell_difference_samples

    @property
    def sc_distance_gap(self) -> float:
        return 2.0 ** (1.5 - self.hierarchy_level)


def plan_choi_sc_samples(
    num_qubits: int,
    hierarchy_level: int,
    *,
    failure_probability: float = 0.01,
) -> ChoiSCSamplePlan:
    """Sufficient budget for exact SC-vs-non-SC under U in C_k, k>=2.

    This is an ideal-model theorem, not a verifier of the promise or noise
    assumptions. GSC-but-not-SC inputs are classified as non-SC, not non-GSC.
    """
    _positive(num_qubits, "num_qubits")
    _positive(hierarchy_level, "hierarchy_level")
    if hierarchy_level < 2:
        raise ValueError("hierarchy_level must be at least 2")
    if not isfinite(failure_probability) or not 0 < failure_probability < 1:
        raise ValueError("failure_probability must lie strictly between zero and one")
    samples = (
        2 * 4 ** (hierarchy_level - 2) * ceil(4 * num_qubits * log(2) - log(failure_probability))
    )
    return ChoiSCSamplePlan(num_qubits, hierarchy_level, failure_probability, samples)
