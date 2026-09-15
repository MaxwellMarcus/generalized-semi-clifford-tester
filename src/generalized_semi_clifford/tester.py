"""Definition-first exhaustive generalized semi-Clifford tester."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .lagrangian import MAX_EXHAUSTIVE_QUBITS, Lagrangian, enumerate_lagrangians, is_lagrangian
from .pauli import all_pauli_labels, infer_num_qubits, pauli_matrix


class GSCStatus(str, Enum):
    """Possible outcomes from the bounded naïve search."""

    GSC = "gsc"
    NOT_GSC = "not_gsc"
    UNKNOWN = "unknown"


@dataclass(frozen=True)
class GSCWitness:
    """Input/output Pauli MASAs witnessing generalized semi-Cliffordness."""

    input_lagrangian: Lagrangian
    output_lagrangian: Lagrangian
    max_leakage: float


@dataclass(frozen=True)
class NaiveGSCResult:
    """Result and diagnostics from :func:`check_gsc_naive`."""

    status: GSCStatus
    num_qubits: int
    arithmetic: str
    tolerance: float
    lagrangians_checked: int
    candidate_pairs_checked: int
    best_leakage: float | None
    witness: GSCWitness | None
    message: str

    @property
    def is_gsc(self) -> bool | None:
        """Return a Boolean result, or ``None`` when the search was not run."""

        if self.status is GSCStatus.UNKNOWN:
            return None
        return self.status is GSCStatus.GSC


def _validated_unitary(unitary: ArrayLike, unitary_tolerance: float) -> NDArray[np.complex128]:
    matrix = np.asarray(unitary, dtype=np.complex128)
    num_qubits = infer_num_qubits(matrix)
    identity = np.eye(1 << num_qubits, dtype=np.complex128)
    residual = float(np.linalg.norm(matrix.conj().T @ matrix - identity, ord=2))
    if residual > unitary_tolerance:
        raise ValueError(
            f"matrix is not unitary within tolerance {unitary_tolerance:g}; residual={residual:g}"
        )
    return matrix


def _lagrangian_leakage(
    unitary: NDArray[np.complex128],
    input_lagrangian: Lagrangian,
    output_lagrangians: tuple[Lagrangian, ...],
) -> NDArray[np.float64]:
    num_qubits = input_lagrangian.num_qubits
    dimension = 1 << num_qubits
    labels = all_pauli_labels(num_qubits)
    label_indices = {label: index for index, label in enumerate(labels)}
    paulis = np.stack([pauli_matrix(label) for label in labels])
    maximum_coefficients = np.zeros(len(labels), dtype=np.float64)
    unitary_dagger = unitary.conj().T

    for label in input_lagrangian.elements:
        image = unitary @ pauli_matrix(label) @ unitary_dagger
        coefficients = np.abs(np.einsum("aij,ij->a", paulis.conj(), image)) / dimension
        maximum_coefficients = np.maximum(maximum_coefficients, coefficients)

    residuals = np.empty(len(output_lagrangians), dtype=np.float64)
    all_indices = np.arange(len(labels))
    for index, output_lagrangian in enumerate(output_lagrangians):
        inside = np.zeros(len(labels), dtype=bool)
        inside[[label_indices[label] for label in output_lagrangian.elements]] = True
        residuals[index] = float(np.max(maximum_coefficients[all_indices[~inside]]))
    return residuals


def check_gsc_naive(
    unitary: ArrayLike,
    *,
    tolerance: float = 1e-9,
    unitary_tolerance: float = 1e-9,
    max_qubits: int = 3,
) -> NaiveGSCResult:
    """Exhaustively test the defining Pauli-MASA condition for small qubit systems.

    For each pair of binary Lagrangians ``(L, S)``, the algorithm checks whether
    every Pauli coefficient of ``U A_L U†`` outside ``A_S`` has magnitude at
    most ``tolerance``. A completed search is exhaustive over the chosen floating-
    point criterion, but it is not an exact symbolic proof.
    """

    if tolerance <= 0 or unitary_tolerance <= 0:
        raise ValueError("tolerances must be positive")
    if max_qubits < 1:
        raise ValueError("max_qubits must be positive")
    matrix = _validated_unitary(unitary, unitary_tolerance)
    num_qubits = infer_num_qubits(matrix)
    search_limit = min(max_qubits, MAX_EXHAUSTIVE_QUBITS)
    if num_qubits > search_limit:
        return NaiveGSCResult(
            status=GSCStatus.UNKNOWN,
            num_qubits=num_qubits,
            arithmetic="complex128 floating point",
            tolerance=tolerance,
            lagrangians_checked=0,
            candidate_pairs_checked=0,
            best_leakage=None,
            witness=None,
            message=f"n={num_qubits} exceeds the exhaustive-search limit {search_limit}",
        )

    lagrangians = enumerate_lagrangians(num_qubits)
    best_leakage = float("inf")
    best_pair: tuple[Lagrangian, Lagrangian] | None = None
    for input_index, input_lagrangian in enumerate(lagrangians):
        residuals = _lagrangian_leakage(matrix, input_lagrangian, lagrangians)
        output_index = int(np.argmin(residuals))
        leakage = float(residuals[output_index])
        if leakage < best_leakage:
            best_leakage = leakage
            best_pair = input_lagrangian, lagrangians[output_index]
        if leakage <= tolerance:
            witness = GSCWitness(input_lagrangian, lagrangians[output_index], leakage)
            return NaiveGSCResult(
                status=GSCStatus.GSC,
                num_qubits=num_qubits,
                arithmetic="complex128 floating point",
                tolerance=tolerance,
                lagrangians_checked=input_index + 1,
                candidate_pairs_checked=(input_index + 1) * len(lagrangians),
                best_leakage=leakage,
                witness=witness,
                message="found Pauli MASAs satisfying U A_L U† = A_S within tolerance",
            )

    assert best_pair is not None
    return NaiveGSCResult(
        status=GSCStatus.NOT_GSC,
        num_qubits=num_qubits,
        arithmetic="complex128 floating point",
        tolerance=tolerance,
        lagrangians_checked=len(lagrangians),
        candidate_pairs_checked=len(lagrangians) ** 2,
        best_leakage=best_leakage,
        witness=None,
        message="no Pauli-MASA pair passed the exhaustive numerical search",
    )


def verify_gsc_witness(
    unitary: ArrayLike,
    witness: GSCWitness,
    *,
    tolerance: float = 1e-9,
    unitary_tolerance: float = 1e-9,
) -> bool:
    """Independently recompute and check a witness from the naïve tester."""

    matrix = _validated_unitary(unitary, unitary_tolerance)
    num_qubits = infer_num_qubits(matrix)
    if (
        num_qubits != witness.input_lagrangian.num_qubits
        or num_qubits != witness.output_lagrangian.num_qubits
        or not is_lagrangian(witness.input_lagrangian.elements, num_qubits)
        or not is_lagrangian(witness.output_lagrangian.elements, num_qubits)
    ):
        return False
    residual = _lagrangian_leakage(
        matrix,
        witness.input_lagrangian,
        (witness.output_lagrangian,),
    )[0]
    return bool(residual <= tolerance)
