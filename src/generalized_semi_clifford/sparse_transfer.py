"""Sparse, exhaustively computed Pauli-transfer representations."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .lagrangian import PauliLabel
from .pauli import all_pauli_labels, infer_num_qubits, pauli_matrix


@dataclass(frozen=True)
class SparsePauliCoefficient:
    """One retained Pauli-transfer coefficient."""

    output_label: PauliLabel
    value: float


@dataclass(frozen=True)
class SparsePauliTransfer:
    """Thresholded Pauli-transfer matrix with numerical provenance.

    Columns are indexed by input Paulis and rows by output Paulis.  Every
    coefficient is computed by a dense trace; sparsity is introduced only
    after extraction.  ``discarded_l2_mass`` records the squared coefficient
    mass omitted from each column, so a thresholded result is not confused
    with an exact zero pattern.
    """

    num_qubits: int
    coefficient_tolerance: float
    unitary_tolerance: float
    unitary_residual: float
    arithmetic: str
    labels: tuple[PauliLabel, ...]
    columns: tuple[tuple[SparsePauliCoefficient, ...], ...]
    discarded_l2_mass: tuple[float, ...]

    def __post_init__(self) -> None:
        expected = 1 << (2 * self.num_qubits)
        if len(self.labels) != expected:
            raise ValueError("labels must contain every phase-free Pauli")
        if len(self.columns) != expected or len(self.discarded_l2_mass) != expected:
            raise ValueError("column metadata must match the Pauli label count")
        label_set = set(self.labels)
        if len(label_set) != expected:
            raise ValueError("Pauli labels must be unique")
        if any(entry.output_label not in label_set for column in self.columns for entry in column):
            raise ValueError("retained coefficients must use known output labels")
        if any(mass < 0 for mass in self.discarded_l2_mass):
            raise ValueError("discarded mass must be nonnegative")

    @property
    def retained_entries(self) -> int:
        """Number of coefficients whose magnitude exceeded the threshold."""

        return sum(len(column) for column in self.columns)

    @property
    def total_entries(self) -> int:
        """Number of coefficients computed before thresholding."""

        return len(self.labels) ** 2

    @property
    def max_discarded_l2_mass(self) -> float:
        """Largest omitted squared coefficient mass in any input column."""

        return max(self.discarded_l2_mass, default=0.0)

    def coefficient(self, input_label: PauliLabel, output_label: PauliLabel) -> float:
        """Return a retained coefficient, or zero when it was thresholded out."""

        try:
            column_index = self.labels.index(input_label)
        except ValueError as error:
            raise KeyError(f"unknown input Pauli label {input_label!r}") from error
        if output_label not in self.labels:
            raise KeyError(f"unknown output Pauli label {output_label!r}")
        for entry in self.columns[column_index]:
            if entry.output_label == output_label:
                return entry.value
        return 0.0

    def to_dense(self) -> NDArray[np.float64]:
        """Reconstruct the retained row-by-column transfer matrix."""

        label_indices = {label: index for index, label in enumerate(self.labels)}
        dense = np.zeros((len(self.labels), len(self.labels)), dtype=np.float64)
        for column_index, column in enumerate(self.columns):
            for entry in column:
                dense[label_indices[entry.output_label], column_index] = entry.value
        dense.setflags(write=False)
        return dense


def dense_pauli_transfer(
    unitary: ArrayLike,
    *,
    unitary_tolerance: float = 1e-9,
    max_qubits: int = 3,
) -> NDArray[np.float64]:
    """Compute every Pauli-transfer coefficient by dense trace evaluation.

    This is an exhaustive floating-point reference calculation, not symbolic
    arithmetic.  The small-qubit cap prevents accidental exponential work.
    """

    matrix, num_qubits, _ = _validated_unitary(unitary, unitary_tolerance, max_qubits)
    labels = all_pauli_labels(num_qubits)
    dimension = 1 << num_qubits
    paulis = np.stack([pauli_matrix(label) for label in labels])
    dagger = matrix.conj().T
    transfer = np.empty((len(labels), len(labels)), dtype=np.float64)
    for column, label in enumerate(labels):
        image = matrix @ pauli_matrix(label) @ dagger
        coefficients = np.einsum("aij,ji->a", paulis, image) / dimension
        transfer[:, column] = np.real_if_close(coefficients, tol=1000).real
    transfer.setflags(write=False)
    return transfer


def extract_sparse_pauli_transfer(
    unitary: ArrayLike,
    *,
    coefficient_tolerance: float = 1e-12,
    unitary_tolerance: float = 1e-9,
    max_qubits: int = 3,
) -> SparsePauliTransfer:
    """Exhaustively compute and threshold a small-qubit Pauli-transfer matrix."""

    if coefficient_tolerance <= 0 or not np.isfinite(coefficient_tolerance):
        raise ValueError("coefficient_tolerance must be finite and positive")
    matrix, num_qubits, residual = _validated_unitary(
        unitary,
        unitary_tolerance,
        max_qubits,
    )
    labels = all_pauli_labels(num_qubits)
    transfer = dense_pauli_transfer(
        matrix,
        unitary_tolerance=unitary_tolerance,
        max_qubits=max_qubits,
    )
    columns: list[tuple[SparsePauliCoefficient, ...]] = []
    discarded: list[float] = []
    for values in transfer.T:
        retained = np.abs(values) > coefficient_tolerance
        columns.append(
            tuple(
                SparsePauliCoefficient(output_label=labels[index], value=float(values[index]))
                for index in np.flatnonzero(retained)
            )
        )
        discarded.append(float(np.sum(np.square(values[~retained]))))
    return SparsePauliTransfer(
        num_qubits=num_qubits,
        coefficient_tolerance=coefficient_tolerance,
        unitary_tolerance=unitary_tolerance,
        unitary_residual=residual,
        arithmetic="complex128 exhaustive traces; float64 real coefficients",
        labels=labels,
        columns=tuple(columns),
        discarded_l2_mass=tuple(discarded),
    )


def _validated_unitary(
    unitary: ArrayLike,
    unitary_tolerance: float,
    max_qubits: int,
) -> tuple[NDArray[np.complex128], int, float]:
    if unitary_tolerance <= 0 or not np.isfinite(unitary_tolerance):
        raise ValueError("unitary_tolerance must be finite and positive")
    if not isinstance(max_qubits, int) or isinstance(max_qubits, bool):
        raise TypeError("max_qubits must be an integer")
    if max_qubits < 1:
        raise ValueError("max_qubits must be positive")
    matrix = np.asarray(unitary, dtype=np.complex128)
    num_qubits = infer_num_qubits(matrix)
    if num_qubits > max_qubits:
        raise ValueError(
            f"n={num_qubits} exceeds the dense Pauli-transfer limit {max_qubits}"
        )
    identity = np.eye(1 << num_qubits, dtype=np.complex128)
    residual = float(np.linalg.norm(matrix.conj().T @ matrix - identity, ord=2))
    if not np.isfinite(residual) or residual > unitary_tolerance:
        raise ValueError(
            f"matrix is not unitary within tolerance {unitary_tolerance:g}; "
            f"residual={residual:g}"
        )
    return matrix, num_qubits, residual
