"""Dense Pauli matrices using binary ``(x | z)`` labels."""

from __future__ import annotations

from functools import cache

import numpy as np
from numpy.typing import NDArray

from .lagrangian import PauliLabel

ComplexMatrix = NDArray[np.complex128]

_I = np.eye(2, dtype=np.complex128)
_X = np.array([[0, 1], [1, 0]], dtype=np.complex128)
_Y = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
_Z = np.array([[1, 0], [0, -1]], dtype=np.complex128)


@cache
def all_pauli_labels(num_qubits: int) -> tuple[PauliLabel, ...]:
    """Return every phase-free Pauli label in integer order."""

    if num_qubits < 1:
        raise ValueError("num_qubits must be positive")
    width = 2 * num_qubits
    return tuple(
        tuple((value >> index) & 1 for index in range(width))
        for value in range(1 << width)
    )


@cache
def pauli_matrix(label: PauliLabel) -> ComplexMatrix:
    """Return the Hermitian Pauli associated with an ``(x | z)`` label."""

    if len(label) == 0 or len(label) % 2 or any(bit not in (0, 1) for bit in label):
        raise ValueError("a Pauli label must contain 2n binary entries")
    num_qubits = len(label) // 2
    matrix = np.array([[1]], dtype=np.complex128)
    for qubit in range(num_qubits):
        x, z = label[qubit], label[num_qubits + qubit]
        local = {(0, 0): _I, (1, 0): _X, (0, 1): _Z, (1, 1): _Y}[(x, z)]
        matrix = np.kron(matrix, local)
    matrix.setflags(write=False)
    return matrix


def infer_num_qubits(matrix: ComplexMatrix) -> int:
    """Infer qubit count from a square power-of-two matrix."""

    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or matrix.shape[0] < 2:
        raise ValueError("unitary must be a square matrix of dimension at least 2")
    dimension = matrix.shape[0]
    if dimension & (dimension - 1):
        raise ValueError("unitary dimension must be a power of two")
    return dimension.bit_length() - 1
