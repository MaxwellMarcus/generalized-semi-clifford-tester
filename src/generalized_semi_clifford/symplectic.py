"""Small exact matrix primitives over the field with two elements."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import TypeAlias

BinaryMatrix: TypeAlias = tuple[tuple[int, ...], ...]


def _normalize(matrix: Iterable[Iterable[int]]) -> BinaryMatrix:
    rows = tuple(tuple(entry for entry in row) for row in matrix)
    if not rows:
        raise ValueError("a matrix must contain at least one row")
    width = len(rows[0])
    if width == 0 or any(len(row) != width for row in rows):
        raise ValueError("matrix rows must have one common nonzero width")
    if any(entry not in (0, 1) for row in rows for entry in row):
        raise ValueError("matrix entries must be 0 or 1")
    return rows


def identity(size: int) -> BinaryMatrix:
    """Return the ``size`` by ``size`` identity matrix over :math:`F_2`."""

    if size < 1:
        raise ValueError("size must be positive")
    return tuple(tuple(int(row == column) for column in range(size)) for row in range(size))


def transpose(matrix: Sequence[Sequence[int]]) -> BinaryMatrix:
    """Return a validated matrix transpose."""

    normalized = _normalize(matrix)
    return tuple(
        tuple(normalized[row][column] for row in range(len(normalized)))
        for column in range(len(normalized[0]))
    )


def matmul(left: Sequence[Sequence[int]], right: Sequence[Sequence[int]]) -> BinaryMatrix:
    """Multiply two binary matrices exactly over :math:`F_2`."""

    left_matrix = _normalize(left)
    right_matrix = _normalize(right)
    if len(left_matrix[0]) != len(right_matrix):
        raise ValueError("inner matrix dimensions must agree")
    right_transpose = transpose(right_matrix)
    return tuple(
        tuple(
            sum(a * b for a, b in zip(row, column, strict=True)) % 2
            for column in right_transpose
        )
        for row in left_matrix
    )


def standard_form(num_qubits: int) -> BinaryMatrix:
    """Return ``[[0, I], [I, 0]]`` for ``(x | z)`` Pauli coordinates."""

    if num_qubits < 1:
        raise ValueError("num_qubits must be positive")
    zero = (0,) * num_qubits
    top = tuple(zero + row for row in identity(num_qubits))
    bottom = tuple(row + zero for row in identity(num_qubits))
    return top + bottom


def is_symplectic(matrix: Sequence[Sequence[int]]) -> bool:
    """Return whether a square binary matrix preserves the standard form."""

    normalized = _normalize(matrix)
    dimension = len(normalized)
    if dimension != len(normalized[0]) or dimension % 2:
        return False
    form = standard_form(dimension // 2)
    return matmul(matmul(transpose(normalized), form), normalized) == form
