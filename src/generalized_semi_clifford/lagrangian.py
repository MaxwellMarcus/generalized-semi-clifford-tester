"""Exact enumeration of small binary Lagrangian subspaces."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from functools import cache
from itertools import combinations
from math import prod
from typing import TypeAlias

PauliLabel: TypeAlias = tuple[int, ...]
MAX_EXHAUSTIVE_QUBITS = 3


@dataclass(frozen=True)
class Lagrangian:
    """A binary Lagrangian represented by a canonical basis and all elements."""

    num_qubits: int
    basis: tuple[PauliLabel, ...]
    elements: tuple[PauliLabel, ...]


def lagrangian_count(num_qubits: int) -> int:
    """Return the number of Lagrangians in ``F_2^(2n)``."""

    if num_qubits < 1:
        raise ValueError("num_qubits must be positive")
    return prod((1 << index) + 1 for index in range(1, num_qubits + 1))


def _label_to_int(label: PauliLabel) -> int:
    return sum(bit << index for index, bit in enumerate(label))


def _int_to_label(value: int, width: int) -> PauliLabel:
    return tuple((value >> index) & 1 for index in range(width))


def _symplectic_pairing(left: int, right: int, num_qubits: int) -> int:
    mask = (1 << num_qubits) - 1
    left_x, left_z = left & mask, left >> num_qubits
    right_x, right_z = right & mask, right >> num_qubits
    return ((left_x & right_z).bit_count() + (left_z & right_x).bit_count()) % 2


def _row_reduce(vectors: tuple[int, ...], width: int) -> tuple[int, ...]:
    rows = list(dict.fromkeys(vector for vector in vectors if vector))
    pivot_row = 0
    for column in reversed(range(width)):
        pivot = next(
            (index for index in range(pivot_row, len(rows)) if rows[index] >> column & 1),
            None,
        )
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        for index in range(len(rows)):
            if index != pivot_row and rows[index] >> column & 1:
                rows[index] ^= rows[pivot_row]
        pivot_row += 1
    return tuple(rows[:pivot_row])


def _span(basis: tuple[int, ...]) -> tuple[int, ...]:
    elements = [0]
    for vector in basis:
        elements += [element ^ vector for element in elements]
    return tuple(sorted(elements))


@cache
def enumerate_lagrangians(num_qubits: int) -> tuple[Lagrangian, ...]:
    """Exhaustively enumerate binary Lagrangians.

    The combination-based algorithm is intentionally transparent and enforces
    the package's three-qubit safety limit.
    """

    if num_qubits < 1:
        raise ValueError("num_qubits must be positive")
    if num_qubits > MAX_EXHAUSTIVE_QUBITS:
        raise ValueError(
            f"naive enumeration is limited to {MAX_EXHAUSTIVE_QUBITS} qubits"
        )
    width = 2 * num_qubits
    subspaces: dict[tuple[int, ...], tuple[int, ...]] = {}
    for candidate in combinations(range(1, 1 << width), num_qubits):
        basis = _row_reduce(candidate, width)
        if len(basis) != num_qubits:
            continue
        if any(
            _symplectic_pairing(left, right, num_qubits)
            for left, right in combinations(basis, 2)
        ):
            continue
        elements = _span(basis)
        subspaces.setdefault(elements, basis)

    expected = lagrangian_count(num_qubits)
    if len(subspaces) != expected:
        raise AssertionError(f"expected {expected} Lagrangians, found {len(subspaces)}")

    return tuple(
        Lagrangian(
            num_qubits=num_qubits,
            basis=tuple(_int_to_label(vector, width) for vector in basis),
            elements=tuple(_int_to_label(vector, width) for vector in elements),
        )
        for elements, basis in sorted(subspaces.items())
    )


def is_lagrangian(labels: tuple[PauliLabel, ...], num_qubits: int) -> bool:
    """Return whether ``labels`` are exactly a binary Lagrangian subspace."""

    width = 2 * num_qubits
    if any(len(label) != width or any(bit not in (0, 1) for bit in label) for label in labels):
        return False
    vectors = tuple(_label_to_int(label) for label in labels)
    vector_set = set(vectors)
    if len(vector_set) != 1 << num_qubits or 0 not in vector_set:
        return False
    if any(left ^ right not in vector_set for left in vector_set for right in vector_set):
        return False
    return not any(
        _symplectic_pairing(left, right, num_qubits)
        for left, right in combinations(vector_set, 2)
    )


def lagrangian_from_basis(
    basis: Iterable[Sequence[int]],
    num_qubits: int,
) -> Lagrangian:
    """Validate and canonicalize a basis of a binary Lagrangian."""

    if num_qubits < 1:
        raise ValueError("num_qubits must be positive")
    normalized = tuple(tuple(label) for label in basis)
    width = 2 * num_qubits
    if len(normalized) != num_qubits:
        raise ValueError(f"a Lagrangian basis must contain exactly {num_qubits} labels")
    if any(
        len(label) != width or any(bit not in (0, 1) for bit in label)
        for label in normalized
    ):
        raise ValueError("a Lagrangian basis must contain binary labels of width 2n")
    vectors = tuple(_label_to_int(label) for label in normalized)
    reduced = _row_reduce(vectors, width)
    if len(reduced) != num_qubits:
        raise ValueError("Lagrangian basis Paulis must be linearly independent")
    if any(
        _symplectic_pairing(left, right, num_qubits)
        for left, right in combinations(reduced, 2)
    ):
        raise ValueError("Lagrangian basis Paulis must commute pairwise")
    return Lagrangian(
        num_qubits=num_qubits,
        basis=tuple(_int_to_label(vector, width) for vector in reduced),
        elements=tuple(_int_to_label(vector, width) for vector in _span(reduced)),
    )
