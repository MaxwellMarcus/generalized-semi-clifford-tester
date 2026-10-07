"""Exact odd-prime Weyl-coordinate and symplectic primitives.

Coordinates use ``(x | z)`` ordering over :math:`F_p`.  With
``X|j> = |j + 1>`` and ``Z|j> = omega**j |j>``, the package convention is

``W(x, z) = omega**((x dot z) / 2) X**x Z**z``.

The symplectic pairing is ``[v, w] = z dot x' - x dot z'``.  Consequently,
``W(v) W(w) = omega**([v, w] / 2) W(v + w)`` and
``W(v) W(w) = omega**[v, w] W(w) W(v)``.  All exponents and coordinates are
computed modulo the declared odd prime.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from functools import cache
from itertools import combinations, product
from math import isqrt, prod
from typing import TypeAlias

PrimeVector: TypeAlias = tuple[int, ...]
PrimeMatrix: TypeAlias = tuple[tuple[int, ...], ...]
MAX_PRIME_LAGRANGIAN_PRIME = 5
MAX_PRIME_LAGRANGIAN_QUDITS = 2


@dataclass(frozen=True)
class PrimeLagrangian:
    """A canonical odd-prime Lagrangian basis and its complete span."""

    prime: int
    num_qudits: int
    basis: tuple[PrimeVector, ...]
    elements: tuple[PrimeVector, ...]


def _validate_odd_prime(prime: int) -> int:
    if isinstance(prime, bool) or not isinstance(prime, int):
        raise TypeError("prime must be an integer")
    if prime < 3 or prime % 2 == 0:
        raise ValueError("prime must be an odd prime")
    if any(prime % divisor == 0 for divisor in range(3, isqrt(prime) + 1, 2)):
        raise ValueError("prime must be an odd prime")
    return prime


def _validate_num_qudits(num_qudits: int) -> int:
    if isinstance(num_qudits, bool) or not isinstance(num_qudits, int):
        raise TypeError("num_qudits must be an integer")
    if num_qudits < 1:
        raise ValueError("num_qudits must be positive")
    return num_qudits


def normalize_prime_label(label: Iterable[int], prime: int) -> PrimeVector:
    """Return a canonical ``(x | z)`` label over the declared odd prime."""

    modulus = _validate_odd_prime(prime)
    entries = tuple(label)
    if not entries or len(entries) % 2:
        raise ValueError("a Weyl label must have positive even width 2n")
    if any(isinstance(entry, bool) or not isinstance(entry, int) for entry in entries):
        raise TypeError("Weyl-label entries must be integers")
    return tuple(entry % modulus for entry in entries)


def prime_symplectic_pairing(
    left: Sequence[int],
    right: Sequence[int],
    prime: int,
) -> int:
    """Return ``z*x' - x*z'`` modulo ``prime`` for two ``(x | z)`` labels."""

    left_label = normalize_prime_label(left, prime)
    right_label = normalize_prime_label(right, prime)
    if len(left_label) != len(right_label):
        raise ValueError("Weyl labels must have the same width")
    num_qudits = len(left_label) // 2
    left_x, left_z = left_label[:num_qudits], left_label[num_qudits:]
    right_x, right_z = right_label[:num_qudits], right_label[num_qudits:]
    return (
        sum(a * b for a, b in zip(left_z, right_x, strict=True))
        - sum(a * b for a, b in zip(left_x, right_z, strict=True))
    ) % prime


def weyl_product_phase_exponent(
    left: Sequence[int],
    right: Sequence[int],
    prime: int,
) -> int:
    """Return ``e`` such that ``W(left) W(right) = omega**e W(left + right)``."""

    modulus = _validate_odd_prime(prime)
    inverse_two = pow(2, -1, modulus)
    return inverse_two * prime_symplectic_pairing(left, right, modulus) % modulus


def prime_standard_form(num_qudits: int, prime: int) -> PrimeMatrix:
    """Return ``[[0, -I], [I, 0]]`` for the package's odd-prime pairing."""

    modulus = _validate_odd_prime(prime)
    _validate_num_qudits(num_qudits)
    zero = (0,) * num_qudits
    identity = tuple(
        tuple(int(row == column) for column in range(num_qudits))
        for row in range(num_qudits)
    )
    top = tuple(zero + tuple((-entry) % modulus for entry in row) for row in identity)
    bottom = tuple(row + zero for row in identity)
    return top + bottom


def _normalize_matrix(matrix: Iterable[Iterable[int]], prime: int) -> PrimeMatrix:
    modulus = _validate_odd_prime(prime)
    rows = tuple(tuple(row) for row in matrix)
    if not rows:
        raise ValueError("a matrix must contain at least one row")
    width = len(rows[0])
    if width == 0 or any(len(row) != width for row in rows):
        raise ValueError("matrix rows must have one common nonzero width")
    if any(
        isinstance(entry, bool) or not isinstance(entry, int)
        for row in rows
        for entry in row
    ):
        raise TypeError("matrix entries must be integers")
    return tuple(tuple(entry % modulus for entry in row) for row in rows)


def _transpose(matrix: PrimeMatrix) -> PrimeMatrix:
    return tuple(
        tuple(matrix[row][column] for row in range(len(matrix)))
        for column in range(len(matrix[0]))
    )


def _matmul(left: PrimeMatrix, right: PrimeMatrix, prime: int) -> PrimeMatrix:
    if len(left[0]) != len(right):
        raise ValueError("inner matrix dimensions must agree")
    columns = _transpose(right)
    return tuple(
        tuple(
            sum(a * b for a, b in zip(row, column, strict=True)) % prime
            for column in columns
        )
        for row in left
    )


def is_prime_symplectic(matrix: Sequence[Sequence[int]], prime: int) -> bool:
    """Return whether a square matrix preserves the odd-prime standard form."""

    normalized = _normalize_matrix(matrix, prime)
    dimension = len(normalized)
    if dimension != len(normalized[0]) or dimension % 2:
        return False
    form = prime_standard_form(dimension // 2, prime)
    return _matmul(_matmul(_transpose(normalized), form, prime), normalized, prime) == form


def prime_lagrangian_count(num_qudits: int, prime: int) -> int:
    """Return the number of Lagrangian subspaces of ``F_prime^(2n)``."""

    dimension = _validate_num_qudits(num_qudits)
    modulus = _validate_odd_prime(prime)
    return prod(modulus**index + 1 for index in range(1, dimension + 1))


def _prime_span(basis: tuple[PrimeVector, ...], prime: int) -> tuple[PrimeVector, ...]:
    width = len(basis[0])
    return tuple(
        sorted(
            {
                tuple(
                    sum(coefficient * vector[column] for coefficient, vector in zip(
                        coefficients, basis, strict=True
                    ))
                    % prime
                    for column in range(width)
                )
                for coefficients in product(range(prime), repeat=len(basis))
            }
        )
    )


def _prime_rref(rows: PrimeMatrix, prime: int) -> PrimeMatrix:
    """Return reduced row echelon form over the declared prime field."""

    matrix = [list(row) for row in rows]
    pivot_row = 0
    for column in range(len(matrix[0])):
        pivot = next(
            (row for row in range(pivot_row, len(matrix)) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        inverse = pow(matrix[pivot_row][column], -1, prime)
        matrix[pivot_row] = [(entry * inverse) % prime for entry in matrix[pivot_row]]
        for row in range(len(matrix)):
            if row == pivot_row or not matrix[row][column]:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                (entry - factor * pivot_entry) % prime
                for entry, pivot_entry in zip(
                    matrix[row], matrix[pivot_row], strict=True
                )
            ]
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    return tuple(tuple(row) for row in matrix if any(row))


def prime_lagrangian_from_basis(
    basis: Iterable[Iterable[int]],
    num_qudits: int,
    prime: int,
) -> PrimeLagrangian:
    """Reconstruct a canonical odd-prime Lagrangian from a supplied basis.

    The basis must have rank ``num_qudits`` and span an isotropic subspace.
    Row operations are performed exactly over ``F_prime``; the returned basis
    is canonical and its full span is recomputed rather than trusted.
    """

    dimension = _validate_num_qudits(num_qudits)
    modulus = _validate_odd_prime(prime)
    rows = tuple(normalize_prime_label(row, modulus) for row in basis)
    if not rows:
        raise ValueError("a Lagrangian basis must contain at least one row")
    if any(len(row) != 2 * dimension for row in rows):
        raise ValueError("Lagrangian basis rows must have width 2n")
    canonical_basis = _prime_rref(rows, modulus)
    if len(canonical_basis) != dimension:
        raise ValueError("a Lagrangian basis must have rank n")
    if any(
        prime_symplectic_pairing(left, right, modulus)
        for left, right in combinations(canonical_basis, 2)
    ):
        raise ValueError("a Lagrangian basis must be isotropic")
    return PrimeLagrangian(
        prime=modulus,
        num_qudits=dimension,
        basis=canonical_basis,
        elements=_prime_span(canonical_basis, modulus),
    )


@cache
def _enumerate_prime_lagrangians(
    num_qudits: int,
    prime: int,
) -> tuple[PrimeLagrangian, ...]:
    width = 2 * num_qudits
    subspaces: list[PrimeLagrangian] = []
    for pivots in combinations(range(width), num_qudits):
        pivot_set = set(pivots)
        free_positions = tuple(
            (row, column)
            for row, pivot in enumerate(pivots)
            for column in range(pivot + 1, width)
            if column not in pivot_set
        )
        for assignment in product(range(prime), repeat=len(free_positions)):
            basis = [
                tuple(int(column == pivot) for column in range(width))
                for pivot in pivots
            ]
            mutable_basis = [list(row) for row in basis]
            for value, (row, column) in zip(assignment, free_positions, strict=True):
                mutable_basis[row][column] = value
            canonical_basis = tuple(tuple(row) for row in mutable_basis)
            if any(
                prime_symplectic_pairing(left, right, prime)
                for left, right in combinations(canonical_basis, 2)
            ):
                continue
            subspaces.append(
                PrimeLagrangian(
                    prime=prime,
                    num_qudits=num_qudits,
                    basis=canonical_basis,
                    elements=_prime_span(canonical_basis, prime),
                )
            )

    expected = prime_lagrangian_count(num_qudits, prime)
    if len(subspaces) != expected:
        raise AssertionError(
            f"expected {expected} odd-prime Lagrangians, found {len(subspaces)}"
        )
    return tuple(subspaces)


def enumerate_prime_lagrangians(
    num_qudits: int,
    prime: int,
) -> tuple[PrimeLagrangian, ...]:
    """Enumerate canonical odd-prime Lagrangians within explicit safety caps.

    Each subspace is represented by its unique reduced-row-echelon basis in
    ``(x | z)`` coordinates.  The intentionally small caps keep this exhaustive
    helper suitable for exact low-dimensional cross-checks rather than implying
    a scalable qudit search.
    """

    dimension = _validate_num_qudits(num_qudits)
    modulus = _validate_odd_prime(prime)
    if modulus > MAX_PRIME_LAGRANGIAN_PRIME:
        raise ValueError(
            "odd-prime Lagrangian enumeration is limited to primes at most "
            f"{MAX_PRIME_LAGRANGIAN_PRIME}"
        )
    if dimension > MAX_PRIME_LAGRANGIAN_QUDITS:
        raise ValueError(
            "odd-prime Lagrangian enumeration is limited to at most "
            f"{MAX_PRIME_LAGRANGIAN_QUDITS} qudits"
        )
    return _enumerate_prime_lagrangians(dimension, modulus)
