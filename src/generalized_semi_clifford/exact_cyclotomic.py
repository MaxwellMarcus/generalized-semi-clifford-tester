"""Exact matrices over a declared cyclotomic field.

This optional module is a representation layer, not an exact GSC tester.  It
stores entries as SymPy ``AlgebraicField`` domain elements and never promotes
floating-point values to exact scalars.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from functools import cached_property
from math import lcm
from numbers import Integral, Rational
from typing import Any

try:
    from sympy import QQ, Float, I, exp, pi, sqrt
    from sympy.polys.polyclasses import ANP
except ImportError as error:  # pragma: no cover - exercised in installations without the extra
    raise ImportError(
        "exact cyclotomic support requires the 'exact' extra: "
        "pip install generalized-semi-clifford-tester[exact]"
    ) from error


@dataclass(frozen=True)
class CyclotomicField:
    """The declared field ``Q(zeta_order)`` for an order divisible by eight."""

    order: int

    def __post_init__(self) -> None:
        if not isinstance(self.order, int) or isinstance(self.order, bool):
            raise TypeError("field order must be an integer")
        if self.order < 8 or self.order % 8:
            raise ValueError("field order must be a positive multiple of eight")

    @cached_property
    def domain(self):
        """SymPy algebraic field used for canonical storage."""

        return QQ.algebraic_field(exp(2 * pi * I / self.order))

    def scalar(self, value: ExactScalarInput) -> CyclotomicScalar:
        """Coerce an exact value into this field.

        Python and SymPy integers/rationals, this field's domain elements, and
        scalars from an embedded cyclotomic subfield are accepted.  Floats and
        complex values are deliberately rejected.
        """

        if isinstance(value, CyclotomicScalar):
            return value.embed(self)
        contains_float = isinstance(value, Float) or (
            hasattr(value, "has") and value.has(Float)
        )
        if isinstance(value, (float, complex)) or isinstance(value, bool) or contains_float:
            raise TypeError("floating-point and complex values are not exact inputs")
        if isinstance(value, ANP):
            if value.mod != self.domain.mod.to_list():
                raise ValueError("a raw domain element must belong to the declared field")
            domain_value = value
        elif isinstance(value, Integral):
            domain_value = self.domain.convert(int(value))
        elif isinstance(value, Rational):
            domain_value = self.domain.convert(QQ(value.numerator, value.denominator))
        else:
            try:
                domain_value = self.domain.from_sympy(value)
            except (TypeError, ValueError) as error:
                raise TypeError("value is not an exact element of the declared field") from error
        return CyclotomicScalar(self, domain_value)

    @property
    def zero(self) -> CyclotomicScalar:
        return CyclotomicScalar(self, self.domain.zero)

    @property
    def one(self) -> CyclotomicScalar:
        return CyclotomicScalar(self, self.domain.one)

    @property
    def zeta(self) -> CyclotomicScalar:
        return self.scalar(exp(2 * pi * I / self.order))

    def root_power(self, exponent: int) -> CyclotomicScalar:
        """Return ``zeta_order**exponent`` exactly."""

        if not isinstance(exponent, int) or isinstance(exponent, bool):
            raise TypeError("root exponent must be an integer")
        return self.zeta**exponent


ExactScalarInput = Any


@dataclass(frozen=True)
class CyclotomicScalar:
    """Canonical scalar stored as a SymPy algebraic-domain element."""

    field: CyclotomicField
    value: ANP

    def __post_init__(self) -> None:
        if (
            not isinstance(value := self.value, ANP)
            or value.mod != self.field.domain.mod.to_list()
        ):
            raise ValueError("scalar storage must belong to its declared field")

    def embed(self, target: CyclotomicField) -> CyclotomicScalar:
        """Embed into a cyclotomic superfield with compatible declared order."""

        if target == self.field:
            return self
        if target.order % self.field.order:
            raise ValueError("target field order must be a multiple of the source order")
        return CyclotomicScalar(
            target,
            target.domain.from_sympy(self.field.domain.to_sympy(self.value)),
        )

    def conjugate(self) -> CyclotomicScalar:
        """Apply the exact automorphism ``zeta -> zeta**-1``."""

        inverse_zeta = self.field.domain.one / self.field.zeta.value
        result = self.field.domain.zero
        for coefficient in self.value.to_list():
            result = result * inverse_zeta + self.field.domain.convert(coefficient)
        return CyclotomicScalar(self.field, result)

    def to_sympy(self):
        """Return an exact SymPy expression for inspection only."""

        return self.field.domain.to_sympy(self.value)

    def _coerce_pair(self, other: ExactScalarInput) -> tuple[CyclotomicScalar, CyclotomicScalar]:
        if isinstance(other, CyclotomicScalar):
            common = common_cyclotomic_field(self.field, other.field)
            return self.embed(common), other.embed(common)
        return self, self.field.scalar(other)

    def __add__(self, other: ExactScalarInput) -> CyclotomicScalar:
        left, right = self._coerce_pair(other)
        return CyclotomicScalar(left.field, left.value + right.value)

    __radd__ = __add__

    def __neg__(self) -> CyclotomicScalar:
        return CyclotomicScalar(self.field, -self.value)

    def __sub__(self, other: ExactScalarInput) -> CyclotomicScalar:
        return self + (-self._coerce_pair(other)[1])

    def __rsub__(self, other: ExactScalarInput) -> CyclotomicScalar:
        return self.field.scalar(other) - self

    def __mul__(self, other: ExactScalarInput) -> CyclotomicScalar:
        left, right = self._coerce_pair(other)
        return CyclotomicScalar(left.field, left.value * right.value)

    __rmul__ = __mul__

    def __truediv__(self, other: ExactScalarInput) -> CyclotomicScalar:
        left, right = self._coerce_pair(other)
        if right.value == right.field.domain.zero:
            raise ZeroDivisionError("division by zero")
        return CyclotomicScalar(left.field, left.value / right.value)

    def __pow__(self, exponent: int) -> CyclotomicScalar:
        if not isinstance(exponent, int) or isinstance(exponent, bool):
            raise TypeError("scalar exponent must be an integer")
        return CyclotomicScalar(self.field, self.value**exponent)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, CyclotomicScalar):
            return False
        common = common_cyclotomic_field(self.field, other.field)
        return self.embed(common).value == other.embed(common).value

    __hash__ = None


def common_cyclotomic_field(*fields: CyclotomicField) -> CyclotomicField:
    """Return the minimal declared cyclotomic field containing every input."""

    if not fields:
        raise ValueError("at least one field is required")
    if any(not isinstance(field, CyclotomicField) for field in fields):
        raise TypeError("common field inputs must be CyclotomicField instances")
    return CyclotomicField(lcm(*(field.order for field in fields)))


@dataclass(frozen=True)
class CyclotomicMatrix:
    """Immutable dense matrix over one declared cyclotomic field."""

    field: CyclotomicField
    rows: tuple[tuple[CyclotomicScalar, ...], ...]

    def __post_init__(self) -> None:
        if not self.rows or not self.rows[0]:
            raise ValueError("matrix must be nonempty")
        width = len(self.rows[0])
        if any(len(row) != width for row in self.rows):
            raise ValueError("matrix rows must have equal length")
        if any(entry.field != self.field for row in self.rows for entry in row):
            raise ValueError("every entry must belong to the matrix field")

    @classmethod
    def from_rows(
        cls,
        field: CyclotomicField,
        rows: Iterable[Iterable[ExactScalarInput]],
    ) -> CyclotomicMatrix:
        """Construct a matrix while enforcing exact scalar input."""

        return cls(field, tuple(tuple(field.scalar(entry) for entry in row) for row in rows))

    @classmethod
    def identity(cls, field: CyclotomicField, dimension: int) -> CyclotomicMatrix:
        if not isinstance(dimension, int) or isinstance(dimension, bool) or dimension < 1:
            raise ValueError("dimension must be a positive integer")
        return cls.from_rows(
            field,
            (
                (1 if row == column else 0 for column in range(dimension))
                for row in range(dimension)
            ),
        )

    @property
    def shape(self) -> tuple[int, int]:
        return len(self.rows), len(self.rows[0])

    def embed(self, target: CyclotomicField) -> CyclotomicMatrix:
        if target == self.field:
            return self
        return CyclotomicMatrix(
            target,
            tuple(tuple(entry.embed(target) for entry in row) for row in self.rows),
        )

    def adjoint(self) -> CyclotomicMatrix:
        row_count, column_count = self.shape
        return CyclotomicMatrix(
            self.field,
            tuple(
                tuple(self.rows[row][column].conjugate() for row in range(row_count))
                for column in range(column_count)
            ),
        )

    def __matmul__(self, other: CyclotomicMatrix) -> CyclotomicMatrix:
        if not isinstance(other, CyclotomicMatrix):
            return NotImplemented
        if self.shape[1] != other.shape[0]:
            raise ValueError("matrix dimensions do not align")
        field = common_cyclotomic_field(self.field, other.field)
        left, right = self.embed(field), other.embed(field)
        return CyclotomicMatrix(
            field,
            tuple(
                tuple(
                    sum(
                        (left.rows[row][index] * right.rows[index][column]
                         for index in range(left.shape[1])),
                        field.zero,
                    )
                    for column in range(right.shape[1])
                )
                for row in range(left.shape[0])
            ),
        )

    def tensor(self, other: CyclotomicMatrix) -> CyclotomicMatrix:
        """Return the exact Kronecker product."""

        if not isinstance(other, CyclotomicMatrix):
            raise TypeError("tensor operand must be a CyclotomicMatrix")
        field = common_cyclotomic_field(self.field, other.field)
        left, right = self.embed(field), other.embed(field)
        return CyclotomicMatrix(
            field,
            tuple(
                tuple(
                    left_entry * right_entry
                    for left_entry in left_row
                    for right_entry in right_row
                )
                for left_row in left.rows
                for right_row in right.rows
            ),
        )

    def is_unitary(self) -> bool:
        """Test unitarity by exact entry equality."""

        if self.shape[0] != self.shape[1]:
            return False
        return self.adjoint() @ self == CyclotomicMatrix.identity(self.field, self.shape[0])


def pauli_x(field: CyclotomicField) -> CyclotomicMatrix:
    return CyclotomicMatrix.from_rows(field, ((0, 1), (1, 0)))


def pauli_y(field: CyclotomicField) -> CyclotomicMatrix:
    imaginary = field.root_power(field.order // 4)
    return CyclotomicMatrix.from_rows(field, ((0, -imaginary), (imaginary, 0)))


def pauli_z(field: CyclotomicField) -> CyclotomicMatrix:
    return CyclotomicMatrix.from_rows(field, ((1, 0), (0, -1)))


def hadamard(field: CyclotomicField) -> CyclotomicMatrix:
    inverse_sqrt_two = field.scalar(1 / sqrt(2))
    return CyclotomicMatrix.from_rows(
        field,
        ((inverse_sqrt_two, inverse_sqrt_two), (inverse_sqrt_two, -inverse_sqrt_two)),
    )


def phase_s(field: CyclotomicField) -> CyclotomicMatrix:
    return CyclotomicMatrix.from_rows(field, ((1, 0), (0, field.root_power(field.order // 4))))


def t_gate(field: CyclotomicField, *, dagger: bool = False) -> CyclotomicMatrix:
    exponent = -field.order // 8 if dagger else field.order // 8
    return CyclotomicMatrix.from_rows(field, ((1, 0), (0, field.root_power(exponent))))


def permutation_gate(
    field: CyclotomicField,
    permutation: Sequence[int],
) -> CyclotomicMatrix:
    """Return the basis permutation mapping column ``j`` to row ``permutation[j]``."""

    permutation = tuple(permutation)
    dimension = len(permutation)
    if dimension < 1 or sorted(permutation) != list(range(dimension)):
        raise ValueError("permutation must contain every index exactly once")
    return CyclotomicMatrix.from_rows(
        field,
        ((1 if permutation[column] == row else 0 for column in range(dimension))
         for row in range(dimension)),
    )


def cnot(field: CyclotomicField) -> CyclotomicMatrix:
    return permutation_gate(field, (0, 1, 3, 2))


def swap(field: CyclotomicField) -> CyclotomicMatrix:
    return permutation_gate(field, (0, 2, 1, 3))


def toffoli(field: CyclotomicField) -> CyclotomicMatrix:
    return permutation_gate(field, (0, 1, 2, 3, 4, 5, 7, 6))


def controlled_phase(
    field: CyclotomicField,
    num_qubits: int,
    exponent: int,
) -> CyclotomicMatrix:
    """Phase the all-ones basis state by ``zeta_order**exponent``."""

    if not isinstance(num_qubits, int) or isinstance(num_qubits, bool) or num_qubits < 2:
        raise ValueError("controlled phase requires at least two qubits")
    dimension = 1 << num_qubits
    diagonal = [field.one] * dimension
    diagonal[-1] = field.root_power(exponent)
    return CyclotomicMatrix.from_rows(
        field,
        ((diagonal[row] if row == column else 0 for column in range(dimension))
         for row in range(dimension)),
    )


def cz(field: CyclotomicField) -> CyclotomicMatrix:
    return controlled_phase(field, 2, field.order // 2)


def ccz(field: CyclotomicField) -> CyclotomicMatrix:
    return controlled_phase(field, 3, field.order // 2)
