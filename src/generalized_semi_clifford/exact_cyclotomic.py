"""Exact matrices and bounded Pauli-MASA search over a cyclotomic field.

This optional module stores entries as SymPy ``AlgebraicField`` domain
elements and never promotes floating-point values to exact scalars. Its exact
GSC classification is limited to completed, explicitly resource-bounded
Pauli-MASA searches.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from enum import Enum
from functools import cached_property
from math import lcm
from numbers import Integral, Rational
from typing import Any

from .lagrangian import (
    MAX_LAGRANGIAN_QUBITS,
    Lagrangian,
    PauliLabel,
    enumerate_lagrangians,
    lagrangian_from_basis,
)

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


@dataclass(frozen=True)
class ExactPauliCoefficient:
    """One exactly computed coefficient in a Pauli expansion."""

    output_label: PauliLabel
    value: CyclotomicScalar


@dataclass(frozen=True)
class ExactPauliExpansion:
    """The complete exact Pauli expansion of one conjugated input Pauli."""

    input_label: PauliLabel
    coefficients: tuple[ExactPauliCoefficient, ...]

    @property
    def support(self) -> tuple[ExactPauliCoefficient, ...]:
        """Return coefficients proved nonzero by exact field equality."""

        return tuple(
            coefficient
            for coefficient in self.coefficients
            if coefficient.value.value != coefficient.value.field.domain.zero
        )


@dataclass(frozen=True)
class ExactWitnessVerification:
    """Exact verification result for one caller-supplied Lagrangian pair.

    This type records no search status: it verifies only the supplied pair and
    makes no claim that a different witness does or does not exist.
    """

    field_order: int
    arithmetic: str
    num_qubits: int
    input_lagrangian: Lagrangian
    output_lagrangian: Lagrangian
    images: tuple[ExactPauliExpansion, ...]
    coefficients_computed: int
    verified: bool


EXACT_SEARCH_SCHEMA = "generalized-semi-clifford.exact-search.v1"


class ExactSearchStatus(str, Enum):
    """Outcomes from the resource-bounded exhaustive exact search."""

    GSC = "gsc"
    NOT_GSC = "not_gsc"
    UNKNOWN = "unknown"


@dataclass(frozen=True)
class ExactSearchResult:
    """Versioned result from :func:`search_exact_lagrangian_witness`.

    ``search_complete`` means every input/output Lagrangian pair was rejected;
    it is therefore true only for an exact ``NOT_GSC`` result. A positive
    result instead carries an independently recomputed witness verification.
    """

    schema_version: str
    status: ExactSearchStatus
    field_order: int
    arithmetic: str
    num_qubits: int
    search_complete: bool
    max_qubits: int
    max_coefficients: int
    max_candidate_pairs: int
    lagrangians_total: int | None
    input_lagrangians_checked: int
    candidate_pairs_checked: int
    coefficients_computed: int
    witness: ExactWitnessVerification | None
    stop_reason: str | None
    message: str

    @property
    def is_gsc(self) -> bool | None:
        """Return a Boolean classification, or ``None`` for ``UNKNOWN``."""

        if self.status is ExactSearchStatus.UNKNOWN:
            return None
        return self.status is ExactSearchStatus.GSC

    def to_dict(self) -> dict[str, Any]:
        """Return the stable, JSON-serializable version-one result schema."""

        witness = None
        if self.witness is not None:
            witness = {
                "verified": self.witness.verified,
                "input_basis": [list(label) for label in self.witness.input_lagrangian.basis],
                "output_basis": [list(label) for label in self.witness.output_lagrangian.basis],
                "verification_coefficients_computed": self.witness.coefficients_computed,
            }
        return {
            "schema_version": self.schema_version,
            "status": self.status.value,
            "field_order": self.field_order,
            "arithmetic": self.arithmetic,
            "num_qubits": self.num_qubits,
            "search_complete": self.search_complete,
            "limits": {
                "max_qubits": self.max_qubits,
                "max_coefficients": self.max_coefficients,
                "max_candidate_pairs": self.max_candidate_pairs,
            },
            "work": {
                "lagrangians_total": self.lagrangians_total,
                "input_lagrangians_checked": self.input_lagrangians_checked,
                "candidate_pairs_checked": self.candidate_pairs_checked,
                "coefficients_computed": self.coefficients_computed,
            },
            "witness": witness,
            "stop_reason": self.stop_reason,
            "message": self.message,
        }


def pauli_x(field: CyclotomicField) -> CyclotomicMatrix:
    return CyclotomicMatrix.from_rows(field, ((0, 1), (1, 0)))


def pauli_y(field: CyclotomicField) -> CyclotomicMatrix:
    imaginary = field.root_power(field.order // 4)
    return CyclotomicMatrix.from_rows(field, ((0, -imaginary), (imaginary, 0)))


def pauli_z(field: CyclotomicField) -> CyclotomicMatrix:
    return CyclotomicMatrix.from_rows(field, ((1, 0), (0, -1)))


def exact_pauli_matrix(
    field: CyclotomicField,
    label: Sequence[int],
) -> CyclotomicMatrix:
    """Return the Hermitian tensor Pauli for an ``(x | z)`` binary label."""

    label = tuple(label)
    if not label or len(label) % 2 or any(bit not in (0, 1) for bit in label):
        raise ValueError("a Pauli label must contain 2n binary entries")
    num_qubits = len(label) // 2
    matrix = CyclotomicMatrix.from_rows(field, ((1,),))
    for qubit in range(num_qubits):
        x, z = label[qubit], label[num_qubits + qubit]
        local = {
            (0, 0): CyclotomicMatrix.identity(field, 2),
            (1, 0): pauli_x(field),
            (0, 1): pauli_z(field),
            (1, 1): pauli_y(field),
        }[(x, z)]
        matrix = matrix.tensor(local)
    return matrix


def _matrix_num_qubits(matrix: CyclotomicMatrix) -> int:
    rows, columns = matrix.shape
    if rows != columns or rows < 2 or rows & (rows - 1):
        raise ValueError("matrix must be square with power-of-two dimension at least two")
    return rows.bit_length() - 1


def _all_pauli_labels(num_qubits: int) -> tuple[PauliLabel, ...]:
    width = 2 * num_qubits
    return tuple(
        tuple((value >> index) & 1 for index in range(width))
        for value in range(1 << width)
    )


def pauli_conjugation_coefficients(
    unitary: CyclotomicMatrix,
    input_label: Sequence[int],
) -> ExactPauliExpansion:
    """Compute every coefficient of ``U P U.adjoint()`` exactly.

    Coefficients are returned in binary-label integer order, including exact
    zeros. The input matrix must be an exactly validated unitary; approximate
    matrices belong to the separate numerical API.
    """

    num_qubits = _matrix_num_qubits(unitary)
    input_label = tuple(input_label)
    if len(input_label) != 2 * num_qubits or any(bit not in (0, 1) for bit in input_label):
        raise ValueError("input Pauli label must be binary and have width 2n")
    if not unitary.is_unitary():
        raise ValueError("exact Pauli coefficients require an exactly unitary matrix")
    return _pauli_conjugation_coefficients_validated(unitary, input_label)


def _pauli_conjugation_coefficients_validated(
    unitary: CyclotomicMatrix,
    input_label: PauliLabel,
) -> ExactPauliExpansion:
    """Compute coefficients after shape, label, and unitarity validation."""

    num_qubits = _matrix_num_qubits(unitary)
    conjugated = unitary @ exact_pauli_matrix(unitary.field, input_label) @ unitary.adjoint()
    dimension = 1 << num_qubits
    coefficients = []
    for output_label in _all_pauli_labels(num_qubits):
        product = exact_pauli_matrix(unitary.field, output_label) @ conjugated
        trace = sum(
            (product.rows[index][index] for index in range(dimension)),
            unitary.field.zero,
        )
        coefficients.append(
            ExactPauliCoefficient(output_label, trace / dimension)
        )
    return ExactPauliExpansion(input_label, tuple(coefficients))


def verify_exact_lagrangian_witness(
    unitary: CyclotomicMatrix,
    input_lagrangian: Lagrangian,
    output_lagrangian: Lagrangian,
) -> ExactWitnessVerification:
    """Independently verify one exact input/output Pauli-MASA witness.

    The Lagrangian bases and spans are revalidated rather than trusted. Only
    the supplied pair is checked; a false result is not an exhaustive negative
    GSC classification.
    """

    num_qubits = _matrix_num_qubits(unitary)
    if not unitary.is_unitary():
        raise ValueError("exact witness verification requires an exactly unitary matrix")
    if not isinstance(input_lagrangian, Lagrangian) or not isinstance(
        output_lagrangian, Lagrangian
    ):
        raise TypeError("input and output witnesses must be Lagrangian instances")
    if (
        input_lagrangian.num_qubits != num_qubits
        or output_lagrangian.num_qubits != num_qubits
    ):
        raise ValueError("witness qubit counts must match the unitary")

    canonical_input = lagrangian_from_basis(input_lagrangian.basis, num_qubits)
    canonical_output = lagrangian_from_basis(output_lagrangian.basis, num_qubits)
    if canonical_input != input_lagrangian or canonical_output != output_lagrangian:
        raise ValueError("witness elements must equal the canonical span of their basis")

    images = tuple(
        _pauli_conjugation_coefficients_validated(unitary, label)
        for label in canonical_input.basis
    )
    output_labels = set(canonical_output.elements)
    verified = all(
        coefficient.output_label in output_labels
        for image in images
        for coefficient in image.support
    )
    return ExactWitnessVerification(
        field_order=unitary.field.order,
        arithmetic="cyclotomic_exact",
        num_qubits=num_qubits,
        input_lagrangian=canonical_input,
        output_lagrangian=canonical_output,
        images=images,
        coefficients_computed=num_qubits * (1 << (2 * num_qubits)),
        verified=verified,
    )


def search_exact_lagrangian_witness(
    unitary: CyclotomicMatrix,
    *,
    max_qubits: int = 3,
    max_coefficients: int = 100_000,
    max_candidate_pairs: int = 100_000,
) -> ExactSearchResult:
    """Search every Pauli-MASA pair subject to explicit exact-work caps.

    A cap reached before a witness is independently verified returns
    ``UNKNOWN``. ``NOT_GSC`` is returned only after every enumerated pair has
    been rejected using exact field-zero tests.
    """

    for name, value in (
        ("max_qubits", max_qubits),
        ("max_coefficients", max_coefficients),
        ("max_candidate_pairs", max_candidate_pairs),
    ):
        if not isinstance(value, int) or isinstance(value, bool) or value < 0:
            raise ValueError(f"{name} must be a nonnegative integer")
    if max_qubits < 1:
        raise ValueError("max_qubits must be positive")

    num_qubits = _matrix_num_qubits(unitary)
    if not unitary.is_unitary():
        raise ValueError("exact search requires an exactly unitary matrix")

    def result(
        status: ExactSearchStatus,
        *,
        search_complete: bool,
        lagrangians_total: int | None,
        input_lagrangians_checked: int,
        candidate_pairs_checked: int,
        coefficients_computed: int,
        witness: ExactWitnessVerification | None = None,
        stop_reason: str | None = None,
        message: str,
    ) -> ExactSearchResult:
        return ExactSearchResult(
            schema_version=EXACT_SEARCH_SCHEMA,
            status=status,
            field_order=unitary.field.order,
            arithmetic="cyclotomic_exact",
            num_qubits=num_qubits,
            search_complete=search_complete,
            max_qubits=max_qubits,
            max_coefficients=max_coefficients,
            max_candidate_pairs=max_candidate_pairs,
            lagrangians_total=lagrangians_total,
            input_lagrangians_checked=input_lagrangians_checked,
            candidate_pairs_checked=candidate_pairs_checked,
            coefficients_computed=coefficients_computed,
            witness=witness,
            stop_reason=stop_reason,
            message=message,
        )

    search_limit = min(max_qubits, MAX_LAGRANGIAN_QUBITS)
    if num_qubits > search_limit:
        return result(
            ExactSearchStatus.UNKNOWN,
            search_complete=False,
            lagrangians_total=None,
            input_lagrangians_checked=0,
            candidate_pairs_checked=0,
            coefficients_computed=0,
            stop_reason="qubit_cap",
            message=f"n={num_qubits} exceeds the exact-search limit {search_limit}",
        )

    lagrangians = enumerate_lagrangians(num_qubits)
    coefficients_per_image = 1 << (2 * num_qubits)
    coefficients_computed = 0
    candidate_pairs_checked = 0
    input_lagrangians_checked = 0

    for input_lagrangian in lagrangians:
        images: list[ExactPauliExpansion] = []
        for label in input_lagrangian.basis:
            if coefficients_computed + coefficients_per_image > max_coefficients:
                return result(
                    ExactSearchStatus.UNKNOWN,
                    search_complete=False,
                    lagrangians_total=len(lagrangians),
                    input_lagrangians_checked=input_lagrangians_checked,
                    candidate_pairs_checked=candidate_pairs_checked,
                    coefficients_computed=coefficients_computed,
                    stop_reason="coefficient_cap",
                    message="exact Pauli-coefficient work cap interrupted the search",
                )
            images.append(_pauli_conjugation_coefficients_validated(unitary, label))
            coefficients_computed += coefficients_per_image

        support = {
            coefficient.output_label
            for image in images
            for coefficient in image.support
        }
        for output_lagrangian in lagrangians:
            if candidate_pairs_checked >= max_candidate_pairs:
                return result(
                    ExactSearchStatus.UNKNOWN,
                    search_complete=False,
                    lagrangians_total=len(lagrangians),
                    input_lagrangians_checked=input_lagrangians_checked,
                    candidate_pairs_checked=candidate_pairs_checked,
                    coefficients_computed=coefficients_computed,
                    stop_reason="candidate_pair_cap",
                    message="input/output candidate-pair cap interrupted the search",
                )
            candidate_pairs_checked += 1
            if not support.issubset(output_lagrangian.elements):
                continue

            verification_work = num_qubits * coefficients_per_image
            if coefficients_computed + verification_work > max_coefficients:
                return result(
                    ExactSearchStatus.UNKNOWN,
                    search_complete=False,
                    lagrangians_total=len(lagrangians),
                    input_lagrangians_checked=input_lagrangians_checked,
                    candidate_pairs_checked=candidate_pairs_checked,
                    coefficients_computed=coefficients_computed,
                    stop_reason="coefficient_cap",
                    message="coefficient cap cannot fund independent witness verification",
                )
            witness = verify_exact_lagrangian_witness(
                unitary,
                input_lagrangian,
                output_lagrangian,
            )
            coefficients_computed += witness.coefficients_computed
            if not witness.verified:
                raise AssertionError("independent exact witness verification disagreed with search")
            return result(
                ExactSearchStatus.GSC,
                search_complete=False,
                lagrangians_total=len(lagrangians),
                input_lagrangians_checked=input_lagrangians_checked + 1,
                candidate_pairs_checked=candidate_pairs_checked,
                coefficients_computed=coefficients_computed,
                witness=witness,
                stop_reason="witness_found",
                message="found and independently verified an exact Pauli-MASA witness",
            )
        input_lagrangians_checked += 1

    return result(
        ExactSearchStatus.NOT_GSC,
        search_complete=True,
        lagrangians_total=len(lagrangians),
        input_lagrangians_checked=input_lagrangians_checked,
        candidate_pairs_checked=candidate_pairs_checked,
        coefficients_computed=coefficients_computed,
        stop_reason="exhausted",
        message="no Pauli-MASA pair passed the exhaustive exact search",
    )


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
