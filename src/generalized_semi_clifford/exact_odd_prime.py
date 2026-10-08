"""Bounded exact witness verification and search for odd-prime Weyl systems."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from enum import Enum
from itertools import product
from typing import Any

from .exact_cyclotomic import (
    CyclotomicField,
    CyclotomicMatrix,
    CyclotomicScalar,
)
from .odd_prime import (
    PrimeLagrangian,
    PrimeVector,
    enumerate_prime_lagrangians,
    normalize_prime_label,
    prime_lagrangian_from_basis,
)

DEFAULT_MAX_EXACT_PRIME_QUDITS = 2
DEFAULT_MAX_EXACT_PRIME_COEFFICIENTS = 100_000
DEFAULT_MAX_EXACT_PRIME_CANDIDATE_PAIRS = 100_000
EXACT_PRIME_SEARCH_SCHEMA = "generalized-semi-clifford.exact-prime-search.v1"


@dataclass(frozen=True)
class ExactPrimeWeylCoefficient:
    """One exactly computed coefficient in an odd-prime Weyl expansion."""

    output_label: PrimeVector
    value: CyclotomicScalar


@dataclass(frozen=True)
class ExactPrimeWeylExpansion:
    """The complete exact Weyl expansion of one conjugated input operator."""

    input_label: PrimeVector
    coefficients: tuple[ExactPrimeWeylCoefficient, ...]

    @property
    def support(self) -> tuple[ExactPrimeWeylCoefficient, ...]:
        """Return coefficients proved nonzero by exact field equality."""

        return tuple(
            coefficient
            for coefficient in self.coefficients
            if coefficient.value.value != coefficient.value.field.domain.zero
        )


class ExactPrimeWitnessStatus(str, Enum):
    """Outcomes for one resource-bounded odd-prime fixed-witness check."""

    VERIFIED = "verified"
    REJECTED = "rejected"
    UNKNOWN = "unknown"


@dataclass(frozen=True)
class ExactPrimeWitnessVerification:
    """Separate result contract for one bounded odd-prime witness check."""

    status: ExactPrimeWitnessStatus
    field_order: int
    arithmetic: str
    prime: int
    num_qudits: int
    input_lagrangian: PrimeLagrangian
    output_lagrangian: PrimeLagrangian
    images: tuple[ExactPrimeWeylExpansion, ...]
    complete: bool
    coefficients_computed: int
    coefficients_required: int
    max_qudits: int
    max_coefficients: int
    stop_reason: str

    @property
    def verified(self) -> bool | None:
        """Return the supplied-pair verdict, or ``None`` when capped."""

        if self.status is ExactPrimeWitnessStatus.UNKNOWN:
            return None
        return self.status is ExactPrimeWitnessStatus.VERIFIED


class ExactPrimeSearchStatus(str, Enum):
    """Outcomes from the resource-bounded odd-prime exact search."""

    GSC = "gsc"
    NOT_GSC = "not_gsc"
    UNKNOWN = "unknown"


@dataclass(frozen=True)
class ExactPrimeSearchResult:
    """Versioned result from :func:`search_exact_prime_lagrangian_witness`.

    ``search_complete`` is true only when every input/output Lagrangian pair
    was rejected exactly. Positive results instead carry an independently
    recomputed fixed-witness verification.
    """

    schema_version: str
    status: ExactPrimeSearchStatus
    field_order: int
    arithmetic: str
    prime: int
    num_qudits: int
    search_complete: bool
    max_qudits: int
    max_coefficients: int
    max_candidate_pairs: int
    lagrangians_total: int | None
    input_lagrangians_checked: int
    candidate_pairs_checked: int
    coefficients_computed: int
    witness: ExactPrimeWitnessVerification | None
    stop_reason: str
    message: str

    @property
    def is_gsc(self) -> bool | None:
        """Return a Boolean classification, or ``None`` for ``UNKNOWN``."""

        if self.status is ExactPrimeSearchStatus.UNKNOWN:
            return None
        return self.status is ExactPrimeSearchStatus.GSC

    def to_dict(self) -> dict[str, Any]:
        """Return the stable, JSON-serializable version-one result schema."""

        witness = None
        if self.witness is not None:
            witness = {
                "status": self.witness.status.value,
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
            "prime": self.prime,
            "num_qudits": self.num_qudits,
            "search_complete": self.search_complete,
            "limits": {
                "max_qudits": self.max_qudits,
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


def exact_prime_weyl_matrix(
    field: CyclotomicField,
    label: Sequence[int],
    prime: int,
) -> CyclotomicMatrix:
    """Return ``W(x,z)`` exactly in the package's odd-prime convention."""

    normalized = normalize_prime_label(label, prime)
    if field.order % prime:
        raise ValueError("field order must be divisible by the Weyl prime")
    num_qudits = len(normalized) // 2
    x_values = normalized[:num_qudits]
    z_values = normalized[num_qudits:]
    omega = field.root_power(field.order // prime)
    inverse_two = pow(2, -1, prime)

    matrix = CyclotomicMatrix.from_rows(field, ((1,),))
    for x_value, z_value in zip(x_values, z_values, strict=True):
        phase_offset = inverse_two * x_value * z_value % prime
        local = CyclotomicMatrix.from_rows(
            field,
            (
                (
                    omega ** ((phase_offset + z_value * column) % prime)
                    if row == (column + x_value) % prime
                    else 0
                    for column in range(prime)
                )
                for row in range(prime)
            ),
        )
        matrix = matrix.tensor(local)
    return matrix


def _prime_power_qudits(dimension: int, prime: int) -> int:
    num_qudits = 0
    residual = dimension
    while residual > 1 and residual % prime == 0:
        residual //= prime
        num_qudits += 1
    if residual != 1 or num_qudits == 0:
        raise ValueError("matrix dimension must be a positive power of the Weyl prime")
    return num_qudits


def _all_prime_labels(num_qudits: int, prime: int) -> tuple[PrimeVector, ...]:
    return tuple(product(range(prime), repeat=2 * num_qudits))


def _prime_conjugation_coefficients_validated(
    unitary: CyclotomicMatrix,
    input_label: PrimeVector,
    prime: int,
) -> ExactPrimeWeylExpansion:
    num_qudits = len(input_label) // 2
    dimension = prime**num_qudits
    conjugated = (
        unitary
        @ exact_prime_weyl_matrix(unitary.field, input_label, prime)
        @ unitary.adjoint()
    )
    coefficients = []
    for output_label in _all_prime_labels(num_qudits, prime):
        output = exact_prime_weyl_matrix(unitary.field, output_label, prime)
        product_matrix = output.adjoint() @ conjugated
        trace = sum(
            (product_matrix.rows[index][index] for index in range(dimension)),
            unitary.field.zero,
        )
        coefficients.append(
            ExactPrimeWeylCoefficient(output_label, trace / dimension)
        )
    return ExactPrimeWeylExpansion(input_label, tuple(coefficients))


def verify_exact_prime_lagrangian_witness(
    unitary: CyclotomicMatrix,
    input_lagrangian: PrimeLagrangian,
    output_lagrangian: PrimeLagrangian,
    *,
    max_qudits: int = DEFAULT_MAX_EXACT_PRIME_QUDITS,
    max_coefficients: int = DEFAULT_MAX_EXACT_PRIME_COEFFICIENTS,
) -> ExactPrimeWitnessVerification:
    """Verify one exact odd-prime witness under explicit work caps.

    Invalid or inconsistent witness objects are rejected with exceptions.
    Resource limits instead produce ``UNKNOWN`` so incomplete work is never
    reported as rejection of the supplied pair.
    """

    for name, value in (("max_qudits", max_qudits), ("max_coefficients", max_coefficients)):
        if not isinstance(value, int) or isinstance(value, bool) or value < 0:
            raise ValueError(f"{name} must be a nonnegative integer")
    if not isinstance(input_lagrangian, PrimeLagrangian) or not isinstance(
        output_lagrangian, PrimeLagrangian
    ):
        raise TypeError("input and output witnesses must be PrimeLagrangian instances")
    if input_lagrangian.prime != output_lagrangian.prime:
        raise ValueError("input and output witnesses must use the same prime")
    prime = input_lagrangian.prime
    if unitary.field.order % prime:
        raise ValueError("field order must be divisible by the witness prime")
    rows, columns = unitary.shape
    if rows != columns:
        raise ValueError("exact witness verification requires a square matrix")
    num_qudits = _prime_power_qudits(rows, prime)
    if not unitary.is_unitary():
        raise ValueError("exact witness verification requires an exactly unitary matrix")
    if (
        input_lagrangian.num_qudits != num_qudits
        or output_lagrangian.num_qudits != num_qudits
    ):
        raise ValueError("witness qudit counts must match the unitary")

    canonical_input = prime_lagrangian_from_basis(
        input_lagrangian.basis, num_qudits, prime
    )
    canonical_output = prime_lagrangian_from_basis(
        output_lagrangian.basis, num_qudits, prime
    )
    if (
        canonical_input.elements != input_lagrangian.elements
        or canonical_output.elements != output_lagrangian.elements
    ):
        raise ValueError("witness elements must equal the canonical span of their basis")

    coefficients_per_image = prime ** (2 * num_qudits)
    coefficients_required = num_qudits * coefficients_per_image

    def result(
        status: ExactPrimeWitnessStatus,
        images: tuple[ExactPrimeWeylExpansion, ...],
        *,
        complete: bool,
        stop_reason: str,
    ) -> ExactPrimeWitnessVerification:
        return ExactPrimeWitnessVerification(
            status=status,
            field_order=unitary.field.order,
            arithmetic="cyclotomic_exact_odd_prime",
            prime=prime,
            num_qudits=num_qudits,
            input_lagrangian=canonical_input,
            output_lagrangian=canonical_output,
            images=images,
            complete=complete,
            coefficients_computed=len(images) * coefficients_per_image,
            coefficients_required=coefficients_required,
            max_qudits=max_qudits,
            max_coefficients=max_coefficients,
            stop_reason=stop_reason,
        )

    if num_qudits > max_qudits:
        return result(
            ExactPrimeWitnessStatus.UNKNOWN,
            (),
            complete=False,
            stop_reason="qudit_cap",
        )

    output_labels = set(canonical_output.elements)
    images: list[ExactPrimeWeylExpansion] = []
    for label in canonical_input.basis:
        if (len(images) + 1) * coefficients_per_image > max_coefficients:
            return result(
                ExactPrimeWitnessStatus.UNKNOWN,
                tuple(images),
                complete=False,
                stop_reason="coefficient_cap",
            )
        image = _prime_conjugation_coefficients_validated(unitary, label, prime)
        images.append(image)
        if any(item.output_label not in output_labels for item in image.support):
            return result(
                ExactPrimeWitnessStatus.REJECTED,
                tuple(images),
                complete=True,
                stop_reason="support_outside_output_lagrangian",
            )

    return result(
        ExactPrimeWitnessStatus.VERIFIED,
        tuple(images),
        complete=True,
        stop_reason="witness_verified",
    )


def search_exact_prime_lagrangian_witness(
    unitary: CyclotomicMatrix,
    prime: int,
    *,
    max_qudits: int = DEFAULT_MAX_EXACT_PRIME_QUDITS,
    max_coefficients: int = DEFAULT_MAX_EXACT_PRIME_COEFFICIENTS,
    max_candidate_pairs: int = DEFAULT_MAX_EXACT_PRIME_CANDIDATE_PAIRS,
) -> ExactPrimeSearchResult:
    """Search every bounded odd-prime input/output Lagrangian pair exactly.

    A resource cap reached before a witness is independently verified returns
    ``UNKNOWN``. ``NOT_GSC`` is returned only after complete rejection of the
    canonical pair enumeration using exact cyclotomic zero tests.
    """

    for name, value in (
        ("max_qudits", max_qudits),
        ("max_coefficients", max_coefficients),
        ("max_candidate_pairs", max_candidate_pairs),
    ):
        if not isinstance(value, int) or isinstance(value, bool) or value < 0:
            raise ValueError(f"{name} must be a nonnegative integer")

    # Validate the declared field before inferring the qudit count from the
    # matrix dimension. The zero label exercises the shared odd-prime checks.
    normalize_prime_label((0, 0), prime)
    if unitary.field.order % prime:
        raise ValueError("field order must be divisible by the Weyl prime")
    rows, columns = unitary.shape
    if rows != columns:
        raise ValueError("exact search requires a square matrix")
    num_qudits = _prime_power_qudits(rows, prime)
    if not unitary.is_unitary():
        raise ValueError("exact search requires an exactly unitary matrix")

    def result(
        status: ExactPrimeSearchStatus,
        *,
        search_complete: bool,
        lagrangians_total: int | None,
        input_lagrangians_checked: int,
        candidate_pairs_checked: int,
        coefficients_computed: int,
        witness: ExactPrimeWitnessVerification | None = None,
        stop_reason: str,
        message: str,
    ) -> ExactPrimeSearchResult:
        return ExactPrimeSearchResult(
            schema_version=EXACT_PRIME_SEARCH_SCHEMA,
            status=status,
            field_order=unitary.field.order,
            arithmetic="cyclotomic_exact_odd_prime",
            prime=prime,
            num_qudits=num_qudits,
            search_complete=search_complete,
            max_qudits=max_qudits,
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

    if num_qudits > max_qudits:
        return result(
            ExactPrimeSearchStatus.UNKNOWN,
            search_complete=False,
            lagrangians_total=None,
            input_lagrangians_checked=0,
            candidate_pairs_checked=0,
            coefficients_computed=0,
            stop_reason="qudit_cap",
            message=f"n={num_qudits} exceeds the odd-prime exact-search limit {max_qudits}",
        )

    lagrangians = enumerate_prime_lagrangians(num_qudits, prime)
    coefficients_per_image = prime ** (2 * num_qudits)
    verification_work = num_qudits * coefficients_per_image
    coefficients_computed = 0
    candidate_pairs_checked = 0
    input_lagrangians_checked = 0

    for input_lagrangian in lagrangians:
        images: list[ExactPrimeWeylExpansion] = []
        for label in input_lagrangian.basis:
            if coefficients_computed + coefficients_per_image > max_coefficients:
                return result(
                    ExactPrimeSearchStatus.UNKNOWN,
                    search_complete=False,
                    lagrangians_total=len(lagrangians),
                    input_lagrangians_checked=input_lagrangians_checked,
                    candidate_pairs_checked=candidate_pairs_checked,
                    coefficients_computed=coefficients_computed,
                    stop_reason="coefficient_cap",
                    message="exact Weyl-coefficient work cap interrupted the search",
                )
            images.append(
                _prime_conjugation_coefficients_validated(unitary, label, prime)
            )
            coefficients_computed += coefficients_per_image

        support = {
            coefficient.output_label
            for image in images
            for coefficient in image.support
        }
        for output_lagrangian in lagrangians:
            if candidate_pairs_checked >= max_candidate_pairs:
                return result(
                    ExactPrimeSearchStatus.UNKNOWN,
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

            if coefficients_computed + verification_work > max_coefficients:
                return result(
                    ExactPrimeSearchStatus.UNKNOWN,
                    search_complete=False,
                    lagrangians_total=len(lagrangians),
                    input_lagrangians_checked=input_lagrangians_checked,
                    candidate_pairs_checked=candidate_pairs_checked,
                    coefficients_computed=coefficients_computed,
                    stop_reason="coefficient_cap",
                    message="coefficient cap cannot fund independent witness verification",
                )
            witness = verify_exact_prime_lagrangian_witness(
                unitary,
                input_lagrangian,
                output_lagrangian,
                max_qudits=max_qudits,
                max_coefficients=verification_work,
            )
            coefficients_computed += witness.coefficients_computed
            if witness.status is not ExactPrimeWitnessStatus.VERIFIED:
                raise AssertionError(
                    "independent odd-prime witness verification disagreed with search"
                )
            return result(
                ExactPrimeSearchStatus.GSC,
                search_complete=False,
                lagrangians_total=len(lagrangians),
                input_lagrangians_checked=input_lagrangians_checked + 1,
                candidate_pairs_checked=candidate_pairs_checked,
                coefficients_computed=coefficients_computed,
                witness=witness,
                stop_reason="witness_found",
                message="found and independently verified an exact odd-prime witness",
            )
        input_lagrangians_checked += 1

    return result(
        ExactPrimeSearchStatus.NOT_GSC,
        search_complete=True,
        lagrangians_total=len(lagrangians),
        input_lagrangians_checked=input_lagrangians_checked,
        candidate_pairs_checked=candidate_pairs_checked,
        coefficients_computed=coefficients_computed,
        stop_reason="exhausted",
        message="no odd-prime Lagrangian pair passed the exhaustive exact search",
    )
