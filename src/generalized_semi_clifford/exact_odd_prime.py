"""Exact fixed-witness verification for bounded odd-prime Weyl systems.

This module deliberately verifies only one caller-supplied input/output
Lagrangian pair. It never turns rejection of that pair into an exhaustive
generalized semi-Clifford classification.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from enum import Enum
from itertools import product

from .exact_cyclotomic import (
    CyclotomicField,
    CyclotomicMatrix,
    CyclotomicScalar,
)
from .odd_prime import (
    PrimeLagrangian,
    PrimeVector,
    normalize_prime_label,
    prime_lagrangian_from_basis,
)

DEFAULT_MAX_EXACT_PRIME_QUDITS = 2
DEFAULT_MAX_EXACT_PRIME_COEFFICIENTS = 100_000


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
