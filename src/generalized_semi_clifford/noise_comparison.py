"""Compare sparse-transfer predictions with finite-shot Bell sampling."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

from .gsc_qiskit import (
    empirical_lagrangian_leakage,
    lagrangian_leakage_upper_bound,
)
from .lagrangian import Lagrangian, PauliLabel, lagrangian_from_basis
from .semi_clifford_qiskit import PauliConjugationObservation
from .sparse_transfer import SparsePauliTransfer


@dataclass(frozen=True)
class TransferSamplingComparison:
    """Ideal numerical prediction beside finite-shot leakage evidence.

    Sparse-transfer prediction bounds include every thresholded coefficient in
    the upper endpoint.  The sampling upper bound is a simultaneous one-sided
    Clopper--Pearson bound and is valid for an independently fixed witness.
    Neither bound turns a noisy simulation into an exact algebraic claim.
    """

    input_lagrangian: Lagrangian
    output_lagrangian: Lagrangian
    predicted_leakage_intervals: tuple[tuple[float, float], ...]
    maximum_predicted_leakage_lower_bound: float
    maximum_predicted_leakage_upper_bound: float
    maximum_empirical_leakage: float
    maximum_sampling_leakage_upper_bound: float
    confidence_level: float
    shots_per_input: tuple[int, ...]


def compare_sparse_transfer_with_sampling(
    transfer: SparsePauliTransfer,
    observations: Iterable[PauliConjugationObservation],
    *,
    output_basis: Iterable[PauliLabel],
    confidence_level: float,
) -> TransferSamplingComparison:
    """Compare an ideal transfer matrix with samples for one fixed witness."""

    samples = tuple(observations)
    if not samples:
        raise ValueError("observations must not be empty")
    input_labels = tuple(sample.input_pauli for sample in samples)
    if len(set(input_labels)) != len(samples):
        raise ValueError("observations must contain one sample per input-basis generator")
    input_lagrangian = lagrangian_from_basis(input_labels, transfer.num_qubits)
    output_lagrangian = lagrangian_from_basis(output_basis, transfer.num_qubits)
    allowed = set(output_lagrangian.elements)
    label_indices = {label: index for index, label in enumerate(transfer.labels)}

    intervals: list[tuple[float, float]] = []
    empirical: list[float] = []
    shots: list[int] = []
    for sample in samples:
        if len(sample.input_pauli) != 2 * transfer.num_qubits:
            raise ValueError("observations and sparse transfer use different qubit counts")
        try:
            column_index = label_indices[sample.input_pauli]
        except KeyError as error:
            raise ValueError("observation input is absent from the sparse transfer") from error
        retained_outside_mass = sum(
            coefficient.value**2
            for coefficient in transfer.columns[column_index]
            if coefficient.output_label not in allowed
        )
        discarded_mass = transfer.discarded_l2_mass[column_index]
        intervals.append(
            (
                retained_outside_mass,
                min(1.0, retained_outside_mass + discarded_mass),
            )
        )
        empirical.append(empirical_lagrangian_leakage(sample, output_lagrangian))
        shots.append(sample.shots)

    return TransferSamplingComparison(
        input_lagrangian=input_lagrangian,
        output_lagrangian=output_lagrangian,
        predicted_leakage_intervals=tuple(intervals),
        maximum_predicted_leakage_lower_bound=max(lower for lower, _ in intervals),
        maximum_predicted_leakage_upper_bound=max(upper for _, upper in intervals),
        maximum_empirical_leakage=max(empirical),
        maximum_sampling_leakage_upper_bound=lagrangian_leakage_upper_bound(
            samples,
            output_lagrangian,
            confidence_level=confidence_level,
        ),
        confidence_level=confidence_level,
        shots_per_input=tuple(shots),
    )
