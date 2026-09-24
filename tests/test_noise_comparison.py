from __future__ import annotations

import numpy as np
import pytest

from generalized_semi_clifford import (
    PauliConjugationObservation,
    compare_sparse_transfer_with_sampling,
    extract_sparse_pauli_transfer,
)


def _observation(
    input_pauli: tuple[int, ...],
    counts: tuple[tuple[str, int], ...],
) -> PauliConjugationObservation:
    shots = sum(count for _, count in counts)
    return PauliConjugationObservation(
        input_pauli=input_pauli,
        dominant_output_pauli=input_pauli,
        dominant_probability=max(count for _, count in counts) / shots,
        shots=shots,
        counts=counts,
    )


def test_comparison_separates_ideal_prediction_from_noisy_sampling_bound() -> None:
    transfer = extract_sparse_pauli_transfer(np.eye(2))
    observation = _observation((1, 0), (("01", 5), ("10", 95)))

    comparison = compare_sparse_transfer_with_sampling(
        transfer,
        (observation,),
        output_basis=((1, 0),),
        confidence_level=0.95,
    )

    assert comparison.predicted_leakage_intervals == ((0.0, 0.0),)
    assert comparison.maximum_predicted_leakage_upper_bound == 0.0
    assert comparison.maximum_empirical_leakage == pytest.approx(0.05)
    assert comparison.maximum_sampling_leakage_upper_bound > 0.05
    assert comparison.shots_per_input == (100,)


def test_prediction_upper_endpoint_includes_thresholded_mass() -> None:
    angle = 1e-3
    x = np.array([[0, 1], [1, 0]], dtype=complex)
    rotation = np.cos(angle) * np.eye(2) - 1j * np.sin(angle) * x
    transfer = extract_sparse_pauli_transfer(rotation, coefficient_tolerance=0.01)
    observation = _observation((0, 1), (("01", 100),))

    comparison = compare_sparse_transfer_with_sampling(
        transfer,
        (observation,),
        output_basis=((0, 1),),
        confidence_level=0.95,
    )

    lower, upper = comparison.predicted_leakage_intervals[0]
    assert lower == 0.0
    assert upper == pytest.approx(np.sin(2 * angle) ** 2)


def test_comparison_rejects_duplicate_or_noncommuting_input_observations() -> None:
    transfer = extract_sparse_pauli_transfer(np.eye(4))
    x0 = _observation((1, 0, 0, 0), (("1000", 10),))
    z0 = _observation((0, 0, 1, 0), (("0100", 10),))

    with pytest.raises(ValueError, match="one sample"):
        compare_sparse_transfer_with_sampling(
            transfer,
            (x0, x0),
            output_basis=((1, 0, 0, 0), (0, 1, 0, 0)),
            confidence_level=0.95,
        )
    with pytest.raises(ValueError, match="commute"):
        compare_sparse_transfer_with_sampling(
            transfer,
            (x0, z0),
            output_basis=((1, 0, 0, 0), (0, 1, 0, 0)),
            confidence_level=0.95,
        )
