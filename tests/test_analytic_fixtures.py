import math

import numpy as np
import pytest

from generalized_semi_clifford import (
    GSCStatus,
    analytic_gsc_fixtures,
    check_gsc_naive,
    one_qubit_equal_axis_rotation_fixture,
    one_qubit_phase_fixture,
    two_qubit_controlled_phase_fixture,
    verify_gsc_witness,
)


@pytest.mark.parametrize(
    "fixture",
    [one_qubit_phase_fixture(0.37), two_qubit_controlled_phase_fixture(-0.61)],
)
def test_positive_families_have_independently_verifiable_witnesses(fixture) -> None:
    assert fixture.expected_status is GSCStatus.GSC
    assert fixture.expected_semi_clifford
    assert fixture.witness is not None
    assert verify_gsc_witness(fixture.unitary, fixture.witness)
    assert check_gsc_naive(fixture.unitary).status is GSCStatus.GSC
    assert all(
        not verify_gsc_witness(fixture.unitary, witness)
        for witness in fixture.rejected_witnesses
    )


def test_equal_axis_rotation_has_the_claimed_exact_pauli_pattern() -> None:
    fixture = one_qubit_equal_axis_rotation_fixture()
    paulis = (
        np.array([[0, 1], [1, 0]], dtype=complex),
        np.array([[0, -1j], [1j, 0]], dtype=complex),
        np.diag([1, -1]).astype(complex),
    )

    for pauli in paulis:
        image = fixture.unitary @ pauli @ fixture.unitary.conj().T
        coefficients = sorted(
            abs(np.trace(output @ image) / 2) for output in paulis
        )
        assert coefficients == pytest.approx([1 / 3, 2 / 3, 2 / 3])

    assert fixture.expected_status is GSCStatus.NOT_GSC
    assert not fixture.expected_semi_clifford
    assert fixture.witness is None
    assert check_gsc_naive(fixture.unitary).status is GSCStatus.NOT_GSC


def test_fixture_catalog_is_deterministic_and_matrices_are_read_only() -> None:
    fixtures = analytic_gsc_fixtures()

    assert tuple(fixture.name for fixture in fixtures) == (
        "one-qubit phase",
        "two-qubit controlled phase",
        "one-qubit equal-axis rotation",
    )
    assert all(not fixture.unitary.flags.writeable for fixture in fixtures)


@pytest.mark.parametrize("angle", [math.inf, -math.inf, math.nan])
def test_fixture_families_reject_nonfinite_angles(angle: float) -> None:
    with pytest.raises(ValueError, match="finite"):
        one_qubit_phase_fixture(angle)
    with pytest.raises(ValueError, match="finite"):
        two_qubit_controlled_phase_fixture(angle)
