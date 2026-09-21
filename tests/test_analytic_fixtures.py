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
from generalized_semi_clifford.analytic_fixtures import (
    three_qubit_cyclic_permutation_fixture,
)
from generalized_semi_clifford.lagrangian import enumerate_lagrangians
from generalized_semi_clifford.pauli import all_pauli_labels


def _permutation_pauli_image(
    permutation: tuple[int, ...],
    label: tuple[int, ...],
) -> tuple[int, ...] | None:
    """Conjugate a phase-free Pauli by a basis permutation over exact bits."""

    num_qubits = len(label) // 2
    dimension = 1 << num_qubits
    inverse = [0] * dimension
    for source, target in enumerate(permutation):
        inverse[target] = source

    def mask(bits: tuple[int, ...]) -> int:
        return sum(bit << (num_qubits - 1 - index) for index, bit in enumerate(bits))

    input_x = mask(label[:num_qubits])
    input_z = mask(label[num_qubits:])
    displacements = []
    signs = []
    for target in range(dimension):
        source = inverse[target]
        displacements.append(permutation[source ^ input_x] ^ target)
        signs.append((input_z & source).bit_count() % 2)
    if len(set(displacements)) != 1:
        return None

    output_x = displacements[0]
    for output_z in range(dimension):
        phases = {
            sign ^ ((output_z & target).bit_count() % 2)
            for target, sign in enumerate(signs)
        }
        if len(phases) == 1:
            x_bits = tuple(
                (output_x >> (num_qubits - 1 - index)) & 1
                for index in range(num_qubits)
            )
            z_bits = tuple(
                (output_z >> (num_qubits - 1 - index)) & 1
                for index in range(num_qubits)
            )
            return x_bits + z_bits
    return None


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
        "three-qubit cyclic permutation",
    )
    assert all(not fixture.unitary.flags.writeable for fixture in fixtures)


def test_cyclic_permutation_is_exactly_gsc_but_not_semi_clifford() -> None:
    fixture = three_qubit_cyclic_permutation_fixture()
    assert fixture.basis_permutation is not None
    assert fixture.witness is not None

    recognized = {
        label
        for label in all_pauli_labels(3)
        if any(label)
        and _permutation_pauli_image(fixture.basis_permutation, label) is not None
    }

    assert recognized == {(0, 0, 0, 1, 0, 0)}  # Z on the leftmost qubit.
    assert not any(
        all(label in recognized for label in lagrangian.elements if any(label))
        for lagrangian in enumerate_lagrangians(3)
    )
    assert fixture.expected_status is GSCStatus.GSC
    assert not fixture.expected_semi_clifford
    assert verify_gsc_witness(fixture.unitary, fixture.witness)
    assert check_gsc_naive(fixture.unitary).status is GSCStatus.GSC
    assert all(
        not verify_gsc_witness(fixture.unitary, witness)
        for witness in fixture.rejected_witnesses
    )


def test_cyclic_permutation_rejects_wrong_semi_clifford_classification() -> None:
    fixture = three_qubit_cyclic_permutation_fixture()

    assert fixture.basis_permutation == (0, 1, 2, 3, 4, 6, 7, 5)
    assert fixture.expected_semi_clifford is not True


@pytest.mark.parametrize("angle", [math.inf, -math.inf, math.nan])
def test_fixture_families_reject_nonfinite_angles(angle: float) -> None:
    with pytest.raises(ValueError, match="finite"):
        one_qubit_phase_fixture(angle)
    with pytest.raises(ValueError, match="finite"):
        two_qubit_controlled_phase_fixture(angle)
