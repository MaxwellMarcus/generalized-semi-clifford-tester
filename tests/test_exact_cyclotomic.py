import json
from fractions import Fraction

import pytest

sympy = pytest.importorskip("sympy")

from generalized_semi_clifford.exact_cyclotomic import (  # noqa: E402
    EXACT_SEARCH_SCHEMA,
    CyclotomicField,
    CyclotomicMatrix,
    ExactSearchStatus,
    ExactWitnessVerification,
    ccz,
    cnot,
    common_cyclotomic_field,
    controlled_phase,
    cz,
    hadamard,
    pauli_conjugation_coefficients,
    pauli_x,
    pauli_y,
    pauli_z,
    permutation_gate,
    phase_s,
    search_exact_lagrangian_witness,
    swap,
    t_gate,
    toffoli,
    verify_exact_lagrangian_witness,
)
from generalized_semi_clifford.lagrangian import lagrangian_from_basis  # noqa: E402


def test_field_rejects_implicit_approximate_inputs() -> None:
    field = CyclotomicField(8)

    with pytest.raises(TypeError, match="floating-point"):
        field.scalar(0.5)
    with pytest.raises(TypeError, match="floating-point"):
        field.scalar(1j)
    with pytest.raises(TypeError, match="floating-point"):
        field.scalar(sympy.Float("0.5"))
    with pytest.raises(ValueError, match="multiple of eight"):
        CyclotomicField(12)

    assert field.scalar(Fraction(1, 2)) + field.scalar(Fraction(1, 2)) == field.one


def test_common_field_embedding_is_exact_and_declared() -> None:
    field8 = CyclotomicField(8)
    field24 = CyclotomicField(24)
    common = common_cyclotomic_field(field8, field24)

    assert common.order == 24
    assert field8.zeta.embed(common) == common.root_power(3)
    assert t_gate(field8).embed(common) == t_gate(common)


@pytest.mark.parametrize(
    "constructor",
    [
        pauli_x,
        pauli_y,
        pauli_z,
        hadamard,
        phase_s,
        t_gate,
        cnot,
        cz,
        swap,
        toffoli,
        ccz,
    ],
)
def test_fixed_gate_constructors_are_exactly_unitary(constructor) -> None:
    matrix = constructor(CyclotomicField(8))

    assert matrix.is_unitary()
    assert matrix.adjoint() @ matrix == CyclotomicMatrix.identity(matrix.field, matrix.shape[0])


def test_controlled_root_phase_and_tensor_products() -> None:
    field = CyclotomicField(16)
    gate = controlled_phase(field, 3, 1)
    product = hadamard(field).tensor(t_gate(field))

    assert gate.rows[-1][-1] == field.zeta
    assert gate.is_unitary()
    assert product.shape == (4, 4)
    assert product.is_unitary()


def test_permutation_validation_and_deliberately_nonunitary_matrix() -> None:
    field = CyclotomicField(8)
    nonunitary = CyclotomicMatrix.from_rows(field, ((1, 1), (0, 1)))

    assert not nonunitary.is_unitary()
    with pytest.raises(ValueError, match="every index"):
        permutation_gate(field, (0, 0))
    with pytest.raises(ValueError, match="equal length"):
        CyclotomicMatrix.from_rows(field, ((1, 0), (1,)))


def test_exact_pauli_coefficients_retain_zeros_and_cyclotomic_support() -> None:
    field = CyclotomicField(8)

    hadamard_expansion = pauli_conjugation_coefficients(hadamard(field), (0, 1))
    t_expansion = pauli_conjugation_coefficients(t_gate(field), (1, 0))

    assert len(hadamard_expansion.coefficients) == 4
    assert tuple(item.output_label for item in hadamard_expansion.support) == ((1, 0),)
    assert hadamard_expansion.support[0].value == field.one
    assert tuple(item.output_label for item in t_expansion.support) == ((1, 0), (1, 1))
    assert all(item.value == field.scalar(1 / sympy.sqrt(2)) for item in t_expansion.support)


def test_supplied_exact_witness_is_revalidated_without_an_exhaustive_claim() -> None:
    field = CyclotomicField(8)
    z_lagrangian = lagrangian_from_basis(((0, 1),), 1)
    x_lagrangian = lagrangian_from_basis(((1, 0),), 1)

    accepted = verify_exact_lagrangian_witness(hadamard(field), z_lagrangian, x_lagrangian)
    rejected = verify_exact_lagrangian_witness(hadamard(field), z_lagrangian, z_lagrangian)

    assert isinstance(accepted, ExactWitnessVerification)
    assert accepted.verified
    assert accepted.field_order == 8
    assert accepted.arithmetic == "cyclotomic_exact"
    assert accepted.coefficients_computed == 4
    assert not hasattr(accepted, "tolerance")
    assert not rejected.verified


def test_exact_witness_rejects_nonunitaries_and_untrusted_spans() -> None:
    field = CyclotomicField(8)
    z_lagrangian = lagrangian_from_basis(((0, 1),), 1)
    malformed = type(z_lagrangian)(1, z_lagrangian.basis, ((0, 0),))

    with pytest.raises(ValueError, match="exactly unitary"):
        verify_exact_lagrangian_witness(
            CyclotomicMatrix.from_rows(field, ((1, 1), (0, 1))),
            z_lagrangian,
            z_lagrangian,
        )
    with pytest.raises(ValueError, match="canonical span"):
        verify_exact_lagrangian_witness(hadamard(field), malformed, z_lagrangian)


def test_bounded_exact_search_finds_and_independently_verifies_witness() -> None:
    field = CyclotomicField(8)

    result = search_exact_lagrangian_witness(hadamard(field))

    assert result.status is ExactSearchStatus.GSC
    assert result.is_gsc is True
    assert not result.search_complete
    assert result.stop_reason == "witness_found"
    assert result.witness is not None and result.witness.verified
    assert result.coefficients_computed == 8
    assert result.schema_version == EXACT_SEARCH_SCHEMA
    payload = result.to_dict()
    assert payload["status"] == "gsc"
    assert payload["witness"]["verified"] is True
    json.dumps(payload)


def test_completed_exact_search_can_prove_no_pauli_masa_pair_exists() -> None:
    field = CyclotomicField(8)
    non_gsc = t_gate(field) @ hadamard(field) @ t_gate(field)

    result = search_exact_lagrangian_witness(non_gsc)

    assert result.status is ExactSearchStatus.NOT_GSC
    assert result.is_gsc is False
    assert result.search_complete
    assert result.stop_reason == "exhausted"
    assert result.input_lagrangians_checked == 3
    assert result.candidate_pairs_checked == 9
    assert result.coefficients_computed == 12
    assert result.witness is None


@pytest.mark.parametrize(
    ("unitary_factory", "limits", "stop_reason"),
    [
        (lambda field: hadamard(field), {"max_coefficients": 3}, "coefficient_cap"),
        (
            lambda field: t_gate(field) @ hadamard(field) @ t_gate(field),
            {"max_candidate_pairs": 2},
            "candidate_pair_cap",
        ),
        (lambda field: cnot(field), {"max_qubits": 1}, "qubit_cap"),
    ],
)
def test_exact_search_caps_preserve_unknown(unitary_factory, limits, stop_reason) -> None:
    result = search_exact_lagrangian_witness(unitary_factory(CyclotomicField(8)), **limits)

    assert result.status is ExactSearchStatus.UNKNOWN
    assert result.is_gsc is None
    assert not result.search_complete
    assert result.stop_reason == stop_reason
    assert result.witness is None


def test_exact_search_rejects_invalid_caps_and_nonunitaries() -> None:
    field = CyclotomicField(8)

    with pytest.raises(ValueError, match="max_coefficients"):
        search_exact_lagrangian_witness(hadamard(field), max_coefficients=-1)
    with pytest.raises(ValueError, match="exactly unitary"):
        search_exact_lagrangian_witness(
            CyclotomicMatrix.from_rows(field, ((1, 1), (0, 1)))
        )
