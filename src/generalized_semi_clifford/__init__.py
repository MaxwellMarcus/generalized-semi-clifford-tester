"""Foundational tools for generalized semi-Clifford analysis."""

from .gsc_qiskit import (
    CircuitGSCWitness,
    GSCSamplingResult,
    empirical_lagrangian_leakage,
    find_gsc_sampling_witness,
    output_pauli_probabilities,
    run_gsc_sampling_test,
    run_gsc_witness_test,
)
from .lagrangian import (
    Lagrangian,
    PauliLabel,
    enumerate_lagrangians,
    lagrangian_containing,
    lagrangian_count,
    lagrangian_from_basis,
)
from .semi_clifford_qiskit import (
    PauliConjugationObservation,
    SemiCliffordSamplingResult,
    SemiCliffordWitness,
    bell_bitstring_to_pauli,
    bell_counts_to_observation,
    build_pauli_conjugation_test_circuit,
    build_semi_clifford_test_circuits,
    find_lagrangian_witness,
    maximum_isotropic_dimension,
    run_semi_clifford_sampling_test,
    run_semi_clifford_witness_test,
    symplectic_pairing,
)
from .symplectic import identity, is_symplectic, matmul, standard_form, transpose
from .tester import GSCStatus, GSCWitness, NaiveGSCResult, check_gsc_naive, verify_gsc_witness

__all__ = [
    "CircuitGSCWitness",
    "GSCStatus",
    "GSCSamplingResult",
    "GSCWitness",
    "Lagrangian",
    "NaiveGSCResult",
    "PauliConjugationObservation",
    "PauliLabel",
    "SemiCliffordSamplingResult",
    "SemiCliffordWitness",
    "bell_bitstring_to_pauli",
    "bell_counts_to_observation",
    "build_pauli_conjugation_test_circuit",
    "build_semi_clifford_test_circuits",
    "check_gsc_naive",
    "empirical_lagrangian_leakage",
    "enumerate_lagrangians",
    "find_gsc_sampling_witness",
    "find_lagrangian_witness",
    "identity",
    "is_symplectic",
    "lagrangian_count",
    "lagrangian_containing",
    "lagrangian_from_basis",
    "matmul",
    "maximum_isotropic_dimension",
    "output_pauli_probabilities",
    "run_gsc_sampling_test",
    "run_gsc_witness_test",
    "run_semi_clifford_sampling_test",
    "run_semi_clifford_witness_test",
    "standard_form",
    "symplectic_pairing",
    "transpose",
    "verify_gsc_witness",
]
__version__ = "0.3.0"
