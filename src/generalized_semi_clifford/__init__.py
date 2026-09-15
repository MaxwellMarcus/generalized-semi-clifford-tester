"""Foundational tools for generalized semi-Clifford analysis."""

from .lagrangian import Lagrangian, PauliLabel, enumerate_lagrangians, lagrangian_count
from .semi_clifford_qiskit import (
    PauliConjugationObservation,
    SemiCliffordSamplingResult,
    SemiCliffordWitness,
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
    "GSCStatus",
    "GSCWitness",
    "Lagrangian",
    "NaiveGSCResult",
    "PauliConjugationObservation",
    "PauliLabel",
    "SemiCliffordSamplingResult",
    "SemiCliffordWitness",
    "bell_counts_to_observation",
    "build_pauli_conjugation_test_circuit",
    "build_semi_clifford_test_circuits",
    "check_gsc_naive",
    "enumerate_lagrangians",
    "find_lagrangian_witness",
    "identity",
    "is_symplectic",
    "lagrangian_count",
    "matmul",
    "maximum_isotropic_dimension",
    "run_semi_clifford_sampling_test",
    "run_semi_clifford_witness_test",
    "standard_form",
    "symplectic_pairing",
    "transpose",
    "verify_gsc_witness",
]
__version__ = "0.2.0"
