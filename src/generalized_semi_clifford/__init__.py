"""Foundational tools for generalized semi-Clifford analysis."""

from .lagrangian import Lagrangian, PauliLabel, enumerate_lagrangians, lagrangian_count
from .symplectic import identity, is_symplectic, matmul, standard_form, transpose
from .tester import GSCStatus, GSCWitness, NaiveGSCResult, check_gsc_naive, verify_gsc_witness

__all__ = [
    "GSCStatus",
    "GSCWitness",
    "Lagrangian",
    "NaiveGSCResult",
    "PauliLabel",
    "check_gsc_naive",
    "enumerate_lagrangians",
    "identity",
    "is_symplectic",
    "lagrangian_count",
    "matmul",
    "standard_form",
    "transpose",
    "verify_gsc_witness",
]
__version__ = "0.2.0"
