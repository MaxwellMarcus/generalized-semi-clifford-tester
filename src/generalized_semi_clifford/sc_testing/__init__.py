"""Standalone promised semi-Clifford testing, separate from GSC testing.

Qiskit is loaded lazily by circuit construction and execution functions.
"""

from .core import (
    ChoiSCAnalysis,
    ChoiSCSamplePlan,
    analyze_choi_bell_differences,
    plan_choi_sc_samples,
)
from .qiskit import (
    ChoiSCResult,
    SCDecision,
    build_choi_bell_sampling_circuit,
    decode_choi_bell_outcome,
    run_choi_sc_test,
)

__all__ = [
    "ChoiSCAnalysis",
    "ChoiSCResult",
    "ChoiSCSamplePlan",
    "SCDecision",
    "analyze_choi_bell_differences",
    "build_choi_bell_sampling_circuit",
    "decode_choi_bell_outcome",
    "plan_choi_sc_samples",
    "run_choi_sc_test",
]
