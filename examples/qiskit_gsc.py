"""Verify a GSC-but-not-semi-Clifford computational-basis permutation."""

import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit.library import UnitaryGate

from generalized_semi_clifford import (
    run_gsc_witness_test,
    run_semi_clifford_sampling_test,
    zero_event_shots_required,
)

permutation = [6, 2, 3, 4, 0, 5, 1, 7]
matrix = np.zeros((8, 8), dtype=np.complex128)
matrix[permutation, np.arange(8)] = 1
unitary = QuantumCircuit(3, name="permutation")
unitary.append(UnitaryGate(matrix), unitary.qubits)

z_basis = (
    (0, 0, 0, 1, 0, 0),
    (0, 0, 0, 0, 1, 0),
    (0, 0, 0, 0, 0, 1),
)
leakage_threshold = 0.01
confidence_level = 0.95
verification_shots = zero_event_shots_required(
    leakage_threshold,
    confidence_level=confidence_level,
    simultaneous_tests=len(z_basis),
)
gsc_result = run_gsc_witness_test(
    unitary,
    z_basis,
    z_basis,
    shots=verification_shots,
    leakage_threshold=leakage_threshold,
    confidence_level=confidence_level,
    seed=8,
)
semi_clifford_result = run_semi_clifford_sampling_test(unitary, shots=512, seed=11)

print("GSC witness found:", gsc_result.has_candidate_witness)
print("semi-Clifford witness found:", semi_clifford_result.has_candidate_witness)
print("empirical GSC leakage:", gsc_result.best_empirical_leakage)
print("GSC leakage upper bound:", gsc_result.maximum_leakage_upper_bound)
print("confidence-certified witness:", gsc_result.has_confidence_certified_witness)
