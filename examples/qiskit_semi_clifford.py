"""Build and run the experimental Qiskit semi-Clifford tester."""

from qiskit import QuantumCircuit

from generalized_semi_clifford.semi_clifford_qiskit import run_semi_clifford_sampling_test

unitary = QuantumCircuit(2, name="U")
unitary.h(0)
unitary.cx(0, 1)
unitary.t(1)

result = run_semi_clifford_sampling_test(unitary, shots=256, seed=7)
print("candidate witness found:", result.has_candidate_witness)
if result.witness is not None:
    print("input Lagrangian basis:", result.witness.input_basis)
    print("output Lagrangian basis:", result.witness.output_basis)
    print("minimum Pauli probability:", result.witness.minimum_dominant_probability)
