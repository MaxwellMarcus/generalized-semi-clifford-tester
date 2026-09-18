"""Emit a reproducible one-qubit GSC benchmark record as JSON."""

from qiskit import QuantumCircuit

from generalized_semi_clifford import run_gsc_sampling_benchmark

benchmark = run_gsc_sampling_benchmark(
    QuantumCircuit(1),
    shots=256,
    leakage_threshold=0.01,
    seed=7,
)
print(benchmark.record.to_json())
