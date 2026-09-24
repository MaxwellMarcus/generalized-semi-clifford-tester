"""Compare ideal sparse-transfer predictions with a noisy Aer experiment."""

from qiskit import QuantumCircuit, transpile
from qiskit.primitives import StatevectorSampler
from qiskit_aer.noise import NoiseModel, depolarizing_error
from qiskit_aer.primitives import SamplerV2 as AerSampler

from generalized_semi_clifford import (
    compare_sparse_transfer_with_sampling,
    extract_sparse_pauli_transfer,
    run_gsc_witness_test,
)

SHOTS = 2_000
CONFIDENCE = 0.95
X_BASIS = ((1, 0),)


class TranspilingSampler:
    """Compile custom unitary instructions before submitting them to Aer."""

    def __init__(self, sampler, basis_gates):
        self.sampler = sampler
        self.basis_gates = basis_gates

    def run(self, circuits, *, shots):
        compiled = [
            transpile(circuit, basis_gates=self.basis_gates, optimization_level=0)
            for circuit in circuits
        ]
        return self.sampler.run(compiled, shots=shots)


def summarize(name, transfer, result):
    comparison = compare_sparse_transfer_with_sampling(
        transfer,
        result.observations,
        output_basis=X_BASIS,
        confidence_level=CONFIDENCE,
    )
    print(name)
    print("  ideal sparse-transfer leakage interval:", comparison.predicted_leakage_intervals[0])
    print("  empirical Bell-sampling leakage:", comparison.maximum_empirical_leakage)
    print(
        f"  {CONFIDENCE:.0%} one-sided leakage upper bound:",
        comparison.maximum_sampling_leakage_upper_bound,
    )


def main():
    identity = QuantumCircuit(1, name="identity")
    transfer = extract_sparse_pauli_transfer([[1, 0], [0, 1]])
    ideal = run_gsc_witness_test(
        identity,
        X_BASIS,
        X_BASIS,
        sampler=StatevectorSampler(seed=17),
        shots=SHOTS,
        confidence_level=CONFIDENCE,
    )

    noise_model = NoiseModel()
    noise_model.add_all_qubit_quantum_error(depolarizing_error(0.01, 1), ["h"])
    noise_model.add_all_qubit_quantum_error(depolarizing_error(0.03, 2), ["cx"])
    noisy_sampler = TranspilingSampler(
        AerSampler(
            seed=17,
            options={"backend_options": {"noise_model": noise_model}},
        ),
        noise_model.basis_gates,
    )
    noisy = run_gsc_witness_test(
        identity,
        X_BASIS,
        X_BASIS,
        sampler=noisy_sampler,
        shots=SHOTS,
        confidence_level=CONFIDENCE,
    )

    summarize("ideal statevector", transfer, ideal)
    summarize("Aer depolarizing noise model", transfer, noisy)
    print(
        "Finite-shot confidence bounds are evidence for this fixed witness, "
        "not exact membership."
    )


if __name__ == "__main__":
    main()
