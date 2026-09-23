"""Standalone SC theorem demo; not a GSC membership tester."""

from qiskit import QuantumCircuit
from qiskit.circuit.library import UnitaryGate

from generalized_semi_clifford.analytic_fixtures import three_qubit_cyclic_permutation_fixture
from generalized_semi_clifford.sc_testing import run_choi_sc_test


def main():
    t_gate = QuantumCircuit(1)
    t_gate.t(0)
    result = run_choi_sc_test(
        t_gate,
        hierarchy_level=3,  # This particular gate is independently known to be in C3.
        assume_hierarchy_membership=True,
        assume_ideal_sampling=True,
        seed=12,
        batch_shots=17,
    )
    print(f"T gate: {result.decision.value}")
    print(f"Forward queries: {result.unitary_queries}; inverse queries: {result.inverse_queries}")
    print(f"Conditional false-acceptance bound: {result.conditional_false_accept_bound}")

    permutation = three_qubit_cyclic_permutation_fixture()
    negative = run_choi_sc_test(
        UnitaryGate(permutation.unitary),
        hierarchy_level=3,  # Planning parameter only: no hierarchy promise asserted here.
        assume_ideal_sampling=True,
        difference_samples=128,
        seed=31,
    )
    print(f"GSC-but-not-SC permutation: {negative.decision.value}")
    print(negative.message)
    print("Local floating-point simulation; all theorem bounds are conditional on ideal premises.")


if __name__ == "__main__":
    main()
