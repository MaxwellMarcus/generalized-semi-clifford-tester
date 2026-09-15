"""Run the exhaustive tester on a positive and a negative one-qubit example."""

import numpy as np

from generalized_semi_clifford import check_gsc_naive


def main() -> None:
    t_gate = np.diag([1, np.exp(1j * np.pi / 4)])
    identity = np.eye(2, dtype=complex)
    x = np.array([[0, 1], [1, 0]], dtype=complex)
    y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    z = np.array([[1, 0], [0, -1]], dtype=complex)
    angle = 0.37
    generic_rotation = np.cos(angle) * identity - 1j * np.sin(angle) * (x + y + z) / np.sqrt(3)

    for name, gate in [("T", t_gate), ("generic rotation", generic_rotation)]:
        result = check_gsc_naive(gate)
        print(f"{name}: {result.status.value}; best leakage={result.best_leakage}")
        if result.witness:
            print("  input basis:", result.witness.input_lagrangian.basis)
            print("  output basis:", result.witness.output_lagrangian.basis)


if __name__ == "__main__":
    main()
