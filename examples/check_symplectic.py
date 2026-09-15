"""Validate familiar one-qubit Clifford actions in binary coordinates."""

from generalized_semi_clifford import is_symplectic


def main() -> None:
    identity_action = ((1, 0), (0, 1))
    hadamard_action = ((0, 1), (1, 0))
    invalid_action = ((1, 0), (1, 0))

    print("identity:", is_symplectic(identity_action))
    print("Hadamard:", is_symplectic(hadamard_action))
    print("singular map:", is_symplectic(invalid_action))


if __name__ == "__main__":
    main()
