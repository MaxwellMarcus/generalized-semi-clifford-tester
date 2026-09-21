# Analytic fixtures

The deterministic fixture catalog separates exact mathematical claims from
the dense tester's floating-point results. The one- and two-qubit positive
families have explicit Pauli-MASA witnesses, while the equal-axis one-qubit
rotation has an exact Bloch-sphere obstruction.

## Minimal GSC-but-not-semi-Clifford permutation

The three-qubit fixture is the computational-basis permutation

\[
  5\mapsto 6,\qquad 6\mapsto 7,\qquad 7\mapsto 5,
\]

with indices 0 through 4 fixed. As a basis permutation it maps the complete
diagonal matrix algebra to itself, so the Z-generated Pauli MASA is a GSC
witness.

Its non-semi-Clifford claim is checked independently with exact bit arithmetic.
For every phase-free Pauli \(X^x Z^z\), the test evaluates the induced basis
displacement and sign on all eight bit strings. A conjugate is Pauli only when
the displacement is constant and the sign is an affine parity. Of the 63
nonidentity Paulis, only \(Z_0\) passes. A three-qubit Pauli Lagrangian needs
three independent generators (seven nonidentity elements), so no
semi-Clifford witness exists.

Three qubits are minimal for this phenomenon. Every permutation on one or two
bits is affine: for two bits, the affine group has
\(4\cdot |\mathrm{GL}(2,2)|=4\cdot6=24\) elements, exactly all permutations of
four basis states. Such permutations are Clifford. After aligning its input
and output Pauli MASAs with Cliffords, any GSC gate has monomial form \(\Pi D\),
where \(\Pi\) is a basis permutation and \(D\) is diagonal. Diagonal gates are
semi-Clifford, and multiplying on either side by Cliffords preserves the
property. Consequently every one- or two-qubit generalized semi-Clifford gate
is semi-Clifford.

The test still runs the dense GSC checker as a regression cross-check, but that
numerical result is not used as the proof of either classification.
