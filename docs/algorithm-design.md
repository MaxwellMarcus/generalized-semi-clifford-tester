# Algorithm design and correctness contract

## Goal

The public naïve API accepts a dense qubit unitary and returns one of three
outcomes:

- `GSC`, together with a certificate that can be verified independently;
- `NOT_GSC`, together with a mathematical obstruction or finite certificate;
- `UNKNOWN`, when the configured exhaustive-search limit is exceeded.

The third outcome is important. An unperformed or incomplete search is never
reported as a negative decision. Current positive and negative decisions are
with respect to a stated complex128 coefficient tolerance, not exact symbolic
arithmetic.

## Naïve algorithm

For each binary Lagrangian (L\subset\mathbb F_2^{2n}), the Paulis labeled by
(L) form a basis of the Pauli MASA \(\mathcal A_L\). The implementation:

1. enumerates every input Lagrangian (L);
2. computes (UP_\ell U^\dagger) for every \(\ell\in L\);
3. expands each image in the orthonormal Pauli basis;
4. enumerates every output Lagrangian (S); and
5. accepts when every coefficient outside (S) is at most the configured
   tolerance.

The inclusion (U\mathcal A_LU^\dagger\subseteq\mathcal A_S) is then equality
because both algebras have dimension (2^n). The returned witness stores (L),
(S), and the measured maximum leakage.

The number of Lagrangians is \(\prod_{j=1}^n(2^j+1)\): 3, 15, 135, and 2,295
for one through four qubits. Because this implementation considers all ordered
input/output pairs and uses dense matrices, it defaults to at most three qubits.

## Correctness requirements before a tester is released

1. State the exact definition of generalized semi-Clifford used by the package.
2. Explain the equivalence between Pauli-support containment and equality of the
   two Pauli MASAs.
3. Expose arithmetic, tolerance, and failure modes in the result.
4. Make positive witnesses independently checkable.
5. Include analytically derived fixtures, edge cases, and known counterexamples.

## Proposed package boundaries

- `symplectic.py`: exact binary linear algebra and symplectic validation.
- `pauli.py`: phase-aware Pauli representations and conjugation actions.
- `lagrangian.py`: isotropic/Lagrangian subspace enumeration and validation.
- `tester.py`: result types, exhaustive orchestration, and witness verification.

## Performance plan

Start with transparent exhaustive algorithms at small qubit counts. Add optimized
backends only after regression fixtures exist. Benchmarks should report qubit
count, representation size, search limits, peak memory, and whether a result was
exact or heuristic.
