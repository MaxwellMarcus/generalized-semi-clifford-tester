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
- `semi_clifford_qiskit.py`: oracle-style conjugation circuits, Bell-count
  decoding, and finite-shot Lagrangian-witness search.

## Circuit-based semi-Clifford prototype

For an exact unitary oracle, define

\[
K_U=\{p\in\mathbb F_2^{2n}:UP_pU^\dagger\text{ is a Pauli}\}.
\]

The set \(K_U\) is a binary subspace: products of two successful Pauli
conjugates are again Pauli conjugates. The unitary \(U\) is semi-Clifford if
and only if \(K_U\) contains an \(n\)-dimensional isotropic subspace. Such a
subspace is automatically Lagrangian in \(\mathbb F_2^{2n}\).

For each input Pauli \(P\), the Qiskit prototype:

1. prepares \(n\) Bell pairs, giving a maximally entangled state \(|\Phi\rangle\);
2. applies \(U^\dagger\), then \(P\), then \(U\) to one half;
3. applies the inverse Bell preparation and measures; and
4. interprets the result as a sampled output Pauli \(Q\).

The probability of output \(Q\) is

\[
\left|2^{-n}\operatorname{tr}(Q^\dagger UPU^\dagger)\right|^2,
\]

so an exact Pauli conjugate gives one deterministic Bell outcome. After
thresholding the empirical dominant probabilities, a backtracking search finds
\(n\) corresponding input/output labels that are independently linearly
independent and pairwise symplectically orthogonal.

The implementation currently checks every nonidentity input Pauli and applies
an empirical probability threshold. A full property tester still requires a
promise formulation, concentration bounds, and completeness/soundness proofs.

## Performance plan

Start with transparent exhaustive algorithms at small qubit counts. Add optimized
backends only after regression fixtures exist. Benchmarks should report qubit
count, representation size, search limits, peak memory, and whether a result was
exact or heuristic.
