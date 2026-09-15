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
- `gsc_qiskit.py`: full Pauli-support decoding, GSC discovery, and candidate
  Pauli-MASA verification.

## Circuit-based semi-Clifford prototype

For an exact unitary oracle, define

\[
K_U=\{p\in\mathbb F_2^{2n}:UP_pU^\dagger\text{ is a Pauli}\}.
\]

The set \(K_U\) is a binary subspace: products of two successful Pauli
conjugates are again Pauli conjugates. The unitary \(U\) is semi-Clifford if
and only if \(K_U\) contains an \(n\)-dimensional isotropic subspace. Such a
subspace is automatically Lagrangian in \(\mathbb F_2^{2n}\).

If \(B\) is a basis of \(K_U\), let \(G=BJB^\mathsf T\) be the restriction of
the standard symplectic form. The maximum isotropic dimension inside \(K_U\)
is

\[
\dim K_U-\frac{1}{2}\operatorname{rank}_{\mathbb F_2}(G).
\]

Thus exact data admits a polynomial-time containment test after \(K_U\) has
been identified. The finite-shot implementation instead uses backtracking on
individually accepted labels, because taking the span of a noisy accepted set
could silently introduce untested Paulis.

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
When a candidate input Lagrangian basis is already known, the witness-testing
entry point checks it using only \(n\) conjugation circuits. Exhaustive discovery
can submit circuits in configurable batches to accommodate sampler limits.

## Circuit-based generalized semi-Clifford prototype

For candidate input and output Lagrangians \(L,S\), it is enough to test a
basis \(p_1,\ldots,p_n\) of \(L\). If every Pauli coefficient of
\(UP_{p_i}U^\dagger\) is supported inside \(S\), then products of those images
show that

\[
U\mathcal A_LU^\dagger\subseteq\mathcal A_S.
\]

Both algebras have dimension \(2^n\), so the inclusion is equality. Bell
sampling gives the relevant Pauli-coefficient probabilities without forming a
dense matrix. The empirical leakage for one generator is the measured
probability mass outside \(S\); a candidate pair is accepted when the maximum
generator leakage is below the configured threshold.

Witness verification therefore needs only \(n\) circuits. Exhaustive discovery
currently enumerates the 3, 15, or 135 possible Lagrangians on each side, but
each required Pauli circuit is executed only once and reused. A future search
should construct the output isotropic subspace directly from sampled supports,
removing the quadratic Lagrangian-pair loop.

## Performance plan

Start with transparent exhaustive algorithms at small qubit counts. Add optimized
backends only after regression fixtures exist. Benchmarks should report qubit
count, representation size, search limits, peak memory, and whether a result was
exact or heuristic.
