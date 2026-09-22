# Conjugation sampling and a first GSC/SC-distance tester

This implementation has two deliberately distinct components:

1. A backend-independent randomized test of `Gamma_j(U) <= C_2`. The
   finite-group argument below supplies a false-acceptance bound **when group
   operations and the Clifford oracle are exact**. It never enumerates an
   earlier conjugation group. It bounds group-operation count, not physical
   oracle-query complexity or bit complexity.
2. A three-qubit-limited numerical research frontend. It explores either
   `Gamma_2(U) <= C_(k-2)` or `Gamma_(k-2)(U) <= C_2` using short words. It can
   return a direct numerical GSC witness, a numerical SC-distance bound, or
   an inconclusive result. It does not attach the first component's probability
   theorem to the short-word sampler.

The new group-theoretic argument is documented here for review, rather than
being treated as an already externally established result.

## Definitions and the logical contract

All groups are projective: global phases are identified. Let `P` be the
projective n-qubit Pauli group, of order `4^n`, and set

```
Gamma_1(U) = U P U†
Gamma_(j+1)(U) = <g P g† : g in Gamma_j(U)>.
```

For j >= 2, Gamma_j contains P and is the smallest subgroup containing P
that is normalized by Gamma_(j-1). This is true even if Gamma_j is infinite.

The research premise supplied by the author is the sufficient condition
`Gamma_(k-2)(U) <= C_2 => U is GSC`, with the stated equivalent shallow
condition `Gamma_2(U) <= C_(k-2)`. Applying that premise requires the paper's
actual hypotheses and indexing; this software does not independently verify
the manuscript or the promise `U in C_k`.

The necessary condition used for SC-distance evidence is

```
U in SC intersect C_k => Gamma_j(U) <= C_(k-j),  1 <= j <= k-2.
```

Here is its justification. Write a semi-Clifford U as C_L D C_R, with D
diagonal and C_L,C_R Clifford. Clifford invariance of hierarchy membership
puts D in the diagonal level-k subgroup D_k. The diagonal level-l gates
form a group, normalized by Paulis, and the finite difference of a diagonal
level-l gate is diagonal level-(l-1). Thus P D_l is a group, is contained in
C_l, and conjugating P by an element of P D_l gives an element of
P D_(l-1). Induction gives Gamma_j(D) <= P D_(k-j). Conjugating back by C_L
gives the assertion. The diagonal hierarchy structure is described by
[Cui, Gottesman and Krishna](https://arxiv.org/abs/1608.06596).

Failure of these conditions therefore excludes SC **within C_k**; it does
not exclude GSC. Nor does it give distance from all SC: an arbitrary
diagonal phase gate is SC but need not satisfy a fixed-level condition.
There is a regression test for precisely this distinction.

## Short words: practical, replayable, no mixing claim

Draw a uniformly random Pauli P and return U P U† to sample Gamma_1 exactly
uniformly. For j >= 2, draw independent samples g_i at depth j-1 and Paulis
Q_i, and return the product of m factors g_i Q_i g_i†. Every output belongs
to Gamma_j; no enumeration or closure calculation occurs.

In particular, a depth-two word is

```
product_i [(U P_i U†) Q_i (U P_i U†)†].
```

Products are essential: testing only the defining conjugates can miss
failure of containment in a higher hierarchy level, since those levels are
not generally groups. Under the hierarchy promise, the individual defining
conjugates at depth two already belong to C_(k-2).

With length m at every recursive level, the expanded number of U/U† uses is

```
q_1 = 2,  q_j = 2 m q_(j-1) = 2 (2m)^(j-1).
```

The sampler records this cost and a replayable expression. For fixed depth
and polynomial m it uses polynomially many seed calls, but there is **no
proven inverse-polynomial detection probability** for this short-word
distribution. It is not claimed to be uniform on Gamma_j, which may be
infinite. Passing all these samples is not an acceptance certificate.

## Finite-subgroup escape without previous-group access

The following construction addresses the missing-conjugator problem for the
Clifford-target route. It is not a rapid-mixing assertion.

### Two elementary facts

**Random subproducts.** If T generates K and J is a proper subgroup of K,
a product of the elements of T, each independently included with probability
1/2 in a fixed order, lies outside J with probability at least 1/2. Pair
inclusion/exclusion of the last generator outside J; the remaining suffix
lies in J, so at most one member of each pair can lie in J.

**Bounded normalizers.** Suppose P <= K and |K| <= 2^b. Its normalizer N(K)
in the full projective unitary group is finite and satisfies

```
|N(K)| <= 4^n |K|^(2n) <= 2^(2n(b+1)).
```

Indeed, a normalizing unitary sends each of the 2n standard Pauli generators
to an element of K, giving at most |K|^(2n) possible projective image tuples.
Two unitaries with the same tuple differ by a unitary that projectively
centralizes every Pauli. Such a unitary is itself projectively a Pauli,
so each tuple has at most 4^n preimages. This proof uses the projective
normalizer, not an unbounded group of scalar phases. Counting only images
of the standard Paulis gives a sharper bound than counting all automorphisms
of K.

### Recursive sampler E(j,b)

For any **fixed finite subgroup** H of order at most 2^b, the sampler has

```
Gamma_j not contained in H => Pr[E(j,b) outside H] >= 1/4.
```

It does not need H, its membership oracle, or its generators as input.

- For j=1, return an exactly uniform element U P U†. Intersection with H is
  a proper subgroup, so escape probability is at least 1/2.
- For j>=2 and b<2n, return a uniform Pauli. H cannot contain P.
- Otherwise initialize a pool T to the 2n standard Pauli generators.
  Repeat R=32(b+1) times:
  1. Draw g = E(j-1, 2n(b+1)), with fresh randomness.
  2. Draw a random subproduct h of the current pool T.
  3. Append g h g† to T.
  Finally return an independent random subproduct of T.

All pool elements stay in Gamma_j. Until the pool escapes H, K=<T> is
contained in H, contains P, and is a proper subgroup of Gamma_j. Hence
Gamma_(j-1) cannot normalize K: otherwise K would contain the entire
normal closure defining Gamma_j. By the normalizer bound and induction,
g lies outside N(K) with probability at least 1/4. Conditional on this,
K intersect g^-1 K g is a proper subgroup of K, so h escapes it with
probability at least 1/2. Thus each round grows K with probability at least
1/8, until an element escapes H.

Each strict growth doubles |K|. There cannot be b+1 strict growth events
while K stays inside H. Coupling the growth events with independent
Bernoulli(1/8) trials gives

```
Pr[the pool stays in H] <= Pr[Binomial(32(b+1), 1/8) <= b]
                         <= exp(-9(b+1)/8) < 1/2.
```

Once the pool contains an element outside H, its final random subproduct
escapes H with probability at least 1/2. Therefore the unconditional
escape probability is at least 1/4. Independence is fresh randomness
conditional on the existing pool; the argument does not assume the pool
is independent of past draws.

### The Clifford-containment test

The projective Clifford group has order

```
|C_2| = 2^(n^2+2n) product_(a=1..n) (4^a-1) < 2^(2n^2+3n).
```

For j>=2 set b=2n^2+3n. Start with Pauli generators T. In each round draw
g=E(j-1,2n(b+1)), draw a random subproduct h of T, and test g h g† for
Clifford membership. A failed test yields a replayable containment
violation; otherwise append the word to T.

Under exact arithmetic and an exact oracle, all accepted pool elements
generate a subgroup of C_2. If Gamma_j is not contained in C_2, each round
either detects this or strictly grows that subgroup with probability at
least p=1/8 (p=1/4 for j=2 because Gamma_1 sampling is uniform). The same
doubling argument gives false acceptance at most delta after

```
R = ceil(max(2(b+1), 8 log(1/delta)) / p)
```

rounds. A depth-one test simply checks the 2n conjugated Pauli generators
and has no sampling error. Unknown oracle answers or exhausted budgets
are INCONCLUSIVE, never successful containment tests.

This is a frequentist false-acceptance bound, **not a posterior probability
that a particular accepted U is GSC**. Applying the author's sufficient
criterion at j=k-2 turns exact-oracle acceptance into a randomized GSC
acceptance statement. The general API intentionally returns a containment
result and does not silently apply a manuscript theorem.

### What is and is not polynomial

For fixed j, b starts polynomial in n and is transformed finitely many
times by b -> 2n(b+1). The number of group operations and outer Clifford
oracle calls is polynomial in n and log(1/delta), with a rapidly growing
degree as j increases. Neither Gamma_j nor the preceding groups need an
order bound.

However, repeated multiplication of previous pool words can make the
**expanded U-query count exponential**. The expression DAG saves storage
and permits reuse of already computed group elements in a matrix/group
backend; it does not save physical oracle calls. Exact representations can
also develop large bit complexity. Therefore this is NOT yet a polynomial-
query quantum algorithm, or a polynomial-bit-time exact algorithm for
arbitrary succinct circuits. An efficient compact representation or a
different short-word escape theorem is still needed. Work caps are especially
important at depth three and above. The recursion is not uniformly
polynomial when k is an input parameter.

The subgroup argument does not apply directly to target C_(k-2) when that
target is not a group. The shallow route remains short-word exploration.

### Explicit complexity of the current construction

The following bounds analyze a full run of the recursive Clifford-target
algorithm, before the prototype's operation and depth caps. They do not
describe the short-word sampler or the exhaustive numerical GSC fallback.
Let d=k-2, with k>=4, and let delta be the desired false-acceptance bound.
Define

```
b_0 = 2n^2 + 3n,
b_(i+1) = 2n(b_i+1),
R = ceil(max(2(b_0+1), 8 log(1/delta)) / p),
p = 1/4 for d=2, and 1/8 for d>=3.
```

In particular, R=O(n^2+log(1/delta)). The full outer test performs R
Clifford-membership checks, unless it terminates early with a violation.

**Abstract group operations.** Let A_s(b) count multiplications and
adjoints used by one escape sample E(s,b). For the b>=2n calls occurring
in this recursion, the pool size grows linearly with the number of rounds,
and a subproduct scans the whole current pool. Consequently,

```
A_1(b) = O(1),
A_s(b) = O(b A_(s-1)(2n(b+1)) + b^2),  s>=2,
G = O(R^2 + R A_(d-1)(b_1)).
```

For fixed k and fixed delta, expanding this recurrence gives

```
G = O_k(n^((k^2-k-4)/2)).
```

For variable delta, a sufficient bound for k>=4 is

```
G = O_k(R^2 + R n^((k^2-k-8)/2)).
```

The k=4 expression is intentionally loose in its second term, which is
already absorbed by R^2. Constants can depend on k; there is no claim of
a polynomial bound uniform in k. The d=1 case (k=3) is separate: it checks
2n conjugated generators with O(n) group operations.

**Expanded seed calls.** Let Q_s(b) bound the number of U/U-dagger
occurrences in one fully expanded escape word. It is not the size of its
shared expression DAG. At a round of E(s,b), if the pool's total seed count
is S and the preceding-level sample has at most Q seed uses, the new word
g h g-dagger has at most 2Q+S seed uses. Appending it raises the pool total
to at most 2S+2Q. With r=32(b+1) rounds this gives

```
Q_1(b) = 2,
Q_s(b) <= 2(2^(32(b+1))-1) Q_(s-1)(2n(b+1)),
Q_words <= 2(2^R-1) Q_(d-1)(b_1).
```

Here Q_words counts executing each outer candidate once. The upper bound
applies to literal word expansion, without resynthesis or cancellation.
Since b_i=O_k(n^(i+2)), it implies

```
Q_words = 2^O_k(n^(k-2) + log(1/delta)).
```

For fixed delta, the comparison is:

| k | Tested group | Abstract group-operation upper bound | Expanded seed-call upper bound |
| --- | --- | --- | --- |
| 4 | Gamma_2(U) | O(n^4) | 2^O(n^2) |
| 5 | Gamma_3(U) | O(n^8) | 2^O(n^3) |
| 6 | Gamma_4(U) | O(n^13) | 2^O(n^4) |

This is a cost of the current construction, not a lower bound on the
GSC-testing problem. A physical membership procedure that makes M calls
to each candidate and/or its adjoint can incur up to M times its expanded
seed cost; executing a candidate once is not a complete membership test.

**Total runtime is a third quantity.** If a faithful compact backend has
group-operation cost g(n) and exact membership cost c(n), the algebraic
part has cost O(G g(n) + R c(n)), in addition to representation, random
sampling, and bookkeeping costs. Polynomial-size integers suffice for
the seed-count metadata at fixed k, but that does not bound the bit size
of a general exact group-element representation. Standard dense matrix
multiplication alone costs O(8^n) arithmetic operations per product.
Thus neither the dense reference backend nor merely declaring U calls
free establishes polynomial total runtime. Re-executing the non-U gates
in a long word can also be expensive.

### Comparison with quantum property-testing literature

Actual oracle applications are a principal resource measure in unitary
property testing. Access to U-dagger, entangled inputs, and collective
measurements is stated explicitly rather than assumed free.

- [Wang, *Property testing of unitary operators*,
  Theorem 8](https://arxiv.org/html/1110.1133) gives an O(1/epsilon^2)-query Clifford
  tester with access to U and U-dagger, independent of n, for the paper's
  normalized distance. The paper also distinguishes query-efficient
  procedures from procedures known to have efficient total implementations.
- [Hinsche et al., *Clifford testing: algorithms and lower bounds*,
  Theorem 1.1](https://arxiv.org/html/2510.07164v1) gives an inverse-free
  four-query experiment that accepts Clifford gates and rejects with
  probability at least min(1/4, epsilon/2) when Clifford fidelity is at
  most 1-epsilon. Repetition is needed for high confidence. This epsilon
  is an infidelity threshold, not the same distance parameter as Wang's.
- [Low, *Learning and Testing Algorithms for the Clifford Group*,
  Theorems 6 and 8](https://arxiv.org/html/0907.2833) learns a promised
  Clifford with O(n) oracle queries and O(n^2) time. The extension learns
  a promised C_k gate with O(n^(k-1)) queries for fixed k. Promised learning
  is not a general membership tester. Learning and resynthesizing accepted
  Clifford pool elements is a possible direction for reducing word growth,
  not an optimization already established for the full recursive tester.

There is a second gap beyond expanded word length: our probability theorem
assumes exact Clifford membership, whereas these physical property testers
distinguish membership from being at least epsilon-far away. Merely finding
a non-Clifford word gives no lower bound on its distance from C_2. We have
not proved a separation bound that would let a finite-precision tester
replace every exact membership call while preserving the probability
argument. Therefore the current construction has neither a complete
physical quantum-query complexity theorem nor a polynomial-query GSC
property-testing guarantee. The next two targets are word compression and
a quantitative distance-gap argument.

## Quantitative distance evidence

Use d(A,B)=min_phi ||A-exp(i phi)B||_F/sqrt(2^n). Conjugation by A is
2-Lipschitz in A for this metric. Descending an operator W through r-1
Pauli conjugations ends at a leaf F(W). Every V in C_r has F(V) Pauli, so

```
d(W,C_r) >= d(F(W),P) / 2^(r-1).
```

For a recorded conjugation word W(U) using q occurrences of U/U†,
telescoping gives d(W(U),W(V)) <= q d(U,V). If all V in SC intersect C_k
have W(V) in the target C_r, then

```
d(U, SC intersect C_k) >= d(F(W(U)),P) / (q 2^(r-1)).
```

This bound depends on the word's **expanded** query count, not DAG size.
Long words can yield extremely weak distance bounds even when they detect
exact nonmembership.

For up to three qubits, a separate exhaustive calculation provides an
unrestricted SC-distance bound:

```
d(U, SC) >= (1/2) min_L max_(P in L) d(U P U†, P_n),
```

where L ranges over all Pauli Lagrangians. Every SC V maps all Paulis of
some L to Paulis, and the 2-Lipschitz inequality proves the bound. This
fallback is not a polynomial-in-n algorithm. A zero bound, or a bound
below the requested epsilon, is inconclusive.

Floating-point implementations subtract a supplied numerical margin before
using these formulas. That margin is NOT a rigorous interval-arithmetic
error enclosure, so all frontend outcomes are explicitly numerical.

## Using the implementation

```python
import numpy as np
from generalized_semi_clifford import (
    DenseConjugationBackend, check_gsc_sandwich, sample_conjugation_word,
    test_clifford_conjugation_containment,
)

u = np.diag([1, np.exp(1j*np.pi/16)])
sample = sample_conjugation_word(
    DenseConjugationBackend(u), depth=2, word_length=4, seed=7,
)
print(sample.expanded_seed_uses)  # 16

result = check_gsc_sandwich(u, hierarchy_level=5, seed=7)
print(result.status, result.gsc_witness)

# The backend-independent bounded test, here with a NUMERICAL backend.
backend = DenseConjugationBackend(np.eye(2))
containment = test_clifford_conjugation_containment(
    backend, backend.clifford_membership, depth=2, seed=7,
)
assert containment.rigorous_false_accept_bound is None
```

A custom exact backend implements the small `ConjugationBackend` protocol
and supplies a genuinely exact Clifford oracle, with `oracle_is_exact=True`.
The tests include exact dihedral coordinates for one-qubit phase/Pauli
groups, compare against exhaustive closures, and check replayed words.
Neither a backend nor an oracle should declare exactness merely because
its numerical tolerance is small.

The numerical sandwich frontend returns a direct Pauli-MASA witness on the
GSC side. Its distance side reports the unrestricted or restricted target
explicitly. GSC and "far from SC" need not be disjoint outcomes: GSC is
larger than SC. The implementation prefers a found GSC witness. A complete
efficient two-outcome tester with no inconclusive region is still a research
objective, not a claim of this prototype.
