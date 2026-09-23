# Standalone promised semi-Clifford tester

Research derivation and circuit implementation, 22 September 2026.
The proof below combines hierarchy separation with standard Bell-difference
sampling. It is supplied for review, not as a novelty claim. The new code
does not infer exact promises or noise assumptions from observed data.
It lives in `generalized_semi_clifford.sc_testing`, separate from the GSC
tester. No existing GSC entry point calls or is replaced by this algorithm.

## The result

For an n-qubit U promised to belong exactly to C_k, k>=2, there is an
ideal quantum algorithm with

```
U queries:           O(4^(k-2) (n + log(1/delta))),
U-dagger queries:    0,
other computation:  O(4^(k-2) (n^3 + n^2 log(1/delta))),
classical workspace: O(n^2) bits,
coherent workspace: 4n qubits using sequential pairs of Choi states.
```

It accepts every SC input with probability one and rejects a non-SC
promised input with probability at least 1-delta. With the
input-distance gap proved in Section 7 below, its outputs
can be interpreted, with failure probability at most delta, as

```
SC, hence GSC;  or  d(U,all SC gates) >= sqrt(2)/2^(k-1).
```

Thus it meets the previously proposed **GSC-or-far-from-SC** goal in
polynomial query complexity AND total non-oracle work for fixed k, without
sampling any Gamma group, finding an input MASA by exhaustive search, or
assuming a mixing rate. If a supplied U circuit has cost T_U, add its
query count times T_U to the total execution cost.

This is NOT a membership tester accepting all GSC gates. GSC-but-not-SC
gates take the non-SC branch, which is permitted by the relaxed output
contract but not by full GSC membership testing. The guaranteed distance
is the displayed k-dependent value; this proves an epsilon-far conclusion
for epsilon no larger than that value, not for arbitrary larger epsilon.

## 1. The Choi stabilizer group encodes the Pauli-preserving subspace

Set D=2^n and prepare the 2n-qubit state

```
|U> = (U tensor I) (sum_x |x,x>/sqrt(D)).
```

Each copy needs one U query and O(n) Clifford gates, with no inverse.
For Hermitian Paulis Q on the output and P on the reference,

```
chi_(Q,P) = <U| Q tensor P |U> = tr(U† Q U P^T)/D.
```

Let S_U be the unsigned stabilizer space, the binary labels (Q,P) for
which |chi_(Q,P)|=1. The equality case in Cauchy-Schwarz implies

```
(Q,P) in S_U  iff  U P^T U† = +/- Q.
```

Transposition changes only a Pauli sign, not its binary label. Therefore
the reference projection of S_U is exactly

```
K_U = {P : U P U† is projectively Pauli}.
```

K_U is a binary linear subspace. U is SC precisely when K_U contains an
n-dimensional isotropic subspace. If K_U has basis B and restricted
symplectic Gram matrix G, its maximum isotropic dimension is

```
dim(K_U) - rank(G)/2.
```

This gives a polynomial postprocessing test once S_U is known. Merely
counting Choi stabilizers is insufficient: their reference projection
can have a nontrivial symplectic form.

## 2. The hierarchy promise supplies a uniform expectation gap

From the hierarchy separation lemma, distinct projective A,B in C_j obey
d(A,B)>=a_j, where a_j=sqrt(2)/2^(j-1).

If (Q,P) is not in S_U, then U P^T U† and Q are distinct elements of
C_(k-1). Since d(A,B)^2=2-2|tr(A†B)/D|, we obtain

```
|chi_(Q,P)| <= 1 - lambda_k,
lambda_k = a_(k-1)^2/2 = 4^(-(k-2)).
```

Every Choi Pauli expectation is consequently either exactly +/-1 or is
bounded away from those values by lambda_k. This is the gap that the
conjugation-group words did not inherit: here only the single conjugate
U P U† appears, so it stays in C_(k-1).

## 3. Bell-difference samples reveal all stabilizer constraints

For each trial, make four independent Choi states. Bell-measure the first
two against each other to get label b_1, do the same to the second pair
to get b_2, and return z=b_1+b_2 over F_2. Each label has 4n bits.
The two pairs may be prepared and measured sequentially, needing at most
4n coherent qubits. No complex-conjugated state or U inverse is needed.

Bell-difference sampling and its stabilizer-support properties are
standard; see [Grewal, Iyer, Kretschmer and Liang, Section 2.3,
Proposition 2.17 and Fact 2.18](https://arxiv.org/html/2305.13409).
For a Pauli label a and symplectic pairing [a,z], the relevant identity is

```
E_z (-1)^[a,z] = chi_a^4,
Pr_z[[a,z]=1] = (1-chi_a^4)/2.
```

One direct derivation: the character of an ordinary Bell outcome has
expectation (-1)^(number of Y factors in a) chi_a^2. XORing two independent
outcomes squares that expression and removes the transpose sign.

For a in S_U, every sample satisfies [a,z]=0. For a outside S_U, the
expectation gap gives

```
Pr[[a,z]=1] >= (1-(1-lambda_k)^4)/2
             >= lambda_k/2
             = beta_k = 1/(2 * 4^(k-2)).
```

Take M independent difference samples, let V be their span, and compute
S_hat=V^perp. Always S_U is contained in S_hat. Any fixed a outside S_U
survives in S_hat with probability at most exp(-beta_k M). There are at
most 2^(4n) labels, so

```
Pr[S_hat != S_U] <= 2^(4n) exp(-beta_k M).
```

It suffices to choose

```
M >= ceil(2 * 4^(k-2) * (4n log(2) + log(1/delta))).
```

The prototype rounds the parenthesized factor upward before multiplying
by the integer prefactor. It uses four U queries per difference sample,
so its sufficient total query budget is 4M.

## 4. One-sided classification even before exact recovery

Project S_hat to the reference Pauli space, obtaining K_hat, and compute
its maximum isotropic dimension. Accept if it equals n; otherwise reject.

Because S_U is always contained in S_hat under ideal sampling, K_U is
always contained in K_hat. Thus if U is SC, K_hat contains an n-dimensional
isotropic space and the algorithm always accepts. Conversely, if U is
non-SC and S_hat=S_U, the algorithm rejects. Its false-acceptance
probability is at most delta under the hierarchy promise.

A negative compatibility result excludes SC even without the hierarchy
promise, assuming ideal samples. The hierarchy promise is needed for the
finite uniform sample bound and the numerical distance-from-SC conclusion.
A positive result is a bounded-error classification, not a deterministic
certificate that the sampled span is complete.

## 5. Cost and limits

Incrementally storing the sample span uses O(n^2) bits and O(M n^2) bit
operations with elementary binary elimination. Taking a symplectic
complement, projecting, and computing the Gram rank costs O(n^3).
Each sample uses O(n) additional quantum gates. These give the resource
bounds at the beginning of this note. Constants are not optimized; the
stronger expression (1-(1-lambda_k)^4)/2 could reduce the sample budget.

The important limits are:

- U must really belong to the supplied C_k; a numerical level parameter
  does not establish this. For an unpromised near-SC gate the stabilizer
  expectation gap can vanish.
- The result is polynomial for fixed k, not uniformly polynomial in k.
- Bell measurements couple two Choi copies; this is not a single-copy,
  ancilla-free, or experimentally noise-tolerant result.
- Exact zero-probability constraints drive one-sided soundness. Noise
  can corrupt the span; ordinary unmodified Gaussian elimination is not
  a robust replacement.
- This solves the relaxed SC/GSC-distance task, not full GSC discovery
  or the missing necessity direction of the Gamma containment criterion.

## 6. Implementation status

`generalized_semi_clifford.sc_testing` provides a sample planner, streaming
binary postprocessor, two-Choi-state Bell circuit, and SamplerV2 runner.
Labels use

```
(x_output, x_reference | z_output, z_reference), each block n bits.
```

They must be Bell-difference outcomes on pairs of 2n-qubit Choi states,
not the 2n-bit outputs of the earlier U P U† Pauli tester. The result's
`semi_clifford_compatible` flag reports finite-data compatibility; the
planner and ideal promise are needed for the theorem's confidence bound.
The postprocessor deliberately does not attach a hierarchy guarantee to
arbitrary input samples.

### Running it

```python
from qiskit import QuantumCircuit
from generalized_semi_clifford.sc_testing import run_choi_sc_test

u = QuantumCircuit(1)
u.t(0)  # Independently known to belong to C_3.
result = run_choi_sc_test(
    u,
    hierarchy_level=3,
    failure_probability=0.01,
    assume_hierarchy_membership=True,
    assume_ideal_sampling=True,
    seed=12,
)
print(result.decision.value, result.unitary_queries)
```

The default local StatevectorSampler is a floating-point simulation,
exponential in the 4n simulated qubits and capped at 12 qubits by default.
The polynomial runtime theorem describes a quantum implementation, not
classical statevector simulation. A supplied `sampler=` follows SamplerV2
and returns ordered `data.bell.get_bitstrings()` observations. Such a
sampler may require external transpilation/configuration for hardware;
no hardware job is automatically submitted by the default local mode.

The runner pairs two independent circuit shots, preserving pairs across
batch boundaries. It uses a persistent random Generator locally so an
integer seed is not restarted for every batch. Custom samplers must supply
fresh independent draws across jobs themselves. Counts alone are rejected:
sorting or expanding aggregated counts into adjacent pairs would not
preserve the independent sampling distribution.

`SCDecision` has only `SC_SUPPORTED`, `NON_SC`, and `INCONCLUSIVE`.
It has no GSC/non-GSC verdict. Both assumptions default to False. Without
them, the analysis is available but no positive theorem-backed decision
is issued. An acknowledged ideal-sampling negative can exclude SC without
the hierarchy promise; the distance bound requires both assumptions.
Acknowledging assumptions is not verifying them: all reported probability
and distance fields are explicitly named `conditional_*`.

The full sample budget is required for positive acceptance. Explicitly
shorter runs remain inconclusive on compatible data. Exceeding a sample
or local-simulation cap returns inconclusive before execution, never a
negative conclusion. `unitary_queries` counts actual executed calls,
including repeated circuit shots; `inverse_queries` is always zero.

Tests compute the actual two-copy Bell distribution of small Choi states,
convolve two independent outcomes, check the Fourier identity, and compare
the recovered stabilizers with direct Pauli expectations. They also check
reference-label indexing, the symplectic-rank criterion, and basis-only
processing beyond the dense cap. These are sanity checks, not a substitute
for the proof or a noise model. Integration tests exercise actual Qiskit
circuits, asymmetrical qubit ordering, odd-sized batches, reproducibility,
resource caps, malformed sampler results, and a GSC-but-not-SC gate whose
negative result must not be interpreted as non-GSC.

Run `python examples/choi_semi_clifford.py` for a local demonstration and
`python -m pytest tests/sc_testing` for the standalone test suite. Install
the optional circuit dependency with `pip install -e ".[qiskit]"`.

## 7. Self-contained hierarchy and SC distance gaps

The distance throughout is

```
d(A,B) = min_phi ||A-exp(i phi)B||_F/sqrt(2^n).
```

Distinct projective A,B in C_j have distance at least
a_j=sqrt(2)/2^(j-1). At j=1, distinct Paulis are orthogonal, proving the
base case. If some Pauli conjugates APA† and BPB† differ projectively,
induction and the triangle inequality give

```
a_(j-1) <= d(APA†,BPB†) <= 2 d(A,B).
```

If all these conjugates agree projectively, B†A projectively centralizes
every Pauli. Distinct Pauli conjugation characters force it to be a
single Pauli up to phase. It is nonidentity, so d(A,B)=sqrt(2).
This proves the separation without requiring C_j to be a group.

Now let U be in C_k and V be any SC gate, without requiring V to belong
to the hierarchy. V maps n independent commuting Paulis P_i to Paulis
Q_i. If d(U,V)<a_k, then

```
d(U P_i U†, Q_i) <= 2 d(U,V) < a_(k-1).
```

Both compared gates lie in C_(k-1); the separation forces equality up to
phase for each generator. Thus U is itself SC. Contrapositively, a non-SC
U in C_k is at least a_k from **all** SC gates. This is the distance
conclusion used by this separate SC tester; no GSC converse is assumed.
