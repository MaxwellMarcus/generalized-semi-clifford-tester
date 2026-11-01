# Conjugation-group theorem and implementation indexing

This note audits the recursive sampler against Theorem 5.5 of
Marcus, Subramanian, and Dall'Agnol,
[*The Generalized Semi-Clifford Conjecture Holds at Level 4*](https://arxiv.org/html/2609.38751v1#S5),
version 1. It records only a contract match; it is not an independent review
of the paper's proof.

## The theorem's contract

The theorem assumes:

- a prime local dimension `d` and any number `n` of qudits;
- a hierarchy level `k >= 4`;
- an exact input promise `U in C_k`; and
- projective conjugation groups

  ```text
  Gamma_1(U) = U P U-dagger,
  Gamma_j(U) = <V P V-dagger : V in Gamma_(j-1)(U)>  for j > 1.
  ```

Under those hypotheses, `Gamma_(k-2)(U) <= C_2` is sufficient for `U` to be
generalized semi-Clifford. The proof also derives
`Gamma_j(U) <= C_(k-j)` for every `1 <= j <= k-2` from that deep containment.
The converse is not stated. In particular, the shallow containment
`Gamma_2(U) <= C_(k-2)` does not imply the theorem's deep premise and must not
be described as equivalent to it.

## Index map

| Paper object | Implementation object | Audit result |
| --- | --- | --- |
| `Gamma_1(U)` | `_Sampler.short_word(depth=1)` and `_Sampler.escape_word(depth=1)` return `U P U-dagger` | Exact index match |
| A generator `V P V-dagger` of `Gamma_j` | the recursive `previous` sample conjugating `random_pauli()` or a pool subproduct | Exact index match: `previous` has depth `j-1` |
| Products generating `Gamma_j` | `short_word(depth=j)` and the escape sampler's growing pool | Every returned word lies in `Gamma_j`; neither distribution is claimed uniform |
| `Gamma_j <= C_2` | `test_clifford_conjugation_containment(..., depth=j)` | Exact depth/target match when the backend and oracle are exact |
| Theorem 5.5 premise | call the containment test with `depth=k-2` | Correct caller mapping; the function intentionally returns containment evidence, not a GSC verdict |
| Shallow numerical route | `probe_gsc_conjugation_condition(..., route="shallow")` uses `depth=2`, target `k-2` | Exploratory necessary-condition check only |
| Deep numerical route | `probe_gsc_conjugation_condition(..., route="deep")` uses `depth=k-2`, target `2` | Index match, but finite short-word sampling and floating-point membership do not satisfy the theorem's exact whole-group premise |

The implementation is currently qubit-only, so its local dimension is the
prime `d=2`. That specializes the paper's dimension hypothesis correctly but
does not implement the full prime-qudit theorem.

## Remaining proof and implementation gaps

1. No public entry point verifies or requires an exact certificate that
   `U in C_k`. The integer `hierarchy_level` is metadata supplied by the caller.
2. `test_clifford_conjugation_containment` has the stated false-acceptance
   guarantee only with exact projective group operations and an exact Clifford
   oracle. The bundled dense backend supplies neither.
3. The numerical deep route uses the short-word sampler, not the recursive
   finite-subgroup escape test, and therefore has no containment
   false-acceptance bound.
4. The recursive exact test bounds abstract group operations, while expanded
   calls to `U` and `U-dagger` can be exponential. The paper's sufficient
   condition does not remove that algorithmic cost gap.
5. Passing the shallow route is not a positive GSC certificate. A shallow
   violation can still provide conditional evidence against semi-Cliffordness
   within `C_k`, which is a different claim.

Consequently, the theorem justifies a GSC conclusion only for an exact
`depth=k-2` containment result combined with an independently established
`U in C_k` promise. Existing numerical results and incomplete runs remain
exploratory or `INCONCLUSIVE`.
