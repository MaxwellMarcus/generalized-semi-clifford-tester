# Fixed-witness tolerant testing contract

This document states the guarantee implemented by `run_gsc_witness_test`.
It is a finite-sample statement about one input/output Pauli-MASA pair fixed
before verification data are collected. It is not a complete property-testing
theorem for generalized semi-Clifford unitaries.

## Experiment and leakage parameter

Fix \(n\)-qubit Lagrangians \(L,S\subset\mathbb F_2^{2n}\), and fix a basis
\(p_1,\ldots,p_n\) of \(L\). For generator \(i\), Bell sampling has output
distribution

\[
 q_i(y)=\left|2^{-n}\operatorname{tr}
 \left(P_y^\dagger U P_{p_i}U^\dagger\right)\right|^2.
\]

The population leakage outside the proposed output MASA is

\[
 \lambda_i=\sum_{y\notin S}q_i(y),\qquad
 \lambda_{\max}=\max_i\lambda_i.
\]

The tested tolerant promise is therefore about the fixed pair \((L,S)\):

- accept when the data certify \(\lambda_{\max}\leq\varepsilon\);
- reject this witness otherwise.

Rejection is not a claim that no other witness exists. At positive
\(\varepsilon\), acceptance is not an exact GSC-membership claim and does not,
without a separate error-propagation result, bound every element of
\(U\mathcal A_LU^\dagger\) merely from the tested basis generators.

## Sampling hypotheses

The guarantee requires all of the following.

1. The input/output witness and threshold are fixed independently of the
   verification samples. A witness selected by discovery must be checked with
   fresh samples.
2. For each generator, decoded outside-\(S\) indicators are identically
   distributed Bernoulli trials with a stationary parameter \(\lambda_i\).
   The implementation's binomial calculation additionally assumes independent
   shots within that generator's job.
3. Counts include every verification shot and Bell outcomes are decoded using
   the same qubit convention as the proposed Lagrangians.
4. The requested familywise confidence \(1-\alpha\) is strictly between zero
   and one. The implementation allocates \(\alpha/n\) to each generator.

The \(\lambda_i\) describe the distribution actually sampled. Noise is folded
into those parameters; interpreting them as properties of an ideal unitary
requires an independently justified noise model.

## Soundness theorem

Let \(K_i\sim\operatorname{Binomial}(N_i,\lambda_i)\) be the outside-support
count and let \(u_i(K_i)\) be its one-sided Clopper--Pearson upper bound at
confidence \(1-\alpha/n\). The implementation accepts exactly when

\[
 \max_i u_i(K_i)\leq\varepsilon.
\]

Clopper--Pearson coverage gives

\[
 \Pr\{\lambda_i>u_i(K_i)\}\leq\alpha/n.
\]

The union bound therefore gives simultaneous coverage

\[
 \Pr\{\lambda_i\leq u_i(K_i)\text{ for every }i\}\geq1-\alpha.
\]

On that simultaneous-coverage event, acceptance implies
\(\lambda_{\max}\leq\varepsilon\). Consequently, for every fixed witness with
\(\lambda_{\max}>\varepsilon\), the probability of false certification is at
most \(\alpha\). This statement is exact up to the numerical evaluation of the
beta quantile; it is not a normal approximation.

## Completeness with an explicit margin

Finite data do not provide uniform high acceptance probability at the boundary
\(\lambda_{\max}=\varepsilon\). Completeness therefore needs a promised margin.
For generator \(i\), define

\[
 k_i^*=\max\{k\in\{0,\ldots,N_i\}:u_i(k)\leq\varepsilon\},
\]

with \(k_i^*=-1\) when the set is empty. Because the upper bound is monotone in
the observed count, the tester accepts if and only if every \(K_i\leq k_i^*\).
If the promise is \(\lambda_i\leq\eta_i<\varepsilon\), binomial stochastic
monotonicity and another union bound give the finite-sample completeness bound

\[
 \Pr\{\text{accept}\}\geq
 1-\sum_{i=1}^n
 \Pr\{\operatorname{Binomial}(N_i,\eta_i)>k_i^*\}.
\]

This formula is an auditable shot-planning criterion: choose the \(N_i\) so
the right-hand side meets the desired power. The existing
`zero_event_shots_required` helper is only the best-case special case ensuring
\(u_i(0)\leq\varepsilon\); it does not promise that zero events will be
observed when \(\lambda_i>0\).

## Exact endpoint

Under the ideal unitary Bell-sampling identity, \(\lambda_i=0\) for every basis
generator if and only if every \(UP_{p_i}U^\dagger\) has Pauli support inside
\(S\). Products of the basis images then give
\(U\mathcal A_LU^\dagger\subseteq\mathcal A_S\), and equal dimensions make the
inclusion equality. This algebraic implication concerns the population
endpoint \(\lambda_{\max}=0\); no finite-shot observation proves that endpoint.
