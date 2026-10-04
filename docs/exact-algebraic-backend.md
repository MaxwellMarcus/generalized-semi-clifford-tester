# Exact algebraic backend: supported domain and API contract

This note fixes the scope of the exact backend. The representation, supplied-
witness verification, and bounded exhaustive-search milestones are implemented
in `generalized_semi_clifford.exact_cyclotomic`. The current dense tester
continues to use `complex128`, a positive tolerance, and the claim boundaries
documented in `algorithm-design.md`.

## Decision

The first backend will represent every matrix entry in one declared
cyclotomic number field

\[
K_m = \mathbb{Q}(\zeta_m), \qquad \zeta_m = e^{2\pi i/m},
\]

with `m` a positive multiple of eight. Elements must be stored canonically as
number-field elements, not as general symbolic expressions. Complex
conjugation is the exact automorphism `zeta_m -> zeta_m**-1`; equality and
zero tests are exact field operations.

Eight is included in the field order so that `sqrt(2)`, `i`, and the phase of
the T gate have exact representations. If a circuit also uses an `r`th root
of unity, its entries can be embedded in `K_lcm(8, r)`. The field order is
input provenance and must be retained in results and serialized witnesses.
The implementation must not infer a field from floating-point samples.

This domain is deliberately narrower than "all algebraic matrices." It is
closed under the additions, products, inverses of nonzero scalars, Kronecker
products, and adjoints used by the tester, while admitting a simple and
auditable input contract.

## Exactly supported inputs

The initial exact constructors should cover the following finite gate
vocabulary:

- Pauli, Hadamard, phase, CNOT, CZ, SWAP, and other explicitly supplied
  Clifford circuits;
- T/T-dagger and therefore arbitrary Clifford+T circuits;
- Toffoli, CCZ, and other permutation or signed-permutation gates;
- diagonal controlled phases whose phase is an explicitly supplied power of
  `zeta_m`;
- tensor products and products of the gates above; and
- an explicit dense matrix whose entries are already elements of the same
  declared `K_m` and which passes exact unitarity validation.

The contract excludes arbitrary Python complex values, NumPy arrays, decimal
approximations, unconstrained symbolic parameters, noisy channels, and gates
defined by a generic real rotation angle. Such an input may still be handled
by the numerical API, but it must never be silently promoted to an exact
claim. An algebraic gate outside a cyclotomic field is also outside version
one even if a more general number-field backend could represent it later.

OpenQASM or Qiskit interoperability should parse only a documented gate
vocabulary with exact parameters such as rational multiples of `pi`. It must
reject or return an unsupported-input result for approximate parameters;
recovering rational angles from floats is not part of the exact interface.

## Representation candidates

| Candidate | Exact equality | Deployment cost | Decision |
| --- | --- | --- | --- |
| SymPy `AlgebraicField`/`cyclotomic_field` domain elements | Canonical within a declared field | Pure-Python package with a familiar optional dependency; adequate for the existing low-qubit cap | Use for the first reference backend |
| General SymPy expressions plus `simplify` | Equality can depend on expression normalization and heuristic simplification | Easy to prototype but difficult to audit | Do not use as the storage or equality layer |
| FLINT/Calcium exact algebraic numbers | Exact algebraic-number predicates and faster compiled arithmetic | Compiled dependency and a broader representation/API surface | Revisit after reference-backend profiling |
| Custom coefficient vectors modulo a cyclotomic polynomial | Fully controlled canonical form | Must implement field embeddings, inversion, conjugation, and substantial validation | Avoid until measurements justify it |
| SageMath or PARI bindings | Mature exact number-field arithmetic | Too heavy for the package's normal optional-install path | Keep as an independent verifier option, not the first backend |

SymPy exposes cyclotomic fields through its polynomial-domain system and
documents exact arithmetic on `AlgebraicField` elements. The implementation
should use those domain elements directly rather than repeatedly converting
through expression trees. A later optimized backend may implement the same
protocol, but cross-backend conversion must preserve the declared field and
must be tested before its results are called exact.

## Exact tester boundary

For a validated exact unitary `U`, the backend can compute every coefficient

\[
2^{-n}\operatorname{tr}(Q U P U^\dagger)
\]

in `K_m`. A coefficient is in or out of an output Pauli algebra according to
exact zero, with no numerical tolerance. The existing binary Lagrangian
enumeration is already exact, so a completed exhaustive search has the
following meaning:

- a positive result contains input and output Lagrangians that an independent
  exact pass verifies;
- a negative result proves that no Pauli-MASA pair passes the defining
  inclusion test within the enumerated qubit system; and
- a qubit, time, memory, or operation cap produces `UNKNOWN`, never a negative
  result.

This is an exact decision for the package's Pauli-MASA definition at the
completed size. It is not a polynomial-time algorithm, a result for noisy
hardware, or a theorem about a matrix supplied only approximately.

Exact and numerical results must be different data types. The exact result
should record at least `field_order`, `arithmetic`, search completeness, and
work counts. It should have no coefficient tolerance and should not encode
exactness as `tolerance=0`. The numerical `NaiveGSCResult` and CLI schema keep
their present meaning. Any explicit conversion from an exact matrix to NumPy
starts a numerical computation and must say so in its result provenance.

## Validation and resource policy

Before classification, an exact matrix must satisfy `U.adjoint() * U == I`
entry by entry. Constructors must test their dimensions and field embeddings.
The first implementation keeps the current default maximum of three qubits;
raising the cap is explicit because dense exact matrices and exhaustive
Lagrangian pairs remain exponential.

Required regression fixtures are `I`, `H`, `S`, `T`, CNOT, SWAP, Toffoli,
CCZ, a controlled root-of-unity phase, and one deliberately nonunitary exact
matrix. Low-qubit exact results must be cross-checked against the numerical
tester away from tolerance boundaries, while the exact assertions themselves
must inspect field zeros rather than floats.

## Implementation milestones

1. **Complete:** add an optional `exact` dependency and a small cyclotomic scalar/matrix
   layer with fixed-gate constructors, common-field embedding, adjoint, and
   exact unitarity tests.
2. **Complete:** add exact Pauli-conjugation coefficients and independent
   verification of a supplied Lagrangian witness.
3. **Complete:** add the bounded exhaustive exact search and a separate
   versioned JSON schema that preserves field order, completeness, and
   `UNKNOWN` states.
4. **Complete:** benchmark one-, two-, and three-qubit fixtures before
   considering FLINT or a custom representation.

Each milestone is independently testable. In particular, completing the
representation layer alone must not be advertised as an exact GSC tester.

## Implemented representation API

Install the optional dependency with
`pip install generalized-semi-clifford-tester[exact]`, then import from
`generalized_semi_clifford.exact_cyclotomic`. `CyclotomicField(m)` requires an
explicit order divisible by eight. `CyclotomicScalar` stores a SymPy
`AlgebraicField` domain element, and `CyclotomicMatrix` implements common-field
matrix products, Kronecker products, adjoints, and entrywise exact unitarity
tests.

The fixed constructors cover `X`, `Y`, `Z`, `H`, `S`, `T`, `T`-dagger, CNOT,
CZ, SWAP, Toffoli, CCZ, arbitrary basis permutations, and an all-controls
root-of-unity phase. Cross-field operations embed into the declared
`Q(zeta_lcm)` field. Python floats, complex values, malformed permutations,
and non-cyclotomic symbolic values are rejected rather than guessed.

`pauli_conjugation_coefficients` returns all `4**n` coefficients, including
exact zeros, in deterministic binary-label order. `verify_exact_lagrangian_witness`
reconstructs both supplied Lagrangians from their bases, validates the unitary
entry by entry, and tests whether the exact support of every input-basis image
is contained in the supplied output algebra. Its `ExactWitnessVerification`
records the field order, arithmetic provenance, and coefficient work count,
but deliberately has no tolerance or exhaustive-search status. Bases are
canonicalized before use, while the stored element tuple must still equal the
reconstructed span; this accepts valid enumerated bases without trusting
caller-supplied elements. A rejected pair means only that this pair is not a
witness; it is not a negative GSC result.

`search_exact_lagrangian_witness` enumerates deterministic input/output
Lagrangian pairs under explicit qubit, coefficient, and candidate-pair caps.
It returns the separate `ExactSearchResult` type and the stable
`generalized-semi-clifford.exact-search.v1` dictionary schema. A cap produces
`UNKNOWN`; `NOT_GSC` is possible only after every pair is rejected. Positive
results contain an independently recomputed `ExactWitnessVerification` and
record both search and verification coefficient work. Run
`python examples/exact_search.py` for positive, negative, and capped results.

`run_exact_search_benchmark` keeps exact coefficient, candidate-pair, and
Lagrangian counts separate from host wall time and Python-managed peak memory.
The checked-in `benchmarks/exact-search-small-qubits.json` records positive
searches through three qubits, a completed one-qubit negative, and a capped
three-qubit `UNKNOWN`. These measurements profile the SymPy reference backend
on one host; they are not portable complexity bounds. Reproduce them with
`python examples/benchmark_exact_search.py`.
