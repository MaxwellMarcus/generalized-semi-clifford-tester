# Two-qutrit exact-search benchmark design

This note fixes the first two-qutrit fixture and its resource budgets in the
checked benchmark suite. It keeps the initial cap-semantics baseline alongside
a separately bounded, independently verified positive classification.

## Exact fixture

Use the qutrit SUM gate over the declared field
`K_24 = Q(zeta_24)`. In the computational basis ordered by
`index(a, b) = 3 * a + b`, define

```text
SUM |a, b> = |a, a + b mod 3>.
```

The corresponding column-to-row permutation is

```text
(0, 1, 2, 4, 5, 3, 8, 6, 7).
```

Construct it with `permutation_gate(CyclotomicField(24), permutation)`. Every
entry is the exact integer zero or one, so no floating-point or inferred
algebraic input enters the fixture.

SUM is a non-product two-qutrit Clifford. With the package convention
`U P U.adjoint()`, its generator action is

```text
X_1 -> X_1 X_2,    X_2 -> X_2,
Z_1 -> Z_1,        Z_2 -> Z_1^-1 Z_2.
```

This analytic action explains why the fixture is expected to have a
Lagrangian witness. The capped benchmark below deliberately does not claim to
have computed or verified that witness.

## Portable workload and checked budget

For two qutrits, canonical enumeration contains 40 Lagrangians. A complete
Weyl expansion of one basis generator contains
`3 ** (2 * 2) = 81` exact coefficients, and one input Lagrangian has two basis
generators. Independent verification of a discovered witness would require a
second 162 coefficients.

The checked record uses these explicit limits:

```text
max_qudits = 2
max_coefficients = 80
max_candidate_pairs = 1
```

Because 80 is one less than the indivisible 81-coefficient expansion, the
expected outcome is `UNKNOWN` with stop reason `coefficient_cap`, zero
coefficients computed, zero candidate pairs checked, and 40 Lagrangians
enumerated. This is an intentional pre-expansion boundary: it exercises exact
9-by-9 unitarity validation, two-qutrit enumeration, result serialization, and
honest incomplete-state handling without making a classification claim.

## Bounded verified run

The optimized exact arithmetic makes a positive run practical without
weakening the cap semantics. The second two-qutrit checked record uses:

```text
max_qudits = 2
max_coefficients = 324
max_candidate_pairs = 1
```

Exactly 162 coefficient slots fund both basis-generator expansions for the
first input Lagrangian. The first deterministic output candidate contains the
resulting support, and another 162 slots fund an independent reconstruction
and verification of that same witness. The run therefore returns `GSC` with a
verified witness after exactly 324 coefficients and one candidate pair.

These limits do not assume that every fixture has a witness in its first
candidate pair. If that prefix does not find a witness, the search still
returns `UNKNOWN` at the candidate-pair cap; it never promotes a truncated
prefix to a global negative result.

## Host-time budget

On 2026-10-13, CPython 3.12.14 on macOS 14.7.3 arm64 measured 1.82 seconds and
929,485 bytes of Python peak memory for the 80-coefficient cap. The checked
example has a practical advisory budget of 10 seconds on a development host;
wall time and peak memory remain observational fields, not pass/fail portable
invariants.

For comparison, before the monomial-support optimization, allowing one
81-coefficient expansion with `max_coefficients = 161` measured 153.00 seconds
and 2,489,436 bytes before returning `UNKNOWN` ahead of the second basis
generator. On 2026-10-15, the same CPython 3.12.14 macOS arm64 host measured
5.54 seconds and 1,437,472 bytes after replacing dense output-Weyl adjoints and
matrix products with the mathematically identical exact trace inner product on
the Weyl matrix's nine nonzero entries. The regression suite independently
cross-checks the optimized coefficient against the dense formula and verifies
the analytic single-Weyl support of SUM's first expansion.

Profiling that fixed expansion then showed that each output coefficient and
each adjoint entry reconstructed the same primitive `K_24` root through SymPy.
On 2026-10-16, caching that immutable field element reduced the same
`max_coefficients = 161` search measurement to 0.22 seconds and 203,686 bytes
of Python peak memory on the same host. Three direct 81-coefficient expansions
measured 0.035--0.038 seconds each. The result remained `UNKNOWN` with exactly
81 coefficients computed, and the independently checked support remained the
single analytic label `(0, 0, 2, 1)` with coefficient one.

The 324-coefficient positive run was measured repeatedly on 2026-10-18 on the
same host at 0.77--0.95 seconds and 273,375--596,899 bytes of Python peak
memory. This bounded positive measurement is not extrapolated into a full
search runtime. Runtime and Python peak memory remain machine-dependent
observations rather than portable pass/fail thresholds.

## Checked-suite contract

The suite and its regression test now:

- construct the exact permutation above rather than loading approximate data;
- record the explicit limits and the expected capped `UNKNOWN` outcome;
- assert the portable counts and stop reason, but only require host runtime and
  peak memory to be positive finite measurements; and
- keep the capped and verified-positive two-qutrit records alongside, but
  semantically separate from, the one-qutrit positive, completed-negative, and
  capped records.
