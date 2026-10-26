# Two-qutrit exact-search benchmark design

This note fixes the first two-qutrit fixture and its resource budgets in the
checked benchmark suite. It keeps the initial cap-semantics baseline alongside
a separately bounded, independently verified positive classification and an
exactly rejected input-Lagrangian prefix for a non-Clifford fixture.

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

## Bounded rejection prefix

For a complementary negative-prefix workload, let `D = diag(1, 1, zeta_24)`,
let `F` be the exact qutrit Fourier transform over `K_24`, and construct

```text
U = (D F D) tensor I_3.
```

The one-qutrit factor has already been exhaustively rejected by the exact
one-qutrit search, so in particular it is not Clifford. The tensor-product
fixture is therefore also non-Clifford, but this benchmark makes no claim
about its global two-qutrit GSC classification.

The first canonical input Lagrangian has basis labels `(1, 0, 0, 0)` and
`(0, 1, 0, 0)`. Their exact conjugation expansions have support sizes seven
and one, respectively. Their eight-label union is contained in none of the 40
canonical output Lagrangians. The checked record uses:

```text
max_qudits = 2
max_coefficients = 162
max_candidate_pairs = 40
```

The run therefore computes both 81-coefficient expansions and rejects all 40
output candidates for exactly one input Lagrangian. At the next input
Lagrangian it reaches the coefficient cap and returns `UNKNOWN`. The exact
prefix rejection is complete; the overall search is deliberately incomplete
and is not a `NOT_GSC` result.

The next canonical input Lagrangian has basis labels `(1, 0, 0, 0)` and
`(0, 1, 0, 1)`. Its two exact expansions again have support sizes seven and
one, and their eight-label union is contained in none of the 40 canonical
output Lagrangians. A second checked record raises the budgets to

```text
max_qudits = 2
max_coefficients = 324
max_candidate_pairs = 80
```

This checks all 40 output candidates for each of the first two input
Lagrangians, then returns `UNKNOWN` at the coefficient cap before expanding
the third input. The longer prefix remains evidence about only those 80 exact
pairs, not a global negative classification.

The third canonical input Lagrangian has basis labels `(1, 0, 0, 0)` and
`(0, 1, 0, 2)`. Its exact support sizes are again seven and one, and no
canonical output Lagrangian contains their eight-label union. A third checked
record uses

```text
max_qudits = 2
max_coefficients = 486
max_candidate_pairs = 120
```

It rejects all output candidates for exactly the first three input
Lagrangians and then returns `UNKNOWN` before expanding the fourth. This is
still a bounded prefix, not a global `NOT_GSC` result.

The fourth canonical input Lagrangian has basis labels `(1, 0, 0, 1)` and
`(0, 1, 1, 0)`. Their exact conjugation expansions have support sizes seven
and three, and no canonical output Lagrangian contains their ten-label union.
A fourth checked record uses

```text
max_qudits = 2
max_coefficients = 648
max_candidate_pairs = 160
```

It rejects all output candidates for exactly the first four input
Lagrangians and then returns `UNKNOWN` before expanding the fifth. The result
remains bounded prefix evidence rather than a global classification.

The fifth canonical input Lagrangian has basis labels `(1, 0, 0, 1)` and
`(0, 1, 1, 1)`. Their exact conjugation expansions have support sizes seven
and three, and no canonical output Lagrangian contains their ten-label union.
A fifth checked record uses

```text
max_qudits = 2
max_coefficients = 810
max_candidate_pairs = 200
```

It rejects all output candidates for exactly the first five input
Lagrangians and then returns `UNKNOWN` before expanding the sixth. This
remains exact evidence about a bounded prefix, not a global classification.

The sixth canonical input Lagrangian has basis labels `(1, 0, 0, 1)` and
`(0, 1, 1, 2)`. Their exact conjugation expansions have support sizes seven
and three, and no canonical output Lagrangian contains their ten-label union.
A sixth checked record uses

```text
max_qudits = 2
max_coefficients = 972
max_candidate_pairs = 240
```

It rejects all output candidates for exactly the first six input Lagrangians
and then returns `UNKNOWN` before expanding the seventh. This is still a
bounded exact prefix and does not establish a global negative result.

The seventh canonical input Lagrangian has basis labels `(1, 0, 0, 2)` and
`(0, 1, 2, 0)`. Their exact conjugation expansions have support sizes seven
and three, and no canonical output Lagrangian contains their ten-label union.
A seventh checked record uses

```text
max_qudits = 2
max_coefficients = 1134
max_candidate_pairs = 280
```

It rejects all output candidates for exactly the first seven input
Lagrangians and then returns `UNKNOWN` before expanding the eighth. This
remains bounded exact prefix evidence rather than a global classification.

The eighth canonical input Lagrangian has basis labels `(1, 0, 0, 2)` and
`(0, 1, 2, 1)`. Their exact conjugation expansions have support sizes seven
and three, and no canonical output Lagrangian contains their ten-label union.
An eighth checked record uses

```text
max_qudits = 2
max_coefficients = 1296
max_candidate_pairs = 320
```

It rejects all output candidates for exactly the first eight input
Lagrangians and then returns `UNKNOWN` before expanding the ninth. This is
still a bounded exact prefix and makes no global GSC classification.

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

On 2026-10-19, the 162-coefficient rejection prefix measured 0.57 seconds and
196,574 bytes of Python peak memory on the development host. These host
measurements remain observational; the checked assertions cover only the
exact work counters, stop reason, and `UNKNOWN` semantics.

On 2026-10-20, the 324-coefficient two-input rejection prefix measured 1.05
seconds and 540,004 bytes of Python peak memory on the development host. The
suite does not turn either host measurement into a portable pass/fail limit.

On 2026-10-21, the 486-coefficient three-input rejection prefix measured 1.58
seconds and 643,128 bytes of Python peak memory on the development host. The
portable regression checks its exact-work counters and incomplete outcome,
not these machine-dependent measurements.

On 2026-10-22, the 648-coefficient four-input rejection prefix measured 1.98
seconds and 712,536 bytes of Python peak memory on the development host. Its
portable regression likewise checks the exact-work counters and `UNKNOWN`
outcome rather than these observational host measurements.

On 2026-10-23, the 810-coefficient five-input rejection prefix measured 5.43
seconds and 582,918 bytes of Python peak memory on the development host. The
portable regression again checks exact work and incomplete semantics only.

On 2026-10-25, the 972-coefficient six-input rejection prefix measured 3.10
seconds and 882,056 bytes of Python peak memory on the development host. The
checked contract continues to exclude these machine-dependent measurements.

On 2026-10-25, the 1,134-coefficient seven-input rejection prefix measured
3.69 seconds and 538,772 bytes of Python peak memory on the development host.
The portable regression checks the exact work counters and incomplete
semantics rather than these machine-dependent observations.

On 2026-10-26, the 1,296-coefficient eight-input rejection prefix measured
3.84 seconds and 761,607 bytes of Python peak memory on the development host.
The checked contract again excludes these observational host measurements.

## Checked-suite contract

The suite and its regression test now:

- construct the exact permutation above rather than loading approximate data;
- record the explicit limits and the expected capped `UNKNOWN` outcome;
- assert the portable counts and stop reason, but only require host runtime and
  peak memory to be positive finite measurements; and
- keep the capped and verified-positive two-qutrit records alongside, but
  semantically separate from, the one-qutrit positive, completed-negative, and
  capped records; and
- record exact one- through eight-input-Lagrangian rejection prefixes without
  promoting any bounded result to a global negative classification.
