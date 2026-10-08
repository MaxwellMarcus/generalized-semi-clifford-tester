# Odd-prime Weyl conventions

This document fixes the coordinate, phase, and symplectic conventions used by
the bounded odd-prime Lagrangian enumeration, exact fixed-witness verifier, and
exact exhaustive search. The implementation is deliberately restricted to
small prime-dimensional systems and does not claim unbounded tractability.

## Coordinates and single-qudit generators

Let `p` be an odd prime, let `omega = exp(2 pi i / p)`, and index the
computational basis by `j` in `F_p`. Define

```text
X |j> = |j + 1>,       Z |j> = omega^j |j>.
```

Thus `Z X = omega X Z`. An `n`-qudit projective Pauli label is
`v = (x | z)` in `F_p^(2n)`, with all coordinates represented canonically by
integers in `0, ..., p - 1`. Tensor powers are ordered by increasing qudit
index in both halves of the label.

Because `p` is odd, `2` has a multiplicative inverse. The phase-fixed Weyl
representative is

```text
W(x, z) = omega^((x dot z) / 2) X^x Z^z,
```

where division by `2`, vector addition, dot products, and exponents are all
computed modulo `p`.

## Pairing and multiplication

For `v = (x | z)` and `w = (x' | z')`, define

```text
[v, w] = z dot x' - x dot z'  (mod p).
```

For column coordinates this is `v^T J w`, where

```text
J = [[0, -I],
     [I,  0]].
```

The chosen Weyl representatives obey the exact identities

```text
W(v) W(w) = omega^([v,w] / 2) W(v + w),
W(v) W(w) = omega^[v,w] W(w) W(v).
```

Therefore two projective Paulis commute exactly when `[v,w] = 0`. A linear
coordinate map represented on column labels by `S` is symplectic exactly when
`S^T J S = J` over `F_p`.

## Independently checkable qutrit examples

For `p = 3`, `2^(-1) = 2` and

```text
X = [[0, 0, 1],       Z = diag(1, omega, omega^2),
     [1, 0, 0],
     [0, 1, 0]].
```

Using labels `e_X = (1 | 0)` and `e_Z = (0 | 1)` gives
`[e_X,e_Z] = -1 = 2 (mod 3)` and `[e_Z,e_X] = 1`. Hence

```text
X Z = omega^2 Z X,
W(1,1) = omega^2 X Z,
X Z = omega W(1,1),
Z X = omega^2 W(1,1).
```

The one-qutrit symplectic form itself,
`[[0, 2], [1, 0]]`, is symplectic over `F_3`; `diag(1,2)` is not. The four
one-dimensional Lagrangians are spanned by `(1|0)`, `(0|1)`, `(1|1)`, and
`(1|2)`. Canonical enumeration returns exactly these four lines; the two-qutrit
qutrit space has 40 canonical Lagrangian planes.

The executable regression in `tests/test_odd_prime.py` constructs the displayed
qutrit matrices and checks both Weyl identities for all 81 ordered label pairs.

## API boundary and result claims

`normalize_prime_label`, `prime_symplectic_pairing`,
`weyl_product_phase_exponent`, `prime_standard_form`, and
`is_prime_symplectic` implement the finite-field conventions.
`enumerate_prime_lagrangians` adds canonical enumeration for `p <= 5` and at
most two qudits. The optional exact module constructs dense cyclotomic Weyl
matrices, verifies supplied pairs, and searches the bounded canonical pair
space with explicit qudit, coefficient, and candidate-pair caps.

The search result is intentionally separate from fixed-witness verification.
`GSC` carries an independently recomputed witness, `NOT_GSC` means every
bounded pair was rejected exactly, and any interrupted search is `UNKNOWN`.
