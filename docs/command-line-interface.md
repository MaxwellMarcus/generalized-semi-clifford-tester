# Command-line interface

The `gsc-test` command emits one deterministic, versioned JSON document. Dense
tests read an arbitrary NumPy `.npy` unitary:

```console
gsc-test dense unitary.npy --tolerance 1e-9 --unitary-tolerance 1e-9
```

The result preserves `"status": "unknown"` and `"is_gsc": null` when the
configured exhaustive-search limit prevents a decision. Its claim-scope field
also records that a completed dense search is a complex128 numerical result,
not an exact symbolic proof.

Circuit tests read an OpenQASM 2 circuit and separate exploratory discovery
from fresh fixed-witness verification:

```console
gsc-test circuit circuit.qasm --discovery-shots 1024 --verification-shots 4096
```

The JSON includes every Bell-count histogram, both shot costs, the empirical
leakage evidence, the requested confidence level, and the simultaneous upper
bound from verification. `has_confidence_certified_witness` is never inferred
from discovery samples. Circuit support requires the optional Qiskit extra.

Exact odd-prime searches read a versioned JSON matrix over an explicitly
declared cyclotomic field:

```console
gsc-test exact-prime qutrit.json --prime 3 --max-qudits 2 \
  --max-coefficients 100000 --max-candidate-pairs 100000
```

The input schema is
`generalized-semi-clifford/exact-cyclotomic-matrix-v1`. `field_order` declares
`Q(zeta_m)`, and every matrix entry in `rows` is a nonempty list of integer or
rational-string coefficients in ascending powers of `zeta_m`. For example,
`["1/2", "0", "-1/3"]` denotes
`1/2 - zeta_m**2 / 3`. Floating-point JSON values and arbitrary symbolic
expressions are rejected rather than promoted to purportedly exact inputs.

```json
{
  "schema": "generalized-semi-clifford/exact-cyclotomic-matrix-v1",
  "field_order": 24,
  "rows": [
    [["1"], ["0"], ["0"]],
    [["0"], ["1"], ["0"]],
    [["0"], ["0"], ["1"]]
  ]
}
```

The command emits the existing
`generalized-semi-clifford.exact-prime-search.v1` result unchanged. Exact work
counters remain under `work`; host runtime and memory are not mixed into that
record. Qudit, coefficient, and candidate-pair caps preserve `"unknown"` with
`"search_complete": false`. Exact CLI support requires the optional `exact`
extra.
