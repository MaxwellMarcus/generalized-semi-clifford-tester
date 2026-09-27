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
