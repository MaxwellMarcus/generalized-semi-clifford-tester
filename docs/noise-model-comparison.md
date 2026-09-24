# Qiskit noise-model comparison

`examples/qiskit_noise_model.py` compares the ideal sparse Pauli-transfer
prediction for a fixed one-qubit witness with Bell-sampling results from both a
statevector sampler and an explicit Aer depolarizing noise model. Install the
noise-model dependencies with:

```console
python -m pip install -e ".[noise]"
python examples/qiskit_noise_model.py
```

The sparse-transfer side reports an interval. Its lower endpoint is retained
leakage mass outside the proposed output Lagrangian; its upper endpoint also
includes all coefficient mass discarded by numerical thresholding. It is a
small-qubit complex128 calculation for the ideal unitary, not a certificate of
exact zeros and not a model of device noise.

The sampling side reports observed leakage separately from a simultaneous
one-sided Clopper--Pearson upper bound. That confidence statement applies to
the fixed witness used by the example. It does not promote finite shots,
simulator output, or a fitted noise model to an exact generalized
semi-Clifford membership claim. A witness selected from exploratory samples
must be verified with fresh samples before using the bound.

The example transpiles custom unitary instructions to the noise model's basis
at optimization level zero. The model applies one-qubit depolarizing error to
every resulting `h` and two-qubit depolarizing error to every resulting `cx`,
including Bell-pair preparation and measurement-basis changes in the test
circuit. It is a documented reproducible simulation, not a calibrated IBM
hardware model.
