"""Circuit execution for the standalone promised semi-Clifford theorem.

The theorem is conditional on exact hierarchy membership and ideal independent
sampling. User acknowledgement is not a machine verification of either fact.
No result here is a full generalized-semi-Clifford decision.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from enum import Enum
from typing import Any

import numpy as np

from .core import (
    ChoiSCAnalysis,
    ChoiSCSamplePlan,
    _positive,
    analyze_choi_bell_differences,
    plan_choi_sc_samples,
)


def _unitary_gate(unitary: Any) -> Any:
    try:
        from qiskit import QuantumCircuit
        from qiskit.circuit import Gate
    except ModuleNotFoundError as exc:  # pragma: no cover - optional dependency
        raise ModuleNotFoundError('Install Qiskit support with pip install -e ".[qiskit]"') from exc
    if isinstance(unitary, QuantumCircuit):
        if unitary.num_clbits or unitary.parameters:
            raise ValueError("unitary must have no classical bits and all parameters bound")
        try:
            unitary = unitary.to_gate(label="U")
        except Exception as exc:
            raise ValueError("unitary circuit must be convertible to a unitary gate") from exc
    if not isinstance(unitary, Gate):
        raise TypeError("unitary must be a Qiskit Gate or QuantumCircuit")
    if unitary.is_parameterized():
        raise ValueError("unitary must have all parameters bound")
    _positive(unitary.num_qubits, "unitary qubit count")
    return unitary


def build_choi_bell_sampling_circuit(unitary: Any, *, measure: bool = True) -> Any:
    """Bell-measure two Choi copies using two U calls on 4n qubits.

    Repeat independently twice and XOR the outcomes for one difference sample.
    Each copy orders its qubits as (output, reference), each block n qubits.
    Classical bits store z in the low 2n positions, then x in the high 2n.
    ``decode_choi_bell_outcome`` converts the displayed bitstring to (x|z).
    """
    gate = _unitary_gate(unitary)
    from qiskit import ClassicalRegister, QuantumCircuit

    n = gate.num_qubits
    circuit = QuantumCircuit(4 * n, name="choi_sc_bell_pair")
    for offset in (0, 2 * n):
        for qubit in range(n):
            circuit.h(offset + qubit)
            circuit.cx(offset + qubit, offset + n + qubit)
        circuit.append(gate, list(range(offset, offset + n)))
    for qubit in range(2 * n):
        circuit.cx(qubit, 2 * n + qubit)
        circuit.h(qubit)
    if measure:
        bell = ClassicalRegister(4 * n, "bell")
        circuit.add_register(bell)
        circuit.measure(range(4 * n), bell)
    circuit.metadata = {
        "experiment": "semi_clifford_choi_bell_pair",
        "unitary_qubits": n,
        "unitary_queries_per_shot": 2,
        "inverse_queries_per_shot": 0,
        "classical_bit_order": "z_output,z_reference,x_output,x_reference",
    }
    return circuit


def decode_choi_bell_outcome(bitstring: str, num_qubits: int) -> tuple[int, ...]:
    """Decode Qiskit's big-endian display into 4n (x|z) label bits."""
    _positive(num_qubits, "num_qubits")
    if not isinstance(bitstring, str):
        raise ValueError("Bell outcome must be a binary bitstring")
    bits = bitstring.replace(" ", "")
    if len(bits) != 4 * num_qubits or any(bit not in "01" for bit in bits):
        raise ValueError("Bell outcome must contain exactly 4n binary bits")
    low_first = tuple(int(bit) for bit in reversed(bits))
    return low_first[2 * num_qubits :] + low_first[: 2 * num_qubits]


class SCDecision(str, Enum):
    SC_SUPPORTED = "semi_clifford_supported"
    NON_SC = "not_semi_clifford"
    INCONCLUSIVE = "inconclusive"


@dataclass(frozen=True)
class ChoiSCResult:
    """Conditional SC result, never a non-GSC verdict or noise certificate."""

    decision: SCDecision
    plan: ChoiSCSamplePlan
    analysis: ChoiSCAnalysis | None
    circuit_shots: int
    unitary_queries: int
    sampler_jobs: int
    hierarchy_promise_acknowledged: bool
    ideal_sampling_acknowledged: bool
    conditional_false_accept_bound: float | None
    conditional_sc_distance_lower_bound: float | None
    message: str

    @property
    def inverse_queries(self) -> int:
        return 0


def run_choi_sc_test(
    unitary: Any,
    *,
    hierarchy_level: int,
    failure_probability: float = 0.01,
    assume_hierarchy_membership: bool = False,
    assume_ideal_sampling: bool = False,
    sampler: Any | None = None,
    seed: int | None = None,
    difference_samples: int | None = None,
    batch_shots: int = 1024,
    max_difference_samples: int = 100_000,
    max_simulation_qubits: int = 12,
) -> ChoiSCResult:
    """Run the separate SC tester with a SamplerV2 backend.

    Default: theorem budget and local StatevectorSampler, capped at 4n=12
    simulation qubits. Custom samplers bypass that cap and must expose ordered
    ``data.bell.get_bitstrings()``; counts alone cannot be paired safely.

    Sampling must be fresh across jobs. The local sampler uses a persistent
    Generator so integer reseeding cannot repeat batches. ``seed`` applies
    only to the local sampler; custom samplers require their own configuration.

    Positive decisions require both assumption acknowledgements and the full
    planned budget. An ideal-sampling negative excludes SC without a hierarchy
    promise, but then has no distance bound. Acknowledgements state premises,
    not verified facts. Bounds remain conditional even in floating simulation.
    Work caps never imply rejection. GSC membership is not decided here.
    """
    for value, name in (
        (batch_shots, "batch_shots"),
        (max_difference_samples, "max_difference_samples"),
        (max_simulation_qubits, "max_simulation_qubits"),
    ):
        _positive(value, name)
    if not isinstance(assume_hierarchy_membership, bool) or not isinstance(
        assume_ideal_sampling, bool
    ):
        raise ValueError("assumption acknowledgements must be Boolean")
    gate = _unitary_gate(unitary)
    n = gate.num_qubits
    plan = plan_choi_sc_samples(n, hierarchy_level, failure_probability=failure_probability)
    requested = plan.bell_difference_samples if difference_samples is None else difference_samples
    _positive(requested, "difference_samples")
    completed_shots = 0
    jobs = 0

    def finish(
        decision: SCDecision,
        analysis: ChoiSCAnalysis | None,
        message: str,
    ) -> ChoiSCResult:
        premises = assume_hierarchy_membership and assume_ideal_sampling
        full = analysis is not None and analysis.samples_used >= plan.bell_difference_samples
        return ChoiSCResult(
            decision,
            plan,
            analysis,
            completed_shots,
            2 * completed_shots,
            jobs,
            assume_hierarchy_membership,
            assume_ideal_sampling,
            failure_probability if premises and full else None,
            plan.sc_distance_gap if premises and decision is SCDecision.NON_SC else None,
            message,
        )

    if requested > max_difference_samples:
        return finish(
            SCDecision.INCONCLUSIVE, None, "Requested sample budget exceeds the work cap."
        )
    if sampler is None:
        if 4 * n > max_simulation_qubits:
            return finish(
                SCDecision.INCONCLUSIVE,
                None,
                "Local statevector size cap exceeded; supply a suitable sampler or raise the cap.",
            )
        from qiskit.primitives import StatevectorSampler

        sampler = StatevectorSampler(seed=np.random.default_rng(seed))
    circuit = build_choi_bell_sampling_circuit(gate)

    def differences() -> Iterator[tuple[int, ...]]:
        nonlocal completed_shots, jobs
        pending = None
        remaining = 2 * requested
        while remaining:
            shots = min(batch_shots, remaining)
            result = sampler.run([circuit], shots=shots).result()
            jobs += 1
            try:
                bitstrings = result[0].data.bell.get_bitstrings()
            except AttributeError as exc:
                raise ValueError(
                    "sampler must return ordered Bell bitstrings, not just counts"
                ) from exc
            if len(bitstrings) != shots:
                raise ValueError("sampler returned a different number of shots than requested")
            completed_shots += shots
            remaining -= shots
            for bitstring in bitstrings:
                label = decode_choi_bell_outcome(bitstring, n)
                if pending is None:
                    pending = label
                else:
                    yield tuple(a ^ b for a, b in zip(pending, label, strict=True))
                    pending = None
        if pending is not None:  # Even total shot count rules this out.
            raise AssertionError("unpaired Bell outcome")

    analysis = analyze_choi_bell_differences(n, differences())
    if not assume_ideal_sampling:
        return finish(
            SCDecision.INCONCLUSIVE,
            analysis,
            "Raw SC compatibility only: ideal independent sampling has not been acknowledged.",
        )
    if not analysis.semi_clifford_compatible:
        return finish(
            SCDecision.NON_SC,
            analysis,
            "Incompatible with SC under ideal sampling; this is not a non-GSC conclusion.",
        )
    if not assume_hierarchy_membership or requested < plan.bell_difference_samples:
        return finish(
            SCDecision.INCONCLUSIVE,
            analysis,
            "SC-compatible data, but the hierarchy promise or full theorem budget is missing.",
        )
    return finish(
        SCDecision.SC_SUPPORTED,
        analysis,
        "SC supported at the requested error bound, conditional on acknowledged ideal premises.",
    )
