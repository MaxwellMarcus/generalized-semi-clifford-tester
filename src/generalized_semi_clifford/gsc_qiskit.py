"""Circuit-sampling prototypes for generalized semi-Clifford testing.

The quantum experiment is the same Pauli-conjugation Bell measurement used by
the semi-Clifford tester.  The post-processing is different: a generalized
semi-Clifford witness only requires the *support* of each conjugated generator
to lie in a common output Pauli MASA, rather than each conjugate being a single
Pauli.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import Any

from .lagrangian import (
    Lagrangian,
    PauliLabel,
    enumerate_lagrangians,
    lagrangian_containing,
    lagrangian_from_basis,
)
from .semi_clifford_qiskit import (
    PauliConjugationObservation,
    bell_bitstring_to_pauli,
    run_semi_clifford_sampling_test,
)


@dataclass(frozen=True)
class CircuitGSCWitness:
    """Candidate Pauli-MASA witness obtained from finite-shot circuits."""

    input_lagrangian: Lagrangian
    output_lagrangian: Lagrangian
    maximum_empirical_leakage: float


@dataclass(frozen=True)
class GSCSamplingResult:
    """Result and diagnostics for a circuit-based GSC experiment."""

    num_qubits: int
    shots_per_pauli: int
    leakage_threshold: float
    observations: tuple[PauliConjugationObservation, ...]
    witness: CircuitGSCWitness | None
    best_empirical_leakage: float
    input_lagrangians_checked: int
    candidate_pairs_checked: int
    search_mode: str

    @property
    def has_candidate_witness(self) -> bool:
        """Return whether the sampled data supports a Pauli-MASA witness."""

        return self.witness is not None


def output_pauli_probabilities(
    observation: PauliConjugationObservation,
) -> dict[PauliLabel, float]:
    """Decode one Bell-count histogram as a Pauli probability distribution."""

    num_qubits = len(observation.input_pauli) // 2
    probabilities: dict[PauliLabel, float] = {}
    for bitstring, count in observation.counts:
        label = bell_bitstring_to_pauli(bitstring, num_qubits)
        probabilities[label] = probabilities.get(label, 0.0) + count / observation.shots
    return probabilities


def empirical_lagrangian_leakage(
    observation: PauliConjugationObservation,
    output_lagrangian: Lagrangian,
) -> float:
    """Return observed probability mass outside an output Lagrangian."""

    if len(observation.input_pauli) != 2 * output_lagrangian.num_qubits:
        raise ValueError("observation and output Lagrangian qubit counts differ")
    allowed = set(output_lagrangian.elements)
    inside_probability = sum(
        probability
        for label, probability in output_pauli_probabilities(observation).items()
        if label in allowed
    )
    return max(0.0, min(1.0, 1.0 - inside_probability))


def find_gsc_sampling_witness(
    observations: Iterable[PauliConjugationObservation],
    *,
    input_lagrangians: Iterable[Lagrangian],
    output_lagrangians: Iterable[Lagrangian],
    leakage_threshold: float,
) -> tuple[CircuitGSCWitness | None, float, int, int]:
    """Search sampled conjugations for an input/output Pauli-MASA pair."""

    if not 0 <= leakage_threshold <= 1:
        raise ValueError("leakage_threshold must lie in [0, 1]")
    by_input = {observation.input_pauli: observation for observation in observations}
    inputs = tuple(input_lagrangians)
    outputs = tuple(output_lagrangians)
    if not inputs or not outputs:
        raise ValueError("candidate Lagrangian collections must not be empty")

    best_leakage = 1.0
    allowed_output_elements = {output.elements for output in outputs}
    candidate_pairs_checked = 0
    eligible_inputs: list[
        tuple[Lagrangian, tuple[PauliConjugationObservation, ...]]
    ] = []
    for input_lagrangian in inputs:
        if any(label not in by_input for label in input_lagrangian.basis):
            continue
        basis_observations = tuple(by_input[label] for label in input_lagrangian.basis)
        eligible_inputs.append((input_lagrangian, basis_observations))

    # First try the output subspace forced by every outcome heavier than the
    # entire allowed leakage budget. This reduces an exact-support search from
    # all LxS pairs to at most one constructed S per input Lagrangian.
    constructed_outputs: dict[tuple[PauliLabel, ...], Lagrangian | None] = {}
    for input_index, (input_lagrangian, basis_observations) in enumerate(
        eligible_inputs,
        start=1,
    ):
        required_output_labels = {
            label
            for observation in basis_observations
            for label, probability in output_pauli_probabilities(observation).items()
            if probability > leakage_threshold
        }
        constructed_output = lagrangian_containing(
            required_output_labels,
            input_lagrangian.num_qubits,
        )
        constructed_outputs[input_lagrangian.elements] = constructed_output
        if (
            constructed_output is not None
            and constructed_output.elements in allowed_output_elements
        ):
            candidate_pairs_checked += 1
            leakage = max(
                empirical_lagrangian_leakage(observation, constructed_output)
                for observation in basis_observations
            )
            best_leakage = min(best_leakage, leakage)
            if leakage <= leakage_threshold:
                return (
                    CircuitGSCWitness(
                        input_lagrangian=input_lagrangian,
                        output_lagrangian=constructed_output,
                        maximum_empirical_leakage=leakage,
                    ),
                    best_leakage,
                    input_index,
                    candidate_pairs_checked,
                )

    # With noisy data, the forced heavy-label span may not determine the best
    # output Lagrangian. Fall back to the complete candidate-pair search.
    for input_index, (input_lagrangian, basis_observations) in enumerate(
        eligible_inputs,
        start=1,
    ):
        constructed_output = constructed_outputs[input_lagrangian.elements]
        for output_lagrangian in outputs:
            if (
                constructed_output is not None
                and output_lagrangian.elements == constructed_output.elements
            ):
                continue
            if output_lagrangian.num_qubits != input_lagrangian.num_qubits:
                raise ValueError("input and output Lagrangians must use the same qubit count")
            candidate_pairs_checked += 1
            leakage = max(
                empirical_lagrangian_leakage(observation, output_lagrangian)
                for observation in basis_observations
            )
            best_leakage = min(best_leakage, leakage)
            if leakage <= leakage_threshold:
                return (
                    CircuitGSCWitness(
                        input_lagrangian=input_lagrangian,
                        output_lagrangian=output_lagrangian,
                        maximum_empirical_leakage=leakage,
                    ),
                    best_leakage,
                    input_index,
                    candidate_pairs_checked,
                )
    return None, best_leakage, len(eligible_inputs), candidate_pairs_checked


def _unitary_num_qubits(unitary: Any) -> int:
    num_qubits = getattr(unitary, "num_qubits", None)
    if not isinstance(num_qubits, int) or num_qubits < 1:
        raise TypeError("unitary must be a Qiskit circuit or instruction acting on qubits")
    return num_qubits


def run_gsc_sampling_test(
    unitary: Any,
    *,
    unitary_dagger: Any | None = None,
    sampler: Any | None = None,
    shots: int = 1024,
    leakage_threshold: float = 0.01,
    seed: int | None = None,
    batch_size: int | None = None,
) -> GSCSamplingResult:
    """Exhaustively search small-qubit Pauli MASAs using circuit samples.

    One canonical basis is fixed for every input Lagrangian.  Only Paulis that
    occur in those bases need circuits; their observations are reused across
    every candidate input/output pair.
    """

    num_qubits = _unitary_num_qubits(unitary)
    lagrangians = enumerate_lagrangians(num_qubits)
    required_inputs = tuple(
        dict.fromkeys(label for lagrangian in lagrangians for label in lagrangian.basis)
    )
    sampling = run_semi_clifford_sampling_test(
        unitary,
        unitary_dagger=unitary_dagger,
        sampler=sampler,
        shots=shots,
        pauli_probability_threshold=1.0,
        seed=seed,
        input_paulis=required_inputs,
        batch_size=batch_size,
    )
    witness, best_leakage, inputs_checked, pairs_checked = find_gsc_sampling_witness(
        sampling.observations,
        input_lagrangians=lagrangians,
        output_lagrangians=lagrangians,
        leakage_threshold=leakage_threshold,
    )
    return GSCSamplingResult(
        num_qubits=num_qubits,
        shots_per_pauli=shots,
        leakage_threshold=leakage_threshold,
        observations=sampling.observations,
        witness=witness,
        best_empirical_leakage=best_leakage,
        input_lagrangians_checked=inputs_checked,
        candidate_pairs_checked=pairs_checked,
        search_mode="exhaustive-lagrangian-discovery",
    )


def run_gsc_witness_test(
    unitary: Any,
    input_basis: Iterable[Sequence[int]],
    output_basis: Iterable[Sequence[int]],
    *,
    unitary_dagger: Any | None = None,
    sampler: Any | None = None,
    shots: int = 1024,
    leakage_threshold: float = 0.01,
    seed: int | None = None,
) -> GSCSamplingResult:
    """Check a proposed GSC Pauli-MASA pair using only ``n`` circuits."""

    num_qubits = _unitary_num_qubits(unitary)
    input_lagrangian = lagrangian_from_basis(input_basis, num_qubits)
    output_lagrangian = lagrangian_from_basis(output_basis, num_qubits)
    sampling = run_semi_clifford_sampling_test(
        unitary,
        unitary_dagger=unitary_dagger,
        sampler=sampler,
        shots=shots,
        pauli_probability_threshold=1.0,
        seed=seed,
        input_paulis=input_lagrangian.basis,
    )
    witness, best_leakage, inputs_checked, pairs_checked = find_gsc_sampling_witness(
        sampling.observations,
        input_lagrangians=(input_lagrangian,),
        output_lagrangians=(output_lagrangian,),
        leakage_threshold=leakage_threshold,
    )
    return GSCSamplingResult(
        num_qubits=num_qubits,
        shots_per_pauli=shots,
        leakage_threshold=leakage_threshold,
        observations=sampling.observations,
        witness=witness,
        best_empirical_leakage=best_leakage,
        input_lagrangians_checked=inputs_checked,
        candidate_pairs_checked=pairs_checked,
        search_mode="candidate-witness-verification",
    )
