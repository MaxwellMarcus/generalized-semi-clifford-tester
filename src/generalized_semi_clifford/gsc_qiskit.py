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
from .statistics import binomial_proportion_upper_bound


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
    confidence_level: float | None = None
    maximum_leakage_upper_bound: float | None = None

    @property
    def has_candidate_witness(self) -> bool:
        """Return whether the sampled data supports a Pauli-MASA witness."""

        return self.witness is not None

    @property
    def has_confidence_certified_witness(self) -> bool:
        """Return whether a fixed witness meets the threshold with confidence.

        Certification is only available from :func:`run_gsc_witness_test`.
        Discovery chooses a witness from the same samples and therefore needs
        independent confirmation before this claim is statistically valid.
        """

        return (
            self.witness is not None
            and self.maximum_leakage_upper_bound is not None
            and self.maximum_leakage_upper_bound <= self.leakage_threshold
        )


@dataclass(frozen=True)
class GSCDiscoveryVerificationResult:
    """Separate exploratory discovery from fixed-witness verification.

    ``discovery`` and ``verification`` always come from distinct sampler runs.
    The verification result is absent when discovery finds no candidate. Shot
    costs count all circuit executions in each phase.
    """

    discovery: GSCSamplingResult
    verification: GSCSamplingResult | None
    discovery_shot_cost: int
    verification_shot_cost: int

    def __post_init__(self) -> None:
        expected_discovery_cost = sum(
            observation.shots for observation in self.discovery.observations
        )
        expected_verification_cost = (
            0
            if self.verification is None
            else sum(observation.shots for observation in self.verification.observations)
        )
        if self.discovery_shot_cost != expected_discovery_cost:
            raise ValueError("discovery_shot_cost must match the discovery observations")
        if self.verification_shot_cost != expected_verification_cost:
            raise ValueError("verification_shot_cost must match the verification observations")

    @property
    def has_discovered_candidate(self) -> bool:
        """Whether the exploratory phase selected a candidate witness."""

        return self.discovery.has_candidate_witness

    @property
    def has_confidence_certified_witness(self) -> bool:
        """Whether fresh verification certified the selected fixed witness."""

        return (
            self.verification is not None
            and self.verification.has_confidence_certified_witness
        )


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


def lagrangian_leakage_upper_bound(
    observations: Iterable[PauliConjugationObservation],
    output_lagrangian: Lagrangian,
    *,
    confidence_level: float,
) -> float:
    """Bound every generator's leakage with familywise confidence.

    Each generator is a binomial experiment whose event is an output label
    outside ``output_lagrangian``.  One-sided exact Clopper--Pearson bounds are
    combined with a Bonferroni correction.  Consequently, all returned
    per-generator bounds hold simultaneously with probability at least
    ``confidence_level`` for a witness fixed independently of these samples.
    """

    samples = tuple(observations)
    if not samples:
        raise ValueError("observations must not be empty")
    if not 0 < confidence_level < 1:
        raise ValueError("confidence_level must lie strictly between zero and one")
    allowed = set(output_lagrangian.elements)
    per_test_failure_probability = (1.0 - confidence_level) / len(samples)
    per_test_confidence = 1.0 - per_test_failure_probability
    bounds = []
    for observation in samples:
        if len(observation.input_pauli) != 2 * output_lagrangian.num_qubits:
            raise ValueError("observation and output Lagrangian qubit counts differ")
        outside_count = sum(
            count
            for bitstring, count in observation.counts
            if bell_bitstring_to_pauli(bitstring, output_lagrangian.num_qubits)
            not in allowed
        )
        bounds.append(
            binomial_proportion_upper_bound(
                outside_count,
                observation.shots,
                confidence_level=per_test_confidence,
            )
        )
    return max(bounds)


def _output_support_index(
    output_lagrangians: tuple[Lagrangian, ...],
) -> dict[PauliLabel, frozenset[int]]:
    indices: dict[PauliLabel, set[int]] = {}
    for output_index, output_lagrangian in enumerate(output_lagrangians):
        for label in output_lagrangian.elements:
            if any(label):
                indices.setdefault(label, set()).add(output_index)
    return {label: frozenset(items) for label, items in indices.items()}


def _support_aware_output_indices(
    observations: tuple[PauliConjugationObservation, ...],
    output_lagrangians: tuple[Lagrangian, ...],
    *,
    leakage_threshold: float,
    support_index: dict[PauliLabel, frozenset[int]],
) -> tuple[int, ...]:
    """Return a threshold-complete set of noisy output candidates.

    Every output Lagrangian contains the identity. If an observation's
    identity mass is smaller than ``1 - leakage_threshold``, any acceptable
    output must also contain at least one observed nonidentity label. The union
    of indexed Lagrangians containing those labels is therefore a necessary
    candidate set. Intersecting that set across observations cannot discard a
    threshold-feasible witness.
    """

    candidates = set(range(len(output_lagrangians)))
    required_inside_probability = 1.0 - leakage_threshold
    for observation in observations:
        probabilities = output_pauli_probabilities(observation)
        identity = (0,) * len(observation.input_pauli)
        if probabilities.get(identity, 0.0) >= required_inside_probability:
            continue
        supported_outputs: set[int] = set()
        for label, probability in probabilities.items():
            if probability > 0.0 and any(label):
                supported_outputs.update(support_index.get(label, ()))
        candidates.intersection_update(supported_outputs)
        if not candidates:
            break
    return tuple(index for index in range(len(output_lagrangians)) if index in candidates)


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
    support_index = _output_support_index(outputs)
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
    # output Lagrangian. Use observed support to retain every threshold-feasible
    # output while avoiding scores for candidates that cannot contain enough
    # mass. A final exhaustive diagnostic pass runs only when no witness exists,
    # preserving the exact historical best-leakage value on negative results.
    for input_index, (input_lagrangian, basis_observations) in enumerate(
        eligible_inputs,
        start=1,
    ):
        constructed_output = constructed_outputs[input_lagrangian.elements]
        support_candidates = _support_aware_output_indices(
            basis_observations,
            outputs,
            leakage_threshold=leakage_threshold,
            support_index=support_index,
        )
        evaluated_indices: set[int] = set()
        for output_index in support_candidates:
            output_lagrangian = outputs[output_index]
            if (
                constructed_output is not None
                and output_lagrangian.elements == constructed_output.elements
            ):
                continue
            evaluated_indices.add(output_index)
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

        for output_index, output_lagrangian in enumerate(outputs):
            if output_index in evaluated_indices or (
                constructed_output is not None
                and output_lagrangian.elements == constructed_output.elements
            ):
                continue
            candidate_pairs_checked += 1
            leakage = max(
                empirical_lagrangian_leakage(observation, output_lagrangian)
                for observation in basis_observations
            )
            best_leakage = min(best_leakage, leakage)
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
        search_mode="support-aware-lagrangian-discovery",
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
    confidence_level: float = 0.95,
    seed: int | None = None,
) -> GSCSamplingResult:
    """Check a fixed GSC Pauli-MASA pair using only ``n`` circuits.

    In addition to the empirical decision, the result contains a simultaneous
    one-sided confidence bound for the maximum generator leakage.  The bound
    is statistically valid when the proposed witness was fixed independently
    of these samples; use fresh samples after exploratory discovery.
    """

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
    leakage_upper_bound = lagrangian_leakage_upper_bound(
        sampling.observations,
        output_lagrangian,
        confidence_level=confidence_level,
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
        confidence_level=confidence_level,
        maximum_leakage_upper_bound=leakage_upper_bound,
    )


def run_gsc_discovery_then_verification(
    unitary: Any,
    *,
    unitary_dagger: Any | None = None,
    discovery_sampler: Any | None = None,
    verification_sampler: Any | None = None,
    discovery_shots: int = 1024,
    verification_shots: int = 1024,
    leakage_threshold: float = 0.01,
    confidence_level: float = 0.95,
    discovery_seed: int | None = None,
    verification_seed: int | None = None,
    discovery_batch_size: int | None = None,
) -> GSCDiscoveryVerificationResult:
    """Discover a GSC candidate, then verify it using a fresh sampler run.

    The exploratory counts are used only to choose a candidate. When one is
    found, its bases are passed to :func:`run_gsc_witness_test`, which obtains
    new counts and applies the fixed-witness confidence bound. Custom samplers
    must likewise return new observations on each ``run`` call.
    """

    discovery = run_gsc_sampling_test(
        unitary,
        unitary_dagger=unitary_dagger,
        sampler=discovery_sampler,
        shots=discovery_shots,
        leakage_threshold=leakage_threshold,
        seed=discovery_seed,
        batch_size=discovery_batch_size,
    )
    verification = None
    if discovery.witness is not None:
        verification = run_gsc_witness_test(
            unitary,
            discovery.witness.input_lagrangian.basis,
            discovery.witness.output_lagrangian.basis,
            unitary_dagger=unitary_dagger,
            sampler=verification_sampler,
            shots=verification_shots,
            leakage_threshold=leakage_threshold,
            confidence_level=confidence_level,
            seed=verification_seed,
        )
    return GSCDiscoveryVerificationResult(
        discovery=discovery,
        verification=verification,
        discovery_shot_cost=sum(
            observation.shots for observation in discovery.observations
        ),
        verification_shot_cost=(
            0
            if verification is None
            else sum(observation.shots for observation in verification.observations)
        ),
    )
