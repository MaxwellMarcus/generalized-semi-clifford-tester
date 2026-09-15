"""Qiskit circuits for experimental semi-Clifford property testing.

This module deliberately lives beside, rather than inside, the dense-matrix
tester.  It treats a unitary as an operation that can be queried in a circuit.
For each phase-free Pauli ``P``, a Bell-basis measurement samples the Pauli
coefficients of ``U P U†``.  Classical binary-symplectic post-processing then
looks for ``n`` independent, commuting successful inputs.

Finite-shot concentration on a Pauli label is evidence, not an exact proof,
that a conjugate is Pauli.  The present implementation is therefore an
experimental prototype; it does not yet prove a tolerant-testing theorem.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any

from .lagrangian import PauliLabel, lagrangian_from_basis
from .pauli import all_pauli_labels


@dataclass(frozen=True)
class PauliConjugationObservation:
    """Finite-shot evidence about one conjugate ``U P U†``."""

    input_pauli: PauliLabel
    dominant_output_pauli: PauliLabel
    dominant_probability: float
    shots: int
    counts: tuple[tuple[str, int], ...]


@dataclass(frozen=True)
class SemiCliffordWitness:
    """Candidate input/output Lagrangian bases found from sampled data."""

    input_basis: tuple[PauliLabel, ...]
    output_basis: tuple[PauliLabel, ...]
    minimum_dominant_probability: float


@dataclass(frozen=True)
class SemiCliffordSamplingResult:
    """Result of the hybrid circuit-sampling and symplectic search."""

    num_qubits: int
    shots_per_pauli: int
    pauli_probability_threshold: float
    observations: tuple[PauliConjugationObservation, ...]
    witness: SemiCliffordWitness | None
    exhaustive_pauli_search: bool

    @property
    def has_candidate_witness(self) -> bool:
        """Return whether the sampled data contains a Lagrangian witness."""

        return self.witness is not None


def _qiskit_types() -> tuple[type[Any], type[Any], type[Any], type[Any]]:
    try:
        from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister
        from qiskit.circuit import Instruction
    except ModuleNotFoundError as exc:  # pragma: no cover - depends on optional install
        raise ModuleNotFoundError(
            'Qiskit support requires: pip install -e ".[qiskit]"'
        ) from exc
    return QuantumCircuit, QuantumRegister, ClassicalRegister, Instruction


def _as_instruction(operation: Any, *, name: str) -> Any:
    QuantumCircuit, _, _, Instruction = _qiskit_types()
    if isinstance(operation, QuantumCircuit):
        if operation.num_clbits:
            raise ValueError(f"{name} must not contain classical bits or measurements")
        if operation.parameters:
            raise ValueError(f"{name} must have all parameters bound")
        return operation.to_gate(label=name)
    if not isinstance(operation, Instruction):
        raise TypeError(f"{name} must be a Qiskit QuantumCircuit or Instruction")
    return operation


def _validate_pauli_label(label: Sequence[int]) -> PauliLabel:
    normalized = tuple(label)
    if not normalized or len(normalized) % 2 or any(bit not in (0, 1) for bit in normalized):
        raise ValueError("a Pauli label must contain 2n binary entries")
    return normalized


def build_pauli_conjugation_test_circuit(
    unitary: Any,
    input_pauli: Sequence[int],
    *,
    unitary_dagger: Any | None = None,
    measure: bool = True,
) -> Any:
    """Build a Bell-sampling circuit for the Pauli expansion of ``U P U†``.

    ``unitary_dagger`` is a separate argument to mirror the oracle access model
    used in the literature.  For an explicitly described circuit it may be
    omitted, in which case Qiskit constructs the inverse instruction.

    Labels use ``(x | z)`` coordinates indexed by Qiskit qubit number.  The
    classical register stores the Bell measurement as ``(z | x)`` so that its
    displayed bit string can be decoded by :func:`bell_counts_to_observation`.
    """

    QuantumCircuit, QuantumRegister, ClassicalRegister, _ = _qiskit_types()
    label = _validate_pauli_label(input_pauli)
    num_qubits = len(label) // 2
    unitary_instruction = _as_instruction(unitary, name="unitary")
    dagger_instruction = (
        unitary_instruction.inverse()
        if unitary_dagger is None
        else _as_instruction(unitary_dagger, name="unitary_dagger")
    )
    if unitary_instruction.num_qubits != num_qubits:
        raise ValueError("unitary qubit count does not match the Pauli label")
    if dagger_instruction.num_qubits != num_qubits:
        raise ValueError("unitary_dagger qubit count does not match the Pauli label")

    system = QuantumRegister(num_qubits, "system")
    reference = QuantumRegister(num_qubits, "reference")
    bell = ClassicalRegister(2 * num_qubits, "bell")
    circuit = QuantumCircuit(system, reference, bell, name="pauli_conjugation_test")

    # Prepare |Phi> = 2^(-n/2) sum_x |x>|x>.
    for qubit in range(num_qubits):
        circuit.h(system[qubit])
        circuit.cx(system[qubit], reference[qubit])

    # Circuit order is U†, P, U, so the net operation is U P U†.
    circuit.append(dagger_instruction, list(system))
    for qubit, (x_bit, z_bit) in enumerate(
        zip(label[:num_qubits], label[num_qubits:], strict=True)
    ):
        if x_bit and z_bit:
            circuit.y(system[qubit])
        elif x_bit:
            circuit.x(system[qubit])
        elif z_bit:
            circuit.z(system[qubit])
    circuit.append(unitary_instruction, list(system))

    # Invert Bell-state preparation.  For each pair, the system result is z
    # and the reference result is x for the corresponding output Pauli.
    for qubit in range(num_qubits):
        circuit.cx(system[qubit], reference[qubit])
        circuit.h(system[qubit])
    if measure:
        for qubit in range(num_qubits):
            circuit.measure(system[qubit], bell[qubit])
            circuit.measure(reference[qubit], bell[num_qubits + qubit])

    circuit.metadata = {
        "experiment": "pauli_conjugation_bell_sampling",
        "input_pauli_xz": list(label),
        "num_qubits": num_qubits,
    }
    return circuit


def build_semi_clifford_test_circuits(
    unitary: Any,
    *,
    unitary_dagger: Any | None = None,
    input_paulis: Iterable[Sequence[int]] | None = None,
) -> tuple[tuple[PauliLabel, Any], ...]:
    """Build conjugation-test circuits for selected nonidentity Paulis.

    If ``input_paulis`` is omitted, all ``4**n - 1`` nonidentity Paulis are
    used. Repeated labels are removed while preserving input order.
    """

    instruction = _as_instruction(unitary, name="unitary")
    num_qubits = instruction.num_qubits
    if num_qubits < 1:
        raise ValueError("unitary must act on at least one qubit")
    if input_paulis is None:
        labels = tuple(label for label in all_pauli_labels(num_qubits) if any(label))
    else:
        unique_labels: dict[PauliLabel, None] = {}
        for candidate in input_paulis:
            label = _validate_pauli_label(candidate)
            if len(label) != 2 * num_qubits:
                raise ValueError("input Pauli qubit count does not match the unitary")
            if not any(label):
                raise ValueError("the identity cannot contribute to a Lagrangian basis")
            unique_labels.setdefault(label, None)
        labels = tuple(unique_labels)
        if not labels:
            raise ValueError("input_paulis must contain at least one nonidentity Pauli")
    experiments = []
    for label in labels:
        experiments.append(
            (
                label,
                build_pauli_conjugation_test_circuit(
                    instruction,
                    label,
                    unitary_dagger=unitary_dagger,
                ),
            )
        )
    return tuple(experiments)


def bell_bitstring_to_pauli(bitstring: str, num_qubits: int) -> PauliLabel:
    """Decode the Qiskit Bell register into an ``(x | z)`` Pauli label."""

    compact = bitstring.replace(" ", "")
    if len(compact) != 2 * num_qubits or any(bit not in "01" for bit in compact):
        raise ValueError(f"expected a {2 * num_qubits}-bit Bell outcome, got {bitstring!r}")
    classical_bits = tuple(int(bit) for bit in reversed(compact))
    z_bits = classical_bits[:num_qubits]
    x_bits = classical_bits[num_qubits:]
    return x_bits + z_bits


def bell_counts_to_observation(
    input_pauli: Sequence[int],
    counts: Mapping[str, int],
) -> PauliConjugationObservation:
    """Decode Qiskit Bell-measurement counts into Pauli evidence."""

    label = _validate_pauli_label(input_pauli)
    num_qubits = len(label) // 2
    if not counts or any(count < 0 for count in counts.values()):
        raise ValueError("counts must contain nonnegative shot counts")
    shots = sum(counts.values())
    if shots <= 0:
        raise ValueError("counts must contain at least one shot")
    dominant_bitstring, dominant_count = max(counts.items(), key=lambda item: item[1])
    return PauliConjugationObservation(
        input_pauli=label,
        dominant_output_pauli=bell_bitstring_to_pauli(dominant_bitstring, num_qubits),
        dominant_probability=dominant_count / shots,
        shots=shots,
        counts=tuple(sorted(counts.items())),
    )


def _label_to_int(label: PauliLabel) -> int:
    return sum(bit << index for index, bit in enumerate(label))


def _binary_rank(labels: Iterable[PauliLabel]) -> int:
    rows = [_label_to_int(label) for label in labels]
    rank = 0
    while rows:
        pivot = max(rows)
        if pivot == 0:
            break
        rank += 1
        pivot_bit = 1 << (pivot.bit_length() - 1)
        rows = [row ^ pivot if row & pivot_bit else row for row in rows if row != pivot]
    return rank


def maximum_isotropic_dimension(
    labels: Iterable[Sequence[int]],
    *,
    num_qubits: int,
) -> int:
    """Return the largest isotropic dimension in the span of ``labels``.

    This criterion is intended for an exact subspace such as ideal ``K_U``.
    It must not be applied to a noisy accepted set whose span has not itself
    been verified: taking that span could manufacture untested conjugations.
    """

    if num_qubits < 1:
        raise ValueError("num_qubits must be positive")
    normalized = tuple(_validate_pauli_label(label) for label in labels)
    if any(len(label) != 2 * num_qubits for label in normalized):
        raise ValueError("Pauli label width does not match num_qubits")

    # Extract an independent basis greedily over F_2.
    basis: list[PauliLabel] = []
    for label in normalized:
        if _binary_rank((*basis, label)) > len(basis):
            basis.append(label)
    dimension = len(basis)
    if dimension == 0:
        return 0

    restricted_form = tuple(
        tuple(symplectic_pairing(left, right) for right in basis) for left in basis
    )
    restricted_rank = _binary_rank(restricted_form)
    if restricted_rank % 2:
        raise AssertionError("an alternating bilinear form must have even rank")
    return dimension - restricted_rank // 2


def symplectic_pairing(left: Sequence[int], right: Sequence[int]) -> int:
    """Return the binary Pauli symplectic pairing of two ``(x | z)`` labels."""

    left_label = _validate_pauli_label(left)
    right_label = _validate_pauli_label(right)
    if len(left_label) != len(right_label):
        raise ValueError("Pauli labels must have the same length")
    num_qubits = len(left_label) // 2
    return sum(
        left_label[index] * right_label[num_qubits + index]
        + left_label[num_qubits + index] * right_label[index]
        for index in range(num_qubits)
    ) % 2


def find_lagrangian_witness(
    observations: Iterable[PauliConjugationObservation],
    *,
    num_qubits: int,
    pauli_probability_threshold: float,
) -> SemiCliffordWitness | None:
    """Find ``n`` independent commuting input/output Pauli pairs.

    Pairwise commutation and independence are checked on both sides.  The
    output checks are redundant in the exact oracle model but guard against an
    inconsistent finite-shot decoding.
    """

    if num_qubits < 1:
        raise ValueError("num_qubits must be positive")
    if not 0 < pauli_probability_threshold <= 1:
        raise ValueError("pauli_probability_threshold must lie in (0, 1]")
    eligible = sorted(
        (
            observation
            for observation in observations
            if observation.dominant_probability >= pauli_probability_threshold
            and any(observation.input_pauli)
        ),
        key=lambda observation: (-observation.dominant_probability, observation.input_pauli),
    )

    def search(
        start: int,
        selected: tuple[PauliConjugationObservation, ...],
    ) -> tuple[PauliConjugationObservation, ...] | None:
        if len(selected) == num_qubits:
            return selected
        if len(selected) + len(eligible) - start < num_qubits:
            return None
        for index in range(start, len(eligible)):
            candidate = eligible[index]
            if any(
                symplectic_pairing(candidate.input_pauli, previous.input_pauli)
                or symplectic_pairing(
                    candidate.dominant_output_pauli, previous.dominant_output_pauli
                )
                for previous in selected
            ):
                continue
            if _binary_rank(
                tuple(item.input_pauli for item in selected) + (candidate.input_pauli,)
            ) != len(selected) + 1:
                continue
            if _binary_rank(
                tuple(item.dominant_output_pauli for item in selected)
                + (candidate.dominant_output_pauli,)
            ) != len(selected) + 1:
                continue
            result = search(index + 1, selected + (candidate,))
            if result is not None:
                return result
        return None

    selected = search(0, ())
    if selected is None:
        return None
    return SemiCliffordWitness(
        input_basis=tuple(observation.input_pauli for observation in selected),
        output_basis=tuple(observation.dominant_output_pauli for observation in selected),
        minimum_dominant_probability=min(
            observation.dominant_probability for observation in selected
        ),
    )


def _validate_input_lagrangian_basis(
    input_basis: Iterable[Sequence[int]],
    *,
    num_qubits: int,
) -> tuple[PauliLabel, ...]:
    return lagrangian_from_basis(input_basis, num_qubits).basis


def run_semi_clifford_sampling_test(
    unitary: Any,
    *,
    unitary_dagger: Any | None = None,
    sampler: Any | None = None,
    shots: int = 1024,
    pauli_probability_threshold: float = 0.99,
    seed: int | None = None,
    input_paulis: Iterable[Sequence[int]] | None = None,
    batch_size: int | None = None,
) -> SemiCliffordSamplingResult:
    """Run Pauli conjugation circuits and search for a candidate witness.

    ``sampler`` may be any Qiskit SamplerV2 implementation.  If omitted, the
    local :class:`qiskit.primitives.StatevectorSampler` is used. Omitting
    ``input_paulis`` performs exhaustive discovery; passing selected labels is
    an incomplete search unless they comprise every nonidentity Pauli.
    """

    if shots < 1:
        raise ValueError("shots must be positive")
    if not 0 < pauli_probability_threshold <= 1:
        raise ValueError("pauli_probability_threshold must lie in (0, 1]")
    if batch_size is not None and batch_size < 1:
        raise ValueError("batch_size must be positive")
    experiments = build_semi_clifford_test_circuits(
        unitary,
        unitary_dagger=unitary_dagger,
        input_paulis=input_paulis,
    )
    num_qubits = len(experiments[0][0]) // 2
    exhaustive_labels = {
        label for label in all_pauli_labels(num_qubits) if any(label)
    }
    exhaustive_pauli_search = {label for label, _ in experiments} == exhaustive_labels
    if sampler is None:
        try:
            from qiskit.primitives import StatevectorSampler
        except ModuleNotFoundError as exc:  # pragma: no cover - optional install
            raise ModuleNotFoundError(
                'Qiskit support requires: pip install -e ".[qiskit]"'
            ) from exc
        sampler = StatevectorSampler(seed=seed)

    effective_batch_size = batch_size or len(experiments)
    observations_list: list[PauliConjugationObservation] = []
    for start in range(0, len(experiments), effective_batch_size):
        batch = experiments[start : start + effective_batch_size]
        job_result = sampler.run([circuit for _, circuit in batch], shots=shots).result()
        observations_list.extend(
            bell_counts_to_observation(label, publication_result.data.bell.get_counts())
            for (label, _), publication_result in zip(batch, job_result, strict=True)
        )
    observations = tuple(observations_list)
    witness = find_lagrangian_witness(
        observations,
        num_qubits=num_qubits,
        pauli_probability_threshold=pauli_probability_threshold,
    )
    return SemiCliffordSamplingResult(
        num_qubits=num_qubits,
        shots_per_pauli=shots,
        pauli_probability_threshold=pauli_probability_threshold,
        observations=observations,
        witness=witness,
        exhaustive_pauli_search=exhaustive_pauli_search,
    )


def run_semi_clifford_witness_test(
    unitary: Any,
    input_basis: Iterable[Sequence[int]],
    *,
    unitary_dagger: Any | None = None,
    sampler: Any | None = None,
    shots: int = 1024,
    pauli_probability_threshold: float = 0.99,
    seed: int | None = None,
) -> SemiCliffordSamplingResult:
    """Test a proposed input Lagrangian using only ``n`` quantum circuits."""

    instruction = _as_instruction(unitary, name="unitary")
    basis = _validate_input_lagrangian_basis(
        input_basis,
        num_qubits=instruction.num_qubits,
    )
    return run_semi_clifford_sampling_test(
        instruction,
        unitary_dagger=unitary_dagger,
        sampler=sampler,
        shots=shots,
        pauli_probability_threshold=pauli_probability_threshold,
        seed=seed,
        input_paulis=basis,
    )
