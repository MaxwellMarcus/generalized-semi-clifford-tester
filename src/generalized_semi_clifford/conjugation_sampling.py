"""Conjugation words without enumerating any conjugation group.

The finite-subgroup escape sampler has a distributional guarantee, not a
uniform-sampling guarantee. See ``docs/conjugation-sampling.md`` for the proof
and, especially, the distinction between group operations and queries to U.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from enum import Enum
from math import ceil, isfinite, log
from typing import Generic, Protocol, TypeVar

import numpy as np

from .lagrangian import PauliLabel

T = TypeVar("T")


class ConjugationBackend(Protocol[T]):
    """Faithful projective-unitary operations; ``exact`` is a caller contract.

    Products have matrix order ``left @ right``. Paulis use this package's
    Hermitian (x | z) convention. No group enumeration or membership is needed
    by the sampler. A Clifford membership oracle is supplied separately.
    """

    num_qubits: int
    seed: T
    exact: bool
    arithmetic: str

    def identity(self) -> T: ...
    def pauli(self, label: PauliLabel) -> T: ...
    def multiply(self, left: T, right: T) -> T: ...
    def adjoint(self, value: T) -> T: ...


@dataclass(frozen=True, eq=False)
class ConjugationWord:
    """A replayable, shared expression DAG, not an expanded circuit.

    ``seed_uses`` counts U and U-dagger occurrences after fully expanding the
    expression. Sharing storage does not make those physical queries free.
    Nodes deliberately use identity equality to avoid expanding shared trees.
    """

    kind: str
    children: tuple[ConjugationWord, ...] = ()
    label: PauliLabel | None = None
    seed_uses: int = 0


@dataclass(frozen=True)
class ConjugationWordNode:
    """One node in a child-first indexed conjugation-word plan."""

    kind: str
    children: tuple[int, ...] = ()
    label: PauliLabel | None = None
    seed_uses: int = 0


@dataclass(frozen=True)
class ConjugationWordPlan:
    """An immutable indexed DAG that stores each shared expression node once.

    The plan is a compact classical representation only. Its
    ``expanded_seed_uses`` property still counts every physical U/U-dagger
    occurrence in a literal execution of the expanded expression.
    """

    nodes: tuple[ConjugationWordNode, ...]
    root: int

    @property
    def stored_nodes(self) -> int:
        return len(self.nodes)

    @property
    def group_operations(self) -> int:
        return sum(node.kind in ("product", "adjoint") for node in self.nodes)

    @property
    def expanded_seed_uses(self) -> int:
        if not 0 <= self.root < len(self.nodes):
            raise ValueError("conjugation-word plan root is out of range")
        return self.nodes[self.root].seed_uses


@dataclass(frozen=True)
class ConjugationSeedFactor:
    """One sequential seed or projective-Pauli factor."""

    kind: str
    adjoint: bool = False
    label: PauliLabel | None = None


@dataclass(frozen=True)
class ConjugationSeedProgram:
    """A freely reduced sequential program for a conjugation word.

    Adjacent inverse seed calls are cancelled and adjacent Paulis are combined
    projectively, so dense replay can differ by a global phase. This is an
    executable resynthesis, not a shared DAG: its length is the number of
    retained sequential factors.
    """

    factors: tuple[ConjugationSeedFactor, ...]

    @property
    def seed_uses(self) -> int:
        return sum(factor.kind == "unitary" for factor in self.factors)

    @property
    def pauli_uses(self) -> int:
        return sum(factor.kind == "pauli" for factor in self.factors)


@dataclass(frozen=True)
class ConjugationSample(Generic[T]):
    value: T
    word: ConjugationWord
    level: int
    group_operations: int

    @property
    def expanded_seed_uses(self) -> int:
        return self.word.seed_uses


class SamplingLimitExceeded(RuntimeError):
    """A work cap was reached; this is not mathematical evidence either way."""


def _integer(value: int, name: str, minimum: int = 1) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")


class _Sampler(Generic[T]):
    def __init__(
        self, backend: ConjugationBackend[T], seed: int | None, max_operations: int
    ) -> None:
        _integer(backend.num_qubits, "num_qubits")
        _integer(max_operations, "max_operations")
        self.backend = backend
        self.rng = np.random.default_rng(seed)
        self.max_operations = max_operations
        self.operations = 0
        self.u = (backend.seed, ConjugationWord("unitary", seed_uses=1))

    def charge(self) -> None:
        if self.operations >= self.max_operations:
            raise SamplingLimitExceeded("group-operation budget exhausted")
        self.operations += 1

    def identity(self) -> tuple[T, ConjugationWord]:
        return self.backend.identity(), ConjugationWord("identity")

    def pauli(self, label: PauliLabel) -> tuple[T, ConjugationWord]:
        return self.backend.pauli(label), ConjugationWord("pauli", label=label)

    def random_pauli(self) -> tuple[T, ConjugationWord]:
        label = tuple(int(bit) for bit in self.rng.integers(2, size=2 * self.backend.num_qubits))
        return self.pauli(label)

    def generators(self) -> list[tuple[T, ConjugationWord]]:
        width = 2 * self.backend.num_qubits
        return [self.pauli(tuple(int(i == j) for i in range(width))) for j in range(width)]

    def multiply(
        self, left: tuple[T, ConjugationWord], right: tuple[T, ConjugationWord]
    ) -> tuple[T, ConjugationWord]:
        self.charge()
        return self.backend.multiply(left[0], right[0]), ConjugationWord(
            "product", (left[1], right[1]), seed_uses=left[1].seed_uses + right[1].seed_uses
        )

    def adjoint(self, item: tuple[T, ConjugationWord]) -> tuple[T, ConjugationWord]:
        self.charge()
        return self.backend.adjoint(item[0]), ConjugationWord(
            "adjoint", (item[1],), seed_uses=item[1].seed_uses
        )

    def conjugate(
        self, conjugator: tuple[T, ConjugationWord], item: tuple[T, ConjugationWord]
    ) -> tuple[T, ConjugationWord]:
        return self.multiply(self.multiply(conjugator, item), self.adjoint(conjugator))

    def subproduct(self, pool: Sequence[tuple[T, ConjugationWord]]) -> tuple[T, ConjugationWord]:
        result = self.identity()
        for item in pool:
            if self.rng.integers(2):
                result = self.multiply(result, item)
        return result

    def short_word(self, depth: int, length: int) -> tuple[T, ConjugationWord]:
        if depth == 1:
            # Gamma_1 = U P U-dagger: exactly uniform projectively.
            return self.conjugate(self.u, self.random_pauli())
        result = self.identity()
        for _ in range(length):
            previous = self.short_word(depth - 1, length)
            factor = self.conjugate(previous, self.random_pauli())
            result = self.multiply(result, factor)
        return result

    def escape_word(self, depth: int, bound: int) -> tuple[T, ConjugationWord]:
        """Escape any fixed H of order <= 2**bound with probability >= 1/4.

        Assumes Gamma_depth is not contained in H. H need not be known, and
        Gamma_depth need not be finite. The statement requires exact operations.
        """
        if depth == 1:
            return self.conjugate(self.u, self.random_pauli())
        if bound < 2 * self.backend.num_qubits:
            # Such an H cannot contain the whole Pauli group.
            return self.random_pauli()
        pool = self.generators()
        normalizer_bound = 2 * self.backend.num_qubits * (bound + 1)
        for _ in range(32 * (bound + 1)):
            previous = self.escape_word(depth - 1, normalizer_bound)
            candidate = self.conjugate(previous, self.subproduct(pool))
            pool.append(candidate)
        return self.subproduct(pool)


def _validate_depth(depth: int) -> None:
    _integer(depth, "depth")
    if depth > 8:
        raise ValueError("depth exceeds the prototype's recursion limit of 8")


def sample_conjugation_word(
    backend: ConjugationBackend[T],
    *,
    depth: int = 2,
    word_length: int = 4,
    seed: int | None = None,
    max_operations: int = 100_000,
) -> ConjugationSample[T]:
    """Sample a short valid Gamma_depth word, without a detection guarantee.

    Every higher-level word is a product of ``word_length`` independently
    generated conjugates. U-query cost is 2*(2*word_length)**(depth-1).
    For depth > 1 this is neither uniform nor a proven rapid-mixing sampler.
    """
    _validate_depth(depth)
    _integer(word_length, "word_length")
    sampler = _Sampler(backend, seed, max_operations)
    value, word = sampler.short_word(depth, word_length)
    return ConjugationSample(value, word, depth, sampler.operations)


def sample_subgroup_escape_word(
    backend: ConjugationBackend[T],
    *,
    depth: int,
    subgroup_log2_bound: int,
    seed: int | None = None,
    max_operations: int = 100_000,
) -> ConjugationSample[T]:
    """Sample with the finite-subgroup escape guarantee described above.

    Budget exhaustion raises SamplingLimitExceeded, never a successful sample.
    Polynomial *group-operation* count for fixed depth does not imply a
    polynomial number of physical U queries or efficient exact arithmetic.
    """
    _validate_depth(depth)
    _integer(subgroup_log2_bound, "subgroup_log2_bound", 0)
    sampler = _Sampler(backend, seed, max_operations)
    value, word = sampler.escape_word(depth, subgroup_log2_bound)
    return ConjugationSample(value, word, depth, sampler.operations)


def evaluate_conjugation_word(
    word: ConjugationWord,
    backend: ConjugationBackend[T],
    *,
    max_operations: int = 100_000,
) -> T:
    """Replay a word on another seed/backend, evaluating shared nodes once."""
    _integer(max_operations, "max_operations")
    values: dict[int, T] = {}
    pending = [(word, False)]
    operations = 0
    while pending:
        node, ready = pending.pop()
        key = id(node)
        if key in values:
            continue
        if not ready:
            pending.append((node, True))
            pending.extend((child, False) for child in reversed(node.children))
            continue
        if node.kind in ("product", "adjoint"):
            if operations >= max_operations:
                raise SamplingLimitExceeded("word-replay operation budget exhausted")
            operations += 1
        if node.kind == "identity":
            value = backend.identity()
        elif node.kind == "unitary":
            value = backend.seed
        elif node.kind == "pauli" and node.label is not None:
            value = backend.pauli(node.label)
        elif node.kind == "product" and len(node.children) == 2:
            value = backend.multiply(*(values[id(child)] for child in node.children))
        elif node.kind == "adjoint" and len(node.children) == 1:
            value = backend.adjoint(values[id(node.children[0])])
        else:
            raise ValueError("malformed conjugation-word node")
        values[key] = value
    return values[id(word)]


def compile_conjugation_word(word: ConjugationWord) -> ConjugationWordPlan:
    """Compile a shared expression DAG to stable child-first integer indices.

    Identity-based indexing preserves intentional sharing without attempting
    algebraic simplification or claiming any reduction in physical seed calls.
    """
    ordered: list[ConjugationWord] = []
    indices: dict[int, int] = {}
    pending = [(word, False)]
    while pending:
        node, ready = pending.pop()
        key = id(node)
        if key in indices:
            continue
        if not ready:
            pending.append((node, True))
            pending.extend((child, False) for child in reversed(node.children))
            continue
        indices[key] = len(ordered)
        ordered.append(node)
    nodes: list[ConjugationWordNode] = []
    for node in ordered:
        children = tuple(indices[id(child)] for child in node.children)
        if node.kind == "unitary" and not children and node.label is None:
            seed_uses = 1
        elif node.kind == "identity" and not children and node.label is None:
            seed_uses = 0
        elif node.kind == "pauli" and not children and node.label is not None:
            seed_uses = 0
        elif node.kind == "product" and len(children) == 2 and node.label is None:
            seed_uses = sum(nodes[child].seed_uses for child in children)
        elif node.kind == "adjoint" and len(children) == 1 and node.label is None:
            seed_uses = nodes[children[0]].seed_uses
        else:
            raise ValueError("malformed conjugation-word node")
        nodes.append(ConjugationWordNode(node.kind, children, node.label, seed_uses))
    return ConjugationWordPlan(tuple(nodes), indices[id(word)])


def evaluate_conjugation_plan(
    plan: ConjugationWordPlan,
    backend: ConjugationBackend[T],
    *,
    max_operations: int = 100_000,
) -> T:
    """Replay an indexed plan, rejecting malformed or non-topological nodes."""
    _integer(max_operations, "max_operations")
    if not 0 <= plan.root < len(plan.nodes):
        raise ValueError("conjugation-word plan root is out of range")
    values: list[T] = []
    operations = 0
    for index, node in enumerate(plan.nodes):
        if any(child < 0 or child >= index for child in node.children):
            raise ValueError("conjugation-word plan is not child-first")
        if node.kind in ("product", "adjoint"):
            if operations >= max_operations:
                raise SamplingLimitExceeded("word-replay operation budget exhausted")
            operations += 1
        if node.kind == "identity" and not node.children and node.label is None:
            value = backend.identity()
            seed_uses = 0
        elif node.kind == "unitary" and not node.children and node.label is None:
            value = backend.seed
            seed_uses = 1
        elif node.kind == "pauli" and not node.children and node.label is not None:
            value = backend.pauli(node.label)
            seed_uses = 0
        elif node.kind == "product" and len(node.children) == 2 and node.label is None:
            value = backend.multiply(*(values[child] for child in node.children))
            seed_uses = sum(plan.nodes[child].seed_uses for child in node.children)
        elif node.kind == "adjoint" and len(node.children) == 1 and node.label is None:
            value = backend.adjoint(values[node.children[0]])
            seed_uses = plan.nodes[node.children[0]].seed_uses
        else:
            raise ValueError("malformed conjugation-word plan node")
        if node.seed_uses != seed_uses:
            raise ValueError("conjugation-word plan has inconsistent seed-use accounting")
        values.append(value)
    return values[plan.root]


def resynthesize_conjugation_word(
    word: ConjugationWord,
    *,
    max_expanded_factors: int = 100_000,
) -> ConjugationSeedProgram:
    """Expand and freely reduce a word to sequential seed/Pauli factors.

    Products are flattened, adjoints are pushed to leaves, consecutive Paulis
    are multiplied by XORing their binary labels (global Pauli phases are
    projectively irrelevant), and adjacent ``U``/``U-dagger`` pairs cancel.
    The expansion cap is charged before cancellation so a compact result never
    conceals excessive resynthesis work.
    """
    _integer(max_expanded_factors, "max_expanded_factors")
    factors: list[ConjugationSeedFactor] = []
    pending = [(word, False)]
    expanded = 0

    def append(factor: ConjugationSeedFactor) -> None:
        if factor.kind == "unitary":
            if factors and factors[-1].kind == "unitary":
                if factors[-1].adjoint != factor.adjoint:
                    factors.pop()
                    return
            factors.append(factor)
            return
        if factor.kind != "pauli" or factor.adjoint or factor.label is None:
            raise ValueError("malformed conjugation seed factor")
        if factors and factors[-1].kind == "pauli":
            previous = factors.pop()
            assert previous.label is not None
            if len(previous.label) != len(factor.label):
                raise ValueError("inconsistent Pauli-label widths")
            label = tuple(
                left ^ right for left, right in zip(previous.label, factor.label, strict=True)
            )
            if any(label):
                factors.append(ConjugationSeedFactor("pauli", label=label))
            return
        if any(factor.label):
            factors.append(factor)

    while pending:
        node, take_adjoint = pending.pop()
        if node.kind == "product" and len(node.children) == 2 and node.label is None:
            left, right = node.children
            if take_adjoint:
                pending.extend(((left, True), (right, True)))
            else:
                pending.extend(((right, False), (left, False)))
            continue
        if node.kind == "adjoint" and len(node.children) == 1 and node.label is None:
            pending.append((node.children[0], not take_adjoint))
            continue
        if node.children:
            raise ValueError("malformed conjugation-word node")
        if node.kind == "identity" and node.label is None:
            continue
        if node.kind == "unitary" and node.label is None:
            factor = ConjugationSeedFactor("unitary", adjoint=take_adjoint)
        elif node.kind == "pauli" and node.label is not None:
            if not node.label or any(bit not in (0, 1) for bit in node.label):
                raise ValueError("Pauli labels must be nonempty binary tuples")
            factor = ConjugationSeedFactor("pauli", label=node.label)
        else:
            raise ValueError("malformed conjugation-word node")
        if expanded >= max_expanded_factors:
            raise SamplingLimitExceeded("seed-program resynthesis factor budget exhausted")
        expanded += 1
        append(factor)
    return ConjugationSeedProgram(tuple(factors))


def evaluate_conjugation_seed_program(
    program: ConjugationSeedProgram,
    backend: ConjugationBackend[T],
    *,
    max_operations: int = 100_000,
) -> T:
    """Execute a resynthesized program from left to right, projectively."""
    _integer(max_operations, "max_operations")
    value = backend.identity()
    operations = 0
    for factor in program.factors:
        if factor.kind == "unitary" and factor.label is None:
            item = backend.adjoint(backend.seed) if factor.adjoint else backend.seed
        elif factor.kind == "pauli" and not factor.adjoint and factor.label is not None:
            if len(factor.label) != 2 * backend.num_qubits:
                raise ValueError("seed-program Pauli width does not match backend")
            item = backend.pauli(factor.label)
        else:
            raise ValueError("malformed conjugation seed factor")
        if operations >= max_operations:
            raise SamplingLimitExceeded("seed-program replay operation budget exhausted")
        value = backend.multiply(value, item)
        operations += 1
    return value


class ContainmentStatus(str, Enum):
    SUPPORTED = "containment_supported"
    VIOLATION = "containment_violation"
    INCONCLUSIVE = "inconclusive"


@dataclass(frozen=True)
class ConjugationContainmentResult(Generic[T]):
    status: ContainmentStatus
    depth: int
    arithmetic: str
    rounds_completed: int
    rounds_required: int
    group_operations: int
    oracle_calls: int
    ideal_false_accept_bound: float | None
    rigorous_false_accept_bound: float | None
    witness: ConjugationSample[T] | None
    message: str


def test_clifford_conjugation_containment(
    backend: ConjugationBackend[T],
    clifford_oracle: Callable[[T], bool | None],
    *,
    depth: int,
    failure_probability: float = 0.01,
    oracle_is_exact: bool = False,
    seed: int | None = None,
    max_operations: int = 100_000,
) -> ConjugationContainmentResult[T]:
    """Test Gamma_depth <= C_2 without accessing any full preceding group.

    In exact arithmetic with an exact Clifford oracle, this has perfect
    completeness and false acceptance at most ``failure_probability``. A
    partial run or an indeterminate oracle returns INCONCLUSIVE, with no
    acceptance probability claim. Numerical backends receive only the
    explicitly conditional ``ideal_false_accept_bound``.

    Applying a GSC theorem, a hierarchy promise, or a distance interpretation
    is deliberately left to the caller; a violation is not a non-GSC verdict.
    """
    _validate_depth(depth)
    if not isfinite(failure_probability) or not 0 < failure_probability < 1:
        raise ValueError("failure_probability must be finite and strictly between 0 and 1")
    sampler = _Sampler(backend, seed, max_operations)
    n = backend.num_qubits
    bound = 2 * n * n + 3 * n  # |projective C_2| < 2**bound
    growth_probability = 0.25 if depth == 2 else 0.125
    rounds = (
        2 * n
        if depth == 1
        else ceil(max(2 * (bound + 1), 8 * -log(failure_probability)) / growth_probability)
    )
    completed = 0
    calls = 0

    def result(
        status: ContainmentStatus,
        message: str,
        witness: ConjugationSample[T] | None = None,
    ) -> ConjugationContainmentResult[T]:
        ideal = (
            (0.0 if depth == 1 else failure_probability)
            if (status is ContainmentStatus.SUPPORTED)
            else None
        )
        rigorous = ideal if backend.exact and oracle_is_exact else None
        return ConjugationContainmentResult(
            status,
            depth,
            backend.arithmetic,
            completed,
            rounds,
            sampler.operations,
            calls,
            ideal,
            rigorous,
            witness,
            message,
        )

    pool = sampler.generators()
    try:
        for index in range(rounds):
            if depth == 1:
                candidate = sampler.conjugate(sampler.u, pool[index])
            else:
                previous = sampler.escape_word(depth - 1, 2 * n * (bound + 1))
                candidate = sampler.conjugate(previous, sampler.subproduct(pool))
            calls += 1
            decision = clifford_oracle(candidate[0])
            if decision is None:
                return result(ContainmentStatus.INCONCLUSIVE, "Clifford oracle was indeterminate")
            if not isinstance(decision, (bool, np.bool_)):
                raise TypeError("Clifford oracle must return bool or None")
            completed += 1
            if not decision:
                witness = ConjugationSample(candidate[0], candidate[1], depth, sampler.operations)
                return result(
                    ContainmentStatus.VIOLATION,
                    "sampled a non-Clifford word; this is not a non-GSC conclusion",
                    witness,
                )
            if depth > 1:
                pool.append(candidate)
    except SamplingLimitExceeded as exc:
        return result(ContainmentStatus.INCONCLUSIVE, str(exc))
    return result(
        ContainmentStatus.SUPPORTED,
        "completed the containment test; probability guarantee requires exact operations "
        "and an exact Clifford oracle",
    )
