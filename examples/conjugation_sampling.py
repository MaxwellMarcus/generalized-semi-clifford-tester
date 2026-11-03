"""Replayable short words and honest small-system sandwich-test outcomes."""

import numpy as np

from generalized_semi_clifford import (
    DenseConjugationBackend,
    check_gsc_sandwich,
    compile_conjugation_word,
    evaluate_conjugation_plan,
    evaluate_conjugation_seed_program,
    evaluate_conjugation_word,
    probe_gsc_conjugation_condition,
    resynthesize_conjugation_word,
    sample_conjugation_word,
)


def main() -> None:
    phase = np.diag([1, np.exp(1j * np.pi / 16)])
    backend = DenseConjugationBackend(phase)
    sample = sample_conjugation_word(backend, depth=2, word_length=4, seed=7)
    assert np.allclose(evaluate_conjugation_word(sample.word, backend), sample.value)
    plan = compile_conjugation_word(sample.word)
    assert np.allclose(evaluate_conjugation_plan(plan, backend), sample.value)
    program = resynthesize_conjugation_word(sample.word)
    resynthesized = evaluate_conjugation_seed_program(program, backend)
    overlap = np.vdot(sample.value, resynthesized)
    assert np.allclose(resynthesized, overlap / abs(overlap) * sample.value)
    print(f"Gamma_2 short word: {sample.expanded_seed_uses} expanded U/U-dagger uses")
    print(f"Indexed DAG: {plan.stored_nodes} stored nodes, {plan.group_operations} operations")
    print(f"Freely reduced program: {program.seed_uses} sequential U/U-dagger uses")
    for route in ("shallow", "deep"):
        probe = probe_gsc_conjugation_condition(
            phase, hierarchy_level=5, route=route, samples=8, seed=7
        )
        print(f"{route}: {probe.words_checked} words, distance evidence: {probe.evidence}")
        print(probe.message)

    x = np.array([[0, 1], [1, 0]])
    y = np.array([[0, -1j], [1j, 0]])
    z = np.diag([1, -1])
    rotation = np.cos(0.37) * np.eye(2) - 1j * np.sin(0.37) * (x + y + z) / np.sqrt(3)
    for name, gate in (("C5 phase", phase), ("generic rotation", rotation)):
        result = check_gsc_sandwich(gate, hierarchy_level=5, samples=8, seed=7)
        print(
            f"{name}: {result.status.value}; "
            f"SC distance bound={result.distance_from_sc_lower_bound}"
        )
        print(result.message)
    print("All dense results are floating-point evidence, not interval-certified proofs.")


if __name__ == "__main__":
    main()
