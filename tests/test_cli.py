from __future__ import annotations

import json

import numpy as np

from generalized_semi_clifford.cli import (
    CIRCUIT_RESULT_SCHEMA,
    DENSE_RESULT_SCHEMA,
    circuit_result_payload,
    dense_result_payload,
    main,
)
from generalized_semi_clifford.gsc_qiskit import (
    GSCDiscoveryVerificationResult,
    GSCSamplingResult,
)
from generalized_semi_clifford.semi_clifford_qiskit import PauliConjugationObservation
from generalized_semi_clifford.tester import check_gsc_naive


def test_dense_payload_preserves_unknown_as_null() -> None:
    result = check_gsc_naive(np.eye(16), max_qubits=3)
    payload = dense_result_payload(result, unitary_tolerance=1e-8, max_qubits=3)

    assert payload["schema"] == DENSE_RESULT_SCHEMA
    assert payload["result"]["status"] == "unknown"
    assert payload["result"]["is_gsc"] is None
    assert payload["result"]["best_leakage"] is None
    assert payload["search"]["complete"] is False
    assert payload["tolerances"] == {
        "pauli_coefficient": 1e-9,
        "unitarity": 1e-8,
    }


def test_dense_cli_emits_deterministic_json(tmp_path, capsys) -> None:
    matrix_path = tmp_path / "identity.npy"
    np.save(matrix_path, np.eye(2))

    assert main(["dense", str(matrix_path)]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["schema"] == DENSE_RESULT_SCHEMA
    assert payload["result"]["status"] == "gsc"
    assert payload["result"]["witness"]["max_leakage"] == 0.0


def test_circuit_payload_keeps_counts_and_absent_verification() -> None:
    observation = PauliConjugationObservation(
        input_pauli=(1, 0),
        dominant_output_pauli=(1, 0),
        dominant_probability=0.75,
        shots=4,
        counts=(("01", 3), ("00", 1)),
    )
    discovery = GSCSamplingResult(
        num_qubits=1,
        shots_per_pauli=4,
        leakage_threshold=0.01,
        observations=(observation,),
        witness=None,
        best_empirical_leakage=0.25,
        input_lagrangians_checked=3,
        candidate_pairs_checked=9,
        search_mode="support-aware-lagrangian-discovery",
    )
    result = GSCDiscoveryVerificationResult(discovery, None, 4, 0)

    payload = circuit_result_payload(result)
    assert payload["schema"] == CIRCUIT_RESULT_SCHEMA
    assert payload["has_discovered_candidate"] is False
    assert payload["has_confidence_certified_witness"] is False
    assert payload["verification"] is None
    assert payload["discovery"]["observations"][0]["counts"] == {"00": 1, "01": 3}
