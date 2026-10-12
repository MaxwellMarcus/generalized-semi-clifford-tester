"""Machine-readable command-line interface for dense and circuit GSC tests."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np

from .gsc_qiskit import (
    GSCDiscoveryVerificationResult,
    GSCSamplingResult,
    run_gsc_discovery_then_verification,
)
from .lagrangian import Lagrangian
from .tester import NaiveGSCResult, check_gsc_naive

DENSE_RESULT_SCHEMA = "generalized-semi-clifford/dense-result-v1"
CIRCUIT_RESULT_SCHEMA = "generalized-semi-clifford/circuit-result-v1"
EXACT_MATRIX_SCHEMA = "generalized-semi-clifford/exact-cyclotomic-matrix-v1"


def _lagrangian_payload(lagrangian: Lagrangian) -> dict[str, object]:
    return {
        "basis": [list(label) for label in lagrangian.basis],
        "elements": [list(label) for label in lagrangian.elements],
        "num_qubits": lagrangian.num_qubits,
    }


def dense_result_payload(
    result: NaiveGSCResult,
    *,
    unitary_tolerance: float,
    max_qubits: int,
) -> dict[str, object]:
    """Return a versioned JSON-compatible dense-test record."""

    witness = None
    if result.witness is not None:
        witness = {
            "input_lagrangian": _lagrangian_payload(result.witness.input_lagrangian),
            "output_lagrangian": _lagrangian_payload(result.witness.output_lagrangian),
            "max_leakage": result.witness.max_leakage,
        }
    return {
        "schema": DENSE_RESULT_SCHEMA,
        "method": "bounded-exhaustive-dense-search",
        "claim_scope": "complex128 numerical criterion; not an exact symbolic proof",
        "result": {
            "status": result.status.value,
            "is_gsc": result.is_gsc,
            "message": result.message,
            "num_qubits": result.num_qubits,
            "arithmetic": result.arithmetic,
            "best_leakage": result.best_leakage,
            "witness": witness,
        },
        "tolerances": {
            "pauli_coefficient": result.tolerance,
            "unitarity": unitary_tolerance,
        },
        "search": {
            "max_qubits": max_qubits,
            "lagrangians_checked": result.lagrangians_checked,
            "candidate_pairs_checked": result.candidate_pairs_checked,
            "complete": result.is_gsc is not None,
        },
    }


def _sampling_payload(result: GSCSamplingResult) -> dict[str, object]:
    witness = None
    if result.witness is not None:
        witness = {
            "input_lagrangian": _lagrangian_payload(result.witness.input_lagrangian),
            "output_lagrangian": _lagrangian_payload(result.witness.output_lagrangian),
            "maximum_empirical_leakage": result.witness.maximum_empirical_leakage,
        }
    return {
        "num_qubits": result.num_qubits,
        "shots_per_pauli": result.shots_per_pauli,
        "leakage_threshold": result.leakage_threshold,
        "search_mode": result.search_mode,
        "has_candidate_witness": result.has_candidate_witness,
        "has_confidence_certified_witness": result.has_confidence_certified_witness,
        "best_empirical_leakage": result.best_empirical_leakage,
        "confidence_level": result.confidence_level,
        "maximum_leakage_upper_bound": result.maximum_leakage_upper_bound,
        "input_lagrangians_checked": result.input_lagrangians_checked,
        "candidate_pairs_checked": result.candidate_pairs_checked,
        "witness": witness,
        "observations": [
            {
                "input_pauli": list(observation.input_pauli),
                "dominant_output_pauli": list(observation.dominant_output_pauli),
                "dominant_probability": observation.dominant_probability,
                "shots": observation.shots,
                "counts": dict(observation.counts),
            }
            for observation in result.observations
        ],
    }


def circuit_result_payload(result: GSCDiscoveryVerificationResult) -> dict[str, object]:
    """Return finite-shot evidence without promoting it to exact membership."""

    return {
        "schema": CIRCUIT_RESULT_SCHEMA,
        "method": "bell-sampling-discovery-then-fresh-verification",
        "claim_scope": (
            "finite-shot candidate evidence; confidence certification applies only to "
            "the independently verified fixed witness"
        ),
        "has_discovered_candidate": result.has_discovered_candidate,
        "has_confidence_certified_witness": result.has_confidence_certified_witness,
        "shot_cost": {
            "discovery": result.discovery_shot_cost,
            "verification": result.verification_shot_cost,
        },
        "discovery": _sampling_payload(result.discovery),
        "verification": (
            None if result.verification is None else _sampling_payload(result.verification)
        ),
    }


def _load_qasm(path: Path) -> Any:
    try:
        from qiskit import qasm2
    except ModuleNotFoundError as exc:  # pragma: no cover - optional dependency guard
        raise ModuleNotFoundError(
            'circuit tests require: pip install -e ".[qiskit]"'
        ) from exc
    return qasm2.load(path)


def _load_exact_cyclotomic_matrix(path: Path) -> Any:
    """Load an exact matrix without parsing arbitrary symbolic expressions."""

    from .exact_cyclotomic import CyclotomicField, CyclotomicMatrix

    payload = json.loads(path.read_text())
    if not isinstance(payload, dict) or payload.get("schema") != EXACT_MATRIX_SCHEMA:
        raise ValueError(f"exact matrix schema must be {EXACT_MATRIX_SCHEMA!r}")
    field_order = payload.get("field_order")
    if not isinstance(field_order, int) or isinstance(field_order, bool):
        raise TypeError("field_order must be an integer")
    rows = payload.get("rows")
    if not isinstance(rows, list) or not rows:
        raise ValueError("rows must be a nonempty list")

    field = CyclotomicField(field_order)

    def scalar(coefficients: object) -> Any:
        if not isinstance(coefficients, list) or not coefficients:
            raise ValueError("each matrix entry must be a nonempty coefficient list")
        result = field.zero
        for power, coefficient in enumerate(coefficients):
            if isinstance(coefficient, bool) or not isinstance(coefficient, (int, str)):
                raise TypeError("cyclotomic coefficients must be integers or rational strings")
            try:
                rational = Fraction(coefficient)
            except (ValueError, ZeroDivisionError) as error:
                raise ValueError(f"invalid rational coefficient {coefficient!r}") from error
            result += field.scalar(rational) * field.root_power(power)
        return result

    parsed_rows = []
    for row in rows:
        if not isinstance(row, list) or not row:
            raise ValueError("each matrix row must be a nonempty list")
        parsed_rows.append(tuple(scalar(entry) for entry in row))
    return CyclotomicMatrix.from_rows(field, parsed_rows)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="gsc-test",
        description="Emit versioned JSON for bounded generalized semi-Clifford tests.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    dense = subparsers.add_parser("dense", help="run the bounded dense numerical tester")
    dense.add_argument("matrix", type=Path, help="NumPy .npy file containing a unitary")
    dense.add_argument("--tolerance", type=float, default=1e-9)
    dense.add_argument("--unitary-tolerance", type=float, default=1e-9)
    dense.add_argument("--max-qubits", type=int, default=3)

    circuit = subparsers.add_parser(
        "circuit",
        help="discover and independently verify a candidate from an OpenQASM 2 circuit",
    )
    circuit.add_argument("qasm", type=Path, help="OpenQASM 2 circuit file")
    circuit.add_argument("--discovery-shots", type=int, default=1024)
    circuit.add_argument("--verification-shots", type=int, default=1024)
    circuit.add_argument("--leakage-threshold", type=float, default=0.01)
    circuit.add_argument("--confidence-level", type=float, default=0.95)
    circuit.add_argument("--discovery-seed", type=int)
    circuit.add_argument("--verification-seed", type=int)
    circuit.add_argument("--batch-size", type=int)

    exact_prime = subparsers.add_parser(
        "exact-prime",
        help="run a bounded exact odd-prime Weyl-Lagrangian search",
    )
    exact_prime.add_argument(
        "matrix",
        type=Path,
        help="JSON exact cyclotomic matrix in the version-one coefficient schema",
    )
    exact_prime.add_argument("--prime", type=int, required=True)
    exact_prime.add_argument("--max-qudits", type=int, default=2)
    exact_prime.add_argument("--max-coefficients", type=int, default=100_000)
    exact_prime.add_argument("--max-candidate-pairs", type=int, default=100_000)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the selected test and print one deterministic JSON document."""

    args = _parser().parse_args(argv)
    if args.command == "dense":
        matrix = np.load(args.matrix, allow_pickle=False)
        result = check_gsc_naive(
            matrix,
            tolerance=args.tolerance,
            unitary_tolerance=args.unitary_tolerance,
            max_qubits=args.max_qubits,
        )
        payload = dense_result_payload(
            result,
            unitary_tolerance=args.unitary_tolerance,
            max_qubits=args.max_qubits,
        )
    elif args.command == "circuit":
        circuit = _load_qasm(args.qasm)
        result = run_gsc_discovery_then_verification(
            circuit,
            discovery_shots=args.discovery_shots,
            verification_shots=args.verification_shots,
            leakage_threshold=args.leakage_threshold,
            confidence_level=args.confidence_level,
            discovery_seed=args.discovery_seed,
            verification_seed=args.verification_seed,
            discovery_batch_size=args.batch_size,
        )
        payload = circuit_result_payload(result)
    else:
        from .exact_odd_prime import search_exact_prime_lagrangian_witness

        matrix = _load_exact_cyclotomic_matrix(args.matrix)
        result = search_exact_prime_lagrangian_witness(
            matrix,
            args.prime,
            max_qudits=args.max_qudits,
            max_coefficients=args.max_coefficients,
            max_candidate_pairs=args.max_candidate_pairs,
        )
        payload = result.to_dict()
    print(json.dumps(payload, allow_nan=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":  # pragma: no cover - console-script entry point
    raise SystemExit(main())
