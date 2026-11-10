import json
import runpy
from pathlib import Path


def test_exact_cyclotomic_containment_example_preserves_claim_boundaries(capsys) -> None:
    example = Path(__file__).parents[1] / "examples" / "exact_cyclotomic_containment.py"

    runpy.run_path(str(example), run_name="__main__")

    records = {item["name"]: item for item in json.loads(capsys.readouterr().out)}
    supported = records["t_gate_supported"]
    violating = records["eighth_turn_violation"]
    inconclusive = records["coefficient_capped_oracle"]
    two_qubit_inconclusive = records["two_qubit_coefficient_capped_oracle"]

    assert supported["status"] == "containment_supported"
    assert supported["rigorous_false_accept_bound"] == 0.0
    assert supported["arithmetic"] == "cyclotomic_exact"
    assert supported["exact_oracle_coefficients_computed"] is None
    assert not supported["has_violation_witness"]

    assert violating["status"] == "containment_violation"
    assert violating["rigorous_false_accept_bound"] is None
    assert violating["has_violation_witness"]

    assert inconclusive["status"] == "inconclusive"
    assert inconclusive["rigorous_false_accept_bound"] is None
    assert inconclusive["oracle_calls"] == 1
    assert inconclusive["exact_oracle_coefficients_computed"] == 4
    assert inconclusive["oracle_stop_reason"] == "coefficient_cap"
    assert not inconclusive["has_violation_witness"]

    assert two_qubit_inconclusive["status"] == "inconclusive"
    assert two_qubit_inconclusive["rigorous_false_accept_bound"] is None
    assert two_qubit_inconclusive["rounds_completed"] == 0
    assert two_qubit_inconclusive["oracle_calls"] == 1
    assert two_qubit_inconclusive["exact_oracle_coefficients_computed"] == 16
    assert two_qubit_inconclusive["oracle_stop_reason"] == "coefficient_cap"
    assert not two_qubit_inconclusive["has_violation_witness"]
