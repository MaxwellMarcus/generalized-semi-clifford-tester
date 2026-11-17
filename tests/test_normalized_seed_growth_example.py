import json
import runpy
from pathlib import Path


def test_normalized_seed_growth_example_records_exact_parameter_grid(capsys) -> None:
    example = Path(__file__).parents[1] / "examples" / "normalized_seed_growth.py"
    baseline = Path(__file__).parents[1] / "benchmarks" / "normalized-seed-growth-grid.json"

    runpy.run_path(str(example), run_name="__main__")

    report = json.loads(capsys.readouterr().out)
    assert report == json.loads(baseline.read_text())
    assert [case["max_normalized_program_length"] for case in report["maxima"]] == [
        3,
        19,
        123,
        556,
        4267,
        27051,
        144983,
    ]
    assert all(
        case["max_constructed_rope_nodes"]
        >= case["constructed_rope_nodes_at_first_length_maximizer"]
        for case in report["maxima"]
    )
    assert all(
        case["max_retained_rope_nodes"] <= case["max_constructed_rope_nodes"]
        and case["retained_rope_nodes_at_first_length_maximizer"]
        <= case["constructed_rope_nodes_at_first_length_maximizer"]
        for case in report["maxima"]
    )
