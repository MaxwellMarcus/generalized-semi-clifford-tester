import json
import runpy
from pathlib import Path


def test_normalized_seed_growth_example_records_exact_depth_prefix(capsys) -> None:
    example = Path(__file__).parents[1] / "examples" / "normalized_seed_growth.py"

    runpy.run_path(str(example), run_name="__main__")

    report = json.loads(capsys.readouterr().out)
    assert report == {
        "arithmetic": "exact dihedral group coordinates",
        "cases": [
            {
                "abstract_operations": operations,
                "depth": depth,
                "expanded_seed_calls": expanded,
                "normalized_program_length": program_length,
                "normalized_seed_calls": normalized_seed_calls,
                "stored_dag_nodes": nodes,
            }
            for depth, nodes, operations, program_length, normalized_seed_calls, expanded in (
                (1, 5, 3, 3, 2, 2),
                (2, 20, 14, 8, 4, 8),
                (3, 50, 36, 15, 8, 32),
                (4, 110, 80, 75, 38, 128),
                (5, 230, 168, 376, 188, 512),
                (6, 470, 344, 751, 376, 2048),
                (7, 950, 696, 895, 448, 8192),
            )
        ],
        "schema_version": 1,
        "seed": 19,
        "word_length": 2,
    }
