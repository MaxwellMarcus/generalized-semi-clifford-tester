# Contributing

This project is at an early research stage. Mathematical corrections, test
cases, and small implementation changes are welcome.

## Setup

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
python -m ruff check src tests examples
python -m pytest
```

Any proposed decision procedure must document its assumptions, cite the result
that establishes correctness, and include positive and negative test cases. If
an implementation uses numerical tolerances, label it as a heuristic rather
than an exact decision procedure.
