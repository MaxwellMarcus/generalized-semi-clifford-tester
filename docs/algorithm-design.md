# Algorithm design and correctness contract

## Goal

The eventual public API should accept a precisely specified representation of a
unitary (or of its induced algebraic data) and return one of three outcomes:

- `GSC`, together with a certificate that can be verified independently;
- `NOT_GSC`, together with a mathematical obstruction or finite certificate;
- `UNKNOWN`, when the input lies outside the proved scope or a configured search
  limit is reached.

The third outcome is important. An incomplete search must never be reported as a
negative decision.

## Correctness requirements before a tester is released

1. State the exact definition of generalized semi-Clifford used by the package.
2. Cite or prove the equivalence between that definition and the implemented
   decision criterion for the supported input class.
3. Specify whether arithmetic is exact. If a numerical front end is added, make
   its tolerance and failure mode visible in the result.
4. Make every positive or negative answer independently checkable.
5. Include analytically derived fixtures, edge cases, and known counterexamples.

## Proposed package boundaries

- `symplectic.py`: exact binary linear algebra and symplectic validation.
- `pauli.py`: phase-aware Pauli representations and conjugation actions.
- `invariants.py`: isotropic/Lagrangian subspace enumeration or search.
- `certificates.py`: serializable witnesses and independent verification.
- `tester.py`: orchestration only after the correctness contract is fixed.

## Performance plan

Start with transparent exhaustive algorithms at small qubit counts. Add optimized
backends only after regression fixtures exist. Benchmarks should report qubit
count, representation size, search limits, peak memory, and whether a result was
exact or heuristic.
