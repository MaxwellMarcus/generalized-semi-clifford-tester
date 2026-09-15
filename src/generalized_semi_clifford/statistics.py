"""Finite-shot confidence bounds for sampled Pauli leakage."""

from __future__ import annotations

from math import ceil, log


def binomial_proportion_upper_bound(
    observed_events: int,
    trials: int,
    *,
    confidence_level: float,
) -> float:
    """Return the one-sided Clopper--Pearson upper confidence bound.

    ``observed_events`` is the number of Bell outcomes outside a proposed
    output Lagrangian.  The returned value is an exact binomial upper bound,
    not a normal approximation.
    """

    if trials <= 0:
        raise ValueError("trials must be positive")
    if not 0 <= observed_events <= trials:
        raise ValueError("observed_events must lie between zero and trials")
    if not 0 < confidence_level < 1:
        raise ValueError("confidence_level must lie strictly between zero and one")
    if observed_events == trials:
        return 1.0

    try:
        from scipy.stats import beta
    except ModuleNotFoundError as exc:  # pragma: no cover - optional dependency guard
        raise ModuleNotFoundError(
            'confidence bounds require: pip install -e ".[qiskit]"'
        ) from exc

    failure_probability = 1.0 - confidence_level
    return float(
        beta.ppf(
            1.0 - failure_probability,
            observed_events + 1,
            trials - observed_events,
        )
    )


def zero_event_shots_required(
    leakage_threshold: float,
    *,
    confidence_level: float,
    simultaneous_tests: int = 1,
) -> int:
    """Return shots needed to certify a threshold after zero leakage events.

    The calculation inverts the exact one-sided Clopper--Pearson bound and
    allocates the failure probability equally across ``simultaneous_tests``.
    It is therefore the best-case shot count: observing any outside-support
    event can require more samples or prevent certification at the target.
    """

    if not 0 < leakage_threshold < 1:
        raise ValueError("leakage_threshold must lie strictly between zero and one")
    if not 0 < confidence_level < 1:
        raise ValueError("confidence_level must lie strictly between zero and one")
    if simultaneous_tests < 1:
        raise ValueError("simultaneous_tests must be positive")
    per_test_failure_probability = (1.0 - confidence_level) / simultaneous_tests
    return ceil(log(per_test_failure_probability) / log(1.0 - leakage_threshold))
