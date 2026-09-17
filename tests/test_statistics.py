import pytest

from generalized_semi_clifford import (
    binomial_proportion_upper_bound,
    zero_event_shots_required,
)


def test_all_observations_leaking_has_unit_upper_bound() -> None:
    assert binomial_proportion_upper_bound(12, 12, confidence_level=0.95) == 1.0


@pytest.mark.parametrize("trials", [0, -1])
def test_binomial_bound_requires_positive_trials(trials: int) -> None:
    with pytest.raises(ValueError, match="trials must be positive"):
        binomial_proportion_upper_bound(0, trials, confidence_level=0.95)


@pytest.mark.parametrize("observed_events,trials", [(-1, 10), (11, 10)])
def test_binomial_bound_rejects_event_counts_outside_trials(
    observed_events: int,
    trials: int,
) -> None:
    with pytest.raises(ValueError, match="between zero and trials"):
        binomial_proportion_upper_bound(
            observed_events,
            trials,
            confidence_level=0.95,
        )


@pytest.mark.parametrize("confidence_level", [0.0, 1.0, float("nan")])
def test_binomial_bound_requires_open_unit_confidence(confidence_level: float) -> None:
    with pytest.raises(ValueError, match="strictly between zero and one"):
        binomial_proportion_upper_bound(0, 10, confidence_level=confidence_level)


@pytest.mark.parametrize("leakage_threshold", [0.0, 1.0, float("nan")])
def test_shot_planner_requires_open_unit_leakage(leakage_threshold: float) -> None:
    with pytest.raises(ValueError, match="leakage_threshold"):
        zero_event_shots_required(leakage_threshold, confidence_level=0.95)


@pytest.mark.parametrize("confidence_level", [0.0, 1.0, float("nan")])
def test_shot_planner_requires_open_unit_confidence(confidence_level: float) -> None:
    with pytest.raises(ValueError, match="confidence_level"):
        zero_event_shots_required(0.01, confidence_level=confidence_level)


@pytest.mark.parametrize("simultaneous_tests", [0, -1])
def test_shot_planner_requires_positive_simultaneous_tests(
    simultaneous_tests: int,
) -> None:
    with pytest.raises(ValueError, match="simultaneous_tests must be positive"):
        zero_event_shots_required(
            0.01,
            confidence_level=0.95,
            simultaneous_tests=simultaneous_tests,
        )
