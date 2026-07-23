"""Shared pytest fixtures."""

import pytest

from two_qubit_transmon_auto_tuning import default_cphase_experiment


@pytest.fixture
def cphase_experiment():
    """The two-transmon CPHASE experiment used across the test suite."""
    return default_cphase_experiment(t_exp=1000)
