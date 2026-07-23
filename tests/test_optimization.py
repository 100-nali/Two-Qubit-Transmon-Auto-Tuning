"""End-to-end test: Bayesian optimization actually finds a high-fidelity X gate."""

import numpy as np
import pytest

from two_qubit_transmon_auto_tuning import (
    Experiment,
    Qubit,
    x_gate_objective,
    run_bayesian_optimization,
)


@pytest.mark.slow
def test_bayesian_optimization_tunes_x_gate():
    q1 = Qubit(i=1, w_q=6 * 2 * np.pi, a_q=-0.3 * 2 * np.pi, r=0.01 * 2 * np.pi,
               w_d=6 * 2 * np.pi, w_ex12=5.3 * 2 * np.pi)
    q2 = Qubit(i=2, w_q=5 * 2 * np.pi, a_q=-0.3 * 2 * np.pi, r=0.01 * 2 * np.pi,
               w_d=5 * 2 * np.pi, w_ex12=5 * 2 * np.pi)
    experiment = Experiment(
        qubits=[q1, q2], g=0.005 * 2 * np.pi, t_exp=81,
        gate_type="X", drive_shape="Gaussian",
    )

    objective, pbounds = x_gate_objective(experiment)
    optimizer = run_bayesian_optimization(
        objective, pbounds,
        init_points=10, n_iter=20,
        acquisition_weight=5.0, alpha=1e-3,
        random_state=1,
    )
    assert optimizer.max["target"] > 0.9
