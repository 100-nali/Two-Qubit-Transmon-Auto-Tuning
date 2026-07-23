"""Auto-tune a single-qubit X gate with Bayesian optimization.

Starting from no knowledge of the control amplitudes, the optimizer reaches
~0.99 process fidelity in ~40 evaluations -- and discovers the physics on its
own: it drives a single quadrature on qubit 1 and leaves the others near zero.

Run:
    python examples/optimize_x_gate.py
"""

import numpy as np

from two_qubit_transmon_auto_tuning import (
    Experiment,
    Qubit,
    x_gate_objective,
    run_bayesian_optimization,
)


def main():
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
        init_points=15, n_iter=25,
        acquisition_weight=5.0, alpha=1e-3,
    )

    print(f"Best X-gate process fidelity: {optimizer.max['target']:.4f}")
    print("Best drive amplitudes:",
          {k: round(v, 2) for k, v in optimizer.max["params"].items()})


if __name__ == "__main__":
    main()
