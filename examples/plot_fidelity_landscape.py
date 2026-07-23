"""Plot the X-gate fidelity landscape over two drive amplitudes.

This is why the problem needs Bayesian optimization rather than a gradient step
or a coarse grid sweep: the landscape is non-convex and periodic, with several
ridges and local optima. Sweeping a 2D slice (in-phase vs quadrature amplitude
on qubit 1) makes that structure visible.

Run:
    python examples/plot_fidelity_landscape.py
"""

import numpy as np
import matplotlib.pyplot as plt

from two_qubit_transmon_auto_tuning import Experiment, Qubit, x_gate_objective

GRID = 20  # points per axis (GRID**2 fidelity evaluations)


def main():
    q1 = Qubit(i=1, w_q=6 * 2 * np.pi, a_q=-0.3 * 2 * np.pi, r=0.01 * 2 * np.pi,
               w_d=6 * 2 * np.pi, w_ex12=5.3 * 2 * np.pi)
    q2 = Qubit(i=2, w_q=5 * 2 * np.pi, a_q=-0.3 * 2 * np.pi, r=0.01 * 2 * np.pi,
               w_d=5 * 2 * np.pi, w_ex12=5 * 2 * np.pi)
    experiment = Experiment(
        qubits=[q1, q2], g=0.005 * 2 * np.pi, t_exp=81,
        gate_type="X", drive_shape="Gaussian",
    )
    objective, _ = x_gate_objective(experiment)

    i1 = np.linspace(0, 200, GRID)
    q1_amp = np.linspace(0, 200, GRID)
    fidelity = np.array([
        [objective(I1p_amp=i, Q1p_amp=q, I2p_amp=0, Q2p_amp=0) for i in i1]
        for q in q1_amp
    ])

    x_mesh, y_mesh = np.meshgrid(i1, q1_amp)
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    ax.plot_surface(x_mesh, y_mesh, fidelity, cmap="viridis", edgecolor="none")
    ax.set_xlabel("I1 drive amplitude")
    ax.set_ylabel("Q1 drive amplitude")
    ax.set_zlabel("X-gate fidelity")
    ax.set_title("X-gate fidelity landscape")
    plt.show()


if __name__ == "__main__":
    main()
