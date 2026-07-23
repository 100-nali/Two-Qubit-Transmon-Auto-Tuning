"""Evaluate the CPHASE gate at its auto-tuned optimum and show the process matrix.

The flux-tuning parameters below are the optimum found by the Bayesian optimizer
(``ratio`` = fraction of the gate spent at the |11>-|20> resonance, ``w12`` = the
tuned qubit-1 frequency). Running this reproduces the headline ~0.9998 CPHASE
process fidelity and plots the reconstructed chi matrix.

Run:
    python examples/run_cphase.py
"""

import matplotlib.pyplot as plt
from qutip import qpt_plot_combined

from two_qubit_transmon_auto_tuning import (
    default_cphase_experiment,
    cphase_objective,
    flatten_dict,
)

# Auto-tuned optimum (see optimization.py / the README results).
BEST = {"ratio": 0.592, "w12": 33.68}


def main():
    experiment = default_cphase_experiment(t_exp=1000)

    # Fidelity at the tuned optimum.
    objective, _ = cphase_objective(experiment)
    fidelity = objective(CPHASE_ratio=BEST["ratio"], CPHASE_w12=BEST["w12"])
    print(f"CPHASE process fidelity at tuned optimum: {fidelity:.5f}")

    # Reconstruct and plot the process (chi) matrix.
    kwargs = flatten_dict({
        "I1p": {"amp": 0, "center": experiment.t_exp / 2, "std": experiment.t_exp / 6},
        "Q1p": {"amp": 0, "center": experiment.t_exp / 2, "std": experiment.t_exp / 6},
        "I2p": {"amp": 0, "center": experiment.t_exp / 2, "std": experiment.t_exp / 6},
        "Q2p": {"amp": 0, "center": experiment.t_exp / 2, "std": experiment.t_exp / 6},
        "CPHASE": {"ratio": BEST["ratio"], "w12": BEST["w12"]},
    })
    chi = experiment.simulate_qpt(**kwargs)["chi"]
    try:
        qpt_plot_combined(chi, [["i", "x", "y", "z"]] * 2, "CPHASE (auto-tuned)")
        plt.show()
    except AttributeError:
        # QuTiP 4.7's 3D tomography plotter uses Axes3D.w_xaxis, removed in
        # matplotlib >= 3.8. The fidelity above is unaffected; pin matplotlib
        # < 3.8 (see pyproject.toml) to render the chi-matrix figure.
        print("(Skipping chi-matrix plot: it requires matplotlib < 3.8.)")


if __name__ == "__main__":
    main()
