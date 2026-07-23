"""Bayesian auto-tuning of two-qubit gate parameters.

The gate-fidelity landscape produced by ``simulation.py`` is non-convex, riddled
with local optima, and expensive to evaluate (each point is a full Schrodinger
propagation followed by quantum process tomography). Gaussian-process Bayesian
optimization is well suited to exactly this regime: it finds high-fidelity
control parameters in a few tens of evaluations instead of an exhaustive sweep.

This module exposes three reusable pieces plus a runnable ``main``:

* :func:`run_bayesian_optimization` -- maximize any fidelity objective over a
  bounded parameter space.
* :func:`run_meta_optimization`     -- a *nested* Bayesian optimization that tunes
  the optimizer's own hyperparameters (acquisition weight, GP noise ``alpha``,
  number of init points and iterations).
* :func:`default_cphase_experiment` / :func:`cphase_objective` -- the CPHASE
  set-up used to produce the results in the README.

Run the full pipeline (meta-optimize, then optimize CPHASE) with::

    python -m two_qubit_transmon_auto_tuning.optimization
"""

from functools import partial

import numpy as np
import matplotlib.pyplot as plt
from bayes_opt import BayesianOptimization
from bayes_opt.util import UtilityFunction

from .simulation import Experiment, Qubit, flatten_dict


def default_cphase_experiment(t_exp: float = 1000.0) -> Experiment:
    """Return the two-transmon CPHASE experiment used for the README results.

    Frequencies are angular (rad/ns): ``6 * 2 * np.pi`` is a 6 GHz transmon.
    """
    q1 = Qubit(i=1, w_q=6 * 2 * np.pi, a_q=-0.3 * 2 * np.pi, r=0.01 * 2 * np.pi,
               w_d=6 * 2 * np.pi, w_ex12=5.3 * 2 * np.pi)
    q2 = Qubit(i=2, w_q=5 * 2 * np.pi, a_q=-0.3 * 2 * np.pi, r=0.01 * 2 * np.pi,
               w_d=5 * 2 * np.pi, w_ex12=5 * 2 * np.pi)
    return Experiment(
        qubits=[q1, q2],
        g=0.005 * 2 * np.pi,
        t_exp=t_exp,
        gate_type="Cphase",
        drive_shape="Gaussian",
    )


def cphase_objective(experiment: Experiment):
    """Build the (objective, bounds) pair for optimizing a CPHASE gate.

    The single-qubit drives are held off (amplitude 0); the CPHASE is activated
    purely by flux-tuning qubit 1 into the ``|11>`` <-> ``|20>`` resonance. The
    optimizer therefore searches the flux-pulse shape:

    * ``ratio`` -- fraction of the gate spent parked at the resonance, and
    * ``w12``   -- the tuned qubit-1 frequency at that resonance.
    """
    t = experiment.t_exp
    drives_off = {
        "I1p": {"amp": 0, "center": t / 2, "std": t / 6},
        "Q1p": {"amp": 0, "center": t / 2, "std": t / 6},
        "I2p": {"amp": 0, "center": t / 2, "std": t / 6},
        "Q2p": {"amp": 0, "center": t / 2, "std": t / 6},
    }
    objective = partial(experiment.fidelity_Cphase, **flatten_dict(drives_off))
    pbounds = flatten_dict({"CPHASE": {"ratio": (0, 1), "w12": (5 * 2 * np.pi, 6 * 2 * np.pi)}})
    return objective, pbounds


def x_gate_objective(experiment: Experiment):
    """Build the (objective, bounds) pair for optimizing a single-qubit X gate.

    The optimizer searches the four drive amplitudes (in-phase/quadrature on each
    qubit); pulse centers and widths are fixed to the gate window. A well-tuned
    X on qubit 1 is expected to drive a single quadrature and leave the others
    near zero -- a nice sanity check that the optimizer finds the physics.
    """
    t = experiment.t_exp
    fixed = flatten_dict({
        "I1p": {"center": t / 2, "std": t / 6},
        "Q1p": {"center": t / 2, "std": t / 6},
        "I2p": {"center": t / 2, "std": t / 6},
        "Q2p": {"center": t / 2, "std": t / 6},
        "CPHASE": {"ratio": 0.5, "w12": 5.5 * 2 * np.pi},  # unused for X; required by simulate_qpt
    })
    objective = partial(experiment.fidelity_X, **fixed)
    pbounds = {
        "I1p_amp": (0, 200), "Q1p_amp": (0, 200),
        "I2p_amp": (0, 50), "Q2p_amp": (0, 50),
    }
    return objective, pbounds


def run_bayesian_optimization(objective, pbounds, *, init_points, n_iter,
                              acquisition_weight, alpha, random_state=1):
    """Maximize ``objective`` over ``pbounds`` with a GP Bayesian optimizer.

    Uses the Expected-Improvement acquisition function; ``acquisition_weight``
    (kappa) trades exploration against exploitation and ``alpha`` is the assumed
    observation noise of the Gaussian process. Returns the fitted optimizer, so
    ``optimizer.max`` and ``optimizer.res`` are available to the caller.
    """
    optimizer = BayesianOptimization(
        f=objective,
        pbounds=pbounds,
        random_state=random_state,
        verbose=0,
        allow_duplicate_points=True,
    )
    optimizer.set_gp_params(alpha=alpha)
    utility_function = UtilityFunction(kind="ei", kappa=acquisition_weight)
    optimizer.maximize(init_points=init_points, n_iter=n_iter,
                       acquisition_function=utility_function)
    return optimizer


def run_meta_optimization(objective, pbounds, *, init_points=5, n_iter=10, random_state=1):
    """Meta-optimize the optimizer: tune its hyperparameters with an outer BO.

    Each evaluation of the outer loop runs a *complete* inner Bayesian
    optimization for a candidate hyperparameter set and reports the best fidelity
    it reached. The outer optimizer then searches for the hyperparameters that
    maximize that fidelity. Returns the outer optimizer's ``max`` dict
    (``{'target': ..., 'params': {...}}``).
    """
    def hyperparameter_objective(acquisition_weight, num_iterations, num_init_points, alpha):
        optimizer = run_bayesian_optimization(
            objective, pbounds,
            init_points=int(num_init_points),
            n_iter=int(num_iterations),
            acquisition_weight=acquisition_weight,
            alpha=alpha,
            random_state=random_state,
        )
        return optimizer.max["target"]

    hyperparameter_bounds = {
        "acquisition_weight": (1, 10),
        "num_iterations": (10, 30),
        "num_init_points": (5, 20),
        "alpha": (1e-5, 1e-1),
    }
    meta_optimizer = BayesianOptimization(
        f=hyperparameter_objective,
        pbounds=hyperparameter_bounds,
        random_state=random_state,
    )
    meta_optimizer.maximize(init_points=init_points, n_iter=n_iter)
    return meta_optimizer.max


def main():
    """Meta-optimize the optimizer, then auto-tune the CPHASE gate, then plot."""
    experiment = default_cphase_experiment()
    objective, pbounds = cphase_objective(experiment)

    print("Meta-optimizing the Bayesian optimizer's hyperparameters ...")
    best = run_meta_optimization(objective, pbounds)
    hp = best["params"]
    print(f"  best hyperparameters : {hp}")
    print(f"  fidelity at meta-optimum : {best['target']:.6f}")

    print("Auto-tuning the CPHASE gate with the tuned hyperparameters ...")
    optimizer = run_bayesian_optimization(
        objective, pbounds,
        init_points=int(hp["num_init_points"]),
        n_iter=int(hp["num_iterations"]),
        acquisition_weight=hp["acquisition_weight"],
        alpha=hp["alpha"],
    )
    print(f"  best CPHASE fidelity : {optimizer.max['target']:.6f}")
    print(f"  best parameters      : {optimizer.max['params']}")

    target_vals = [res["target"] for res in optimizer.res]
    plt.figure()
    plt.plot(target_vals)
    plt.xlabel("Evaluation")
    plt.ylabel("CPHASE process fidelity")
    plt.title("Bayesian optimization convergence")
    plt.show()


if __name__ == "__main__":
    main()
