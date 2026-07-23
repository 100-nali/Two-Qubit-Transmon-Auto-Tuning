"""Two-qubit transmon gate simulation and Bayesian auto-tuning.

Public API
----------
Simulation (``two_qubit_transmon_auto_tuning.simulation``):
    ``Experiment``, ``Qubit``, ``Drive``, the pulse/flux classes
    (``GaussianPulse``, ``DRAG``, ``LinearFluxTuning``, ...) and the helpers
    ``flatten_dict`` / ``unflatten_dict`` / ``process_fidelity``.

Optimization (``two_qubit_transmon_auto_tuning.optimization``):
    ``run_bayesian_optimization``, ``run_nested_optimization``,
    ``default_cphase_experiment`` and ``cphase_objective``.
"""

from .simulation import (
    Experiment,
    Qubit,
    Drive,
    GaussianPulse,
    DRAG,
    GRAPE,
    LinearFluxTuning,
    SquareFluxTuning,
    drive_shape,
    flatten_dict,
    unflatten_dict,
    virtual_Z_cphase,
    process_fidelity,
    operational_basis,
)
from .optimization import (
    run_bayesian_optimization,
    run_nested_optimization,
    default_cphase_experiment,
    cphase_objective,
    x_gate_objective,
)

__version__ = "0.1.0"

__all__ = [
    "Experiment",
    "Qubit",
    "Drive",
    "GaussianPulse",
    "DRAG",
    "GRAPE",
    "LinearFluxTuning",
    "SquareFluxTuning",
    "drive_shape",
    "flatten_dict",
    "unflatten_dict",
    "virtual_Z_cphase",
    "process_fidelity",
    "operational_basis",
    "run_bayesian_optimization",
    "run_nested_optimization",
    "default_cphase_experiment",
    "cphase_objective",
    "x_gate_objective",
]
