"""Regression test for the headline result: the auto-tuned CPHASE gate.

Reproduces the ~0.9998 process fidelity obtained by flux-tuning qubit 1 into the
|11>-|20> resonance and applying the virtual-Z correction.
"""

from two_qubit_transmon_auto_tuning import cphase_objective

# Optimum found by the Bayesian optimizer.
BEST_RATIO = 0.592
BEST_W12 = 33.68


def test_cphase_reaches_high_fidelity_at_optimum(cphase_experiment):
    objective, _ = cphase_objective(cphase_experiment)
    fidelity = objective(CPHASE_ratio=BEST_RATIO, CPHASE_w12=BEST_W12)
    assert fidelity > 0.999


def test_virtual_z_correction_helps(cphase_experiment):
    """Without flux tuning (ratio 0), the CPHASE fidelity should be far from 1."""
    objective, _ = cphase_objective(cphase_experiment)
    off_resonance = objective(CPHASE_ratio=0.0, CPHASE_w12=BEST_W12)
    on_resonance = objective(CPHASE_ratio=BEST_RATIO, CPHASE_w12=BEST_W12)
    assert on_resonance > off_resonance
