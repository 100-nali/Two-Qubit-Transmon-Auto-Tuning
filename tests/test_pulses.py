"""Fast, sim-free tests of the control-pulse building blocks and dict helpers."""

from two_qubit_transmon_auto_tuning import (
    GaussianPulse,
    DRAG,
    LinearFluxTuning,
    flatten_dict,
    unflatten_dict,
)


def test_flatten_unflatten_roundtrip():
    nested = {"I1p": {"amp": 1.0, "center": 2.0, "std": 3.0}}
    assert unflatten_dict(flatten_dict(nested)) == nested
    assert flatten_dict(nested) == {"I1p_amp": 1.0, "I1p_center": 2.0, "I1p_std": 3.0}


def test_gaussian_is_symmetric_about_its_center():
    pulse = GaussianPulse(amp=10.0, center=5.0, std=1.0)
    assert pulse(5.0) > pulse(4.0)          # peaks at the center
    assert abs(pulse(4.0) - pulse(6.0)) < 1e-12  # symmetric


def test_drag_returns_amplitude_and_derivative_components():
    drag = DRAG(GaussianParams={"amp": 50, "center": 500, "std": 100}, lamb=0.75, alpha=1)
    amplitude, derivative = drag(500)
    assert derivative == 0.0  # the DRAG correction vanishes at the pulse center
    assert amplitude > 0.0


def test_linear_flux_tuning_ramps_and_parks():
    # t0=0 -> ramp (0,100) -> park [100,200) -> ramp back [200,300) -> baseline
    flux = LinearFluxTuning(t0=0, y0=6.0, t_tune=100, y_max=5.0, tau=100)
    assert flux(-1) == 6.0                    # before the pulse: baseline
    assert abs(flux(150) - 5.0) < 1e-12       # parked at the tuned frequency
    assert 5.0 < flux(50) < 6.0               # mid-ramp: between the two levels
