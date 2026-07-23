"""Sanity checks on the process-tomography / fidelity machinery.

The process fidelity of an ideal gate with itself must be exactly 1, which both
validates the metric's normalization and confirms the operator basis is correct.
"""

from qutip import Qobj, sigmax, sigmay, qeye, tensor, spre, spost, qpt

from two_qubit_transmon_auto_tuning import operational_basis, process_fidelity

OP_BASIS = operational_basis(2)

IDEAL_GATES = {
    "X on q1": tensor(sigmax(), qeye(2)),
    "Y on q1": tensor(sigmay(), qeye(2)),
    "CPHASE": Qobj([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]]),
}


def _ideal_chi(unitary):
    rho = spre(unitary) * spost(unitary.dag())
    return Qobj(qpt(rho, OP_BASIS))


def test_operational_basis_shape():
    basis = operational_basis(2)
    assert len(basis) == 2
    assert all(len(row) == 4 for row in basis)


def test_ideal_gate_self_fidelity_is_one():
    for name, unitary in IDEAL_GATES.items():
        chi = _ideal_chi(unitary)
        assert abs(process_fidelity(chi, chi) - 1.0) < 1e-6, name


def test_distinct_gates_are_distinguishable():
    chi_x = _ideal_chi(IDEAL_GATES["X on q1"])
    chi_cphase = _ideal_chi(IDEAL_GATES["CPHASE"])
    # Two different gates should not look identical to the fidelity metric.
    assert process_fidelity(chi_x, chi_cphase) < 0.5
