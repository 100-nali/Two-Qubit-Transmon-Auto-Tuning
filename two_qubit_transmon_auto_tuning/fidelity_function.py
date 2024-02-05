#Imports
import numpy as np
from qutip import *
import matplotlib.pyplot as plt
import itertools
from dataclasses import dataclass
import math
from functools import partial

def operational_basis(dim):
    if dim == 2:
        return [[sigmax(), sigmay(), sigmaz(), qeye(2)]]*2
    if dim == 3:
        return [[qeye(3), qeye(3), qeye(3), qeye(3), qeye(3), qeye(3), qeye(3), qeye(3), qeye(3)]] * 2


# %%  Define global variables (constants)
# Define dim
dim = 3

# Define a1, a2
a1 = tensor(destroy(dim), qeye(dim))
a2 = tensor(qeye(dim), destroy(dim))

#Define n1 and n2
n1 = a1.dag() * a1
n2 = a2.dag() * a2

#Define op_basis
op_basis = operational_basis(2)

#%% define universal Pulse class

class Pulse:
    def __init__(self,  tmeas, **kwargs):
        self.tmeas = tmeas
        self.params = kwargs

# %%
class GaussianPulse(Pulse):
    def define_pulse(self):
        self.amp = self.params.get('amp', 50.0)
        self.center = self.params.get('center', 0.0)
        self.std = self.params.get('std', 1.0)

    def __call__(self, t, *args, **kwargs):
        self.define_pulse()
        return self.amp * np.exp(-(t - self.center) ** 2 / (2 * self.std ** 2)) / (
                self.std * np.sqrt(2 * np.pi))

#%%
class GRAPE(Pulse):
    def define_pulse(self):
        self.N = self.params.get('N', 5)
        self.amps = self.params.get('amps', [1,1,1,1,1])

    def __call__(self, t, *args, **kwargs):
        self.define_pulse()
        return

def drive_shape(shape, tmeas, **kwargs):
    if shape == 'Gaussian':
        return GaussianPulse(tmeas, **kwargs)
    if shape == 'GRAPE':
        return GRAPE(tmeas, **kwargs)

#%%
class LowPassStepValue:
    def __init__(self, value, step_time, tau, default_value):
        self.value =  value
        self.step_time = step_time
        self.tau = tau
        self.default_value = default_value


    def __call__(self, t, *args, **kwargs):
        return (self.value - self.default_value) * (1 - np.exp(
            -(t - self.step_time) / self.tau)) + self.default_value if t > self.step_time else self.default_value

## Redefine this class.
# class Linear:
#     def __init__(self, slope, start_time, stop_value, default_value=0, current_value=0):
#         self.slope = slope
#         self.start_time = start_time
#         self.stop_value = stop_value
#         self.default_value = default_value
#         self.current_value = current_value
#
#     def store_last(self, current_value):
#         self.current_value = current_value
#     def __call__(self, t, *args, **kwargs):
#         next_value = self.slope*(t-self.start_time) if t > self.start_time and self.current_value < self.stop_value else self.default_value
#         if self.current_value == self.stop_value:
#             next_value = self.stop_value
#         self.store_last(next_value)
#         return next_value


class Linear:
    def __init__(self, m, t0, t_stop):
        self.m = m
        self.t0 = t0
        self.t_stop = t_stop

    def __call__(self, t, *args, **kwargs):
        if t <= self.t_stop:
            y = self.m*(t-self.t0)
            return y

# %% Define a Drive class
class Drive:
    def __init__(self, I, Q):
        self.I = I  # in-phase drive
        self.Q = Q  # quadrature drive

#%%
class Qubit:
    def __init__(self, i, w_q, a_q, r, w_d):
        self.i = i              #index
        self.w_q = w_q          #qubit freq
        self.a_q = a_q          #anharmonicity
        self.r = r              #drive coupling
        self.w_d = w_d          #drive freq

    def set_wq(self, wq_update):
        self.w_q = wq_update
#%%
#%% Class that defines experimental parameters
class Experiment:
    def __init__(self, qubits, g, t_exp, gate_type, drive_shape):
        self.qubits = qubits
        self.g = g
        self.t_exp = t_exp
        self.gate_type = gate_type
        self.drive_shape = drive_shape

    def create_H(self, drives):
        qubits = self.qubits
        g = self.g
        q1, q2 = qubits
        d1, d2 = drives
        # Autonomous
        if type(q1.w_q) != float:
            H_01 = [n1, q1.w_q]
        else:
            H_01 = n1 * q1.w_q
        if type(q2.w_q) != float:
            H_02 = [n2, q2.w_q]
        else:
            H_02 = n2 * q2.w_q
        if type(q1.w_d) != float:
            H_0d1 = [-n1, q1.w_d]
        else:
            H_0d1 = -n1 * q1.w_d
        if type(q2.w_d) != float:
            H_0d2 = [-n2, q2.w_d]
        else:
            H_0d2 = -n2 * q2.w_d

        if dim == 3:
            H_0 = 0.5 * ((q1.a_q * a1.dag() * a1.dag() * a1 * a1) + (q2.a_q * a2.dag() * a2.dag() * a2 * a2))
        else:
            H_0 = np.zeros(n1.shape)

        # Drive terms
        if type(d1.I) != float:
            H_d1_0b = [-0.5 * q1.r * 1j * (a1 - a1.dag()), d1.I]
        else:
            H_d1_0b = -0.5 * q1.r * 1j * (a1 - a1.dag()) * d1.I

        if type(d1.Q) != float:
            H_d1_0a = [- 0.5 * q1.r * (a1 + a1.dag()), d1.Q]
        else:
            H_d1_0a = - 0.5 * q1.r * (a1 + a1.dag()) * d1.Q

        if type(d2.I) != float:
            H_d2_0b = [1j * (a2 - a2.dag()), d2.I]
        else:
            H_d2_0b = 1j * (a2 - a2.dag()) * d2.I

        if type(d2.Q) != float:
            H_d2_0a = [- 0.5 * q2.r * (a2 + a2.dag()), d2.Q]
        else:
            H_d2_0a = - 0.5 * q2.r * (a2 + a2.dag()) * d2.Q

        # delta = q1.w_q - q2.w_q
        if type(q1.w_q) != float:
            H_d1_1 = [g * a1 * a2.dag(), lambda t, *args: (np.exp(-1j * (q1.w_q(t)) - q2.w_q) * t)]
            H_d2_1 = [g * a1.dag() * a2, lambda t, *args: (np.exp(1j * (q1.w_q(t)) - q2.w_q) * t)]
        else:
            H_d1_1 = [g * a1 * a2.dag(), lambda t, *args: np.exp(-1j * (q1.w_q - q2.w_q) * t)]
            H_d2_1 = [g * a1.dag() * a2, lambda t, *args: np.exp(1j * (q1.w_q - q2.w_q) * t)]

        # Total H
        H = [H_0, H_01, H_02, H_0d1, H_0d2, H_d1_0a, H_d1_0b, H_d2_0a, H_d2_0b, H_d1_1, H_d2_1]
        return H

    def simulate_qpt(self, **kwargs):
        params = kwargs
        I1_p = params['I1_p']
        Q1_p = params['Q1_p']
        I2_p = params['I2_p']
        Q2_p = params['Q2_p']

        drive_1 = Drive(
            I=drive_shape(self.drive_shape, self.t_exp, **I1_p),
            Q=drive_shape(self.drive_shape, self.t_exp, **Q1_p)
        )

        drive_2 = Drive(
            I=drive_shape(self.drive_shape, self.t_exp, **I2_p),
            Q=drive_shape(self.drive_shape, self.t_exp, **Q2_p)
        )

        drives = [drive_1, drive_2]

        wq1 = 5.30 * 2 * np.pi  # exchange oscillations for these parameters occur at this freq.
        on_Cphase = False
        if self.gate_type == 'Cphase':
            on_Cphase = True

        qubits = CPhase(on_Cphase, self.qubits, wq1)
        qubit_1, qubit_2 = qubits

        nT: int = 100
        times = np.linspace(0, self.t_exp, nT)
        H = self.create_H(drives)  # drives will change to piecewise constant
        U_psi_real = qutip.propagator(H, times)
        U_psi_real_T = U_psi_real[nT - 1]

        rows = [0,1,3,4]
        cols = [0,1,3,4]

        U_psi_real_T = Qobj(U_psi_real_T[rows][:,cols])
        U_rho_real = spre(U_psi_real_T) * spost(U_psi_real_T.dag())
        U_rho_real = Qobj(U_rho_real)
        chi_real = qpt(Qobj(U_rho_real), op_basis)

        return chi_real

    def fidelity_CNOT(self, **kwargs):
        U_psi_CNOT = Qobj([[1, 0, 0, 0],
                           [0, 1, 0, 0],
                           [0, 0, 0, 1],
                           [0, 0, 1, 0]])
        U_rho_CNOT = spre(U_psi_CNOT) * spost(U_psi_CNOT.dag())
        chi_ideal_CNOT = qpt(U_rho_CNOT, op_basis)
        chi_real = self.simulate_qpt(**kwargs)

        # Evaluate process fidelity
        fidelity = process_fidelity(Qobj(chi_ideal_CNOT), Qobj(chi_real))
        return fidelity

    def fidelity_X(self,**kwargs):
        # sigma_x = basis(3, 0) * basis(3, 1).dag() + basis(3, 1) * basis(3, 0).dag()
        sigma_x = sigmax()
        U_psi_X = tensor(sigmax(), qeye(2))
        U_rho_X = spre(U_psi_X) * spost(U_psi_X.dag())
        chi_ideal_X = qpt(U_rho_X, op_basis)
        chi_real = self.simulate_qpt(**kwargs)
        # %% Evaluate process fidelity  ###########################
        fidelity = process_fidelity(Qobj(chi_ideal_X), Qobj(chi_real))
        return fidelity
    def fidelity_Y(self, **kwargs):
        sigma_y = sigmay()
        U_psi_Y = tensor(sigma_y, qeye(2))
        U_rho_Y = spre(U_psi_Y) * spost(U_psi_Y.dag())
        chi_ideal_Y = qpt(U_rho_Y, op_basis)
        chi_real = self.simulate_qpt(**kwargs)

        # %% Evaluate process fidelity  ###########################
        fidelity = process_fidelity(Qobj(chi_ideal_Y), Qobj(chi_real))
        return fidelity

    def fidelity_Y_90(self, **kwargs):
        U_psi_Y_90 = tensor(Qobj([[1 / np.sqrt(2), 1 / np.sqrt(2)], [-1 / np.sqrt(2), 1 / np.sqrt(2)]]), qeye(2))
        U_rho_Y_90 = spre(U_psi_Y_90) * spost(U_psi_Y_90.dag())
        chi_ideal_Y_90 = qpt(U_rho_Y_90, op_basis)
        chi_real = self.simulate_qpt(**kwargs)

        # %% Evaluate process fidelity
        fidelity = process_fidelity(Qobj(chi_ideal_Y_90), Qobj(chi_real))
        return fidelity

    def fidelity_H(self, **kwargs):
        U_psi_H = tensor(Qobj([[1 / np.sqrt(2), 1 / np.sqrt(2)], [1 / np.sqrt(2), -1 / np.sqrt(2)]]), qeye(2))
        U_rho_H = spre(U_psi_H) * spost(U_psi_H.dag())
        chi_ideal_H = qpt(U_rho_H, op_basis)
        chi_real = self.simulate_qpt( **kwargs)

        # %% Evaluate process fidelity
        fidelity = process_fidelity(Qobj(chi_ideal_H), Qobj(chi_real))
        return fidelity

    def fidelity_iSWAP(self, **kwargs):
        U_psi_SWAP = Qobj([[1, 0, 0, 0],
                           [0, 0, 1, 0],
                           [0, 1, 0, 0],
                           [0, 0, 0, 1]])

        U_psi_i = Qobj([[1, 0, 0, 0],
                        [0, 1j, 1, 0],
                        [0, 1, 1j, 0],
                        [0, 0, 0, 1]])

        U_psi_iSWAP = U_psi_SWAP * U_psi_i
        U_rho_iSWAP = spre(U_psi_iSWAP) * spost(U_psi_iSWAP.dag())
        chi_ideal_iSWAP = qpt(U_rho_iSWAP, op_basis)
        chi_real = self.simulate_qpt(**kwargs)

        # Evaluate process fidelity
        fidelity = process_fidelity(Qobj(chi_ideal_iSWAP), Qobj(chi_real))
        return fidelity

    def fidelity_CZ(self, **kwargs):
        U_psi_CZ = Qobj([[1, 0, 0, 0],
                         [0, 1, 0, 0],
                         [0, 0, 1, 0],
                         [0, 0, 0, -1]])
        U_rho_CZ = spre(U_psi_CZ) * spost(U_psi_CZ.dag())
        chi_ideal_CZ = qpt(U_rho_CZ, op_basis)
        chi_real = self.simulate_qpt( **kwargs)

        # Evaluate process fidelity
        fidelity = process_fidelity(Qobj(chi_ideal_CZ), Qobj(chi_real))
        return fidelity


#%% %%%%%%%%%%%%%%%%%%%%%%%%%
#%% TODO: change this (implement cphase fr).
def CPhase(on:bool, qubits, wq1):
    qubit_1,qubit_2 = qubits
    if on == True:
        #### Make wq1 set to a linear ramp!!
        # wq1 = LowPassStepValue(
        #     value=wq1,
        #     step_time=0,
        #     tau=8,
        #     default_value=qubit_1.w_q)
        qubit_1.set_wq(wq1)
    return qubits

#%%
def process_fidelity(chi_ideal: Qobj, chi_real: Qobj):
    fid = ((chi_ideal * chi_real).tr()).real
    return fid
