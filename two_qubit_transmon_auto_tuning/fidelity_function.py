#Imports
import numpy as np
from qutip import *
import scipy.linalg as la
import matplotlib.pyplot as plt
import itertools
from dataclasses import dataclass
import math
from functools import partial


"""
TODO: figure out how to define is_11. Make some connection with the qpt code to the main function-code : is_11 is
dependent on the op in op_basis_list.
"""

def operational_basis(dim):
    if dim == 2:
        # return [[(sigmax()+ sigmay()), (sigmax()- sigmay()), (qeye(2)+sigmaz()), (qeye(2)-sigmaz())]]*2
        return [[sigmax(), sigmay(), sigmaz(), qeye(2)]]*2


# %%  Define global variables (constants)
# Define dim
dim = 3

# Define a1, a2
a1 = tensor(destroy(dim), qeye(dim))
a2 = tensor(qeye(dim), destroy(dim))

# Define n1 and n2
n1 = a1.dag() * a1
n2 = a2.dag() * a2

# Define op_basis
op_basis = operational_basis(2)
# %%
class GaussianPulse:
    def __init__(self, **kwargs):
        # super().__init__(**kwargs)
        self.params = kwargs
        self.amp = self.params.get('amp', 50.0)
        self.center = self.params.get('center', 0.0)
        self.std = self.params.get('std', 1.0)

    def __call__(self, t, **kwargs):
        # self.define_pulse()
        return self.amp * np.exp(-(t - self.center) ** 2 / (2 * self.std ** 2)) / (
                self.std * np.sqrt(2 * np.pi))

# %%

class DRAG:
    def __init__(self, **kwargs):
        # super().__init__(**kwargs)
        self.params = kwargs
        self.GaussianParams = self.params.get('GaussianParams')
        self.lamb = self.params.get('lamb', 0.75)
        self.alpha = self.params.get('alpha', 1)

    def __call__(self, t, **kwargs):
        # self.define_pulse()
        print(t)
        Gaussian = GaussianPulse(**self.GaussianParams)
        center = self.GaussianParams.get('center')
        std = self.GaussianParams.get('std')
        B_t_coeff = ((t - center)**3)/(2 * std**4)
        A_drag = Gaussian(t)
        B_drag = (self.lamb / self.alpha) * Gaussian(t) * B_t_coeff
        return [A_drag, B_drag]

#%%
class GRAPE:
    def __init__(self, **kwargs):
        self.params = kwargs
        self.N = self.params.get('N', 5)
        self.amps = self.params.get('amps', [1, 1, 1, 1, 1])

    def __call__(self, *args, **kwargs):
        # self.define_pulse()
        return 0

# %%
def drive_shape(shape, **kwargs):

    if shape == 'Gaussian':
        return GaussianPulse(**kwargs)
    if shape == 'DRAG':
        return DRAG(**kwargs)
    if shape == 'GRAPE':
        return GRAPE(**kwargs)

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

# %%
class Linear:
    def __init__(self, t0, y0, t_stop, y_stop):
        self.t0 = t0
        self.y0 = y0
        self.t_stop = t_stop
        self.y_stop = y_stop

    def __call__(self, t, *args, **kwargs):
        if t <= self.t_stop:
            m = (self.y_stop - self.y0)/(self.t_stop - self.t0)
            b = self.y0 - m*self.t0
            y = m*t + b
            return y

# %% Define the piecewise linear pulse class

class LinearFluxTuning:
    def __init__(self, t0, y0, t_tune, y_max, tau):
        self.t0 = t0
        self.y0 = y0
        self.t1 = t0 + t_tune
        self.t2 = t0 + t_tune + tau
        self.t_stop = t0 + 2*t_tune + tau
        self.y_max = y_max

    def __call__(self, t, *args, **kwargs):
        if t>self.t0 and t < self.t1:
            m = (self.y_max - self.y0)/(self.t1 - self.t0)
            b = self.y0 + m*self.t0

        elif t>= self.t1 and t < self.t2:
            m = 0
            b = self.y_max
        elif t>= self.t2 and t< self.t_stop:
            m = -(self.y_max - self.y0)/(self.t_stop - self.t2)
            b = self.y_max - m*self.t2
        else:
            m = 0
            b = self.y0
        y = m*t + b
        return y

#%%
class SquareFluxTuning:
    def __init__(self, t0, y0, t_stop, y_max):
        self.t0 = t0
        self.y0 = y0
        self.t_stop = t_stop
        self.y_max = y_max

    def __call__(self, t, *args, **kwargs):
        if t> self.t0 and t< self.t_stop:
            y = self.y_max
        else:
            y = self.y0
        return y

# %% Define a Drive class
class Drive:
    def __init__(self, I, Q):
        self.I = I  # in-phase drive
        self.Q = Q  # quadrature drive

#%%
class Qubit:
    def __init__(self, i, w_q, a_q, r, w_d, w_ex12):
        self.i = i              #index
        self.w_q = w_q          #qubit freq
        self.w_q_def = w_q      #default qubit frequency
        self.a_q = a_q          #anharmonicity
        self.r = r              #drive coupling
        self.w_d = w_d          #drive freq
        self.w_ex12 = w_ex12    # |11> and |20> exchange oscillation freq

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

    def create_H(self, drives, args):
        qubits = self.qubits
        g = self.g
        q1, q2 = qubits
        d1, d2 = drives
        #
        w1_set = self.Cphase(args)
        if self.gate_type == 'Cphase':
            q1.set_wq(w1_set)

        # Autonomous
        if type(q1.w_q) != float:
            H_01 = [n1, q1.w_q]
        else:
            H_01 = n1 * q1.w_q
        if type(q2.w_q) != float:
            H_02 = [n2, q2.w_q]
        else:
            H_02 = n2 * q2.w_q

        H_0d1 = -n1 * q1.w_d
        H_0d2 = -n2 * q2.w_d

        if dim == 3:
            H_0 = 0.5 * ((q1.a_q * a1.dag() * a1.dag() * a1 * a1) + (q2.a_q * a2.dag() * a2.dag() * a2 * a2))
        else:
            H_0 = np.zeros(n1.shape)

        # Drive terms
        if type(d1.I) != float:
            H_d1_0b = [-0.5 * q1.r * 1j * (a1 - a1.dag()), lambda t, args: d1.I(t)]
        else:
            H_d1_0b = -0.5 * q1.r * 1j * (a1 - a1.dag()) * d1.I

        if type(d1.Q) != float:
            H_d1_0a = [-0.5 * q1.r * (a1 + a1.dag()), lambda t, args: d1.Q(t)]
        else:
            H_d1_0a = -0.5 * q1.r * (a1 + a1.dag()) * d1.Q

        if type(d2.I) != float:
            H_d2_0b = [-0.5*1j * q2.r * (a2 - a2.dag()), lambda t, args: d2.I(t)]
        else:
            H_d2_0b = -0.5*1j * q2.r * (a2 - a2.dag()) * d2.I

        if type(d2.Q) != float:
            H_d2_0a = [-0.5 * q2.r * (a2 + a2.dag()), lambda t, args: d2.Q(t)]
        else:
            H_d2_0a = -0.5 * q2.r * (a2 + a2.dag()) * d2.Q

        delta_d = q1.w_d - q2.w_d
        H_d1_1 = [a1 * a2.dag(), lambda t, args: self.g * np.exp(-1j * delta_d * t)]
        H_d2_1 = [a1.dag() * a2, lambda t, args: self.g * np.exp(1j * delta_d * t)]

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
            I=drive_shape(self.drive_shape, **I1_p),
            Q=drive_shape(self.drive_shape,**Q1_p)
        )

        drive_2 = Drive(
            I=drive_shape(self.drive_shape,  **I2_p),
            Q=drive_shape(self.drive_shape, **Q2_p)
        )

        drives = [drive_1, drive_2]
#%%
        nT: int = 100
        times = np.linspace(0, self.t_exp, nT)

        # tau_ratios = np.linspace(0,1,100)
        tau_ratios = [0.7676767676767677]
        U33_angles = []
        for tau_ratio in tau_ratios:
            H = self.create_H(drives, args = {'tau_ratio': tau_ratio, 't0': 0}) #lambda function of time.
            U_psi_real = qutip.propagator(H, times)

    #%%
            np.set_printoptions(precision=1)
            U_psi_real_T = U_psi_real[nT - 1]
            rows = [0,1,3,4]
            cols = [0,1,3,4]
            U_psi_real_T = (Qobj((U_psi_real_T[rows][:,cols])))
            phis = np.angle([ U_psi_real_T[1,1],  U_psi_real_T[2,2],  U_psi_real_T[3,3]], deg = True)
            U33_angle = phis[2] - phis[1] - phis[0]
            U33_angles.append(U33_angle)
            print(U33_angle)
        print(U33_angles)
        # print(U_psi_real_T)
        U_rho_real = spre(U_psi_real_T) * spost(U_psi_real_T.dag())
        U_rho_real = Qobj(U_rho_real)

        chi_real = qpt(Qobj(U_rho_real), op_basis)

        # %% Plotting
        self.plot_qpt(chi_real, 'Unknown Process')
        return [chi_real,U_psi_real_T]  # the *16 term comes from the newly defined op_basis

    #%% QPT Plotting function
    def plot_qpt(self, chi, title):
        fig = qpt_plot_combined(chi, [["i", "x", "y", "z"]] * 2, title)
        plt.show()

    #%%
    def power_rabi(self, **kwargs):
        amps = np.linspace(0, 100, 100)
        times = np.linspace(0,81,100)
        outputs = np.zeros([len(amps), len(times), 8])

        params = kwargs
        I1_p = params['I1_p']
        Q1_p = params['Q1_p']
        I2_p = params['I2_p']
        Q2_p = params['Q2_p']

        drive_1 = Drive(
            I=drive_shape(self.drive_shape, **I1_p),
            Q=drive_shape(self.drive_shape, **Q1_p)
        )

        drive_2 = Drive(
            I=drive_shape(self.drive_shape, **I2_p),
            Q=drive_shape(self.drive_shape, **Q2_p)
        )

        drives = [drive_1, drive_2]

        I = basis(3, 0) * basis(3, 0).dag() + basis(3, 1) * basis(3, 1).dag()
        sigma_x = basis(3, 0) * basis(3, 1).dag() + basis(3, 1) * basis(3, 0).dag()
        sigma_y = -1j * basis(3, 0) * basis(3, 1).dag() + 1j * basis(3, 1) * basis(3, 0).dag()
        sigma_z = basis(3, 0) * basis(3, 0).dag() - basis(3, 1) * basis(3, 1).dag()

        sigma_x1 = tensor(sigma_x, I)
        sigma_x2 = tensor(I, sigma_x)

        sigma_y1 = tensor(sigma_y, I)
        sigma_y2 = tensor(I, sigma_y)

        sigma_z1 = tensor(sigma_z, I)
        sigma_z2 = tensor(I, sigma_z)

        for i, amp in enumerate(amps):
            # choose which drive to drive based on keyword
            if params['drive'] == "I":
                drive_1.I.amp = amp
            if params['drive'] == "Q":
                drive_1.Q.amp = amp
            # Define starting state
            psi0 = tensor(basis(3,0), basis(3,0))

            # create H
            H = self.create_H(drives, args = {'tau_ratio': 0.7676, 't0': 0})

            # input to mesolve
            res = mesolve(H, psi0, times, [], [sigma_x1, sigma_y1, sigma_z1, sigma_x2, sigma_y2, sigma_z2, n1, n2])
            outputs[i, ...] = np.stack(res.expect, axis=-1)

        fig, ax = plt.subplots(1, 1)
        ax.plot(amps, outputs[:, -1, 3 * (self.qubits[0].i - 1)], label='X', linestyle='--')
        ax.plot(amps, outputs[:, -1, 3 * (self.qubits[0].i - 1) + 1], label='Y', linestyle='--')
        ax.plot(amps, outputs[:, -1, 3 * (self.qubits[0].i - 1) + 2], label='Z')

        ax.set_ylim([-1, 1])
        ax.set_xlabel('Amplitude')
        ax.set_ylabel('Expectation value')
        ax.set_title('Power Rabi Oscillations')
        ax.legend()
        plt.show()
#%%
    def Cphase(self, args):
        exp = self
        on = exp.gate_type == 'Cphase'
        q1, q2 = exp.qubits
        t0 = args['t0']
        tau = args['tau_ratio'] * exp.t_exp
        t_tune = 0.5 * (exp.t_exp - tau)
        y0 = q1.w_q_def
        y_max = q1.w_ex12

        if on is True:
            w1_set = LinearFluxTuning(
                t0 = t0,
                t_tune = t_tune,
                y0 = y0,
                y_max = y_max,
                tau = tau
            )
            # w1_set = SquareFluxTuning(
            #     t0 = t0,
            #     y0 = y0,
            #     t_stop =exp.t_exp,
            #     y_max = y_max
            # )
            return w1_set
#%%
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

    def fidelity_Cphase(self, **kwargs):
        U_psi_Cphase = Qobj([[1, 0, 0, 0],
                           [0, 1, 0, 0],
                           [0, 0, 1, 0],
                           [0, 0, 0, -1]])
        U_rho_Cphase = spre(U_psi_Cphase) * spost(U_psi_Cphase.dag())
        chi_ideal_Cphase = qpt(U_rho_Cphase, op_basis)
        vars= self.simulate_qpt(**kwargs)
        chi_real, U = vars

        # Evaluate process fidelity
        fidelity = process_fidelity(Qobj(chi_ideal_Cphase), Qobj(chi_real))
        return [fidelity, U]

    def fidelity_X(self,**kwargs):
        # sigma_x = basis(3, 0) * basis(3, 1).dag() + basis(3, 1) * basis(3, 0).dag()
        sigma_x = sigmax()
        U_psi_X = tensor(sigmax(), qeye(2))
        U_rho_X = spre(U_psi_X) * spost(U_psi_X.dag())
        chi_ideal_X = qpt(U_rho_X, op_basis)
        vars = self.simulate_qpt(**kwargs)
        chi_real, U = vars
        # %% Evaluate process fidelity  ###########################
        fidelity = process_fidelity(Qobj(chi_ideal_X), Qobj(chi_real))
        return [fidelity,U]
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

#%%
def process_fidelity(chi_ideal: Qobj, chi_real: Qobj):
    fid = ((chi_ideal * chi_real).tr()).real
    return fid
