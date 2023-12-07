#Imports
import numpy as np
from qutip import *
import matplotlib.pyplot as plt
import itertools
from dataclasses import dataclass
import math
from functools import partial

def process_fidelity_paper(chi_ideal: Qobj, chi_real: Qobj):
        fid = ((chi_ideal * chi_real).tr()).real
        return fid

def process_fidelity_t(chi_ideal: Qobj, chi_real: Qobj):
 # Calculate the square root of the ideal process matrix
    sqrt_chi_ideal = chi_ideal.sqrtm()

    # Calculate the matrix square root of sqrt_chi_ideal * chi_real * sqrt_chi_ideal
    sqrt_term = sqrt_chi_ideal * chi_real * sqrt_chi_ideal

    # Calculate the trace of the square root term
    trace_sqrt_term = sqrt_term.tr()

    # Calculate the process fidelity using the formula
    fidelity = (np.abs(trace_sqrt_term) ** 2).real
    return fidelity


def process_fidelity_f(chi_ideal: Qobj, chi_real: Qobj):
 # Calculate the square root of the ideal process matrix
    sqrt_chi_ideal = chi_ideal.sqrtm()

    # Calculate the matrix square root of sqrt_chi_ideal * chi_real * sqrt_chi_ideal
    sqrt_term = sqrt_chi_ideal * chi_real

    # Calculate the trace of the square root term
    another_sqrt = sqrt_term.sqrtm()

    trace_sqrt_term = another_sqrt.tr()

    # Calculate the process fidelity using the formula
    fidelity = (np.abs(trace_sqrt_term) ** 2)

    return fidelity


def fidelity_fn_internal(**kwargs):
    I1_p = kwargs['I1_p']
    I2_p = kwargs['I2_p']
    Q1_p = kwargs['Q1_p']
    Q2_p = kwargs['Q2_p']

    #Define number of energy levels
    dim:int         = 2

    #Define time-stamps
    nT:int          = 100                                 #Number of time steps to propagate over
    tmeas = 64
    times = np.linspace(0,tmeas,nT)

    #Define inter-qubit coupling strength
    g: float        = 0.005 * 2 * np.pi

    # %% Define a Gaussian Pulse class
    class Pulse:
        def __init__(self, amp, tmeas):
            self.amp = amp
            self.center = tmeas/2
            self.std = tmeas/6

        def __call__(self, t, *args, **kwargs):
            return self.amp * np.exp(-(t - self.center) ** 2 / (2 * self.std ** 2)) / (
                    self.std * np.sqrt(2 * np.pi))
        
    # %% Define a Drive class
    class Drive:
        def __init__(self, I, Q):
            self.I = I          #in-phase drive
            self.Q = Q          #quadrature drive

    drive_1 = Drive(
        I = Pulse(I1_p, tmeas),
        Q = Pulse(Q1_p, tmeas)
    )
    
    drive_2 = Drive(
        I = Pulse(I2_p,tmeas),
        Q = Pulse(Q2_p,tmeas)
    )

    # %% Define a Qubit class
    class Qubit:
        def __init__(self, i, w_q, a_q, r, w_d):
            self.i = i              #index
            self.w_q = w_q          #qubit freq
            self.a_q = a_q          #anharmonicity
            self.r = r              #drive coupling
            self.w_d = w_d          #drive freq

    #%% Define the qubits
    qubit_1 = Qubit(
        i = 1,
        w_q = 5 * 2 * np.pi ,
        a_q = -0.3 * 2 * np.pi,
        r = 0.01 * 2 * np.pi,
        w_d = 5 * 2 * np.pi
    )

    qubit_2 = Qubit(
        i = 2,
        w_q = 6 * 2 * np.pi,
        a_q = -0.3 * 2 * np.pi,
        r = 0.01 * 2 * np.pi,
        w_d = 6 * 2 * np.pi
    )

    #%% Define Hamiltonian of 2-qubit system - Obtained from Barnaby's Notes

    def create_H(qubits, drives):
        q1, q2 = qubits
        d1, d2 = drives

        Delta_1 = q1.w_q - q1.w_d
        Delta_2 = q2.w_q - q2.w_d

        delta = q1.w_q - q2.w_q

        #Autonomous
        H_0 = Delta_1 * n1 + Delta_2 * n2

        #Drive terms
        if type(d1.I) != float:
            H_d1_0b = [-0.5*q1.r * 1j*(a1-a1.dag()), d1.I]
        else:
            H_d1_0b = -0.5*q1.r * 1j*(a1-a1.dag())* d1.I

        if type(d1.Q) !=float:
            H_d1_0a = [- 0.5* q1.r*(a1 + a1.dag()),d1.Q]
        else:
            H_d1_0a = - 0.5* q1.r*(a1 + a1.dag())*d1.Q

        if type(d2.I) !=float:
            H_d2_0b = [1j*(a2-a2.dag()), d2.I]
        else:
            H_d2_0b = 1j*(a2-a2.dag())* d2.I

        if type(d2.Q) !=float:
            H_d2_0a = [- 0.5*q2.r*(a2 + a2.dag()),d2.Q]
        else:
            H_d2_0a = - 0.5*q2.r*(a2 + a2.dag())*d2.Q
    
        H_d1_1 = [g*a1*a2.dag(), lambda t, *args: np.exp(-1j*delta*t)]
        H_d2_1 = [g*a1.dag()*a2, lambda t, *args: np.exp(1j*delta*t)]

        #Total H
        H = [H_0, H_d1_0a, H_d1_0b, H_d2_0a, H_d2_0b, H_d1_1, H_d2_1]
        return H

    # %% Define a1, a2
    a1 = tensor(destroy(dim), qeye(dim))
    a2 = tensor(qeye(dim), destroy(dim))

    #Define n1 and n2
    n1 = a1.dag() * a1
    n2 = a2.dag() * a2

    #Define op_basis
    op_basis = [[qeye(2), sigmax(), sigmay(), sigmaz()]] * 2

    #Define plot labels
    op_label = [["i", "x", "y", "z"]] * 2

    #Define sigma matrices for the qubits
    sigma_x1 = tensor(sigmax(), qeye(2))
    sigma_y1 = tensor(sigmay(), qeye(2))
    sigma_z1 = tensor(sigmaz(), qeye(2))

    sigma_x2 = tensor(qeye(2), sigmax())
    sigma_y2 = tensor(qeye(2), sigmay())
    sigma_z2 = tensor(qeye(2), sigmaz())

    II       = tensor(qeye(2), qeye(2))
    XX       = tensor(sigmax(), sigmax())
    YY       = tensor(sigmay(), sigmay())
    ZZ       = tensor(sigmaz(), sigmaz())

    # %% QPT over unknown quantum process  ###########################
    H = create_H([qubit_1, qubit_2], [drive_1, drive_2])
    U_psi_real = qutip.propagator(H, times)                   #List of matrices due to time dependence.
    U_psi_real_T = U_psi_real[nT-1]                           #Take entry from last time step
    U_rho_real = spre(U_psi_real_T) * spost(U_psi_real_T.dag())
    chi_real = qpt(U_rho_real, op_basis)

    return chi_real


#%% Single Qubit: X Gate  ###########################
def fidelity_X(**kwargs):
    op_basis = [[qeye(2), sigmax(), sigmay(), sigmaz()]] * 2

    U_psi_X = tensor(sigmax(),qeye(2))
    U_rho_X = spre(U_psi_X) * spost(U_psi_X.dag())
    chi_ideal_X = qpt(U_rho_X, op_basis)
    chi_real = fidelity_fn_internal(**kwargs)
    #%% Evaluate process fidelity  ###########################
    fidelity = process_fidelity_paper(Qobj(chi_ideal_X), Qobj(chi_real))
    return fidelity


#%% Single Qubit: Y Gate

def fidelity_Y(**kwargs):
    op_basis = [[qeye(2), sigmax(), sigmay(), sigmaz()]] * 2

    U_psi_Y = tensor(sigmay(), qeye(2))
    U_rho_Y = spre(U_psi_Y) * spost(U_psi_Y.dag())
    chi_ideal_Y = qpt(U_rho_Y, op_basis)
    chi_real = fidelity_fn_internal(**kwargs)

    #%% Evaluate process fidelity  ###########################
    fidelity = process_fidelity_paper(Qobj(chi_ideal_Y), Qobj(chi_real))
    print('CNOT Process Fidelity at T = ', fidelity)
    return fidelity

#%% Single Qubit: Y pi/2

def fidelity_Y_90(**kwargs):
    op_basis = [[qeye(2), sigmax(), sigmay(), sigmaz()]] * 2

    U_psi_Y_90 = tensor(Qobj([[1/np.sqrt(2),1/np.sqrt(2)],[-1/np.sqrt(2), 1/np.sqrt(2)]]), qeye(2))
    U_rho_Y_90 = spre(U_psi_Y_90) * spost(U_psi_Y_90.dag())
    chi_ideal_Y_90 = qpt(U_rho_Y_90, op_basis)
    chi_real = fidelity_fn_internal(**kwargs)

    #%% Evaluate process fidelity
    fidelity = process_fidelity_paper(Qobj(chi_ideal_Y_90), Qobj(chi_real))
    return fidelity

#%% Single Qubit: Hadamard

def fidelity_H(**kwargs):
    op_basis = [[qeye(2), sigmax(), sigmay(), sigmaz()]] * 2

    U_psi_H = tensor(Qobj([[1/np.sqrt(2),1/np.sqrt(2)],[1/np.sqrt(2), -1/np.sqrt(2)]]), qeye(2))
    U_rho_H = spre(U_psi_H) * spost(U_psi_H.dag())
    chi_ideal_H = qpt(U_rho_H, op_basis)
    chi_real = fidelity_fn_internal(**kwargs)

    #%% Evaluate process fidelity
    fidelity = process_fidelity_paper(Qobj(chi_ideal_H), Qobj(chi_real))
    return fidelity

#%% iSWAP Gate
def fidelity_iSWAP(**kwargs):
    op_basis = [[qeye(2), sigmax(), sigmay(), sigmaz()]] * 2

    U_psi_SWAP = Qobj([[1, 0, 0, 0],
                       [0, 0, 1, 0],
                        [0, 1, 0, 0],
                        [0, 0, 0, 1]])

    U_psi_i = Qobj([[1, 0, 0, 0],
                        [0, 1j, 1, 0],
                        [0, 1, 1j, 0],
                        [0, 0, 0, 1]])

    U_psi_CNOT = U_psi_SWAP * U_psi_i
    U_rho_CNOT = spre(U_psi_CNOT) * spost(U_psi_CNOT.dag())
    chi_ideal_CNOT = qpt(U_rho_CNOT, op_basis)
    chi_real = fidelity_fn_internal(**kwargs)

    # Evaluate process fidelity
    fidelity = process_fidelity_paper(Qobj(chi_ideal_CNOT), Qobj(chi_real))
    return fidelity

#%% CZ Gate
def fidelity_CZ(**kwargs):
    op_basis = [[qeye(2), sigmax(), sigmay(), sigmaz()]] * 2

    U_psi_CNOT = Qobj([[1, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, 1, 0],
                        [0, 0, 0, -1]])
    U_rho_CNOT = spre(U_psi_CNOT) * spost(U_psi_CNOT.dag())
    chi_ideal_CNOT = qpt(U_rho_CNOT, op_basis)
    chi_real = fidelity_fn_internal(**kwargs)

    # Evaluate process fidelity
    fidelity = process_fidelity_paper(Qobj(chi_ideal_CNOT), Qobj(chi_real))
    return fidelity



#%% CNOT Gate
def fidelity_CNOT(**kwargs):
    op_basis = [[qeye(2), sigmax(), sigmay(), sigmaz()]] * 2

    U_psi_CNOT = Qobj([[1, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, 0, 1],
                        [0, 0, 1, 0]])
    U_rho_CNOT = spre(U_psi_CNOT) * spost(U_psi_CNOT.dag())
    chi_ideal_CNOT = qpt(U_rho_CNOT, op_basis)
    chi_real = fidelity_fn_internal(**kwargs)

    #Evaluate process fidelity
    fidelity = process_fidelity_paper(Qobj(chi_ideal_CNOT), Qobj(chi_real))
    print('CNOT Process Fidelity at T = ', fidelity)
    return fidelity
