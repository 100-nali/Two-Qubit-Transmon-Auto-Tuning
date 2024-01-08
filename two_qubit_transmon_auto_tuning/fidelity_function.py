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

#TODO: define op_basis function
def operational_basis(dim):
    if dim == 2:
        return [[sigmax(), sigmay(), sigmaz(), qeye(2)]]*2
    if dim == 3:
        I = basis(3, 0) * basis(3, 0).dag() + basis(3, 1) * basis(3, 1).dag()

        sigma_x = basis(3, 0) * basis(3, 1).dag() + basis(3, 1) * basis(3, 0).dag()
        sigma_y = -1j * basis(3, 0) * basis(3, 1).dag() + 1j * basis(3, 1) * basis(3, 0).dag()
        sigma_z = basis(3, 0) * basis(3, 0).dag() - basis(3, 1) * basis(3, 1).dag()
        l4 = Qobj([[0, 0, 1], [0, 0, 0], [1, 0, 0]])
        l5 = Qobj([[0, 0, -1j], [0, 0, 0], [1j, 0, 0]])
        l6 = Qobj([[0, 0, 0], [0, 0, 1], [0, 1, 0]])
        l7 = Qobj([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]])
        l8 = Qobj([[1 / np.sqrt(3), 0, 0], [0, 1 / np.sqrt(3), 0], [0, 0, -2 / np.sqrt(3)]])

        testing = Qobj([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

        # op_basis = [[I, sigma_x, sigma_y, sigma_z]] * 2
        return [[sigma_x, sigma_y, sigma_z, l4, l5, l6, l7, l8, testing]] * 2

def fidelity_fn_internal(dim=3, **kwargs):
    I1_p = kwargs['I1_p']
    I2_p = kwargs['I2_p']
    Q1_p = kwargs['Q1_p']
    Q2_p = kwargs['Q2_p']

    #Define time-stamps
    nT:int = 100                                 #Number of time steps to propagate over
    tmeas = 64
    times = np.linspace(0,tmeas,nT)

    #Define inter-qubit coupling strength
    g: float = 0.005 * 2 * np.pi

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

    # %% Define a1, a2
    a1 = tensor(destroy(dim), qeye(dim))
    a2 = tensor(qeye(dim), destroy(dim))

    #Define n1 and n2
    n1 = a1.dag() * a1
    n2 = a2.dag() * a2

    #Define op_basis
    op_basis = operational_basis(dim-1)

    #%% Define Hamiltonian of 2-qubit system - Obtained from Barnaby's Notes
    def create_H(qubits, drives):
        q1, q2 = qubits
        d1, d2 = drives

        Delta_1 = q1.w_q - q1.w_d
        Delta_2 = q2.w_q - q2.w_d

        delta = q1.w_q - q2.w_q

        #Autonomous
        H_0 = Delta_1 * n1 + Delta_2 * n2
        if dim == 3:
            H_0 += 0.5 * ( (q1.a_q * a1.dag() * a1.dag() * a1 * a1) + (q2.a_q * a2.dag() * a2.dag() * a2 * a2) )

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



    # %% QPT over unknown quantum process  ###########################
    H = create_H([qubit_1, qubit_2], [drive_1, drive_2])
    U_psi_real = qutip.propagator(H, times)                   #List of matrices due to time dependence.
    U_psi_real_T = U_psi_real[nT-1]                           #Take entry from last time step
    U_rho_real = spre(U_psi_real_T) * spost(U_psi_real_T.dag())
    U_rho_real = Qobj(U_rho_real[0:16, 0:16])
    chi_real = qpt(Qobj(U_rho_real), op_basis)

    return chi_real


#%% Single Qubit: X Gate  ###########################
def fidelity_X(dim = 3, **kwargs):
    op_basis = operational_basis(dim)
    sigma_x = basis(3, 0) * basis(3, 1).dag() + basis(3, 1) * basis(3, 0).dag()
    U_psi_X = tensor(sigma_x, qeye(3))
    U_rho_X = spre(U_psi_X) * spost(U_psi_X.dag())
    chi_ideal_X = qpt(U_rho_X, op_basis)
    chi_real = fidelity_fn_internal(**kwargs)
    #%% Evaluate process fidelity  ###########################
    fidelity = process_fidelity_paper(Qobj(chi_ideal_X), Qobj(chi_real))
    return fidelity


#%% Single Qubit: Y Gate

def fidelity_Y(dim = 2, **kwargs):
    op_basis = operational_basis(dim)
    sigma_y = -1j * basis(3, 0) * basis(3, 1).dag() + 1j * basis(3, 1) * basis(3, 0).dag()
    U_psi_Y = tensor(sigma_y, qeye(3))
    U_rho_Y = spre(U_psi_Y) * spost(U_psi_Y.dag())
    chi_ideal_Y = qpt(U_rho_Y, op_basis)
    chi_real = fidelity_fn_internal(**kwargs)

    #%% Evaluate process fidelity  ###########################
    fidelity = process_fidelity_paper(Qobj(chi_ideal_Y), Qobj(chi_real))
    print('CNOT Process Fidelity at T = ', fidelity)
    return fidelity

#%% Single Qubit: Y pi/2

def fidelity_Y_90(dim = 2, **kwargs):
    op_basis = operational_basis(dim)
    if kwargs['dim'] == 2:
        U_psi_Y_90 = tensor(
            Qobj([[1 / np.sqrt(2), 1 / np.sqrt(2)], [-1 / np.sqrt(2), 1 / np.sqrt(2)]]), qeye(2))
    if kwargs['dim'] == 3:
        U_psi_Y_90 = tensor(Qobj([[1/np.sqrt(2),1/np.sqrt(2),0],[-1/np.sqrt(2), 1/np.sqrt(2),0], [0,0,0]]), qeye(3))

    U_rho_Y_90 = spre(U_psi_Y_90) * spost(U_psi_Y_90.dag())
    chi_ideal_Y_90 = qpt(U_rho_Y_90, op_basis)
    chi_real = fidelity_fn_internal(**kwargs)

    #%% Evaluate process fidelity
    fidelity = process_fidelity_paper(Qobj(chi_ideal_Y_90), Qobj(chi_real))
    return fidelity

#%% Single Qubit: Hadamard

def fidelity_H(dim = 2, **kwargs):
    op_basis = operational_basis(dim)

    if kwargs['dim'] == 2:
        U_psi_H =  tensor(Qobj([[1/np.sqrt(2),1/np.sqrt(2)],[1/np.sqrt(2), -1/np.sqrt(2)]]), qeye(2))
    if kwargs['dim'] == 3:
        U_psi_H = tensor(Qobj([[1/np.sqrt(2),1/np.sqrt(2),0],[1/np.sqrt(2), -1/np.sqrt(2),0],[0,0,0]]), qeye(3))
    U_rho_H = spre(U_psi_H) * spost(U_psi_H.dag())
    chi_ideal_H = qpt(U_rho_H, op_basis)
    chi_real = fidelity_fn_internal(**kwargs)

    #%% Evaluate process fidelity
    fidelity = process_fidelity_paper(Qobj(chi_ideal_H), Qobj(chi_real))
    return fidelity

#%% iSWAP Gate
def fidelity_iSWAP(dim = 2, **kwargs):
    op_basis = operational_basis(dim)

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
    chi_real = fidelity_fn_internal(**kwargs)

    # Evaluate process fidelity
    fidelity = process_fidelity_paper(Qobj(chi_ideal_iSWAP), Qobj(chi_real))
    return fidelity

#%% CZ Gate
def fidelity_CZ(dim = 2, **kwargs):
    op_basis = operational_basis(dim)
    #TODO: change to 3 dim!
    U_psi_CZ = Qobj([[1, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, 1, 0],
                        [0, 0, 0, -1]])
    U_rho_CZ = spre(U_psi_CZ) * spost(U_psi_CZ.dag())
    chi_ideal_CZ = qpt(U_rho_CZ, op_basis)
    chi_real = fidelity_fn_internal(**kwargs)

    # Evaluate process fidelity
    fidelity = process_fidelity_paper(Qobj(chi_ideal_CZ), Qobj(chi_real))
    return fidelity



#%% CNOT Gate
def fidelity_CNOT(dim = 2, **kwargs):
    op_basis = operational_basis(dim)

    #TODO: change to 3dim
    U_psi_CNOT = Qobj([[1, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, 0, 1],
                        [0, 0, 1, 0]])
    U_rho_CNOT = spre(U_psi_CNOT) * spost(U_psi_CNOT.dag())
    chi_ideal_CNOT = qpt(U_rho_CNOT, op_basis)
    chi_real = fidelity_fn_internal(**kwargs)

    #Evaluate process fidelity
    fidelity = process_fidelity_paper(Qobj(chi_ideal_CNOT), Qobj(chi_real))
    return fidelity
