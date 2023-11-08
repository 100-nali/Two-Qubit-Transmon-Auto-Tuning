#Imports
import numpy as np
from qutip import *
import matplotlib.pyplot as plt
import itertools
from dataclasses import dataclass
import math

#Define number of energy levels
dim:int         = 2

#Define time-stampss
T:float         = 71                                    #Indicator of pulse time
t_meas:float    = 10                                    #Time to take measurement - fixed.
total_T:float   = T + t_meas
nT:int          = 100                                   #Number of time steps to propagate over 
times           = np.linspace(0, total_T, nT)           #Time vector


#Define inter-qubit coupling strength
g: float        = 0.005 * 2 * np.pi

# %% Define a Qubit class
class Qubit:
    def __init__(self, i, w_q, a_q, r, w_d):
        self.i = i              #index
        self.w_q = w_q          #qubit freq
        self.a_q = a_q          #anharmonicity
        self.r = r              #drive coupling
        self.w_d = w_d          #drive freq

# %% Define a Pulse class
class Pulse:
    def __init__(self, amp, center, std):
        self.amp = amp
        self.center = center
        self.std = std

    def __call__(self, t, *args, **kwargs):
        return self.amp * np.exp(-(t - self.center) ** 2 / (2 * self.std ** 2)) / (
                self.std * np.sqrt(2 * np.pi))
    
# %% Define a Drive class
class Drive:
    def __init__(self, I, Q):
        self.I = I          #in-phase drive
        self.Q = Q          #quadrature drive

#%% Define the qubits
qubit_1 = Qubit(
    i = 1,
    w_q = 5 * 2 * np.pi ,
    a_q = 0.3 * 2 * np.pi,
    r = 0.01 * 2 * np.pi,
    w_d = 5 * 2 * np.pi
)

qubit_2 = Qubit(
    i = 2,
    w_q = 5 * 2 * np.pi,
    a_q = 0.3 * 2 * np.pi,
    r = 0.01 * 2 * np.pi,
    w_d = 5 * 2 * np.pi
)

#%% Define the pulses at qubits 1 and 2
drive_1 =  Drive(
    I = Pulse(
    amp = np.pi *100/ 2,
    center= T/2, 
    std = T/6
),
Q = 0
)

drive_2 =  Drive(
    I = 0,
    Q = 0
)

#%% Define Hamiltonian of 2-qubit system - Obtained from Barnaby's Notes

def create_H(qubits, drives):
    q1, q2 = qubits
    d1, d2 = drives

    delta = q1.w_d - q2.w_d

    #Autonomous
    H_0 = n1 + n2 +\
            0.5*(q1.a_q*(a1.dag() * a1.dag() * a1 * a1) +\
            q2.a_q*(a2.dag() * a2.dag() * a2 * a2))

    #Drive terms
    if type(d1.I) == Pulse:
        H_d1_0b = [-0.5*q1.r * 1j*(a1-a1.dag()), d1.I]
    else:
        H_d1_0b = -0.5*q1.r * 1j*(a1-a1.dag())* d1.I

    if type(d1.Q) == Pulse:
        H_d1_0a = [- 0.5* q1.r*(a1 + a1.dag()),d1.Q]
    else:
        H_d1_0a = - 0.5* q1.r*(a1 + a1.dag())*d1.Q

    if type(d2.I) == Pulse:
        H_d2_0b = [1j*(a2-a2.dag()), d2.I]
    else:
        H_d2_0b = 1j*(a2-a2.dag())* d2.I

    if type(d2.Q) == Pulse:
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

# %% TODO: Implement Power Rabi ##################################
###Comment this out if you want the rest of the code to work/run
amps = np.linspace(0, 100, 100)
outputs = np.zeros([len(amps), len(times), 8])
for i,amp in enumerate(amps):
    #change I1 amp to amp
    drive_1.I.amp = amp

    #Define starting states
    psi0 = tensor(basis(2, 0), basis(2, 0))

    #create H
    H = create_H([qubit_1, qubit_2], [drive_1, drive_2])

    #input to mesolve
    res = mesolve(H, psi0, times, [], [sigma_x1, sigma_y1, sigma_z1, sigma_x2, sigma_y2, sigma_z2, n1, n2])
    outputs[i, ...] = np.stack(res.expect, axis=-1)
print(shape(outputs))

fig, ax = plt.subplots(1, 1)
ax.plot(amps, outputs[:, -1, 3 * (qubit_1.i - 1)], label='X', linestyle='--')
ax.plot(amps, outputs[:, -1, 3 * (qubit_1.i - 1) + 1], label='Y', linestyle='--')
ax.plot(amps, outputs[:, -1, 3 * (qubit_1.i - 1) + 2], label='Z')

ax.set_ylim([-1, 1])
ax.set_xlabel('Amplitude')
ax.set_ylabel('Expectation value')
ax.set_title('Power Rabi Oscillations')
ax.legend()
plt.show()


# %% QPT over unknown quantum process  ###########################
H = create_H([qubit_1, qubit_2], [drive_1, drive_2])
U_psi_real = qutip.propagator(H, times)                   #List of matrices due to time dependence.
U_psi_real_T = U_psi_real[nT-1]                           #Take entry from last time step

# TODO: Plot what's going on over the bloch sphere.


U_rho_real = spre(U_psi_real_T) * spost(U_psi_real_T.dag())
chi_real = qpt(U_rho_real, op_basis)

fig = qpt_plot_combined(chi_real, op_label, r'$Actual Process$')
plt.show()

#%% iSWAP ###########################

#  Define the iSWAP gate
U_psi_iSWAP = II * (0.5+0j) + XX * 0.5j + YY * 0.5j + ZZ * (0.5+0j)             #for state vector
U_rho_iSWAP = spre(U_psi_iSWAP) * spost(U_psi_iSWAP.dag())                      #for density matrix

#Find ideal process matrix
chi_ideal_iSWAP = qpt(U_rho_iSWAP, op_basis)
# fig = qpt_plot_combined(chi_ideal_iSWAP, op_label, r'$i$SWAP')
# plt.show()


#%% Cphase #### #######################

#Define U: TODO: convert to tensor product form
U_psi_cphase = Qobj([[1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, -1]])


U_rho_cphase = spre(U_psi_cphase) * spost(U_psi_cphase.dag())
chi_ideal_cphase = qpt(U_rho_cphase, op_basis)

# fig = qpt_plot_combined(chi_ideal_cphase, op_label, r'$CPHASE')
# plt.show()


#%% X Gate  ###########################
U_psi_CNOT = Qobj([[1, 0, 0, 0],
       [0, 1, 0, 0],
        [0, 0, 0, 1],
        [0, 0, 1, 0]])

U_rho_CNOT = spre(U_psi_CNOT) * spost(U_psi_CNOT.dag())
chi_ideal_CNOT = qpt(U_rho_CNOT, op_basis)
fig = qpt_plot_combined(chi_ideal_CNOT, op_label, r'$CPHASE')
plt.show()


#%% Evaluate process fidelity  ###########################

def process_fidelity_s(chi_ideal, chi_real):
    # Calculate the square root of the ideal process matrix
    sqrt_chi_ideal = chi_ideal.sqrtm()

    # Calculate the matrix square root of sqrt_chi_ideal * chi_real * sqrt_chi_ideal
    sqrt_term = sqrt_chi_ideal * chi_real * sqrt_chi_ideal

    # Calculate the trace of the square root term
    trace_sqrt_term = sqrt_term.tr()

    # Calculate the process fidelity using the formula
    fidelity = (np.abs(trace_sqrt_term) ** 2).real

    return fidelity

chi_ideal = Qobj(chi_ideal_iSWAP)
chi_real = Qobj(chi_real)
fidelity = process_fidelity_s(chi_ideal, chi_real)

print('iSWAP Process Fidelity at T = ', fidelity)
