#Imports
import numpy as np
from qutip import *
import matplotlib.pyplot as plt
import itertools
from dataclasses import dataclass
import math

#Define number of energy levels
dim = 2

#Define qubit 1 parameters
w_q1: float     = 1.0
a_q1: float     = 1.0           #anharmonicity for qubit 1
r1:float        = 1.0           #drive coupling for qubit 1

#Define qubit 2 parameters
w_q2: float     = 1.0
a_q2: float     = 1.0           #anharmonicity for qubit 2
r2:float        = 1.0           #drive coupling for qubit 2

#Define inter-qubit coupling strength
g: float        = 1.0

#Define pulse parameters - q1
w_d1            = 1.0           #drive frequency for qubit 1
amp_1           = 1.0           #pulse amplitude for qubit 1
mean_1          = 0             #mean of gaussian pulse for qubit 1
std_1           = 1.0           #std of gaussian pulse for qubit 1


#Define pulse parameters - q2
w_d2            = 1.0
amp_2           = 1.0
mean_2          = 0
std_2           = 1.0


#phi
phi1            = np.pi/4       #Phase angle for I1 and Q1
phi2            = np.pi/4       #Phase angle for I2 and Q2


#Define time-stamps
T               = g / (2*np.pi)                 #Final time step to propagate over
nT              = 100                           #Number of time steps to propagate over 
t               = np.linspace(0, T, nT)         #time vector

#Define driving functions
v_1 = lambda t, *args: amp_1 * np.exp(-(t - mean_1) ** 2 / (2 * std_1 ** 2)) / (
                std_1 * np.sqrt(2 * np.pi)) * np.exp(1j*(w_d1*t + phi1))
v_2 = lambda t, *args: amp_2 * np.exp(-(t - mean_2) ** 2 / (2 * std_2 ** 2)) / (
                std_2 * np.sqrt(2 * np.pi)) * np.exp(1j*(w_d2*t + phi2))

#Define a1, a2
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

#Define Hamiltonian of 2-qubit system
H_0 = 0.5 * w_q1 * sigma_z1 + 0.5*w_q2*sigma_z2 + 0.5* g*(sigma_y1*sigma_y2)
H_d1 = [r1*sigma_y1 ,v_1]
H_d2 = [r2*sigma_y2 ,v_2]

H = [H_0, H_d1, H_d2]

#%% iSWAP

#  Define the iSWAP gate 
U_psi_iSWAP = II * (0.5+0j) + XX * 0.5j + YY * 0.5j + ZZ * (0.5+0j)             #for state vector
U_rho_iSWAP = spre(U_psi_iSWAP) * spost(U_psi_iSWAP.dag())                      #for density matrix

#Find ideal process matrix
chi_ideal_iSWAP = qpt(U_rho_iSWAP, op_basis)
fig = qpt_plot_combined(chi_ideal_iSWAP, op_label, r'$i$SWAP')
plt.show()


#%% Cphase

#Define U: TODO: convert to tensor product form
U_psi_cphase = Qobj([[1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, -1]])

U_rho_cphase = spre(U_psi_cphase) * spost(U_psi_cphase.dag())
chi_ideal_cphase = qpt(U_rho_cphase, op_basis)

fig = qpt_plot_combined(chi_ideal_cphase, op_label, r'$CPHASE')
plt.show()

# %% QPT over unknown quantum process
U_psi_real = qutip.propagator(H,t)                      #List of lists due to time dependence.
U_psi_real = U_psi_real[nT-1]                           #Take entry from last time step
U_rho_real = spre(U_psi_real) * spost(U_psi_real.dag())
chi_real = qpt(U_rho_real, op_basis)

fig = qpt_plot_combined(chi_real, op_label, r'$Actual Process$')
plt.show()

#%% Evaluate process fidelity
U1 = U_rho_real
U2 = U_rho_iSWAP

gate_fidelity = process_fidelity(U_rho_real, U_rho_iSWAP, normalize=True)
print('Process Fidelity at T = ', gate_fidelity)