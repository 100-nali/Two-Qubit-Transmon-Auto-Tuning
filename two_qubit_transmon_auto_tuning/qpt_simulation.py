#Imports
import numpy as np
from qutip import *
import matplotlib.pyplot as plt
import itertools
from dataclasses import dataclass
import math


#Define number of energy levels
dim:int         = 2

#Define qubit 1 parameters
w_q1: float     = 1.0 
a_q1: float     = 2.0           #anharmonicity for qubit 1
r1:float        = 1.0           #drive coupling for qubit 1

#Define qubit 2 parameters
w_q2: float     = 1.1
a_q2: float     = 2.0           #anharmonicity for qubit 2
r2:float        = 1.0           #drive coupling for qubit 2

#Define inter-qubit coupling strength
g: float        = 0.001

#Define pulse parameters - q1
w_d1:float      = 1.0                   #drive frequency for qubit 1
amp_1:float     = 0                   #pulse amplitude for qubit 1
mean_1:float    = 0                     #mean of gaussian pulse for qubit 1
std_1:float     = 1.0 *  np.pi / (2* g)                   #std of gaussian pulse for qubit 1


#Define pulse parameters - q2
w_d2:float      = 1.1
amp_2:float     = 0
mean_2:float    = 0
std_2:float     = 1.0


#phi
phi1            = np.pi/2              #Phase angle for I1 and Q1
phi2            = np.pi/2               #Phase angle for I2 and Q2


#Define args
args = {
    'w_d1' : w_d1,
    'std_1' : std_1,
    'mean_1': mean_1,
    'phi1' : phi1,
    'w_d2' : w_d2,
    'std_2' : std_2,
    'mean_2': mean_2,
    'phi2': phi2
}

#Define alternative H definition's drive term
I1 = np.sin(np.pi/4)
Q1 = np.cos(np.pi/4)
I2 = np.sin(np.pi/4)
Q2 = np.cos(np.pi/4)

#Define time-stamps
T               = 505                        #Final time step to propagate over
nT              = 100                           #Number of time steps to propagate over (doesnt affect t_real)
t               = np.linspace(0, T, nT)         #time vector



#Define driving pulses
def wave_shape1(t, args):
    x = np.exp(1j*(args['w_d1']*t + args['phi1']))
    return x

def wave_shape2(t, args):
    x = np.exp(1j*(args['w_d2']*t + args['phi2']))
    return x

def pulse_shape1(t, args):
    x = 1.0/ (np.sqrt(2.0 * np.pi) * args['std_1']) * np.exp(-np.power((t - args['mean_1']) / args['std_2'], 2.0) / 2)
    return x

def pulse_shape2(t, args):
    x = 1.0/ (np.sqrt(2.0 * np.pi) * args['std_2']) * np.exp(-np.power((t - args['mean_2']) / args['std_2'], 2.0) / 2)
    return x

def v_t1(t, args):
    x = wave_shape1(t, args)
    y = pulse_shape1(t,args)
    v = x*y
    return v

def v_t2(t, args):
    x = wave_shape2(t, args)
    y = pulse_shape2(t,args)
    v = x*y
    return v

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

print(sigma_x1)
print(sigma_y1)
print(sigma_z1)

sigma_x2 = tensor(qeye(2), sigmax())
sigma_y2 = tensor(qeye(2), sigmay())
sigma_z2 = tensor(qeye(2), sigmaz())

II       = tensor(qeye(2), qeye(2))
XX       = tensor(sigmax(), sigmax())
YY       = tensor(sigmay(), sigmay())
ZZ       = tensor(sigmaz(), sigmaz())

#Define Hamiltonian of 2-qubit system


#   OPTION A ---- SELF-INFERRED SUPERPOSITION OF SINGLE-QUBIT DRIVES
# H_0 = Qobj(0.5 * w_q1 * sigma_z1 + 0.5*w_q2*sigma_z2 + 0.5* g*(sigma_y1*sigma_y2))
# H_d1 = [r1*sigma_y1 *amp_1,v_t1]
# H_d2 = [r2*sigma_y2 *amp_2,v_t2]

#   OPTION B ---- BARNABY'S NOTES
delta = w_d1 - w_d2
H_0 = n1 + n2 +\
        0.5*(a_q1*(a1.dag() * a1.dag() * a1 * a1) +\
        a_q2*(a2.dag() * a2.dag() * a2 * a2)) -\
        0.5* r1*(Q1*(a1 + a1.dag()) + 1j*I1*(a1-a1.dag())) -\
        0.5*r2*(Q2*(a2 + a2.dag()) + 1j*I2*(a2-a2.dag()))

H_d1 = [g*a1*a2.dag(), lambda t, *args: np.exp(-1j*delta*t)]
H_d2 = [g*a1.dag()*a2, lambda t, *args: np.exp(1j*delta*t)]

H = [H_0, H_d1, H_d2]
print(H[0][0].shape)

# H = H_0


# %% QPT over unknown quantum process  ###########################
U_psi_real = qutip.propagator(H, t, args=args)                    #List of lists due to time dependence.
# print(U_psi_real)
U_psi_real_T = U_psi_real[nT-1]                           #Take entry from last time step
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
fig = qpt_plot_combined(chi_ideal_iSWAP, op_label, r'$i$SWAP')
# plt.show()


#%% Cphase #### #######################

#Define U: TODO: convert to tensor product form
U_psi_cphase = Qobj([[1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, -1]])


U_rho_cphase = spre(U_psi_cphase) * spost(U_psi_cphase.dag())
chi_ideal_cphase = qpt(U_rho_cphase, op_basis)

fig = qpt_plot_combined(chi_ideal_cphase, op_label, r'$CPHASE')
# plt.show()


#%% X Gate  ###########################
U_psi_CNOT = Qobj([[1, 0, 0, 0],
       [0, 1, 0, 0],
        [0, 0, 0, 1],
        [0, 0, 1, 0]])

U_rho_CNOT = spre(U_psi_CNOT) * spost(U_psi_CNOT.dag())
chi_ideal_CNOT = qpt(U_rho_CNOT, op_basis)
fig = qpt_plot_combined(chi_ideal_CNOT, op_label, r'$CPHASE')
# plt.show()


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
