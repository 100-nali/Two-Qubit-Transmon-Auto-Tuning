#Imports
import numpy as np
from qutip import *
import matplotlib.pyplot as plt
import itertools
from dataclasses import dataclass
import math


@dataclass
class GaussianPulse:
    amplitude: float
    std: float
    center: float

    def __call__(self, t, *args, **kwargs):
        return self.amplitude * np.exp(-(t - self.center) ** 2 / (2 * self.std ** 2)) / (
                self.std * np.sqrt(2 * np.pi))

#Define number of energy levels
dim = 2

#Define time-stamps
t = np.linspace(0, 10, 100)

#Define qubit 1 parameters
w_q1: float     = 1.0
a_q1: float     = 1.0
r1:float  = 1.0

#Define qubit 2 parameters
w_q2: float     = 1.0
a_q2: float     = 1.0
r2:float  = 1.0

#Define inter-qubit coupling strength
g: float        = 1.0

#Define pulse parameters - q1
w_d1 = 1.0
amp_1 = 1.0
shape_1 = lambda t, *args: math.sin(t)


#Define pulse parameters - q2
w_d2 = 1.0
amp_2 = 1.0
shape_2 = lambda t, *args: math.sin(t)

#phi
phi1 = np.pi/4
phi2 = np.pi/4


#Define a1, a2.
a1 = tensor(destroy(dim), qeye(dim))
a2 = tensor(qeye(dim), destroy(dim))

#Define op_basis
op_basis = [[qeye(2), sigmax(), sigmay(), sigmaz()]] * 2

#Define plot labels
op_label = [["i", "x", "y", "z"]] * 2


# %% Define Hamiltonian of 2-qubit system.
# H_q1 = (w_q1 * a1.dag()* a1) + (0.5 * a_q1 * a1.dag() * a1.dag() * a1 * a1)
# H_q2 = (w_q2 * a2.dag()* a2) + (0.5 * a_q2 * a2.dag() * a2.dag() * a2 * a2)
# H_q_int = -1 * g * (a1 - a1.dag())*(a2 - a2.dag())
# H_d1 = [ 
#     [0, lambda t, *args: -0.5 * r1 * amp_1 * shape_1 * np.exp(1j*(w_q1-w_d1)*t + phi1) ],
#     [lambda t, *args: -0.5 * r1 * amp_2 * shape_2 *np.exp(-1j*((w_q1-w_d1)*t + phi1)), 0]
#     ]
# H_d2 = [ 
#     [0, lambda t, *args: -0.5 * r2 * np.exp(1j*(w_q2-w_d2)*t + phi2) ],
#     [lambda t, *args: -0.5 * r2 * np.exp(-1j*((w_q2-w_d2)*t + phi2)), 0]]
# # H_d_int = ?

# H = [*H_q1, *H_q2 , *H_q_int,  *H_d1, *H_d2]





# a1 = tensor(destroy(2), qeye(2)) 
# a2 = tensor(qeye(2), destroy(2)) 
#  # the qubit 1 drift hamiltonian
# H1 = [
#     [(a_q1 / 2) * (a1.dag() * a1.dag() * a1 * a1) - w_d1 * (a1.dag() * a1), 0],
#     [a1.dag() * a1,                                                      w_q1]
# ]

# # the qubit 1 drive hamiltonian
# H1d = [
#     [-r1 * 1j * (a1 - a1.dag()) / 2.,       math.cos(phi1)],
#     [-r1 * (a1 + a1.dag()) / 2.,            math.sin(phi1)]
# ]

# # the qubit 2 drift hamiltonian
# H2 = [
#     [(a_q2 / 2) * (a2.dag() * a2.dag() * a2 * a2) - w_d2 * (a2.dag() * a2), 0],
#     [(a2.dag() * a2),                                                    w_q2]
# ]

# H2d = [
#     [-r2 * 1j * (a2 - a2.dag()) / 2.,            math.cos(phi2)],
#     [-r2 * (a2 + a2.dag()) / 2.,                 math.sin(phi2)]
# ]

# delta_d = w_d1 - w_d2

# # the interaction hamiltonian (includes interaction with drive, i think!)
# H_int = [
#     [a1 * a2.dag(),         lambda t, *args: g * np.exp(-1j * delta_d * t)],
#     [a1.dag() * a2,         lambda t, *args: g * np.exp(1j * delta_d * t)]
# ]

# H = [*H1, *H1d, *H2, *H2d, *H_int]
# print(H)
# print(shape(H))
# H = [H1, H1d, H2, H2d, H_int]
# print(shape(H))

# Create identity operators for the specified dimension
identity = qutip.qeye(dim)
a1 = qutip.tensor(qutip.destroy(dim), identity)
a2 = qutip.tensor(identity, qutip.destroy(dim))

# Define the Hamiltonian terms (adjust the operators and dimensions as needed)
H1 = (a_q1 / 2) * (a1.dag() * a1.dag() * a1 * a1) - w_d1 * (a1.dag() * a1)
H1d = -r1 * (a1 + a1.dag()) * (math.cos(phi1) - 1j * math.sin(phi1))
H2 = (a_q2 / 2) * (a2.dag() * a2.dag() * a2 * a2) - w_d2 * (a2.dag() * a2)
H2d = -r2 * (a2 + a2.dag()) * (math.cos(phi2) - 1j * math.sin(phi2))
delta_d = w_d1 - w_d2
H_int = a1 * a2.dag() + a1.dag() * a2

# Combine Hamiltonian terms into a list (if needed)
H = [H1, H1d, H2, H2d, [H_int, lambda t, *args: g * np.exp(-1j * delta_d * t)]]


# %% QPT over unknown quantum process

# Find chi_real ###################
U_psi = qutip.propagator(H,t) #List of lists due to time dependence.
# U_psi = ### Some transformative operation that converts to 'time-independent' form. (single matrix)
U_rho = spre(U_psi) * spost(U_psi.dag())
chi_real = qpt(U_rho, op_basis)
# for U_psi_i in U_psi:
#     U_rho = spre(U_psi_i) * spost(U_psi_i.dag())
#     chi_real = qpt(U_rho, op_basis)
fig = qpt_plot_combined(chi_real, op_label, r'$Actual Process$')
plt.show()
#%% iSWAP

# Define Pauli matrices and the identity operator for a qutrit
sigma_x = basis(3, 0) * basis(3, 1).dag() + basis(3, 1) * basis(3, 0).dag()
sigma_y = -1j * basis(3, 0) * basis(3, 1).dag() + 1j * basis(3, 1) * basis(3, 0).dag()
sigma_z = basis(3, 0) * basis(3, 0).dag() - basis(3, 1) * basis(3, 1).dag()
identity = qeye(3)

#  Define the iSWAP gate 
U_psi =  Qobj([[1, 0, 0, 0],
                [0, 0, 1j, 0],
                [0, 1j, 0, 0],
                [0, 0, 0, 1]])        
U_rho = spre(U_psi) * spost(U_psi.dag())


#Find ideal process matrix
chi_ideal = qpt(U_rho, op_basis)
fig = qpt_plot_combined(chi_ideal, op_label, r'$i$SWAP')
plt.show()


#%% Cphase

#Define U:
U_psi = Qobj([[1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, -1]])

U_rho = spre(U_psi) * spost(U_psi.dag())
print('U_rho for Cphase below!')
print(U_rho)
chi_ideal = qpt(U_rho, op_basis)
fig = qpt_plot_combined(chi_ideal, op_label, r'$CPHASE')
plt.show()
