#Imports
import numpy as np
from qutip import *
import matplotlib.pyplot as plt
import itertools

#Define number of energy levels
dim = 3

#Define time-stamps
t = np.linspace(0, 10, 100)

#Define basis states
# g = basis(3, 0)
# e = basis(3, 1)
# f = basis(3, 2)
# op_basis_list = [g, e, f]

# op_basis_list = [
# basis(3, 0),
# basis(3, 1),
# (basis(3, 0) + basis(3, 1)).unit(),
# (basis(3, 0) - basis(3, 1)).unit(),
# (basis(3, 0) + 1j * basis(3, 1)).unit(),
# (basis(3, 0) - 1j * basis(3, 1)).unit()
# ]



#Define qubit 1 parameters
w_q1: float     = 1.0
a_q1: float     = 1.0

#Define qubit 2 parameters
w_q2: float     = 1.0
a_q2: float     = 1.0

#Define inter-qubit coupling strength
g: float        = 1.0

#Define a1, a2.
a1 = tensor(destroy(dim), qeye(dim))
a2 = tensor(qeye(dim), destroy(dim)) 

#Define Hamiltonian of 2-qubit system.
H_q1 = (w_q1 * a1.dag()* a1) + (0.5 * a_q1 * a1.dag() * a1.dag() * a1 * a1)
H_q2 = (w_q2 * a2.dag()* a2) + (0.5 * a_q2 * a2.dag() * a2.dag() * a2 * a2)
H_int = -1 * g * (a1 - a1.dag())*(a2 - a2.dag())
H = (H_q1 + H_q2 + H_int)

#%% iSWAP

# #Define U: (Is lacking a parameter that defines gate type...)
# # U = (qutip.propagator(H,t))

# # Define Pauli matrices and the identity operator for a qutrit
sigma_x = basis(3, 0) * basis(3, 1).dag() + basis(3, 1) * basis(3, 0).dag()
sigma_y = -1j * basis(3, 0) * basis(3, 1).dag() + 1j * basis(3, 1) * basis(3, 0).dag()
sigma_z = basis(3, 0) * basis(3, 0).dag() - basis(3, 1) * basis(3, 1).dag()
identity = qeye(3)  

# # Define the iSWAP gate for a qutrit
# iSWAP_qutrit = tensor(identity, identity, identity) + 1j * tensor(sigma_x, sigma_x, identity)
# U_psi = iSWAP_qutrit.full()
# print(U)

U_psi =  Qobj([[1, 0, 0, 0],
                [0, 0, 1j, 0],
                [0, 1j, 0, 0],
                [0, 0, 0, 1]])
U_rho = spre(U_psi) * spost(U_psi.dag())
op_basis = [[qeye(2), sigmax(), sigmay(), sigmaz()]] * 2
op_label = [["i", "x", "y", "z"]] * 2
chi = qpt(U_rho, op_basis)
print(chi)
fig = qpt_plot_combined(chi, op_label, r'$i$SWAP')
plt.show()



#%% Cphase

#Define U:

#Find ideal chi: