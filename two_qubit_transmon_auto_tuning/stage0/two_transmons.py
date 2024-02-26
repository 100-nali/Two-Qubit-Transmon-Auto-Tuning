#Import relevant modules
import numpy as np
import matplotlib.pyplot as plt
from qutip import *

#Define number of dimensions
dim: int = 2

#Define parameters for qubit 1
anharm_1: float = 0.2
omega_q1: float = 1.0
r1: float       = 0.06 #Drive coupling strength for qubit 1
I1              = 0
Q1              = 0
omega_d1:float  = 1.0


#Define parameters for qubit 2
anharm_2: float = 0.2
omega_q2: float = 1.0
r2: float       = 0.06 #Drive coupling strength for qubit 2
I2              = 0
Q2              = 0
omega_d2:float  = 1.0

#Define Inter-coupling strength
g = 0.03

#Define operators
a1: Qobj = tensor(destroy(dim), qeye(dim)) ## meant for a 2-qubit system.
a2: Qobj = tensor(qeye(dim), destroy(dim)) ##

n1: Qobj = a1.dag() * a1
n2: Qobj = a2.dag() * a2

#Define Pauli Matrices
I = basis(dim, 0) * basis(dim, 0).dag() + basis(dim, 1) * basis(dim, 1).dag()
sigma_x = basis(dim, 0) * basis(dim, 1).dag() + basis(dim, 1) * basis(dim, 0).dag()
sigma_y = -1j * basis(dim, 0) * basis(dim, 1).dag() + 1j * basis(dim, 1) * basis(dim, 0).dag()
sigma_z = basis(dim, 0) * basis(dim, 0).dag() - basis(dim, 1) * basis(dim, 1).dag()

sigma_x1 = tensor(sigma_x, I)
sigma_x2 = tensor(I, sigma_x)

sigma_y1 = tensor(sigma_y, I)
sigma_y2 = tensor(I, sigma_y)

sigma_z1 = tensor(sigma_z, I)
sigma_z2 = tensor(I, sigma_z)

#Define drift Hamiltonian for qubit 1
H1_drift = [
        (anharm_1 / 2) * (a1.dag() * a1.dag() * a1 * a1) - omega_q1 * (a1.dag() * a1),
        [a1.dag() * a1, omega_q1]
    ]

#Define drive Hamiltonian for qubit 1
H1_drive = [
    [-r1 * 1j * (a1 - a1.dag()) / 2., I1],
    [-r1 * (a1 + a1.dag()) / 2., Q1]
]

#Define drift Hamiltonian for qubit 2
H2_drift = [
        (anharm_2 / 2) * (a2.dag() * a2.dag() * a2 * a2) - omega_q2 * (a2.dag() * a2),
        [a2.dag() * a2, omega_q2]
    ]


#Define drive Hamiltonian for qubit 2
H2_drive = [
    [-r2 * 1j * (a2 - a2.dag()) / 2., I2],
    [-r2 * (a2 + a2.dag()) / 2., Q2]
]

#Define interaction Hamiltonian
delta_d = omega_d1 - omega_d2

#The interaction hamiltonian (gate specific?)
H_int = [
    [a1 * a2.dag(), g * np.exp(-1j * delta_d * t)],
    [a1.dag() * a2, g * np.exp(1j * delta_d * t)]
]

#Overal H
H = [*H1_drift, *H1_drive, *H2_drift, *H2_drive, *H_int]


#Define initial states
psi_0 = [
    basis(3, 0),
    basis(3, 1),
    (basis(3, 0) + basis(3, 1)).unit(),
    (basis(3, 0) - basis(3, 1)).unit(),
    (basis(3, 0) + 1j * basis(3, 1)).unit(),
    (basis(3, 0) - 1j * basis(3, 1)).unit()
    ]

#Define t
t = np.linspace(0, 10, 100)

#Bloch sphere
b = Bloch()

#Evolve states using Hamiltonians
for psi0 in psi_0:
    res = mesolve(H, psi0, t, [], [sigma_x1, sigma_y1, sigma_z1, sigma_x2, sigma_y2, sigma_z2, n1, n2])
    expectations = np.stack(res.expect, axis=-1)

    #Plot on Bloch Sphere
    b.add_points(expectations.T) ####  color is function of initial state? Plots all times together.
    b.show()