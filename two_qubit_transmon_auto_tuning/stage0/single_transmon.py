#Import relevant modules
import numpy as np
import matplotlib.pyplot as plt
from qutip import *

# Define the qubit's frequency (transition frequency) and anharmonicity
omega_q = 1.0
anharm = 0.2

#Define number of energy levels (dim = 3 if accounting for leakage to 3rd level)
dim = 2

# Create the qubit basis states; 
g = basis(dim, 0)  # Ground state |0⟩
e = basis(dim, 1)  # Excited state |1⟩

# Define operators
c1: Qobj = destroy(dim) 
c2: Qobj = c1.dag()
c2c1 = c1*c2

# Create the qubit Hamiltonian
H_drift = (omega_q*c2c1) + 0.5*anharm*c2c1*(c2c1-qeye(dim))

# Initial state (e.g., qubit in the ground state)
psi0 = g
rho0 = ket2dm(psi0)

print(psi0)

# Time points for simulation
tlist = np.linspace(0, 10, 100)

# Time evolution of the qubit
result = mesolve(H_drift, rho0, tlist, [], [sigmax(), sigmay(), sigmaz()])

# Plot the results
plt.figure()
plt.plot(tlist, result.expect[0], label="X-expectation")
plt.plot(tlist, result.expect[1], label="Y-expectation")
plt.plot(tlist, result.expect[2], label="Z-expectation")
plt.xlabel("Time")
plt.ylabel("Expectation Values")
plt.legend()
plt.show()
