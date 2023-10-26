import numpy as np
import matplotlib.pyplot as plt
from qutip import *

# Define the qubit's frequency (transition frequency) and anharmonicity
omega_q = 1.0
anharm = 0.2

#Define number of energy levels
dim = 3

# Create the qubit basis states; 3 accounts for the 3rd energy level that can be leaked into
g = basis(3, 0)  # Ground state |0⟩
e = basis(3, 1)  # Excited state |1⟩

# Define operators
c1 = destroy(dim)
c2 = c1.dag
c2c1 = tensor(c2,c1)

# Create the qubit Hamiltonian
# H_qubit = (omega_qubit * (sigmaz() / 2)) + (anharmonicity * sigmaz() * sigmaz() / 2)
H_drift = (omega_q*c2c1) + 0.5*anharm*c2c1*(c2c1-1)

# Initial state (e.g., qubit in the ground state)
psi0 = g

# Time points for simulation
tlist = np.linspace(0, 10, 100)

# Time evolution of the qubit
result = mesolve(H_drift, psi0, tlist, [], [sigmax(), sigmay(), sigmaz()])

# Plot the results
plt.figure()
plt.plot(tlist, result.expect[0], label="X-expectation")
plt.plot(tlist, result.expect[1], label="Y-expectation")
plt.plot(tlist, result.expect[2], label="Z-expectation")
plt.xlabel("Time")
plt.ylabel("Expectation Values")
plt.legend()
plt.show()
