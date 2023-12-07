"""
File created to visualize fidelity landscape, as a function of
drive amplitudes.

Things observed (for H/X):
- periodicity: every 100 in Q1_p
- max F: H has it at 0.5 what, and X at ~1.
- local maxima vs global maxima: global-ish in Q1_p, but I1_p has local
"""

import numpy as np
from matplotlib import pyplot as plt
from fidelity_function import *
from mpl_toolkits import mplot3d

objective = fidelity_iSWAP

I1 = np.linspace(0,250,20)
Q1 = np.linspace(0,250,20)
I2 = np.linspace(0,250,20)
Q2 = np.linspace(0,250,20)

I1_mesh, Q1_mesh = np.meshgrid(I1, Q1)
I2_mesh, Q2_mesh = np.meshgrid(I2, Q2)
# %% Plot F vs I1, Q1
fig = plt.figure()
ax = plt.axes(projection='3d')
zdata = ([[objective(**{'I1_p': i1, 'Q1_p': q1, 'I2_p': 0, 'Q2_p':0}) for i1 in I1] for q1 in Q1])
xdata = I1_mesh
ydata = Q1_mesh
ax.plot_surface(np.array(xdata), np.array(ydata), np.array(zdata), cstride=1, cmap='viridis', edgecolor='none')

ax.set_xlabel('I1_p')
ax.set_ylabel('Q1_p')
ax.set_zlabel('Gate Fidelity')

plt.show()


# %% Plot F vs I2, Q2
fig = plt.figure()
ax = plt.axes(projection='3d')
zdata = ([[objective(**{'I1_p': 0, 'Q1_p': 0, 'I2_p': i2, 'Q2_p':q2}) for i2 in I2] for q2 in Q2])
xdata = I2_mesh
ydata = Q2_mesh
ax.plot_surface(np.array(xdata), np.array(ydata), np.array(zdata), cstride=1, cmap='viridis', edgecolor='none')

ax.set_xlabel('I2_p')
ax.set_ylabel('Q2_p')
ax.set_zlabel('Gate Fidelity')

plt.show()












