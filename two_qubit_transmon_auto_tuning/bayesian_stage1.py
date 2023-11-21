#%% Imports
import numpy as np
from matplotlib import pyplot as plt
# import scipy
# from sklearn.gaussian_process import GaussianProcessRegressor
# from sklearn.gaussian_process.kernels import RBF
from bayes_opt import BayesianOptimization
from bayes_opt.util import UtilityFunction
import fidelity_function as fid
from bayes_opt.event import Events, DEFAULT_EVENTS


#%% Notes
"""
Choose a single-pulse gate.
Assume Gaussian envelope.
Assume time of measurement is fixed.
Optimizing over one general parameter: Amplitude.
i.e., total number of parameters to optimize over = 4.

"""

# #%% Define seed
# np.random.seed(10)

iterations = 300
init_points = 300

#%% Define objective function -> gate fidelity as a function of the drives.
objective = fid.fidelity_iSWAP

#%% Bayesian Optimization

"""""
these bounds here are for the case that nothing is known about what the gate
should look like. However, maybe prior knowledge can be taken about
gate amplitudes (in an ideal case), and bounds taken around that.
"""""
pbounds = {'I1_p': (0, 100), 'Q1_p': (0, 100), 'I2_p': (0, 100), 'Q2_p': (0, 100)}

optimizer = BayesianOptimization(
    f=objective,
    pbounds=pbounds,
    random_state=1,
)

# Step 4: Access the scores and iteration index
iteration_index = []
scores = []

# optimizer.probe(
#     params = {'I1_p': 78, 'Q1_p': 0, 'I2_p': 0, 'Q2_p': 0},
#     lazy = True
#     )

optimizer.maximize(
    init_points=init_points, #300 was better for both of these
    n_iter=iterations,
    acquisition_function= UtilityFunction(kind = 'ei')
)

target_vals = [res.get('target') for res in optimizer.res]
plt.figure()
plt.plot(target_vals)
plt.show()