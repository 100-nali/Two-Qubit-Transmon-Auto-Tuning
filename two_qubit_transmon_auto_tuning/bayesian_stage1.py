#%% Imports
import numpy as np
from matplotlib import pyplot as plt
import scipy
from sklearn.gaussian_process import GaussianProcessRegressor
import fidelity_function as fid

#%% Notes
"""
Choose a single-pulse gate.
Assume Gaussian envelope.
Assume time of measurement is fixed.
Optimizing over three general parameters: Amplitude, Mean, and Width (std).
and over 4 drive components.
i.e., total number of parameters to optimize over = 4*3 = 12.

"""


#%% Define objective function -> gate fidelity as a function of the drives.
objective = fid.fidelity_fn

#%% Define surrogate function -> (RBF GP)
"""
X is also 12D since surrogate is an approximation of objective which depends on 12 parameters.
"""
def surrogate(model, X):
    return model.predict(X, return_std = True)

#%% Define acqusition function -> EI
