#%% Imports

from bayes_opt import BayesianOptimization
from bayes_opt import SequentialDomainReductionTransformer
from bayes_opt.util import UtilityFunction
from fidelity_function import *
from functools import partial



#%% Notes
"""
Choose a single-pulse gate.
Assume Gaussian envelope.
Assume time of measurement is fixed.
Optimizing over one general parameter: Amplitude.
i.e., total number of parameters to optimize over = 4.

"""

#%% Define parameters
exp = Experiment(
    [Qubit(
    i = 1,
    w_q = 6 * 2 * np.pi,
    a_q = -0.3 * 2 * np.pi,
    r = 0.01 * 2 * np.pi,
    w_d = 6 * 2 * np.pi,
    w_ex12 = 5.3*2*np.pi),

    Qubit(
    i = 2,
    w_q = 5 * 2 * np.pi,
    a_q = -0.3 * 2 * np.pi,
    r = 0.01 * 2 * np.pi,
    w_d = 5 * 2 * np.pi,
    w_ex12=5.3 * 2 * np.pi)
    ],

    g = 0.005 * 2 * np.pi, ### ? 20MHz for optimal trajectory?

    # t_exp = 64,
    t_exp  = 1000, #CPHASE
    # t_exp = 81, #X

    gate_type = 'Cphase',

    drive_shape = 'Gaussian')


#%% Define Bayes opt parameters

# iterations = 300
# init_points = 300
iterations = 20
init_points = 20

#%% Define objective function -> gate fidelity as a function of the drives.
objective = exp.fidelity_Cphase

#%% input the parameters that are not to be iterated over into objective

kwargs = {'I1p': {'amp': 0, 'center': exp.t_exp/2, 'std': exp.t_exp/6},
    'I2p': {'amp': 0, 'center': exp.t_exp/2, 'std': exp.t_exp/6},
    'Q1p': {'amp': 0, 'center': exp.t_exp/2, 'std': exp.t_exp/6},
    'Q2p': {'amp': 0, 'center': exp.t_exp/2, 'std': exp.t_exp / 6}}
flattened_kwargs = flatten_dict(kwargs)
objective = partial(objective, **flattened_kwargs)


#%% Bayesian Optimization

nested_bounds = {
    # 'I1p': {'amp': (0, 1), 'center': (exp.t_exp/2.5,exp.t_exp/1.5) , 'std': (exp.t_exp/7.5, exp.t_exp/4.5)},
    # 'I2p': {'amp': (0, 1), 'center': (exp.t_exp/2.5,exp.t_exp/1.5) , 'std': (exp.t_exp/7.5, exp.t_exp/4.5)},
    # 'Q1p': {'amp': (0, 1), 'center': (exp.t_exp/2.5,exp.t_exp/1.5) , 'std': (exp.t_exp/7.5, exp.t_exp/4.5)},
    # 'Q2p': {'amp': (0, 1), 'center': (exp.t_exp/2.5,exp.t_exp/1.5) , 'std': (exp.t_exp/7.5, exp.t_exp/4.5)},
    'CPHASE': {'ratio': (0, 1), 'w12': (5*np.pi*2, 6*np.pi*2)}
}


flattened_bounds = flatten_dict(nested_bounds)


bounds_transformer = SequentialDomainReductionTransformer(minimum_window=0.5)

optimizer = BayesianOptimization(
    f=objective,
    pbounds=flattened_bounds,
    random_state=1,
    bounds_transformer=bounds_transformer
)

# Step 4: Access the scores and iteration index
iteration_index = []
scores = []

optimizer.maximize(
    init_points=init_points, #300 was better for both of these
    n_iter=iterations,
    # acquisition_function= UtilityFunction(kind = 'ei')
)

target_vals = [res.get('target') for res in optimizer.res]
plt.figure()
plt.plot(target_vals)
plt.show()