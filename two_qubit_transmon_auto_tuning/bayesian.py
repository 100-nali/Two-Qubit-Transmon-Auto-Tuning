#%% Imports

from bayes_opt import BayesianOptimization
from bayes_opt import SequentialDomainReductionTransformer
from bayes_opt.util import UtilityFunction
from fidelity_function import *
from functools import partial

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
    t_exp  = 800, #CPHASE
    # t_exp = 81, #X

    gate_type = 'Cphase',

    drive_shape = 'Gaussian')


#%% Define Bayes opt parameters

# iterations = 300
# init_points = 300
# iterations = 30
# init_points = 10

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

# bounds_transformer = SequentialDomainReductionTransformer(minimum_window=0.5)

def hyperparameter_objective(acquisition_weight, num_iterations, num_init_points, alpha):

    utility_function = UtilityFunction(kind='ei', kappa=acquisition_weight)

    optimizer = BayesianOptimization(
        f=objective,
        pbounds=flattened_bounds,
        random_state=1,
        verbose = 0,
        # bounds_transformer=bounds_transformer,
        allow_duplicate_points = True
    )

    gp_params = {
        'alpha': alpha
    }

    optimizer.set_gp_params(**gp_params)

    optimizer.maximize(
        init_points= int(num_init_points),
        n_iter= int(num_iterations),
        acquisition_function= utility_function,
    )

    return optimizer.max['target']

# Hyperparameter tuning ranges
hyperparameter_bounds = {
    'acquisition_weight': (1, 10),
    'num_iterations': (10, 30),
    'num_init_points': (5, 20),
    'alpha': (1e-5, 1e-1)
}

# Create the hyperparameter tuning optimizer
hyperparameter_optimizer = BayesianOptimization(
    f=hyperparameter_objective,
    pbounds=hyperparameter_bounds,
    random_state=1
)

# Run hyperparameter tuning
hyperparameter_optimizer.maximize(
    init_points=5,
    n_iter=10
)

# Get the best hyperparameters
best_hyperparameters = hyperparameter_optimizer.max

# Use the best hyperparameters for the main Bayesian optimization
acquisition_weight = best_hyperparameters['params']['acquisition_weight']
iterations = int(best_hyperparameters['params']['num_iterations'])
init_points = int(best_hyperparameters['params']['num_init_points'])
alpha = best_hyperparameters['params']['alpha']


optimizer = BayesianOptimization(
    f = objective,
    pbounds=flattened_bounds,
    random_state=1,
    # bounds_transformer=bounds_transformer,
    allow_duplicate_points=True
)

# standard_optimizer = BayesianOptimization(
#     f = objective,
#     pbounds=flattened_bounds,
#     # verbose=0,
#     random_state=1,
# )

# Step 4: Access the scores and iteration index

gp_params = {'alpha': alpha}

optimizer.set_gp_params(**gp_params)

utility_function = UtilityFunction(kind = 'ei', kappa = acquisition_weight)

optimizer.maximize(
    init_points=init_points,
    n_iter=iterations,
    # kappa=acquisition_weight
    acquisition_function= utility_function
)

# standard_optimizer.maximize(
#     init_points=init_points,
#     n_iter=iterations,
# )

plt.plot(optimizer.space.target, label='Mutated Optimizer')
# plt.plot(standard_optimizer.space.target, label='Standard Optimizer')
# plt.legend()


iteration_index = []
scores = []
target_vals = [res.get('target') for res in optimizer.res]
plt.figure()
plt.plot(target_vals)
plt.show()


iteration_index = []
scores = []
target_vals = [res.get('target') for res in hyperparameter_optimizer.res]
plt.figure()
plt.plot(target_vals)
plt.show()