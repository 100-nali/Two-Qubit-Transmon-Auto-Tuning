from fidelity_function import *
from qutip import *
import numpy as np
from matplotlib import pyplot as plt

exp = Experiment(
    [Qubit(
    i = 1,
    w_q = 5 * 2 * np.pi,
    a_q = -0.3 * 2 * np.pi,
    r = 0.01 * 2 * np.pi,
    w_d = 5 * 2 * np.pi),

    Qubit(
    i = 2,
    w_q = 6 * 2 * np.pi,
    a_q = -0.3 * 2 * np.pi,
    r = 0.01 * 2 * np.pi,
    w_d = 6 * 2 * np.pi)],

    g = 0.005 * 2 * np.pi,

    t_exp = 64,

    gate_type = 'OFFCPhase',

    drive_shape = 'Gaussian')

kwargs = {'I1_p':
              { 'amp': 7.296*0,
                'center': exp.t_exp/2,
                'std': exp.t_exp/6},
          'Q1_p':
              {'amp': 151.2*0 +50,
               'center': exp.t_exp / 2,
               'std': exp.t_exp / 6},
          'I2_p':
              { 'amp':15.76*0,
                'center': exp.t_exp/2,
                'std': exp.t_exp/6},
          'Q2_p': {'amp': 12.61*0,
                'center': exp.t_exp/2,
                'std': exp.t_exp/6}}


x = exp.fidelity_X(**kwargs)