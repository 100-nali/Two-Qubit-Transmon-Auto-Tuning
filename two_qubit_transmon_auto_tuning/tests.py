from fidelity_function import *
from qutip import *
import numpy as np
from matplotlib import pyplot as plt

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
    w_ex12=5 * 2 * np.pi)
    ],

    g = 0.005 * 2 * np.pi, ### ? 20MHz for optimal trajectory?

    # t_exp = 64,
    t_exp  = 1000,
    # t_exp = 81,

    gate_type = 'Cphase',

    drive_shape = 'Gaussian')

kwargs = {'I1_p':
              # { 'amp': 7.296, #X
                { 'amp': 0,
                'center': exp.t_exp/2,
                'std': exp.t_exp/10},
          'Q1_p':
              # {'amp': 151.2, #X
                { 'amp':    0,
               'center': exp.t_exp / 2,
               'std': exp.t_exp / 10},
          'I2_p':
              # { 'amp':15.76, #X
                { 'amp': 0,
                'center': exp.t_exp/2,
                'std': exp.t_exp/60},
          'Q2_p':
          # {'amp': 12.61, #X
                { 'amp': 0,
                'center': exp.t_exp/2,
                'std': exp.t_exp/60}}

#%%

#
# x = exp.fidelity_Cphase(**kwargs)
# np.set_printoptions(precision=3, suppress=True)
# fid, U = x

#%% Virtual Z gates to bring Cphase to expected form

# Z2 = Qobj([[1.   +0.j  ,  0.   +0.j  ,  0.   +0.j  ,  0.   +0.j   ],
#  [0.   +0.j  ,   np.exp(-1j*0.46454206314509705), 0.   +0.j   , 0.   +0.j   ],
#  [0.   +0.j  ,  0.   +0.j   , 1.   +0.j  ,  0.   +0.j   ],
#  [0.   +0.j  ,  0.   +0.j   , 0.   +0.j  ,  np.exp(-1j*0.46454206314509705)]]  )
#
# Z1 = Qobj([[1.   +0.j ,   0.   +0.j  ,  0.   +0.j  ,  0.   +0.j   ],
#  [0.   +0.j  ,  1.   +0.j ,   0.   +0.j  ,  0.   +0.j   ],
#  [0.   +0.j  ,  0.   +0.j  ,   np.exp(-1j*-2.432640966569607) ,0.   +0.j   ],
#  [0.   +0.j  ,  0.   +0.j ,   0.   +0.j ,    np.exp(-1j*-2.432640966569607)]])
#
#
# #%% After applying virtual Z gates
# U_c = Z2*Z1*U


#%% Testing DRAG


kwargs_drag = {
    'GaussianParams': {
    'amp': 50,
    'center': 0,
    'std': 100
    },
    'lamb' : 0.75,
    'alpha'  : 1
}

t = 50
drag = DRAG(t, **kwargs_drag)
gauss = GaussianPulse(t, **kwargs_drag['GaussianParams'])
print(drag)
print(gauss)



#%%
# q1 = exp.qubits[0]
# q2 = exp.qubits[1]
# q1.w_q = LinearFluxTuning(
#     t0= 0,
#     y0 = q1.w_q_def,
#     t_tune = 90/3,
#     y_max = 5*np.pi*2,
#     tau = 90/3
# )

# H =  lambda t, *args: (np.exp(-1j * (q1.w_q(t)) - q2.w_q) * t)


# omegas = np.linspace(5*np.pi*2,6*np.pi*2,300)
# # omegas = [5.30*2*np.pi]
# off_diag = []
# last_entry = []
# middle = []
# for w in omegas:
#     exp.qubits[0].w_ex12 = w
#     x = exp.fidelity_X(**kwargs)
#     fid, U = x
#     U12 = U[1,2]
#     U33 = U[3,3]
#     U22 = U[2,2]
#     off_diag.append(np.abs(U12))
#     last_entry.append(np.abs(U33))
#     middle.append(np.abs(U22))
#
# plt.figure(1)
# plt.plot(omegas, off_diag)
# plt.show()
#
# plt.figure(2)
# plt.plot(omegas, last_entry)
# plt.show()
#
# plt.figure(3)
# plt.plot(omegas, middle)
# plt.show()
#
# w_tune = LinearFluxTuning(
#     t0= 0,
#     y0 = 6,
#     t_tune = 495,
#     y_max = 5,
#     tau = 10
# )
#

#
# omegas = []
# times = np.linspace(1,90,100)
#
# for t in times:
#     w = w_tune(t)
#     omegas.append(w)
#
# plt.plot(times, omegas)
# plt.show()
#



# w_tune = SquareFluxTuning(
#     t0= 0,
#     y0 = 6,
#     t_stop = 90,
#     y_max = 5
# )

# omegas = []
# times = np.linspace(1,1000,500)
#
# for t in times:
#     w = w_tune(t)
#     omegas.append(w)
#
# plt.plot(times, omegas)
# plt.show()
