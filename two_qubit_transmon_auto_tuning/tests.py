from fidelity_function import *
from qutip import *
import numpy as np
from matplotlib import pyplot as plt

"""
fine for cphase. x gate questionable since redefining GaussianPulse. look into that next -- do power rabi again.
"""
# %% Virtual Z gates to bring Cphase to expected form
# U is accepted as operating on the state in vector form.
def virtual_Z_cphase(U):
    diags = np.diag(U, k = 0)
    phi =  np.angle(diags)

    Z1 = Qobj([[1, 0, 0, 0],
              [0, 1, 0, 0],
              [0, 0, np.exp(-1j * phi[2]), 0],
              [0, 0, 0, np.exp(-1j * phi[2])]])

    Z2 =  Qobj([[1, 0, 0, 0],
              [0, np.exp(-1j * phi[1]), 0, 0],
              [0, 0, 1, 0],
              [0, 0, 0, np.exp(-1j * phi[1])]])

    #After applying virtual Z gates
    U_c = Z2*Z1*U

    return U_c


#%%
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
    t_exp  = 1000, #CPHASE
    # t_exp = 81, #X

    gate_type = 'Cphase',

    drive_shape = 'Gaussian')

kwargs = {'I1p':
              # { 'amp': 7.296, #X
                { 'amp': 0,
                'center': exp.t_exp/2,
                'std': exp.t_exp/6},
          'Q1p':
              # {'amp': 151.2, #X
                { 'amp': 0,
               'center': exp.t_exp / 2,
               'std': exp.t_exp / 6},
          'I2p':
              # { 'amp':15.76, #X
                { 'amp': 0,
                'center': exp.t_exp/2,
                'std': exp.t_exp/6},
          'Q2p':
          # {'amp': 12.61, #X
                { 'amp': 0,
                'center': exp.t_exp/2,
                'std': exp.t_exp/6},

          # 'drive': 'Q',

          'CPHASE': {
          # 'ratio': 0.7676,
          # 'ratio': 0.9979,
          'ratio':0.592,                #GIVES 0.9999 FIDELITY
          # 'w12': exp.qubits[0].w_ex12
          # 'w12': 32.47
          'w12':33.68                   #GIVES 0.9999 FIDELITY
          }

          }


#%% Do the flattening
kwargs = flatten_dict(kwargs)
#%%

if exp.gate_type == 'Cphase':
    x = exp.fidelity_Cphase(**kwargs)
    np.set_printoptions(precision=3, suppress=True)
    U_psi_Cphase = Qobj([[1, 0, 0, 0],
                         [0, 1, 0, 0],
                         [0, 0, 1, 0],
                         [0, 0, 0, -1]])
    U_rho_Cphase = spre(U_psi_Cphase) * spost(U_psi_Cphase.dag())
    chi_ideal_Cphase = qpt(U_rho_Cphase, op_basis)
    # fig2 = qpt_plot_combined(qpt(U_rho_Cphase, op_basis), [["i", "x", "y", "z"]] * 2, "ideal")
    # plt.show()
else:
    x = exp.fidelity_X(**kwargs)
    U_psi_X = tensor(sigmax(), qeye(2))
    U_rho_X = spre(U_psi_X) * spost(U_psi_X.dag())
    chi_ideal_X = qpt(U_rho_X, op_basis)
    fig2 = qpt_plot_combined(qpt(U_rho_X, op_basis), [["i", "x", "y", "z"]] * 2, "ideal")
    plt.show()

#%% Power rabi

# exp.power_rabi(**kwargs)

#%% Testing DRAG

kwargs_drag = {
    'GaussianParams': {
    'amp': 50,
    'center': 500,
    'std': 100
    },
    'lamb': 0.75,
    'alpha': 1
}

# %% Checking Gaussian shape
# times = np.linspace(0,1000,400)
# g = []
# for t in times:
#     g.append(GaussianPulse(**kwargs_drag['GaussianParams'])(t))
# plt.plot(times, g)
# plt.show


# t = 50
# drag = DRAG(**kwargs_drag)
# gauss = GaussianPulse(**kwargs_drag['GaussianParams'])

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
