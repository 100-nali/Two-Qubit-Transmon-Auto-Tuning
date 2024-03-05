#
# # # %% Define Hamiltonian of 2-qubit system.
# # # H_q1 = (w_q1 * a1.dag()* a1) + (0.5 * a_q1 * a1.dag() * a1.dag() * a1 * a1)
# # # H_q2 = (w_q2 * a2.dag()* a2) + (0.5 * a_q2 * a2.dag() * a2.dag() * a2 * a2)
# # # H_q_int = -1 * g * (a1 - a1.dag())*(a2 - a2.dag())
# # # H_d1 = [
# # #     [0, lambda t, *args: -0.5 * r1 * amp_1 * shape_1 * np.exp(1j*(w_q1-w_d1)*t + phi1) ],
# # #     [0, lambda t, *args: -0.5 * r1 * amp_1 * shape_1 *np.exp(-1j*((w_q1-w_d1)*t + phi1))]
# # #     ]
# # # H_d1 = r1 * amp_1 * shape_1 *  * sigmay() #
# # # H_d2
# # # H_d2 = [
# # #     [a1*a2.dag(), lambda t, *args: -0.5 * r2 * np.exp(1j*(w_q2-w_d2)*t + phi2) ],
# # #     [a1.dag() * a2], lambda t, *args: -0.5 * r2 * np.exp(-1j*((w_q2-w_d2)*t + phi2))]
# # # H_d1_d2 = tensor(H_d1, H_d2)
#
# # # H = [*H_q1, *H_q2 , *H_q_int,  *H_d1, *H_d2]
# # # H =  [*H_q1, *H_q2 , *H_q_int,  *H_d1_d2]
#
# # # print(H)
#
#
# # #%%
# # # Create identity operators for the specified dimension
# # # identity = qutip.qeye(dim)
# # # a1 = qutip.tensor(qutip.destroy(dim), identity)
# # # a2 = qutip.tensor(identity, qutip.destroy(dim))
#
# # # Define the Hamiltonian terms (adjust the operators and dimensions as needed)
# # # H1 = (a_q1 / 2) * (a1.dag() * a1.dag() * a1 * a1) - w_d1 * (a1.dag() * a1)
# # # H1d = -r1 * (a1 + a1.dag()) * (math.cos(phi1) - 1j * math.sin(phi1))
# # # H2 = (a_q2 / 2) * (a2.dag() * a2.dag() * a2 * a2) - w_d2 * (a2.dag() * a2)
# # # H2d = -r2 * (a2 + a2.dag()) * (math.cos(phi2) - 1j * math.sin(phi2))
# # # delta_d = w_d1 - w_d2
# # # H_int = a1 * a2.dag() + a1.dag() * a2
#
# # # Combine Hamiltonian terms into a list (if needed)
# # # H = [H1, H1d, H2, H2d, [H_int, lambda t, *args: g * np.exp(-1j * delta_d * t)]]
#
# # # %%
#
# # # U_rho = []
# # # for U_psi_i in U_psi:
# # #     U_rho_i = spre(U_psi_i) * spost(U_psi_i.dag())
# # #     U_rho.append(U_rho_i)
#
#
# # # print(U_rho)
# # # chi_real = qpt(U_rho, op_basis)
#
# # # %%
#
# # # @dataclass
# # # class GaussianPulse:
# # #     amplitude: float
# # #     std: float
# # #     center: float
#
# # #     def __call__(self, t, *args, **kwargs):
# # #         return self.amplitude * np.exp(-(t - self.center) ** 2 / (2 * self.std ** 2)) / (
# # #                 self.std * np.sqrt(2 * np.pi))
#
#
# # # %%
# # # %% Plot on Bloch Sphere  ###########################
# # # b = Bloch()
# # # b.make_sphere()
# # # starting_states = [
# # #             tensor(basis(2, 0), basis(2, 0)),
# # #             tensor(basis(2, 0),basis(2, 1)),
# # #             tensor(basis(2, 1),basis(2, 0)),
# # #             tensor(basis(2, 1),basis(2, 1)),
# # #         ]
# # # states = starting_states
#
# # # for U in U_psi_real:
# # #     # b.add_points()                        #add expectation values of states not the states...? w.r.t. x,y,z
# # #     for k,state in enumerate(states):
# # #         states[k] = U* states[k]
# # # b.show()
# # # plt.show()
#
# # #%%
# # #Define driving functions
# # v_1 = lambda t, *args: amp_1 * np.exp(-(t - mean_1) ** 2 / (2 * std_1 ** 2)) / (
# #                 std_1 * np.sqrt(2 * np.pi)) * np.exp(1j*(w_d1*t + phi1))
# # v_2 = lambda t, *args: amp_2 * np.exp(-(t - mean_2) ** 2 / (2 * std_2 ** 2)) / (
# #                 std_2 * np.sqrt(2 * np.pi)) * np.exp(1j*(w_d2*t + phi2))
#
#
# # (amp/ (std * np.sqrt(2 * np.pi))) * np.exp(-(t - mean) ** 2 / (2 * std ** 2))
# #     return x
#
#
#
# # # %% Plot v_1 and v_2 for a sanity check
# # """
# # plt.plot(v_1)
# # plt.show()
# # plt.plot(v_2)
# # plt.show()
#
#
# # gate_fidelity = abs(((U2.dag() * U1).tr()))**2
# # sqrt(sqrt(chi_ideal_iSWAP) * chi_real * sqrt(chi_ideal_iSWAP)
#
#
#
# # """
#
# # # v_1 = wave_shape(w_d1, phi1, t) * pulse_shape(amp_1, mean_1, std_1, t)
# # # v_2 = wave_shape(w_d2, phi2, t) * pulse_shape(amp_2, mean_2, std_2, t)
#
# # v_1 = amp_1 * wave_shape(args = args_1) * pulse_shape(args = args_1)
# # v_2 = amp_2 * wave_shape(args = args_2) * pulse_shape(args =args_2)
#
#
# # #Define constants in the driving expressions (page 37...)
# # del_12 = w_q1 - w_q2
# # print(del_12)
# # mu1_plus = g * a_q1 / ((del_12)* (a_q1 - del_12))
# # mu1_minus = -g * a_q1 / ((del_12)* (a_q1 + del_12))
# # v1_plus = -g * del_12/ ((del_12)* (a_q1 - del_12))
# # v1_minus = -g * del_12/ ((del_12)* (a_q1 + del_12))
#
# # mu2_plus = g * a_q2 / ((del_12)* (a_q2 - del_12))
# # mu2_minus = -g * a_q2 / ((del_12)* (a_q2 + del_12))
# # v2_plus = -g * del_12/ ((del_12)* (a_q2 - del_12))
# # v2_minus = -g * del_12/ ((del_12)* (a_q2 + del_12))
#
#
# # #   OPTION B ---- Drive term descriptions for 2 qubits from page 37 ...
# # H_d1 = [amp_1 * r1*(sigma_x1 + v1_minus * sigma_x2 + mu1_minus*(tensor(sigmaz(), sigmax()))), v_t1]
# # H_d2 = [amp_2 * r2*(sigma_x2 + v2_plus*sigma_x1 + mu2_plus*(tensor(sigmax(),sigmaz()))), v_t2]
#
#
#
# # def wave_shape1(t, **kwargs):
# #     x = np.exp(1j*(args['w_d1']*t + args['phi1']))
# #     return x
#
# # def wave_shape2(t, args):
# #     x = np.exp(1j*(args['w_d2']*t + args['phi2']))
# #     return x
#
#
# # def pulse_shape1(t, args):
# #     x = 1.0/ (np.sqrt(2.0 * np.pi) * args['std_1']) * np.exp(-np.power((t - args['mean_1']) / args['std_2'], 2.0) / 2)
# #     return x
#
# # def pulse_shape2(t, args):
# #     x = 1.0/ (np.sqrt(2.0 * np.pi) * args['std_2']) * np.exp(-np.power((t - args['mean_2']) / args['std_2'], 2.0) / 2)
# #     return x
#
# # def v_t1(t, args):
# #     x = wave_shape1(t, args)
# #     y = pulse_shape1(t,args)
# #     v = x*y
# #     return v
#
# # def v_t2(t, args):
# #     x = wave_shape2(t, args)
# #     y = pulse_shape2(t,args)
# #     v = x*y
# #     return v
#
#
#
# # #Define args
# # args = {
# #     'w_d1' : w_d1,
# #     'std_1' : std_1,
# #     'mean_1': mean_1,
# #     'phi1' : phi1,
# #     'w_d2' : w_d2,
# #     'std_2' : std_2,
# #     'mean_2': mean_2,
# #     'phi2': phi2
# # }
#
#
#
#
# # #phi
# # phi1            = np.pi/2              #Phase angle for I1 and Q1
# # phi2            = np.pi/2               #Phase angle for I2 and Q2
#
#
#
# # # H_d1_0a = [- 0.5* r1*(a1 + a1.dag()), Q1]
# # H_d1_0a = - 0.5* r1*(a1 + a1.dag())*Q1
# # H_d1_0b = [-0.5*r1 * 1j*(a1-a1.dag()), I1]
# # # H_d2_0a = [- 0.5*r2*(a2 + a2.dag()) , Q2]
# # H_d2_0a = - 0.5*r2*(a2 + a2.dag())*Q2
# # # H_d2_0b = [1j*(a2-a2.dag()), I2]
# # H_d2_0b = 1j*(a2-a2.dag())* I2
#
#
#
# # #   OPTION A ---- SELF-INFERRED SUPERPOSITION OF SINGLE-QUBIT DRIVES
# # # H_0 = Qobj(0.5 * w_q1 * sigma_z1 + 0.5*w_q2*sigma_z2 + 0.5* g*(sigma_y1*sigma_y2))
# # # H_d1 = [r1*sigma_y1 *amp_1,v_t1]
# # # H_d2 = [r2*sigma_y2 *amp_2,v_t2]
#
#
#
#
#
# # #Define qubit 1 parameters
# # w_q1: float     = 5 * 2 * np.pi
# # a_q1: float     = 0.3 * 2 * np.pi           #anharmonicity for qubit 1
# # r1:float        = 0.01 * 2 * np.pi           #drive coupling for qubit 1
#
# # #Define qubit 2 parameters
# # w_q2: float     = 5 * 2 * np.pi
# # a_q2: float     = 0.3 * 2 * np.pi            #anharmonicity for qubit 2
# # r2:float        = 0.01 * 2 * np.pi           #drive coupling for qubit 2
#
# # w_d1:float      = 5 * 2 * np.pi                    #drive frequency for qubit 1
# # w_d2:float      = 5 * 2 * np.pi                    #drive frequency for qubit 2
#
#
# # #Define pulse parameters - q1
# # amp_1:float     = np.pi *100/ 2                   #pulse amplitude for qubit 1
# # mean_1:float    = T/2                     #mean of gaussian pulse for qubit 1
# # std_1:float     = T/6                   #std of gaussian pulse for qubit 1
#
# # #Define pulse parameters - q2
# # amp_2:float     = 0
# # mean_2:float    = 0
# # std_2:float     = 1.0
#
#
# # #%% Define driving pulse(s)
# # def pulse_shape(t, args) :
# #     amp = args['amp']
# #     std = args['std']
# #     mean = args['mean']
#
# #     return (amp/ (np.sqrt(2.0 * np.pi* std**2))) * np.exp(-np.power((t - mean) /std, 2.0) / 2)
#
#
#
# # args = {'amp': amp_1, 'std': std_1, 'mean': mean_1}
#
#
#
# # plt.plot(times, I1(times, args))
# # plt.title("Applied pulse of %s ns" % T)
# # plt.show()
#
# # def __init__(self, amp, center, std):
# #         self.amp = amp
# #         self.center = center
#
# # I1 = Pulse(
# #     amp = np.pi *100/ 2,
# #     center= T/2 ,
# #     std = T/6
# # )
#
# # Q1 = 0
# # I2 = 0
# # Q2 = 0
#
#
# #  #Drive Components
# #     H_d1_0a = - 0.5* q1.r*(a1 + a1.dag())*d1.Q
# #     H_d1_0b = [-0.5*q1.r * 1j*(a1-a1.dag()), d1.I]
# #     H_d2_0a = - 0.5*q2.r*(a2 + a2.dag())*d2.Q
# #     H_d2_0b = 1j*(a2-a2.dag())* d2.I
#
#
# #     i: int
# #     w_q: float
# #     a_q: float
# #     r: float
# #     w_d: float
#
# #     amp: float
# #     center: float
# #     std: float
#
#
#
# # psi_0 = [
# #             basis(2, 0),
# #             basis(2, 1),
# #             (basis(2, 0) + basis(2, 1)).unit(),
# #             (basis(2, 0) - basis(2, 1)).unit(),
# #             (basis(2, 0) + 1j * basis(2, 1)).unit(),
# #             (basis(2, 0) - 1j * basis(2, 1)).unit()
# #         ]
# # b = Bloch()
# # for psi0 in psi_0:
# #     res = U_psi_real_T[] * psi0
# #     expectations = np.stack(res.expect, axis=-1)
# #     b.add_points(expectations.T) ####  color is function of initial state? Plots all times together.
# # b.show()
#
#
# # #%% Find first peak amplitude - autocorrelation
#
# # #### ft of autocorrelation = PSD
#
# # # Nearest size with power of 2
# # data = outputs[:, -1, 1]
# # size = 2 ** np.ceil(np.log2(2*len(data) - 1)).astype('int')
#
# # # Variance
# # var = np.var(data)
#
# # # Normalized data
# # ndata = data - np.mean(data)
#
# # # Compute the FFT
# # fft = np.fft.fft(ndata, size)
#
# # # Get the power spectrum
# # pwr = np.abs(fft) ** 2
#
# # # Calculate the autocorrelation from inverse FFT of the power spectrum
# # acorr = np.fft.ifft(pwr).real / var / len(data)
#
# # #Find first peak of acorr
# # peaks = argrelextrema(acorr, np.greater)
# # print(peaks)
#
# # Q =  GaussianPulse(
# #     amp = np.pi *100/ 2,
# #     center= T/2,
# #     std = T/6
# # )
#
#
# # if type(d1.I) == GaussianPulse:
# #         H_d1_0b = [-0.5*q1.r * 1j*(a1-a1.dag()), d1.I]
# #     else:
# #         H_d1_0b = -0.5*q1.r * 1j*(a1-a1.dag())* d1.I
#
# #     if type(d1.Q) == GaussianPulse:
# #         H_d1_0a = [- 0.5* q1.r*(a1 + a1.dag()),d1.Q]
# #     else:
# #         H_d1_0a = - 0.5* q1.r*(a1 + a1.dag())*d1.Q
#
# #     if type(d2.I) == GaussianPulse:
# #         H_d2_0b = [1j*(a2-a2.dag()), d2.I]
# #     else:
# #         H_d2_0b = 1j*(a2-a2.dag())* d2.I
#
# #     if type(d2.Q) == GaussianPulse:
# #         H_d2_0a = [- 0.5*q2.r*(a2 + a2.dag()),d2.Q]
# #     else:
# #         H_d2_0a = - 0.5*q2.r*(a2 + a2.dag())*d2.Q
#
#
#
#
# # class DRAGPulse:
# #     def __init__(self, s_t, q ,l):
# #         self.s_t:Pulse = s_t
# #         self.q: Qubit = q
# #         self.l:float = l
#
# #     def Q(self):
# #         return s_t
#
# #     def I(self):
# #         return [self.l/self.q.a_q,  (self.s_t.deriv)]
#
#
#
#
# #     def deriv(self,t):
# #         return (-2 * t) * self.amp * np.exp(-(t - self.center) ** 2 / (2 * self.std ** 2)) / (
# #                 self.std * np.sqrt(2 * np.pi))
#
#
#
# #Imports
# import numpy as np
# from qutip import *
# import matplotlib.pyplot as plt
# import itertools
# from dataclasses import dataclass
# import math
# from scipy.signal import argrelextrema
#
# #Define number of energy levels
# dim:int         = 2
#
# #Define time-stamps
# T:float         = 80                                    #Indicator of pulse time
# t_meas:float    = 0                                     #Time to take measurement - fixed.
# total_T:float   = T + t_meas
# nT:int          = 100                                   #Number of time steps to propagate over
# times           = np.linspace(0, total_T, nT)           #Time vector
#
#
# #Define inter-qubit coupling strength
# g: float        = 0.005 * 2 * np.pi
#
# # %% Define a Qubit class
# class Qubit:
#     def __init__(self, i, w_q, a_q, r, w_d):
#         self.i = i              #index
#         self.w_q = w_q          #qubit freq
#         self.a_q = a_q          #anharmonicity
#         self.r = r              #drive coupling
#         self.w_d = w_d          #drive freq
#
# # %% Define Pulse classes
# class Pulse:
#     def __init__(self, amp, center, std):
#         self.amp = amp
#         self.center = center
#         self.std = std
#
#     def __call__(self, t, *args, **kwargs):
#         return self.amp * np.exp(-(t - self.center) ** 2 / (2 * self.std ** 2)) / (
#                 self.std * np.sqrt(2 * np.pi))
#
#
# class DRAG:
#     def __init__(self, amp, center, std, l, q):
#         self.amp = amp
#         self.center = center
#         self.std = std
#         self. l = l
#         self.q = q
#
#
#     def __call__(self, t, *args, **kwargs):
#         return (-2 * t * self.l/self.q.a_q) * self.amp * np.exp(-(t - self.center) ** 2 / (2 * self.std ** 2)) / (
#             self.std * np.sqrt(2 * np.pi))
#
#
# # %% Define a Drive class
# class Drive:
#     def __init__(self, I, Q):
#         self.I = I          #in-phase drive
#         self.Q = Q          #quadrature drive
#
# #%% Define the qubits
# qubit_1 = Qubit(
#     i = 1,
#     w_q = 5 * 2 * np.pi ,
#     a_q = -0.3 * 2 * np.pi,
#     r = 0.01 * 2 * np.pi,
#     w_d = 5 * 2 * np.pi
# )
#
# qubit_2 = Qubit(
#     i = 2,
#     w_q = 10 * 2 * np.pi,
#     a_q = -0.3 * 2 * np.pi,
#     r = 0.01 * 2 * np.pi,
#     w_d = 6 * 2 * np.pi
# )
#
# #%% Define the pulses at qubits 1 and 2
#
#
# drive_1 =  Drive(
#     # I = DRAG(
#     # amp = np.pi *100/ 2,
#     # center= T/2,
#     # std = T/6,
#     # l = 0.5,
#     # q = qubit_1
#     # ),
#
#     # I = Pulse(
#     # amp = np.pi *100/ 2,
#     # center= T/2,
#     # std = T/6,
#     # ),
#
#     I = 0,
#
#     Q = Pulse(
#     amp = np.pi *100/ 2,
#     center= T/2,
#     std = T/6,
# )
# )
#
# drive_2 =  Drive(
#     I = 0,
#     Q = 0
# )
#
# #%% Define Hamiltonian of 2-qubit system - Obtained from Barnaby's Notes
#
# def create_H(qubits, drives):
#     q1, q2 = qubits
#     d1, d2 = drives
#
#     Delta_1 = q1.w_q - q1.w_d
#     Delta_2 = q2.w_q - q2.w_d
#
#     delta = q1.w_q - q2.w_q
#
#     #Autonomous
#     H_0 = Delta_1 * n1 + Delta_2 * n2
#
#     if type(d1.I) == Pulse or DRAG:
#         H_d1_0b = [-0.5*q1.r * 1j*(a1-a1.dag()), d1.I]
#     else:
#         H_d1_0b = -0.5*q1.r * 1j*(a1-a1.dag())* d1.I
#
#     if type(d1.Q)== Pulse or DRAG:
#         H_d1_0a = [- 0.5* q1.r*(a1 + a1.dag()),d1.Q]
#     else:
#         H_d1_0a = - 0.5* q1.r*(a1 + a1.dag())*d1.Q
#
#     if type(d2.I) == Pulse or DRAG:
#         H_d2_0b = [1j*(a2-a2.dag()), d2.I]
#     else:
#         H_d2_0b = 1j*(a2-a2.dag())* d2.I
#
#     if type(d2.Q) == Pulse or DRAG:
#         H_d2_0a = [- 0.5*q2.r*(a2 + a2.dag()),d2.Q]
#     else:
#         H_d2_0a = - 0.5*q2.r*(a2 + a2.dag())*d2.Q
#
#     H_d1_1 = [g*a1*a2.dag(), lambda t, *args: np.exp(-1j*delta*t)]
#     H_d2_1 = [g*a1.dag()*a2, lambda t, *args: np.exp(1j*delta*t)]
#
#     #Total H
#     H = [H_0, H_d1_0a, H_d1_0b, H_d2_0a, H_d2_0b, H_d1_1, H_d2_1]
#     return H
#
# # %% Define a1, a2
# a1 = tensor(destroy(dim), qeye(dim))
# a2 = tensor(qeye(dim), destroy(dim))
#
# #Define n1 and n2
# n1 = a1.dag() * a1
# n2 = a2.dag() * a2
#
# #Define op_basis
# op_basis = [[qeye(2), sigmax(), sigmay(), sigmaz()]] * 2
#
# #Define plot labels
# op_label = [["i", "x", "y", "z"]] * 2
#
# #Define sigma matrices for the qubits
# sigma_x1 = tensor(sigmax(), qeye(2))
# sigma_y1 = tensor(sigmay(), qeye(2))
# sigma_z1 = tensor(sigmaz(), qeye(2))
#
# sigma_x2 = tensor(qeye(2), sigmax())
# sigma_y2 = tensor(qeye(2), sigmay())
# sigma_z2 = tensor(qeye(2), sigmaz())
#
# II       = tensor(qeye(2), qeye(2))
# XX       = tensor(sigmax(), sigmax())
# YY       = tensor(sigmay(), sigmay())
# ZZ       = tensor(sigmaz(), sigmaz())
#
# #%%Power Rabi
# amps = np.linspace(0, 500, 100)
# outputs = np.zeros([len(amps), len(times), 8])
# for i,amp in enumerate(amps):
#     #change I1 amp to amp
#     drive_1.Q.amp = amp
#
#     #Define starting states
#     psi0 = tensor(basis(2, 0), basis(2, 0))
#
#     #create H
#     H = create_H([qubit_1, qubit_2], [drive_1, drive_2])
#     print(H)
#
#     #input to mesolve
#     res = mesolve(H, psi0, times, [], [sigma_x1, sigma_y1, sigma_z1, sigma_x2, sigma_y2, sigma_z2, n1, n2])
#     outputs[i, ...] = np.stack(res.expect, axis=-1)
# print(shape(outputs))
#
# fig, ax = plt.subplots(1, 1)
# ax.plot(amps, outputs[:, -1, 3 * (qubit_1.i - 1)], label='X', linestyle='--')
# ax.plot(amps, outputs[:, -1, 3 * (qubit_1.i - 1) + 1], label='Y', linestyle='--')
# ax.plot(amps, outputs[:, -1, 3 * (qubit_1.i - 1) + 2], label='Z')
#
# ax.set_ylim([-1, 1])
# ax.set_xlabel('Amplitude')
# ax.set_ylabel('Expectation value')
# ax.set_title('Power Rabi Oscillations')
# ax.legend()
# plt.show()
#
#
# # %% QPT over unknown quantum process  ###########################
# drive_1.Q.amp = 50
# H = create_H([qubit_1, qubit_2], [drive_1, drive_2])
# U_psi_real = qutip.propagator(H, times)                   #List of matrices due to time dependence.
# U_psi_real_T = U_psi_real[nT-1]                           #Take entry from last time step
#
# # TODO: Plot what's going on over the bloch sphere.
#
#
# U_rho_real = spre(U_psi_real_T) * spost(U_psi_real_T.dag())
# chi_real = qpt(U_rho_real, op_basis)
#
# fig = qpt_plot_combined(chi_real, op_label, r'$Actual Process$')
# # plt.show()
#
# #%% iSWAP ###########################
#
# #  Define the iSWAP gate
# U_psi_iSWAP = II * (0.5+0j) + XX * 0.5j + YY * 0.5j + ZZ * (0.5+0j)             #for state vector
# U_rho_iSWAP = spre(U_psi_iSWAP) * spost(U_psi_iSWAP.dag())                      #for density matrix
#
# #Find ideal process matrix
# chi_ideal_iSWAP = qpt(U_rho_iSWAP, op_basis)
# # fig = qpt_plot_combined(chi_ideal_iSWAP, op_label, r'$i$SWAP')
# # plt.show()
#
#
# #%% Cphase #### #######################
#
# #Define U: TODO: convert to tensor product form
# U_psi_cphase = Qobj([[1, 0, 0, 0],
#                 [0, 1, 0, 0],
#                 [0, 0, 1, 0],
#                 [0, 0, 0, -1]])
#
#
# U_rho_cphase = spre(U_psi_cphase) * spost(U_psi_cphase.dag())
# chi_ideal_cphase = qpt(U_rho_cphase, op_basis)
#
# # fig = qpt_plot_combined(chi_ideal_cphase, op_label, r'$CPHASE')
# # plt.show()
#
#
# #%% X Gate  ###########################
# # U_psi_CNOT = Qobj([[1, 0, 0, 0],
# #        [0, 1, 0, 0],
# #         [0, 0, 0, 1],
# #         [0, 0, 1, 0]])
#
# U_psi_CNOT = tensor(sigmax(),qeye(2))
#
# U_rho_CNOT = spre(U_psi_CNOT) * spost(U_psi_CNOT.dag())
# chi_ideal_CNOT = qpt(U_rho_CNOT, op_basis)
# fig = qpt_plot_combined(chi_ideal_CNOT, op_label, r'$CPHASE')
# plt.show()
#
#
# #%% Evaluate process fidelity  ###########################
#
# def process_fidelity_s(chi_ideal, chi_real):
#     # Calculate the square root of the ideal process matrix
#     sqrt_chi_ideal = chi_ideal.sqrtm()
#
#     # Calculate the matrix square root of sqrt_chi_ideal * chi_real * sqrt_chi_ideal
#     sqrt_term = sqrt_chi_ideal * chi_real * sqrt_chi_ideal
#
#     # Calculate the trace of the square root term
#     trace_sqrt_term = sqrt_term.tr()
#
#     # Calculate the process fidelity using the formula
#     fidelity = (np.abs(trace_sqrt_term) ** 2).real
#
#     return fidelity
#
# chi_ideal = Qobj(chi_ideal_CNOT)
# chi_real = Qobj(chi_real)
# fidelity = process_fidelity_s(chi_ideal, chi_real)
#
# print('X-gate Process Fidelity at T = ', fidelity)
#
#
# #%%
# def CPhase(on:bool, qubits, wq1):
# #     qubit_1,qubit_2 = qubits
# #     if on == True:
# #         #### Make wq1 set to a linear ramp!!
# #         # wq1 = LowPassStepValue(
# #         #     value=wq1,
# #         #     step_time=0,
# #         #     tau=8,
# #         #     default_value=qubit_1.w_q)
# #         qubit_1.set_wq(wq1)
# #     return qubits
#
#
# def Cphase(t, **kwargs):
#     exp = kwargs['exp']
#     on = exp.gate_type == 'Cphase'
#     q1, q2 = exp.qubits
#     t_forward = exp.t_exp / 3
#     tau = 2 * exp.t_exp / 3
#     t_backward = exp.t_exp
#
#     if on == True:
#         if t < t_forward:
#             wq_forward = Linear(
#                 t0=0,
#                 t_stop=exp.t_exp / 3,
#                 y0=q1.w_q_def,
#                 y_stop=q1.w_ex12
#             )
#             q1.set_wq(wq_forward)
#
#         elif t >= t_forward and t < tau:
#             wq_tau = Linear(
#                 t0=exp.t_exp / 3,
#                 t_stop=2 * exp.t_exp / 3,
#                 y0=q1.w_ex12,
#                 y_stop=q1.w_ex12
#             )
#             q1.set_wq(wq_tau)
#         elif t >= tau and t < t_backward:
#             wq_backward = Linear(
#                 t0=2 * exp.t_exp / 3,
#                 t_stop=exp.t_exp,
#                 y0=q1.w_ex12,
#                 y_stop=q1.w_q_def
#             )
#             q1.set_wq(wq_backward)
#         else:
#             q1.set_wq(q1.w_q_def)


def Cphase(self):
    exp = self
    on = exp.gate_type == 'Cphase'
    q1, q2 = exp.qubits
    t_forward = exp.t_exp / 3
    tau = 2 * exp.t_exp / 3
    t_backward = exp.t_exp

    if on == True:
        w1_set = lambda t: Linear(
            t0=0,
            t_stop=exp.t_exp / 3,
            y0=q1.w_q_def,
            y_stop=q1.w_ex12
        ) \
            if t < t_forward else Linear(
            t0=exp.t_exp / 3,
            t_stop=2 * exp.t_exp / 3,
            y0=q1.w_ex12,
            y_stop=q1.w_ex12
        ) if t >= t_forward and t < tau else Linear(
            t0=exp.t_exp / 3,
            t_stop=2 * exp.t_exp / 3,
            y0=q1.w_ex12,
            y_stop=q1.w_ex12
        ) if t >= tau and t < t_backward else q1.w_q_def

        # q1.set_wq(w1_set(t))
        # qubits = [q1, q2]
        return w1_set

    exp = self
    on = exp.gate_type == 'Cphase'
    q1, q2 = exp.qubits
    t_forward = exp.t_exp / 3
    tau = 2 * exp.t_exp / 3
    t_backward = exp.t_exp

    if on == True:
        linear_params = lambda t: [0, exp.t_exp / 3, q1.w_q_def, q1.w_ex12] if t < t_forward \
            else [exp.t_exp / 3, 2 * exp.t_exp / 3, q1.w_ex12, q1.w_ex12] if t >= t_forward and t < tau \
            else [2 * exp.t_exp / 3, exp.t_exp, q1.w_ex12, q1.w_q_def]

        t0, t_stop, y0, y_stop = linear_params

        w1_set = Linear(
            t0=t0,
            t_stop=t_stop,
            y0=y0,
            y_stop=y_stop
        )
