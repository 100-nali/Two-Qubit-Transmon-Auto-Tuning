# def define_pulse(self):
#     # self.Gaussian = GaussianPulse(self.params)
#     self.GaussianParams = self.params.get('GaussianParams')
#     self.lamb = self.params.get('lamb', 0.75)
#     self.alpha = self.params.get('alpha', 1)

# %%
# def define_pulse(self):
#     self.amp = self.params.get('amp', 50.0)
#     self.center = self.params.get('center', 0.0)
#     self.std = self.params.get('std', 1.0)


# %% define universal Pulse class

# class Pulse:
#     def __init__(self, **kwargs):
#         self.params = kwargs

# %%
# def define_pulse(self, **kwargs):
#     self.params = kwargs
#     self.N = self.params.get('N', 5)
#     self.amps = self.params.get('amps', [1,1,1,1,1])
#
#

# %%


# delta = q1.w_q - q2.w_q
# %%
# if type(q1.w_q) != float:
#     H_d1_1 = [g * a1 * a2.dag(), lambda t, *args: (np.exp(-1j * (q1.w_q(t)) - q2.w_q) * t)]
#     H_d2_1 = [g * a1.dag() * a2, lambda t, *args: (np.exp(1j * (q1.w_q(t)) - q2.w_q) * t)]
# #%%
# else:
#     H_d1_1 = [g * a1 * a2.dag(), lambda t, *args: np.exp(-1j * (q1.w_q - q2.w_q) * t)]
#     H_d2_1 = [g * a1.dag() * a2, lambda t, *args: np.exp(1j * (q1.w_q - q2.w_q) * t)]


# # the qubit 1 drift hamiltonian
        # H1 = [
        #     (q1.a_q / 2) * (a1.dag() * a1.dag() * a1 * a1) - q1.w_d * (
        #                 a1.dag() * a1),
        #     [a1.dag() * a1, q1.w_q]
        # ]
        #
        # # the qubit 1 drive hamiltonian
        # H1d = [
        #     [-q1.r * 1j * (a1 - a1.dag()) / 2., d1.I],
        #     [-q1.r * (a1 + a1.dag()) / 2., d1.Q]
        # ]
        #
        # # the qubit 2 drift hamiltonian
        # H2 = [
        #     (q1.a_q / 2) * (a2.dag() * a2.dag() * a2 * a2) - q2.w_d * (
        #                 a2.dag() * a2),
        #     [(a2.dag() * a2), q2.w_q]
        # ]
        #
        # H2d = [
        #     [-q2.r * 1j * (a2 - a2.dag()) / 2., d2.I],
        #     [-q2.r * (a2 + a2.dag()) / 2., d2.Q]
        # ]
        #
        # delta_d = q1.w_d - q2.w_d
        #
        # # the interaction hamiltonian
        # H_int = [
        #     [a1 * a2.dag(), lambda t, *args: self.g * np.exp(-1j * delta_d * t)],
        #     [a1.dag() * a2, lambda t, *args: self.g * np.exp(1j * delta_d * t)]
        # ]
        #
        # H= [*H1, *H1d, *H2, *H2d, *H_int]
