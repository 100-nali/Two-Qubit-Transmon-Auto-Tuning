
# %% Define Hamiltonian of 2-qubit system.
# H_q1 = (w_q1 * a1.dag()* a1) + (0.5 * a_q1 * a1.dag() * a1.dag() * a1 * a1)
# H_q2 = (w_q2 * a2.dag()* a2) + (0.5 * a_q2 * a2.dag() * a2.dag() * a2 * a2)
# H_q_int = -1 * g * (a1 - a1.dag())*(a2 - a2.dag())
# H_d1 = [ 
#     [0, lambda t, *args: -0.5 * r1 * amp_1 * shape_1 * np.exp(1j*(w_q1-w_d1)*t + phi1) ],
#     [0, lambda t, *args: -0.5 * r1 * amp_1 * shape_1 *np.exp(-1j*((w_q1-w_d1)*t + phi1))]
#     ]
# H_d1 = r1 * amp_1 * shape_1 *  * sigmay() #
# H_d2
# H_d2 = [ 
#     [a1*a2.dag(), lambda t, *args: -0.5 * r2 * np.exp(1j*(w_q2-w_d2)*t + phi2) ],
#     [a1.dag() * a2], lambda t, *args: -0.5 * r2 * np.exp(-1j*((w_q2-w_d2)*t + phi2))]
# H_d1_d2 = tensor(H_d1, H_d2)

# H = [*H_q1, *H_q2 , *H_q_int,  *H_d1, *H_d2]
# H =  [*H_q1, *H_q2 , *H_q_int,  *H_d1_d2]

# print(H)


#%%
# Create identity operators for the specified dimension
# identity = qutip.qeye(dim)
# a1 = qutip.tensor(qutip.destroy(dim), identity)
# a2 = qutip.tensor(identity, qutip.destroy(dim))

# Define the Hamiltonian terms (adjust the operators and dimensions as needed)
# H1 = (a_q1 / 2) * (a1.dag() * a1.dag() * a1 * a1) - w_d1 * (a1.dag() * a1)
# H1d = -r1 * (a1 + a1.dag()) * (math.cos(phi1) - 1j * math.sin(phi1))
# H2 = (a_q2 / 2) * (a2.dag() * a2.dag() * a2 * a2) - w_d2 * (a2.dag() * a2)
# H2d = -r2 * (a2 + a2.dag()) * (math.cos(phi2) - 1j * math.sin(phi2))
# delta_d = w_d1 - w_d2
# H_int = a1 * a2.dag() + a1.dag() * a2

# Combine Hamiltonian terms into a list (if needed)
# H = [H1, H1d, H2, H2d, [H_int, lambda t, *args: g * np.exp(-1j * delta_d * t)]]

# %% 

# U_rho = []
# for U_psi_i in U_psi:
#     U_rho_i = spre(U_psi_i) * spost(U_psi_i.dag())
#     U_rho.append(U_rho_i)


# print(U_rho)
# chi_real = qpt(U_rho, op_basis) 

# %%

# @dataclass
# class GaussianPulse:
#     amplitude: float
#     std: float
#     center: float

#     def __call__(self, t, *args, **kwargs):
#         return self.amplitude * np.exp(-(t - self.center) ** 2 / (2 * self.std ** 2)) / (
#                 self.std * np.sqrt(2 * np.pi))

