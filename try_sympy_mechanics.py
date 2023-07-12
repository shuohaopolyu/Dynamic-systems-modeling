from sympy import symbols, diff
from sympy.physics.mechanics import *

m1, c1, k1 = symbols("m1 c1 k1")
m2, c2, k2 = symbols("m2 c2 k2")

u1, u2 = dynamicsymbols("u1, u2")
u1_d, u2_d = dynamicsymbols("u1, u2", 1)
u1_dd, u2_dd = dynamicsymbols("u1, u2", 2)


# define the kinetic energy
T = 0.5 * m1 * u1_d**2 + 0.5 * m2 * u2_d**2

# define the potential energy
V = 0.5 * k1 * (u1 - u2) ** 2 + 0.5 * k2 * u2**2

# define the lagrangian
L = T - V

# define the lagrange equation
LM = LagrangesMethod(L, [u1, u2])
LE = LM.form_lagranges_equations()
m_mtx = LM.mass_matrix
k_mtx = LM.stiffness_matrix
print(m_mtx)
print(k_mtx)
# k_mtx = LM.forcing
# print(m_mtx.subs({m: 2, u_dd: 1}))
# print(k_mtx.subs({k: 3, u: 1}))
