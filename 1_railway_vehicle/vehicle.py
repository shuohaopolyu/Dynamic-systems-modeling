from sympy import symbols, diff
from sympy.physics.mechanics import *
import numpy as np
from lin_struct_dyn import MultiDOF

# define the system parameters
# m_c: mass of the vehicle body
# I_cx, I_cy, I_cz: moments of inertia of the vehicle body
# m_b: mass of the bogie
# I_bx, I_by, I_bz: moments of inertia of the bogie
# m_w: mass of the wheel
# I_wx, I_wz: moments of inertia of the wheel
# k_sx, k_sy, k_sz: stiffness of the second suspension
# k_px, k_py, k_pz: stiffness of the primary suspension
# k_gx, k_gy, k_gz: equivalent stiffness of the ground
# h_1: height from the bogie center to the primary suspension
# h_2: height from the vehicle body center to the second suspension
# h_3: height from the bogie center to the second suspension
# b_1: distance from the wheel center to the primary suspension
# b_2: distance from the bogie center to the second suspension
# l_1: distance from the wheel center to the bogie center
# l_2: distance from the bogie center to the vehicle body center
# a: distance from the wheel center to the ground
# h_w: height from the wheel center to the ground
m_c, I_cx, I_cy, I_cz = symbols("m_c I_cx I_cy I_cz")
m_b, I_bx, I_by, I_bz = symbols("m_b I_bx I_by I_bz")
m_w, I_wx, I_wz = symbols("m_w I_wx I_wz")
k_sx, k_sy, k_sz = symbols("k_sx k_sy k_sz")
k_px, k_py, k_pz = symbols("k_px k_py k_pz")
k_gx, k_gy, k_gz = symbols("k_gx k_gy k_gz")
h_1, h_2, h_3, b_1, b_2, l_1, l_2 = symbols("h_1 h_2 h_3 b_1 b_2 l_1 l_2")
a, h_w = symbols("a h_w")

# set parameters for the system
params = {
    m_c: 2.9e4,
    I_cx: 1.0401e5,
    I_cy: 1.42887e6,
    I_cz: 1.33816e6,
    m_b: 2.43e3,
    I_bx: 1.97e3,
    I_by: 1.32e3,
    I_bz: 2.42e3,
    m_w: 1.744e3,
    I_wx: 0.855e3,
    I_wz: 0.855e3,
    k_sx: 1.12e5,
    k_sy: 1.12e5,
    k_sz: 2.78e5,
    k_px: 9.49e6,
    k_py: 3.99e6,
    k_pz: 6.26e5,
    k_gx: 6e5,
    k_gy: 6e5,
    k_gz: 1e6,
    h_1: 0.2,
    h_2: 0.75,
    h_3: 0.42,
    b_1: 1,
    b_2: 1.23,
    l_1: 1.25,
    l_2: 7.85,
    a: 0.7175,
    h_w: 0.25,
}

# define the generalized coordinates
y_c, z_c, phi_c, theta_c, psi_c = dynamicsymbols("y_c z_c phi_c theta_c psi_c")
y_b1, z_b1, phi_b1, theta_b1, psi_b1 = dynamicsymbols(
    "y_b1 z_b1 phi_b1 theta_b1 psi_b1"
)
y_b2, z_b2, phi_b2, theta_b2, psi_b2 = dynamicsymbols(
    "y_b2 z_b2 phi_b2 theta_b2 psi_b2"
)
y_w1, z_w1, phi_w1, psi_w1 = dynamicsymbols("y_w1 z_w1 phi_w1 psi_w1")
y_w2, z_w2, phi_w2, psi_w2 = dynamicsymbols("y_w2 z_w2 phi_w2 psi_w2")
y_w3, z_w3, phi_w3, psi_w3 = dynamicsymbols("y_w3 z_w3 phi_w3 psi_w3")
y_w4, z_w4, phi_w4, psi_w4 = dynamicsymbols("y_w4 z_w4 phi_w4 psi_w4")

# define the generalized speeds
y_c_d, z_c_d, phi_c_d, theta_c_d, psi_c_d = dynamicsymbols(
    "y_c z_c phi_c theta_c psi_c", 1
)
y_b1_d, z_b1_d, phi_b1_d, theta_b1_d, psi_b1_d = dynamicsymbols(
    "y_b1 z_b1 phi_b1 theta_b1 psi_b1", 1
)
y_b2_d, z_b2_d, phi_b2_d, theta_b2_d, psi_b2_d = dynamicsymbols(
    "y_b2 z_b2 phi_b2 theta_b2 psi_b2", 1
)
y_w1_d, z_w1_d, phi_w1_d, psi_w1_d = dynamicsymbols("y_w1 z_w1 phi_w1 psi_w1", 1)
y_w2_d, z_w2_d, phi_w2_d, psi_w2_d = dynamicsymbols("y_w2 z_w2 phi_w2 psi_w2", 1)
y_w3_d, z_w3_d, phi_w3_d, psi_w3_d = dynamicsymbols("y_w3 z_w3 phi_w3 psi_w3", 1)
y_w4_d, z_w4_d, phi_w4_d, psi_w4_d = dynamicsymbols("y_w4 z_w4 phi_w4 psi_w4", 1)

# define the generalized accelerations
y_c_dd, z_c_dd, phi_c_dd, theta_c_dd, psi_c_dd = dynamicsymbols(
    "y_c z_c phi_c theta_c psi_c", 2
)
y_b1_dd, z_b1_dd, phi_b1_dd, theta_b1_dd, psi_b1_dd = dynamicsymbols(
    "y_b1 z_b1 phi_b1 theta_b1 psi_b1", 2
)
y_b2_dd, z_b2_dd, phi_b2_dd, theta_b2_dd, psi_b2_dd = dynamicsymbols(
    "y_b2 z_b2 phi_b2 theta_b2 psi_b2", 2
)
y_w1_dd, z_w1_dd, phi_w1_dd, psi_w1_dd = dynamicsymbols("y_w1 z_w1 phi_w1 psi_w1", 2)
y_w2_dd, z_w2_dd, phi_w2_dd, psi_w2_dd = dynamicsymbols("y_w2 z_w2 phi_w2 psi_w2", 2)
y_w3_dd, z_w3_dd, phi_w3_dd, psi_w3_dd = dynamicsymbols("y_w3 z_w3 phi_w3 psi_w3", 2)
y_w4_dd, z_w4_dd, phi_w4_dd, psi_w4_dd = dynamicsymbols("y_w4 z_w4 phi_w4 psi_w4", 2)

# define the kinetic energy
T_c = (
    0.5 * m_c * (y_c_d**2 + z_c_d**2)
    + 0.5 * I_cx * phi_c_d**2
    + 0.5 * I_cy * theta_c_d**2
    + 0.5 * I_cz * psi_c_d**2
)
T_b1 = (
    0.5 * m_b * (y_b1_d**2 + z_b1_d**2)
    + 0.5 * I_bx * phi_b1_d**2
    + 0.5 * I_by * theta_b1_d**2
    + 0.5 * I_bz * psi_b1_d**2
)
T_b2 = (
    0.5 * m_b * (y_b2_d**2 + z_b2_d**2)
    + 0.5 * I_bx * phi_b2_d**2
    + 0.5 * I_by * theta_b2_d**2
    + 0.5 * I_bz * psi_b2_d**2
)
T_w1 = (
    0.5 * m_w * (y_w1_d**2 + z_w1_d**2)
    + 0.5 * I_wx * phi_w1_d**2
    + 0.5 * I_wz * psi_w1_d**2
)
T_w2 = (
    0.5 * m_w * (y_w2_d**2 + z_w2_d**2)
    + 0.5 * I_wx * phi_w2_d**2
    + 0.5 * I_wz * psi_w2_d**2
)
T_w3 = (
    0.5 * m_w * (y_w3_d**2 + z_w3_d**2)
    + 0.5 * I_wx * phi_w3_d**2
    + 0.5 * I_wz * psi_w3_d**2
)
T_w4 = (
    0.5 * m_w * (y_w4_d**2 + z_w4_d**2)
    + 0.5 * I_wx * phi_w4_d**2
    + 0.5 * I_wz * psi_w4_d**2
)

T = T_c + T_b1 + T_b2 + T_w1 + T_w2 + T_w3 + T_w4

# define the potential energy
# V_s: potential energy of the second suspension
V_s = (
    0.5
    * k_sx
    * (
        (theta_c * h_2 + theta_b1 * h_3 - psi_c * b_2 + psi_b1 * b_2) ** 2
        + (theta_c * h_2 + theta_b1 * h_3 + psi_c * b_2 - psi_b1 * b_2) ** 2
        + (theta_c * h_2 + theta_b2 * h_3 - psi_c * b_2 + psi_b2 * b_2) ** 2
        + (theta_c * h_2 + theta_b2 * h_3 + psi_c * b_2 - psi_b2 * b_2) ** 2
    )
    + 0.5
    * k_sy
    * (
        (y_c - phi_c * h_2 + psi_c * l_2 - y_b1 - phi_b1 * h_3) ** 2
        + (y_c - phi_c * h_2 + psi_c * l_2 - y_b1 - phi_b1 * h_3) ** 2
        + (y_c - phi_c * h_2 - psi_c * l_2 - y_b2 - phi_b2 * h_3) ** 2
        + (y_c - phi_c * h_2 - psi_c * l_2 - y_b2 - phi_b2 * h_3) ** 2
    )
    + 0.5
    * k_sz
    * (
        (z_c - phi_c * b_2 + theta_c * l_2 - z_b1 + phi_b1 * b_2) ** 2
        + (z_c + phi_c * b_2 + theta_c * l_2 - z_b1 - phi_b1 * b_2) ** 2
        + (z_c - phi_c * b_2 - theta_c * l_2 - z_b2 + phi_b2 * b_2) ** 2
        + (z_c + phi_c * b_2 - theta_c * l_2 - z_b2 - phi_b2 * b_2) ** 2
    )
)

# V_p: potential energy of the primary suspension
V_p = (
    0.5
    * k_px
    * (
        (psi_b1 * b_1 - theta_b1 * h_1 - psi_w1 * b_1) ** 2
        + (psi_b1 * b_1 + theta_b1 * h_1 - psi_w1 * b_1) ** 2
        + (psi_b1 * b_1 - theta_b1 * h_1 - psi_w2 * b_1) ** 2
        + (psi_b1 * b_1 + theta_b1 * h_1 - psi_w2 * b_1) ** 2
        + (psi_b2 * b_1 - theta_b2 * h_1 - psi_w3 * b_1) ** 2
        + (psi_b2 * b_1 + theta_b2 * h_1 - psi_w3 * b_1) ** 2
        + (psi_b2 * b_1 - theta_b2 * h_1 - psi_w4 * b_1) ** 2
        + (psi_b2 * b_1 + theta_b2 * h_1 - psi_w4 * b_1) ** 2
    )
    + 0.5
    * k_py
    * (
        (y_b1 - phi_b1 * h_1 + psi_b1 * l_1 - y_w1) ** 2
        + (y_b1 - phi_b1 * h_1 + psi_b1 * l_1 - y_w1) ** 2
        + (y_b1 - phi_b1 * h_1 - psi_b1 * l_1 - y_w2) ** 2
        + (y_b1 - phi_b1 * h_1 - psi_b1 * l_1 - y_w2) ** 2
        + (y_b2 - phi_b2 * h_1 + psi_b2 * l_1 - y_w3) ** 2
        + (y_b2 - phi_b2 * h_1 + psi_b2 * l_1 - y_w3) ** 2
        + (y_b2 - phi_b2 * h_1 - psi_b2 * l_1 - y_w4) ** 2
        + (y_b2 - phi_b2 * h_1 - psi_b2 * l_1 - y_w4) ** 2
    )
    + 0.5
    * k_pz
    * (
        (z_b1 + phi_b1 * b_1 - theta_b1 * l_1 - z_w1 - phi_w1 * b_1) ** 2
        + (z_b1 - phi_b1 * b_1 - theta_b1 * l_1 - z_w1 + phi_w1 * b_1) ** 2
        + (z_b1 + phi_b1 * b_1 + theta_b1 * l_1 - z_w2 - phi_w2 * b_1) ** 2
        + (z_b1 - phi_b1 * b_1 + theta_b1 * l_1 - z_w2 + phi_w2 * b_1) ** 2
        + (z_b2 + phi_b2 * b_1 - theta_b2 * l_1 - z_w3 - phi_w3 * b_1) ** 2
        + (z_b2 - phi_b2 * b_1 - theta_b2 * l_1 - z_w3 + phi_w3 * b_1) ** 2
        + (z_b2 + phi_b2 * b_1 + theta_b2 * l_1 - z_w4 - phi_w4 * b_1) ** 2
        + (z_b2 - phi_b2 * b_1 + theta_b2 * l_1 - z_w4 + phi_w4 * b_1) ** 2
    )
)

# V_g: potential energy of the ground
V_g = (
    2
    * 0.5
    * k_gx
    * ((psi_w1 * a) ** 2 + (psi_w2 * a) ** 2 + (psi_w3 * a) ** 2 + (psi_w4 * a) ** 2)
    + 0.5
    * k_gy
    * (
        (y_w1 + phi_w1 * h_w) ** 2
        + (y_w1 - phi_w1 * h_w) ** 2
        + (y_w2 + phi_w2 * h_w) ** 2
        + (y_w2 - phi_w2 * h_w) ** 2
        + (y_w3 + phi_w3 * h_w) ** 2
        + (y_w3 - phi_w3 * h_w) ** 2
        + (y_w4 + phi_w4 * h_w) ** 2
        + (y_w4 - phi_w4 * h_w) ** 2
    )
    + 0.5
    * k_gz
    * (
        (z_w1 + phi_w1 * a) ** 2
        + (z_w1 - phi_w1 * a) ** 2
        + (z_w2 + phi_w2 * a) ** 2
        + (z_w2 - phi_w2 * a) ** 2
        + (z_w3 + phi_w3 * a) ** 2
        + (z_w3 - phi_w3 * a) ** 2
        + (z_w4 + phi_w4 * a) ** 2
        + (z_w4 - phi_w4 * a) ** 2
    )
)

V = V_s + V_p + V_g

L = T - V
q = [
    y_c,
    z_c,
    phi_c,
    theta_c,
    psi_c,
    y_b1,
    z_b1,
    phi_b1,
    theta_b1,
    psi_b1,
    y_b2,
    z_b2,
    phi_b2,
    theta_b2,
    psi_b2,
    y_w1,
    z_w1,
    phi_w1,
    psi_w1,
    y_w2,
    z_w2,
    phi_w2,
    psi_w2,
    y_w3,
    z_w3,
    phi_w3,
    psi_w3,
    y_w4,
    z_w4,
    phi_w4,
    psi_w4,
]

# compute the mass matrix and the stiffness matrix
LM = LagrangesMethod(L, q)
mechanics_printing(pretty_print=False)
LM.form_lagranges_equations()
m_mtx = LM.mass_matrix
m_mtx = np.matrix(m_mtx.subs(params))
k_mtx = LM.stiffness_matrix
k_mtx = np.matrix(k_mtx.subs(params))

with open("./sys_matrices/m_mtx.txt", "w") as f:
    for line in m_mtx:
        np.savetxt(f, line, fmt="%.3f")

with open("./sys_matrices/k_mtx.txt", "w") as f:
    for line in k_mtx:
        np.savetxt(f, line, fmt="%.3f")

veh = MultiDOF(np.float64(m_mtx), np.float64(k_mtx))
print(veh)

# compute the damping matrix
params["k_sx"] = 3.5e5
params["k_sy"] = 5e4
params["k_sz"] = 4.5e4
params["k_px"] = 0
params["k_py"] = 0
params["k_pz"] = 3e4
params["k_gx"] = 0
params["k_gy"] = 0
params["k_gz"] = 0
c_mtx = LM.stiffness_matrix
c_mtx = np.matrix(c_mtx.subs(params))
with open("./sys_matrices/c_mtx.txt", "w") as f:
    for line in c_mtx:
        np.savetxt(f, line, fmt="%.3f")
