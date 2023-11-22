"""
A simple demonstration of forced vibration of the 3-span continuous bridge.
A periodic force is applied at the mid-span of the left side span.
"""

import numpy as np
from lsds import MultiDOF
from fembeam import beam_fem
import matplotlib.pyplot as plt

fem_model = beam_fem()
m_mtx = fem_model.assemble_glb_mtx(type="mass")
k_mtx = fem_model.assemble_glb_mtx(type="stiff")

# Define the system
sys = MultiDOF(
    mass_mtx=m_mtx,
    stiff_mtx=k_mtx,
    damp_type="Rayleigh",
    damp_params=(0, 4, 0.03),
    f_dof=[11],
    resp_dof=[i for i in range(86)],
    t_eval=np.linspace(0, 1, 5000),
    f_t=[
        lambda t: 2e5
        * (
            np.sin(2 * np.pi * 8 * t)
            + np.sin(2 * np.pi * 38 * t)
            + np.sin(2 * np.pi * 218 * t)
            + 0.5 * np.sin(2 * np.pi * 100 * t)
        )
    ],
    init_cond=np.zeros(172),
)
print(sys)

# Solve the system
full_resp = sys.response(method="Radau", type="full")
f_mtx = sys.f_mtx()

states = np.vstack((full_resp["displacement"], full_resp["velocity"], f_mtx)).T
acc = full_resp["acceleration"].T

# Save the results
# np.savetxt("states.txt", states)
# np.savetxt("acc.txt", acc)

# plot the force history
plt.plot(sys.t_eval, f_mtx.reshape(-1))
plt.show()

# Plot the results
plt.plot(sys.t_eval, full_resp["displacement"][5, :], label="Node 1")
plt.plot(sys.t_eval, full_resp["displacement"][29, :], label="Node 2")
plt.show()
