"""
A simple demonstration of forced vibration of the 3-span continuous bridge.
Four stochastic forces are applied to the bridge.
"""

import numpy as np
import matplotlib.pyplot as plt
from system import ContinuousBridge
from excitation.excitgen import PSDExcitationGenerator
from excitation.psd import WhiteNoisePSD
import pickle


# Define the force
whitenoisepsd = WhiteNoisePSD(a_v=1e8)
whitenoise = PSDExcitationGenerator(psd=whitenoisepsd, tmax=100, fmax=50)
whitenoise_fun_1 = whitenoise()
whitenoise_fun_2 = whitenoise()
whitenoise_fun_3 = whitenoise()
whitenoise_fun_4 = whitenoise()

# Define the system
cb = ContinuousBridge(
    rho=1e4,
    E=2e10,
    A=1,
    I=1 / 12,
    L=0.5,
    damp_params=(0, 3, 0.03),
    f_dof=[13, 30, 44, 69],
    resp_dof="full",
    t_eval=np.arange(0, 100, 1/100),
    f_t=[whitenoise_fun_1, whitenoise_fun_2, whitenoise_fun_3, whitenoise_fun_4],
    init_cond=np.zeros(172),
)

force_mtx = cb.f_mtx()
time = cb.t_eval

# # generate the data
# full_data = cb.response(type="full", method="Radau")
# full_data["displacement"] = full_data["displacement"].T
# full_data["velocity"] = full_data["velocity"].T
# full_data["acceleration"] = full_data["acceleration"].T
# full_data["force"] = force_mtx.T
# full_data["time"] = time.reshape(-1,1)

# with open('continous_bridge_dataset/full_data.pkl', 'wb') as f:
#     pickle.dump(full_data, f)

# read the data
with open('continous_bridge_dataset/full_data.pkl', 'rb') as f:
    full_data = pickle.load(f)

# # plot the data
fig, ax = plt.subplots(3,1, figsize=(12,6))
ax[0].plot(time, full_data["displacement"][:,30])
ax[0].set_ylabel("Displacement [m]")
ax[1].plot(time, full_data["velocity"][:,31])
ax[1].set_ylabel("Velocity [m/s]")
ax[2].plot(time, full_data["acceleration"][:,32])
ax[2].set_ylabel("Acceleration [m/s^2]")
ax[2].set_xlabel("Time [s]")
plt.show()
