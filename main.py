from b_truss_bridge.truss_bridge import truss_bridge

ele_p_1 = {"A": 45.6e-3, "J": 2.35e-3, "Iy": 1.10e-3, "Iz": 2.48e-3}
ele_p_2 = {"A": 20.1e-3, "J": 552.07e-6, "Iy": 460.31e-6, "Iz": 293.11e-6}
ele_p_3 = {"A": 10.42e-3, "J": 93.63e-6, "Iy": 60.40e-6, "Iz": 60.40e-6}
ele_p_4 = {"A": 34.94e-3, "J": 8.88e-6, "Iy": 320.38e-6, "Iz": 2.22e-3}
ele_p_5 = {"A": 26.84e-3, "J": 6.75e-6, "Iy": 214.45e-6, "Iz": 957.43e-6}
ele_p_6 = {"A": 2.94e-3, "J": 6.52e-6, "Iy": 4.18e-6, "Iz": 4.18e-6}

ele_properties = [ele_p_1, ele_p_2, ele_p_3, ele_p_4, ele_p_5, ele_p_6]
for i, ele in enumerate(ele_properties):
    ele["E"] = 210e9
    ele["G"] = 81e9
    if i != 3 or i != 4:
        ele["rho"] = 7850
    else:
        ele["rho"] = 7850 * 4
tb01 = truss_bridge(5, 7, 5, ele_properties)
tb01.plot_mode_shape()

