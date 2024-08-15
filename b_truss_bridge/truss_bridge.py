import numpy as np
import matplotlib.pyplot as plt
from common_modules import MultiDOF
import pickle

class truss_bridge(MultiDOF):

    def __init__(
        self,
        a,
        b,
        h,
        ele_properties,
        t_eval=np.linspace(0, 6, 12000),
        resp_dof=[41 * 6 + 1],
    ):
        self.a = a
        self.b = b
        self.h = h
        self.ele_properties = ele_properties
        self.mass_mtx = self._assemble_global_matrix(type="mass")
        self.stiff_mtx = self._assemble_global_matrix(type="stiff")
        self.mass_mtx = self.apply_boundary_condition(self.mass_mtx)
        self.stiff_mtx = self.apply_boundary_condition(self.stiff_mtx)
        super().__init__(
            self.mass_mtx,
            self.stiff_mtx,
            f_dof=[31 * 6 + 1],
            t_eval=t_eval,
            f_t=[self.transient_impact_force(factor=1e8)],
            damp_params=(0, 1, 0.02),
            resp_dof=resp_dof,
        )

    def _generate_nodes(self):
        nodes_coord = []
        for i in range(9):
            nodes_coord.append([self.a * i, 0, 0])
        for i in range(9):
            nodes_coord.append([self.a * i, self.b / 3, 0])
        for i in range(9):
            nodes_coord.append([self.a * i, 2 * self.b / 3, 0])
        for i in range(9):
            nodes_coord.append([self.a * i, self.b, 0])
        for i in range(1, 8):
            nodes_coord.append([self.a * i, 0, self.h])
        for i in range(1, 8):
            nodes_coord.append([self.a * i, self.b, self.h])
        return nodes_coord

    def _dof_list(self):
        dof_list = []
        for i in range(50):
            dof_list.append(
                [6 * i, 6 * i + 1, 6 * i + 2, 6 * i + 3, 6 * i + 4, 6 * i + 5]
            )
        return dof_list

    def element_lists(self):
        eles_1_to_6 = [[] for i in range(6)]
        for i in range(6):
            file_name = f"./b_truss_bridge/elements/ele_{i+1}.csv"
            with open(file_name, "r") as f:
                for line in f:
                    eles_1_to_6[i].append([int(i) for i in line.split(sep=",")])
        return eles_1_to_6

    def _element_length(self, node_1_coord, node_2_coord):
        return np.sqrt(
            (node_2_coord[0] - node_1_coord[0]) ** 2
            + (node_2_coord[1] - node_1_coord[1]) ** 2
            + (node_2_coord[2] - node_1_coord[2]) ** 2
        )

    def _beam_element_cos(self, node_1_coord, node_2_coord, node_3_coord):
        L = self._element_length(node_1_coord, node_2_coord)
        X_vec = np.array([node_1_coord[0], node_2_coord[0], node_3_coord[0]])
        Y_vec = np.array([node_1_coord[1], node_2_coord[1], node_3_coord[1]])
        Z_vec = np.array([node_1_coord[2], node_2_coord[2], node_3_coord[2]])
        X_k, X_l = np.meshgrid(X_vec, X_vec)
        Y_k, Y_l = np.meshgrid(Y_vec, Y_vec)
        Z_k, Z_l = np.meshgrid(Z_vec, Z_vec)
        X_kl = X_k - X_l
        Y_kl = Y_k - Y_l
        Z_kl = Z_k - Z_l
        l_x = X_kl[1, 0] / L
        m_x = Y_kl[1, 0] / L
        n_x = Z_kl[1, 0] / L
        const_1 = Y_kl[1, 0] - Z_kl[2, 0] - Y_kl[2, 0] + Z_kl[1, 0]
        const_2 = Z_kl[1, 0] - X_kl[2, 0] - Z_kl[2, 0] + X_kl[1, 0]
        const_3 = X_kl[1, 0] - Y_kl[2, 0] - X_kl[2, 0] + Y_kl[1, 0]
        A_123 = 2 * np.sqrt(const_1**2 + const_2**2 + const_3**2)
        l_z = const_1 / A_123
        m_z = const_2 / A_123
        n_z = const_3 / A_123
        l_y = m_z * n_x - m_x * n_z
        m_y = n_z * l_x - n_x * l_z
        n_y = l_z * m_x - l_x * m_z
        T3 = np.array(
            [
                [l_x, m_x, n_x],
                [l_y, m_y, n_y],
                [l_z, m_z, n_z],
            ]
        )
        return T3

    def _beam_T_matrix(self, node_1_coord, node_2_coord, node_3_coord):
        T3 = self._beam_element_cos(node_1_coord, node_2_coord, node_3_coord)
        T = np.zeros((12, 12))
        T[0:3, 0:3] = T3
        T[3:6, 3:6] = T3
        T[6:9, 6:9] = T3
        T[9:12, 9:12] = T3
        return T

    def bar_ele_matrix(self, type, ele_properties, node_1_coord, node_2_coord):
        if type == "stiff":
            E = ele_properties["E"]
            A = ele_properties["A"]
            K = self._bar_stiffness_matrix(E, A, node_1_coord, node_2_coord)
            return K
        elif type == "mass":
            rho = ele_properties["rho"]
            A = ele_properties["A"]
            M = self._bar_mass_matrix(rho, A, node_1_coord, node_2_coord)
            return M

    def _beam_stiffness_matrix(
        self, G, E, J, Iy, Iz, A, node_1_coord, node_2_coord, node_3_coord
    ):
        L = self._element_length(node_1_coord, node_2_coord)
        Ke = np.array(
            [
                [A * E / L, 0, 0, 0, 0, 0, -A * E / L, 0, 0, 0, 0, 0],
                [
                    0,
                    12 * E * Iz / L**3,
                    0,
                    0,
                    0,
                    6 * E * Iz / L**2,
                    0,
                    -12 * E * Iz / L**3,
                    0,
                    0,
                    0,
                    6 * E * Iz / L**2,
                ],
                [
                    0,
                    0,
                    12 * E * Iy / L**3,
                    0,
                    -6 * E * Iy / L**2,
                    0,
                    0,
                    0,
                    -12 * E * Iy / L**3,
                    0,
                    -6 * E * Iy / L**2,
                    0,
                ],
                [0, 0, 0, G * J / L, 0, 0, 0, 0, 0, -G * J / L, 0, 0],
                [
                    0,
                    0,
                    -6 * E * Iy / L**2,
                    0,
                    4 * E * Iy / L,
                    0,
                    0,
                    0,
                    6 * E * Iy / L**2,
                    0,
                    2 * E * Iy / L,
                    0,
                ],
                [
                    0,
                    6 * E * Iz / L**2,
                    0,
                    0,
                    0,
                    4 * E * Iz / L,
                    0,
                    -6 * E * Iz / L**2,
                    0,
                    0,
                    0,
                    2 * E * Iz / L,
                ],
                [-A * E / L, 0, 0, 0, 0, 0, A * E / L, 0, 0, 0, 0, 0],
                [
                    0,
                    -12 * E * Iz / L**3,
                    0,
                    0,
                    0,
                    -6 * E * Iz / L**2,
                    0,
                    12 * E * Iz / L**3,
                    0,
                    0,
                    0,
                    -6 * E * Iz / L**2,
                ],
                [
                    0,
                    0,
                    -12 * E * Iy / L**3,
                    0,
                    6 * E * Iy / L**2,
                    0,
                    0,
                    0,
                    12 * E * Iy / L**3,
                    0,
                    6 * E * Iy / L**2,
                    0,
                ],
                [0, 0, 0, -G * J / L, 0, 0, 0, 0, 0, G * J / L, 0, 0],
                [
                    0,
                    0,
                    -6 * E * Iy / L**2,
                    0,
                    2 * E * Iy / L,
                    0,
                    0,
                    0,
                    6 * E * Iy / L**2,
                    0,
                    4 * E * Iy / L,
                    0,
                ],
                [
                    0,
                    6 * E * Iz / L**2,
                    0,
                    0,
                    0,
                    2 * E * Iz / L,
                    0,
                    -6 * E * Iz / L**2,
                    0,
                    0,
                    0,
                    4 * E * Iz / L,
                ],
            ]
        )
        T = self._beam_T_matrix(node_1_coord, node_2_coord, node_3_coord)
        K = np.dot(np.dot(T.T, Ke), T)
        return K

    def _beam_mass_matrix(
        self, rho, Iy, Iz, A, L, node_1_coord, node_2_coord, node_3_coord
    ):
        Ix = Iy + Iz
        Me = (
            rho
            * A
            * L
            / 420
            * np.array(
                [
                    [140, 0, 0, 0, 0, 0, 70, 0, 0, 0, 0, 0],
                    [0, 156, 0, 0, 0, 22 * L, 0, 54, 0, 0, 0, -13 * L],
                    [0, 0, 156, 0, -22 * L, 0, 0, 0, 54, 0, 13 * L, 0],
                    [0, 0, 0, 140 * Ix / A, 0, 0, 0, 0, 0, -70 * Ix / A, 0, 0],
                    [0, 0, -22 * L, 0, 4 * L**2, 0, 0, 0, -13 * L, 0, -3 * L**2, 0],
                    [0, 22 * L, 0, 0, 0, 4 * L**2, 0, 13 * L, 0, 0, 0, -3 * L**2],
                    [70, 0, 0, 0, 0, 0, 140, 0, 0, 0, 0, 0],
                    [0, 54, 0, 0, 0, 13 * L, 0, 156, 0, 0, 0, -22 * L],
                    [0, 0, 54, 0, -13 * L, 0, 0, 0, 156, 0, 22 * L, 0],
                    [0, 0, 0, -70 * Ix / A, 0, 0, 0, 0, 0, 140 * Ix / A, 0, 0],
                    [0, 0, 13 * L, 0, -3 * L**2, 0, 0, 0, 22 * L, 0, 4 * L**2, 0],
                    [0, -13 * L, 0, 0, 0, -3 * L**2, 0, 22 * L, 0, 0, 0, 4 * L**2],
                ]
            )
        )
        T = self._beam_T_matrix(node_1_coord, node_2_coord, node_3_coord)
        M = np.dot(np.dot(T.T, Me), T)
        return M

    def beam_ele_matrix(self, type, ele_properties, node_1_coord, node_2_coord):
        if node_1_coord[1] - node_2_coord[1] == 0:
            node_3_coord = [0, node_1_coord[1], 1]
        elif node_1_coord[0] - node_2_coord[0] == 0:
            node_3_coord = [node_1_coord[0], 0, 1]
        else:
            node_3_coord = [node_2_coord[0], node_2_coord[1], 1]
        if type == "stiff":
            G = ele_properties["G"]
            E = ele_properties["E"]
            J = ele_properties["J"]
            Iy = ele_properties["Iy"]
            Iz = ele_properties["Iz"]
            A = ele_properties["A"]
            K = self._beam_stiffness_matrix(
                G, E, J, Iy, Iz, A, node_1_coord, node_2_coord, node_3_coord
            )
            return K
        elif type == "mass":
            rho = ele_properties["rho"]
            Iy = ele_properties["Iy"]
            Iz = ele_properties["Iz"]
            A = ele_properties["A"]
            L = self._element_length(node_1_coord, node_2_coord)
            M = self._beam_mass_matrix(
                rho, Iy, Iz, A, L, node_1_coord, node_2_coord, node_3_coord
            )
            return M

    def _assemble_matrix(
        self, type="stiff", ele_type="bar", ele_list=None, ele_properties=None
    ):
        nodes_coord = self._generate_nodes()
        dof_list = self._dof_list()
        max_dof_num = dof_list[-1][-1] + 1
        ele_matrix = np.zeros((max_dof_num, max_dof_num))
        for element in ele_list:
            node_1_num = element[0]
            node_2_num = element[1]
            node_1_coord = nodes_coord[node_1_num]
            node_2_coord = nodes_coord[node_2_num]
            node_1_dof = dof_list[node_1_num]
            node_2_dof = dof_list[node_2_num]

            if ele_type == "bar":
                ele_dof = np.array([node_1_dof[0:3], node_2_dof[0:3]])
                ele_dof_x, ele_dof_y = np.meshgrid(ele_dof, ele_dof)
                ele_matrix[ele_dof_x, ele_dof_y] += self.bar_ele_matrix(
                    type, ele_properties, node_1_coord, node_2_coord
                )
            elif ele_type == "beam":
                ele_dof = np.array([node_1_dof, node_2_dof])
                ele_dof_x, ele_dof_y = np.meshgrid(ele_dof, ele_dof)
                ele_matrix[ele_dof_x, ele_dof_y] += self.beam_ele_matrix(
                    type, ele_properties, node_1_coord, node_2_coord
                )
        return ele_matrix

    def _assemble_global_matrix(self, type="stiff"):
        ele_1, ele_2, ele_3, ele_4, ele_5, ele_6 = self.element_lists()
        global_matrix = self._assemble_matrix(
            type=type,
            ele_type="beam",
            ele_list=ele_1,
            ele_properties=self.ele_properties[0],
        )
        global_matrix += self._assemble_matrix(
            type=type,
            ele_type="beam",
            ele_list=ele_2,
            ele_properties=self.ele_properties[1],
        )
        global_matrix += self._assemble_matrix(
            type=type,
            ele_type="beam",
            ele_list=ele_3,
            ele_properties=self.ele_properties[2],
        )
        global_matrix += self._assemble_matrix(
            type=type,
            ele_type="beam",
            ele_list=ele_4,
            ele_properties=self.ele_properties[3],
        )
        global_matrix += self._assemble_matrix(
            type=type,
            ele_type="beam",
            ele_list=ele_5,
            ele_properties=self.ele_properties[4],
        )
        global_matrix += self._assemble_matrix(
            type=type,
            ele_type="beam",
            ele_list=ele_6,
            ele_properties=self.ele_properties[5],
        )
        return global_matrix

    def apply_boundary_condition(self, global_matrix):
        node_set = [0, 8, 27, 35]
        contsrained_dof = [[0, 1, 2, 3, 5], [1, 2, 3, 5], [0, 2, 3, 5], [2, 3, 5]]
        dof_list = self._dof_list()
        std_glb_mtx = np.std(np.abs(np.diag(global_matrix)))

        for i, node in enumerate(node_set):
            node_dof = [dof_list[node][j] for j in contsrained_dof[i]]
            for dof in node_dof:
                global_matrix[dof, :] = 0
                global_matrix[:, dof] = 0
                global_matrix[dof, dof] = std_glb_mtx

        return global_matrix

    def plot_mode_shape(self, order=0, factor=40.0):
        node_coord = self._generate_nodes()
        ele_lists = self.element_lists()

        # 3D plot
        fig = plt.figure(figsize=(14, 10))
        ax = fig.add_subplot(111, projection="3d")
        ax.scatter(
            [i[0] for i in node_coord],
            [i[1] for i in node_coord],
            [i[2] for i in node_coord],
            c="gray",
            marker="o",
        )
        for ele_list in ele_lists:
            for ele in ele_list:
                node_1 = ele[0]
                node_2 = ele[1]
                ax.plot(
                    [node_coord[node_1][0], node_coord[node_2][0]],
                    [node_coord[node_1][1], node_coord[node_2][1]],
                    [node_coord[node_1][2], node_coord[node_2][2]],
                    color="gray",
                    alpha=0.4,
                )
        ax.set_axis_off()
        ax.set_aspect("equal")
        v, w = self.freqs_modes()
        print(v)
        # keep to two decimal places of numbers for all elemenbts in v
        v = np.round(v, 3)
        # there are multiple repeated numbers at the start of v, find the first index of the first un-repeated number
        first_index = 0
        # for i in range(len(v)):
        #     if v[i] != v[0]:
        #         first_index = i
        #         break
        order_index = first_index + order
        mode_shape = w[:, order_index]
        for i in range(len(node_coord)):
            node_coord[i][0] += mode_shape[6 * i] * factor
            node_coord[i][1] += mode_shape[6 * i + 1] * factor
            node_coord[i][2] += mode_shape[6 * i + 2] * factor
        ax.scatter(
            [i[0] for i in node_coord],
            [i[1] for i in node_coord],
            [i[2] for i in node_coord],
            c="r",
            marker="o",
        )
        for ele_list in ele_lists:
            for ele in ele_list:
                node_1 = ele[0]
                node_2 = ele[1]
                ax.plot(
                    [node_coord[node_1][0], node_coord[node_2][0]],
                    [node_coord[node_1][1], node_coord[node_2][1]],
                    [node_coord[node_1][2], node_coord[node_2][2]],
                    color="k",
                )
        ax.text2D(
            0.05,
            0.8,
            f"Mode {order+1} Frequency: {v[order_index]} Hz",
            transform=ax.transAxes,
        )
        plt.show()

    def plot_force(self):
        f_func = self.force_func_sine(factor=1e5)
        time = self.t_eval
        force = f_func(time)
        plt.plot(time, force)
        plt.show()

    def force_func(self, factor=1.0):
        import scipy.interpolate as spi

        force = []
        with open("./data/force.csv", "r") as f:
            for line in f:
                force.append([float(i) for i in line.split(sep=",")])
        force = np.array(force)
        xp = force[:, 0]
        yp = force[:, 1] * factor
        func = spi.interp1d(xp, yp, kind="cubic", fill_value=(0, 0), bounds_error=False)
        return func

    def force_func_sine(self, factor=1.0):
        def func(t):
            # if isinstance(t, (int, float)):
            #     if t < 9:
            #         q = t
            #         return (
            #             factor * np.sin(2 * 2 * np.pi * q)
            #             + factor * np.sin(7 * 2 * np.pi * q) * 0.8
            #             + factor * np.sin(15 * 2 * np.pi * q) * 0.6
            #             + factor * np.sin(17 * 2 * np.pi * q) * 0.3
            #         )
            #     else:
            #         return 0
            # else:
            #     # t is an array
            #     value = np.zeros_like(t)
            #     for i in range(len(t)):
            #         if t[i] < 9:
            #             q = t[i]
            #             value[i] = (
            #                 factor * np.sin(2 * 2 * np.pi * q)
            #                 + factor * np.sin(7 * 2 * np.pi * q) * 0.8
            #                 + factor * np.sin(15 * 2 * np.pi * q) * 0.6
            #                 + factor * np.sin(17 * 2 * np.pi * q) * 0.3
            #             )
            #     return value
            return (
                factor * np.sin(2 * 2 * np.pi * t)
                + factor * np.sin(7 * 2 * np.pi * t) * 0.8
                + factor * np.sin(15 * 2 * np.pi * t) * 0.6
                + factor * np.sin(17 * 2 * np.pi * t) * 0.3
            )

        return func

    def transient_impact_force(self, factor=1.0):
        # give the formula of the impact force using gamma function

        def func(t):
            if isinstance(t, (int, float)):
                q = t - 0.3
                if q < 0:
                    return 0
                else:
                    return factor * q**2 * np.exp(-12 * q)
            else:
                # t is an array
                value = np.zeros_like(t)
                for i in range(len(t)):
                    q = t[i] - 0.3
                    if q < 0:
                        value[i] = 0
                    else:
                        value[i] = factor * q**2 * np.exp(-12 * q)
                return value

        return func

    def transient_impact_force_2(self, factor=1.0):
        # give the formula of the impact force using gaussian distribution function

        def func(t):
            if isinstance(t, (int, float)):
                q = t - 0.6
                return factor * np.exp(-40 * (q**2))
            else:
                # t is an array
                value = np.zeros_like(t)
                for i in range(len(t)):
                    q = t[i] - 0.6
                    value[i] = factor * np.exp(-40 * (q**2))
                return value

        return func

    def nrmse(self, y_true, y_pred):
        return np.linalg.norm(y_true - y_pred) / np.linalg.norm(
            y_true - np.mean(y_true)
        )

    def compute_response(self, save_path="./data/"):
        if save_path[-1] == "/":
            file_name = save_path + "collision_resp.pkl"
        else:
            file_name = save_path
        resp_mtx = self.response(method="Radau").T
        # add noise, snr is signal to noise ratio
        with open(file_name, "wb") as f:
            pickle.dump(resp_mtx, f)

        print("Collision response data has been saved to " + file_name)
        return resp_mtx

    def plot_force_responses(self, computed=True):
        # set the fonttype to be Arial
        plt.rcParams["font.family"] = "Times New Roman"
        # set the font size's default value
        plt.rcParams.update({"font.size": 8})
        ts = {"fontname": "Times New Roman"}
        cm = 1 / 2.54  # centimeters in inches
        if computed:
            with open("./data/collision_resp.pkl", "rb") as f:
                resp_mtx = pickle.load(f)
        else:
            resp_mtx = self.compute_response()

        fig, axs = plt.subplots(1, 2, figsize=(16.4 * cm, 7 * cm))
        time = self.t_eval
        axs[0].plot(time, self.f_t[0](time), "k", label="Input force", linewidth=1.2)
        axs[0].set_ylabel("Impact force on node 31 along y axis (N)", **ts)
        axs[0].set_xlabel("Time (s)", **ts)
        axs[0].tick_params(axis="both", which="major", direction="in")
        axs[0].set_ylim(-0.5 * 1e5, 4e5)
        axs[0].set_yticks(np.arange(0, 4e5 + 1, 1e5))
        axs[0].set_xticks(np.arange(0, 6.1, 1))
        axs[0].set_xlim(0, 6)
        # set y-axis to be in scientific notation
        axs[0].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        axs[0].text(-0.1, -0.1, "(a)", transform=axs[0].transAxes, **ts)
        axs[0].grid(True)
        snr = 30
        exp_sq = np.mean(resp_mtx[:, 0] ** 2, axis=0)
        std_noise_sq = exp_sq / snr
        std = np.sqrt(std_noise_sq)
        noise = np.random.normal(0, std, resp_mtx.shape[0])
        axs[1].plot(time, (resp_mtx[:, 0] + noise), "b", linewidth=1.0)
        axs[1].set_ylabel("Acceleration on node 41 along y axis (m/s$^2$)", **ts)
        axs[1].set_xlabel("Time (s)", **ts)
        axs[1].set_ylim(-20, 10.0)
        axs[1].set_xlim(0, 6)
        axs[1].set_yticks(np.arange(-20.0, 11, 10.0))
        axs[1].set_xticks(np.arange(0, 6.1, 1))
        axs[1].tick_params(axis="both", which="major", direction="in")
        axs[1].arrow(0.3, -17, 1.0, 0, head_width=0.8, head_length=0.2, fc="k", ec="k")
        axs[1].arrow(1.3, -17, -1.0, 0, head_width=0.8, head_length=0.2, fc="k", ec="k")
        axs[1].text(0.1, -19, "Impact force dominated", **ts)

        axs[1].arrow(1.5, -12, 4.3, 0, head_width=0.8, head_length=0.2, fc="k", ec="k")
        axs[1].arrow(5.7, -12, -4.1, 0, head_width=0.8, head_length=0.2, fc="k", ec="k")
        axs[1].text(2.8, -14, "Free vibration dominated", **ts)
        # axs[1].legend(
        #     fontsize=8,
        #     facecolor="white",
        #     edgecolor="black",
        #     ncol=2,
        #     loc="upper right",
        #     prop={"family": "Times New Roman"},
        # )
        axs[1].text(-0.1, -0.1, "(b)", transform=axs[1].transAxes, **ts)
        axs[1].grid(True)

        plt.tight_layout(pad=0.2)
        plt.savefig("./manuscript/figures/force_responses.png", dpi=300)
        plt.savefig("./manuscript/figures/force_responses.svg", dpi=300)
        plt.show()