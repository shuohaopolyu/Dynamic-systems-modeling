import numpy as np
from numpy import linalg as LA
from scipy.linalg import expm
from scipy.integrate import solve_ivp


class MultiDOF:
    """
    Establishing the multi-dof model using the mechanical properties of the structural dynamical system.

    Parameters
    ----------
    mass_mtx: array_like, with shape (n, n)
            Input mass matrix
    stiff: array_like, with shape (n, n)
            Input stiffness matrix
    damp_type: {"Rayleigh", "Proportional"} (default: "Rayleigh")
            The type of method to generate damping matrix
    damp_params: tuple, (default: (0, 4, 0.03))
                The parameters for the damping matrix
    f_dof: list, with length n_f
            Degrees of freedom that the force applied to the system
    resp_dof: list, with length n_s
            Degrees of freedom that the responses are measured
    t_eval: array_like, with shape (n_t, 1)
        Time points to evaluate the response
    f_t: list, with length n_f
        Callable functions of time that return the forces applied to the system
    """

    def __init__(
        self,
        mass_mtx,
        stiff_mtx,
        damp_type="Rayleigh",
        damp_params=(0, 4, 0.03),
        f_dof=[0],
        resp_dof=[0],
        t_eval=np.linspace(0, 1, 100).reshape(-1, 1),
        f_t=None,
    ):
        f""" """
        self.mass_mtx = mass_mtx
        self.stiff_mtx = stiff_mtx
        self.damp_mtx = np.zeros_like(mass_mtx)
        self.DOF = mass_mtx.shape[0]
        self.damp_type = damp_type
        self._update_damp_mtx(damp_params)
        self.f_dof = f_dof
        self.resp_dof = resp_dof
        self.t_eval = t_eval
        self.f_t = f_t
        self.n_s = len(resp_dof)
        self.n_f = len(f_dof)
        self.n_t = t_eval.size

    def __str__(self):
        return (
            "A multi-dof dynamical system with %d DOFs" % self.DOF
            + "\n"
            + "The first five natural frequencies are (in Hertz):"
            + "\n"
            + str(np.sort(self.freqs_modes()[0])[0:5])
        )

    def _update_damp_mtx(self, args):
        """
        Update the damping matrix.
        For Rayleigh damping, args = (lower_index, upper_index, damping_ratio)
        For d_mtx damping, args = d_mtx
        """
        if self.damp_type == "Rayleigh":
            freqs_rad_by_s = np.sort(self.freqs_modes()[0]) * 2 * np.pi
            lower_index, upper_index, damping_ratio = args
            alpha = (
                2
                * damping_ratio
                * freqs_rad_by_s[lower_index]
                * freqs_rad_by_s[upper_index]
                / (freqs_rad_by_s[lower_index] + freqs_rad_by_s[upper_index])
            )
            beta = (
                2
                * damping_ratio
                / (freqs_rad_by_s[lower_index] + freqs_rad_by_s[upper_index])
            )
            self.damp_mtx = alpha * self.mass_mtx + beta * self.stiff_mtx
        elif self.damp_type == "d_mtx":
            self.damp_mtx = args
        else:
            raise ValueError("Damping type not supported")

    def freqs_modes(self):
        """
        Return the natural frequencies and mode shapes of the system
        w is a 1*DOF numpy array of natural frequencies in Hertz
        v is a DOF*DOF numpy array of mode shapes
        """
        w, v = LA.eig(LA.inv(self.mass_mtx) @ self.stiff_mtx)
        return np.sqrt(w) / (2 * np.pi), v

    def state_space_mtx(self, type="continuous"):
        """
        type: 'continuous' or 'discrete'
        Return the continuous/discrete time state space matrices of the system A, B, C, D
        Denote the subscript of the state space matrices as 'c' for continuous time and 'd' for discrete time
        A_c, A_d: 2DOF*2DOF numpy arrays
        B_c, B_d: 2DOF*n_f numpy arrays
        C_c, C_d: n_s*2DOF numpy arrays
        D_c, D_d: n_s*n_f numpy arrays
        """
        L = np.zeros((self.DOF, len(self.f_dof)))
        for i in range(self.n_f):
            L[self.f_dof[i], i] = 1
        invM = LA.inv(self.mass_mtx)
        A = np.vstack(
            (
                np.hstack((np.zeros((self.DOF, self.DOF)), np.eye(self.DOF))),
                np.hstack((-invM @ self.stiff_mtx, -invM @ self.damp_mtx)),
            )
        )
        B = np.vstack(
            (
                np.zeros((self.DOF, len(self.f_dof))),
                invM @ L,
            )
        )
        J = np.zeros((self.n_s, self.DOF))
        for i in range(self.n_s):
            J[i, self.resp_dof[i]] = 1
        C = J @ np.hstack((-invM @ self.stiff_mtx, -invM @ self.damp_mtx))
        D = J @ invM @ L
        if type == "continuous":
            return A, B, C, D
        elif type == "discrete":
            dt = self.t_eval[1] - self.t_eval[0]
            A_d = expm(A * dt)
            B_d = (A_d - np.eye(2 * self.DOF)) @ LA.inv(A) @ B
            C_d = C
            D_d = D
            return A_d, B_d, C_d, D_d

    def markov_params(self):
        """
        Return
        -------
        Markov parameters of the system mp
        m_p is a list contains n_t elements, each element is a n_s*n_f numpy array
        """
        A_d, B_d, C_d, D_d = self.state_space_mtx(type="discrete")
        m_p = []
        A_acum = np.eye(2 * self.DOF)
        for i in range(len(self.t_eval)):
            if i == 0:
                m_p.append(D_d)
            else:
                m_p.append(C_d @ A_acum @ B_d)
                A_acum = A_acum @ A_d
        return m_p

    def block_mtx_assemble(self, m_p):
        """
        Assemble the block Hankel matrix
        m_p is a list contains n_t elements, each element is a n_s*n_f numpy array, e.g., [X1, X2, X3, X4]
        H is a (n_s*n_t)*(n_f*n_t) numpy array, i.e., H = [ X1  0  0  0;
                                                            X2 X1  0  0;
                                                            X3 X2 X1  0;
                                                            X4 X3 X2 X1 ].
        """
        H = np.zeros((self.n_s * self.n_t, self.n_f * self.n_t))
        for i, X_i in enumerate(m_p):
            for j in range(self.n_t - i):
                H[
                    (j + i) * self.n_s : (j + i + 1) * self.n_s,
                    j * self.n_f : (j + 1) * self.n_f,
                ] = X_i
        return H

    def transfer_mtx(self, abs_triangular=True):
        """
        abs_triangular: True or False
        Return the time-domain force-response transfer matrix of the system.
        If abs_triangular is True, return the absolute triangular transfer matrix,
        in this case, the force and response vectors are of n_t*n_f and n_t*n_s dimensions, respectively.
        If abs_triangular is False, return the original block-wise transfer matrix,
        in this case, the force and response vectors are of n_f*n_t and n_s*n_t dimensions, respectively.
        """
        m_p = self.markov_params()
        H = self.block_mtx_assemble(m_p)
        if abs_triangular:
            Delta_1 = np.zeros((self.n_t * self.n_s, self.n_t * self.n_s))
            Delta_2 = np.zeros((self.n_t * self.n_f, self.n_t * self.n_f))
            for i in range(self.n_t):
                for j in range(self.n_s):
                    Delta_1[j * self.n_t + i, i * self.n_s + j] = 1
                for k in range(self.n_f):
                    Delta_2[k * self.n_t + i, i * self.n_f + k] = 1
            return Delta_1 @ H @ LA.inv(Delta_2)
        else:
            return H

    def frf_mtx(self, resp_dof, f_dof, omega):
        """
        Return the frequency response function matrix of the system
        The dimension of the matrix is len(omega)*len(resp_dof)*len(f_dof)
        If len(resp_dof) == 1 and len(f_dof) == 1, return a 1D array
        """
        frf_mtx = []
        for i in range(len(omega)):
            inv_H_d = LA.inv(
                -(omega[i] ** 2) * self.mass_mtx
                + 1j * omega[i] * self.damp_mtx
                + self.stiff_mtx
            )
            h_d = inv_H_d[resp_dof, f_dof]
            frf_mtx.append(h_d)
        if len(resp_dof) == 1 and len(f_dof) == 1:
            return np.array(frf_mtx).reshape(-1)
        else:
            return np.array(frf_mtx)

    def f_vec(self, t):
        """
        Return the force vector with respect to time
        """
        return np.array([i_th_force(t) for i_th_force in self.f_t]).reshape(-1, 1)

    def response(self, type="acceleration"):
        """
        type: 'acceleration', 'velocity', 'displacement'
        Return the response of the system to the input forces
        The response is a n_s*n_t numpy array
        Computational method used: Runge-Kutta 4th order
        """
        A_c, B_c, C_c, D_c = self.state_space_mtx(type="continuous")
        print((A_c).shape)
        sol = solve_ivp(
            fun=lambda t, y: A_c @ y + B_c @ self.f_vec(t),
            t_span=(self.t_eval[0], self.t_eval[-1]),
            y0=np.zeros(2 * self.DOF),
            method="RK45",
            t_eval=self.t_eval.reshape(-1),
            vectorized=True,
        )
        if type == "acceleration":
            return C_c @ sol.y + D_c @ self.f_vec(self.t_eval).reshape(self.n_f, -1)
        elif type == "velocity":
            return sol.y[self.DOF + self.resp_dof, :]
        elif type == "displacement":
            return sol.y[self.resp_dof, :]
