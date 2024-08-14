import numpy as np

class Bar3D:
    """
    3D bar element class
    node1coord: 3D coordinates of node 1, 1d numpy array
    node2coord: 3D coordinates of node 2, 1d numpy array
    E: Young's modulus, float
    A: cross-sectional area, float
    """
    def __init__(self, node1coord, node2coord, E, A, rho):
        self.node1coord = node1coord
        self.node2coord = node2coord
        self.E = E
        self.A = A
        self.rho = rho
        self.T = self.transformation_matrix()

    def length(self):
        return np.linalg.norm(self.node2coord - self.node1coord)
    
    def direction_cosines(self):
        L = self.length()
        return (self.node2coord - self.node1coord) / L
    
    def transformation_matrix(self):
        l = self.direction_cosines()
        T = np.array([[l[0], l[1], l[2], 0, 0, 0],
                      [0, 0, 0, l[0], l[1], l[2]]])
        return T
    
    def stiffness_matrix(self):
        L = self.length()
        k = self.A * self.E / L
        K = k * np.array([[1, -1],
                            [-1, 1]])
        K = np.dot(np.dot(self.T.T, K), self.T)
        return K

    def mass_matrix(self, rho):
        L = self.length()
        m = self.A * rho * L / 6
        M = m * np.array([[2, 1],
                            [1, 2]])
        M = np.dot(np.dot(self.T.T, M), self.T)
        return M

class Beam3D:
    """
    3D beam element class
    node1coord: 3D coordinates of node 1, 1d numpy array
    node2coord: 3D coordinates of node 2, 1d numpy array
    E: Young's modulus, float
    A: cross-sectional area, float
    I: moment of inertia, float
    """
    def __init__(self, node1coord, node2coord, node3coord, E, G, A, J, Iy, Iz, rho):
        self.node1coord = node1coord
        self.node2coord = node2coord
        self.node3coord = node3coord
        self.E = E
        self.A = A
        self.G = G
        self.J = J
        self.Iy = Iy
        self.Iz = Iz
        self.rho = rho
        self.T = self.transformation_matrix()

    def length(self):
        return np.linalg.norm(self.node2coord - self.node1coord)
    
    def direction_cosine_mtx(self):
        L = self.length()
        # Create vectors for coordinates
        coords = np.array([self.node1coord, self.node2coord, self.node3coord]).T

        # Calculate differences
        diffs = np.diff(coords, axis=1)

        # Calculate direction cosines
        l_x, m_x, n_x = diffs[:, 0] / L

        # Calculate constants for area calculation
        consts = diffs[:, 0] - diffs[:, 1][::-1]

        # Calculate area
        A_123 = 2 * np.linalg.norm(consts)

        # Calculate direction cosines for z
        l_z, m_z, n_z = consts / A_123

        # Calculate direction cosines for y
        l_y = m_z * n_x - m_x * n_z
        m_y = n_z * l_x - n_x * l_z
        n_y = l_z * m_x - l_x * m_z

        # Transformation matrix
        T3= np.array([
            [l_x, m_x, n_x],
            [l_y, m_y, n_y],
            [l_z, m_z, n_z],
        ])
        return T3
    
    def transformation_matrix(self):
        T3 = self.direction_cosine_mtx()
        T = np.zeros((12, 12))
        T[:3, :3] = T3
        T[3:6, 3:6] = T3
        T[6:9, 6:9] = T3
        T[9:, 9:] = T3
        return T

    def stiffness_matrix(self):
        L = self.length()
        EAL = self.E * self.A / L
        EIzL2 = self.E * self.Iz / L**2
        EIyL2 = self.E * self.Iy / L**2
        GJL = self.G * self.J / L
        K = np.array(
            [
                [EAL, 0, 0, 0, 0, 0, -EAL, 0, 0, 0, 0, 0],
                [0,12 * EIzL2 / L,0,0,0,6 * EIzL2,0,-12 * EIzL2 / L,0,0,0,6 * EIzL2],
                [0,0,12 * EIyL2/L,0,-6 * EIyL2,0,0,0,-12 * EIyL2/L,0,-6 * EIyL2,0],
                [0, 0, 0, GJL, 0, 0, 0, 0, 0, -GJL, 0, 0],
                [0,0,-6 * EIyL2,0,4 * EIyL2*L,0,0,0,6 * EIyL2,0,2 * EIyL2*L,0],
                [0,6 * EIzL2,0,0,0,4 * EIzL2 * L,0,-6 * EIzL2,0,0,0,2 * EIzL2 * L],
                [-EAL, 0, 0, 0, 0, 0, EAL, 0, 0, 0, 0, 0],
                [0,-12 * EIzL2 / L,0,0,0,-6 * EIzL2,0,12 * EIzL2 / L,0,0,0,-6 * EIzL2],
                [0,0,-12 * EIyL2/L,0,6 * EIyL2,0,0,0,12 * EIyL2/L,0,6 * EIyL2,0],
                [0, 0, 0, -GJL, 0, 0, 0, 0, 0, GJL, 0, 0],
                [0,0,-6 * EIyL2,0,2 * EIyL2*L,0,0,0,6 * EIyL2,0,4 * EIyL2*L,0],
                [0,6 * EIzL2,0,0,0,2 * EIzL2*L,0,-6 * EIzL2,0,0,0,4 * EIzL2*L],
            ]
        )
        K = np.dot(np.dot(self.T.T, K), self.T)
        return K

    def mass_matrix(self):
        L = self.length()
        m = self.A * self.rho * L / 420
        Ix = self.Iy + self.Iz
        A = self.A
        M = m * np.array(
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
        M = np.dot(np.dot(self.T.T, M), self.T)
        return M
 
    
        
        

    