import time
import numpy as np
import EvilPlotting as plot

N3 = 160  # p3 precision
N4 = 80  # p4 precision


class solver:
    def __init__(self, v_data=None):
        self.g = None
        self.x0 = None
        self.straight_fac = None
        self.tf_ = None
        self.r2 = None
        self.r1 = None
        self.m_wet_log = None
        self.m_wet = None
        self.p_cs_cos = None
        self.y_gs_cot = None
        self.V_max = None
        self.G_max = None
        self.alpha = None
        self.Isp_inv = None
        if v_data is not None:
            self.set_params(v_data)

    def set_params(self, v_data):
        self.Isp_inv = 1 / v_data['Isp']
        self.alpha = 1 / 9.80665 / v_data['Isp']
        self.G_max = v_data['G_max']
        self.V_max = v_data['V_max']
        self.y_gs_cot = 1 / np.tan(v_data['y_gs'])
        self.p_cs_cos = np.cos(v_data['p_cs'])
        self.m_wet = v_data['m_wet']
        self.m_wet_log = np.log(v_data['m_wet'])
        self.r1 = v_data['T_max'] * v_data['throt'][0]
        self.r2 = v_data['T_max'] * v_data['throt'][1]
        self.tf_ = v_data['tf']
        self.straight_fac = v_data['straight_fac']
        self.x0 = v_data['x0']
        self.g = v_data['g']

    def pack_data(self, N):
        dt = self.tf_ / N
        alpha_dt = self.alpha * dt
        t = np.linspace(0, (N - 1) * dt, N)
        z0_term = self.m_wet - self.alpha * self.r2 * t
        z0_term_inv = (1 / z0_term)
        z0_term_log = np.log(z0_term)
        x0 = self.x0.reshape(6)
        g = self.g.reshape(3)
        sparse_params = np.array((alpha_dt, self.G_max, self.V_max, self.y_gs_cot, self.p_cs_cos, self.m_wet_log, self.r1, self.r2, self.tf_, self.straight_fac))
        sparse_params = sparse_params.reshape(len(sparse_params), 1)
        return x0, z0_term_inv, z0_term_log, g, sparse_params

    def solve_direct(self):
        print("------solve_direct-------")
        import GFOLD_direct_exec as solver_direct
        start = time.time()
        packed_data = self.pack_data(N3)
        obj_opt, x, u, m, s, z = solver_direct.GFOLD_direct(N3, 'p3', packed_data)
        if obj_opt is None:
            print('p3 failed')
            return None
        tf_m = self.tf_
        for i in range(x.shape[1]):
            if (np.linalg.norm(x[0:3, i]) + np.linalg.norm(x[3:6, i])) < 0.1:
                tf_m = i / x.shape[1] * self.tf_
                break
        print('tf_m:' + str(tf_m))
        self.tf_ = tf_m + 0.1 * self.straight_fac
        packed_data = self.pack_data(N4)
        obj_opt, x, u, m, s, z = solver_direct.GFOLD_direct(N4, 'p4', packed_data)
        if obj_opt is None:
            print('p4 failed')
            return None
        print("------solved in %fs-------" % (time.time() - start))
        return tf_m, x, u, m, s, z


if __name__ == '__main__':
    test_vessel = {
        'Isp': 250,
        'G_max': 100,
        'V_max': 200,
        'y_gs': np.radians(45),
        'p_cs': np.radians(45),
        'm_wet': 5.5e3,
        'T_max': 168e3,
        'throt': [0.1, 0.8],
        'x0': np.array([1500, 150, 200, -50, 30, 20]),
        'g': np.array([-9.8, 0, 0]),
        'tf': 40,
        'straight_fac': 5,
    }
    try:
        plot.plot_run3D(*solver(test_vessel).solve_direct(), test_vessel)
    except TypeError:
        print("solve failed")
