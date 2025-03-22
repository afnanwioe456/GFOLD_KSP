import numpy as np
import cvxpy as cp
import warnings

warnings.filterwarnings("ignore", module="cvxpy")

''' As defined in the paper...

 PROBLEM 3: Minimum Landing Error (tf roughly solved)
 MINIMIZE : norm of landing error vector
 SUBJ TO  :
            0) initial conditions satisfied (position, velocity)
            1) final conditions satisfied (altitude, velocity)
            2) dynamics always satisfied
            3) x stays in cone at all times
            4) relaxed convexified mass and thrust constraints
            5) thrust pointing constraint
            6) sub-surface flight constraint

 PROBLEM 4: Minimum Fuel Use
 MAXIMIZE : landing mass, opt variables are dynamical and
 SUBJ TO  :
            0) same constraints as p1, plus:
            1) landing point must be equal or better than that found by p1

'''

class GFOLDSolver:
    def __init__(self, N):
        self.N = N
        self._params: dict[str, cp.Parameter]
        self._vars: dict[str, cp.Variable]
        self._params, self._vars = self._init_vp()
        self._p3 = self._build(3)
        self._p4 = self._build(4)

    def _init_vp(self):
        params = {}
        params['x0'] = cp.Parameter((6,), name='x0')
        params['xf'] = cp.Parameter((6,), name='xf')
        params['g'] = cp.Parameter((3,), name='g')
        params['m_wet'] = cp.Parameter(name='m_wet')
        params['alpha'] = cp.Parameter(name='alpha')
        params['z0_inv'] = cp.Parameter((self.N,), name='z0_inv', pos=True)
        params['z0_log'] = cp.Parameter((self.N,), name='z0_log', pos=True)
        params['tf'] = cp.Parameter(name='tf')
        params['r1'] = cp.Parameter(name='r1', pos=True)
        params['r2'] = cp.Parameter(name='r2', pos=True)  # make DCP checker happy
        params['V_max'] = cp.Parameter(name='V_max')
        params['G_max'] = cp.Parameter(name='G_max')
        params['y_gs_cot'] = cp.Parameter(name='y_gs_cot')
        params['p_cs_cos'] = cp.Parameter(name='p_cs_cos')
        params['straight_fac'] = cp.Parameter(name='straight_fac', pos=True)
        
        vars = {}
        vars['x'] = cp.Variable((6, self.N), name='var_x')      # state vector (3position,3velocity)
        vars['u'] = cp.Variable((3, self.N), name='var_u')      # u = Tc/mass because Tc[:,n]/m[n] is not allowed by DCP
        vars['z'] = cp.Variable(self.N, name='var_z')           # z = ln(mass)
        vars['s'] = cp.Variable(self.N, name='var_s')           # thrust slack parameter
        return params, vars

    def _build(self, pmark):
        x0 = self._params['x0']
        xf = self._params['xf']
        g = self._params['g']
        m_wet = self._params['m_wet']
        alpha = self._params['alpha']
        tf = self._params['tf']
        r1 = self._params['r1']
        r2 = self._params['r2']
        V_max = self._params['V_max']
        G_max = self._params['G_max']
        y_gs_cot = self._params['y_gs_cot']
        p_cs_cos = self._params['p_cs_cos']
        straight_fac = self._params['straight_fac']
        z0_inv = self._params['z0_inv']
        z0_log = self._params['z0_log']

        dt = tf / self.N
        x = self._vars['x']
        u = self._vars['u']
        z = self._vars['z']
        s = self._vars['s']

        con = []                                                            # CONSTRAINTS LIST
        con += [x[:, 0] == x0]                                              # initial pos and vel
        con += [x[3:6, self.N - 1] == xf[3:6]]                              # final vel
        con += [s[self.N - 1] == 0]                                         # thrust at the end must be zero
        con += [u[:, 0] == s[0] * np.array([1, 0, 0])]                      # thrust direction starts straight
        con += [u[:, self.N - 1] == s[self.N - 1] * np.array([1, 0, 0])]    # and ends straight
        con += [z[0] == cp.log(m_wet)]                                      # convexified (7)

        # final pos
        if pmark == 3:
            con += [x[0, self.N - 1] == xf[0]]
        elif pmark == 4:
            con += [x[0:3, self.N - 1] == xf[0:3]]

        for n in range(0, self.N - 1):
            
            con += [x[3:6, n + 1] == x[3:6, n] + (dt * 0.5) * ((u[:, n] + g) + (u[:, n + 1] + g))]
            con += [x[0:3, n + 1] == x[0:3, n] + (dt * 0.5) * (x[3:6, n + 1] + x[3:6, n])]
            con += [cp.norm((x[1:3, n] - xf[1:3])) - y_gs_cot * (x[0, n] - xf[0]) <= 0]  # glideslope cone

            # con += [cp.norm(x[3:6, n]) <= V_max]  # velocity
            con += [z[n + 1] == z[n] - (alpha * dt * 0.5) * (s[n] + s[n + 1])]  # mass decreases
            con += [cp.norm(u[:, n]) <= s[n]]  # limit thrust magnitude & also therefore, mass
            con += [u[0, n] >= p_cs_cos * s[n]]  # Thrust pointing

            if n > 0:
                mu_1 = r1 * z0_inv[n]
                mu_2 = r2 * z0_inv[n]

                # https://www.desmos.com/calculator/wtcfgnepe1
                con += [s[n] >= mu_1 * (1 - (z[n] - z0_log[n]) + cp.square(z[n] - z0_log[n]) * 0.5)]  # lower thrust bound
                con += [s[n] <= mu_2 * (1 - (z[n] - z0_log[n]))]  # upper thrust bound

            # con += [x[0,0:N-1] >= 0] # no

        if pmark == 3:
            expression = 0
            for i in range(self.N):
                expression += cp.norm(x[0:3, i]) * (i / self.N)
            objective = cp.Minimize(expression)
            # objective = cp.Minimize(cp.norm(x[0:3, self.N - 1] - xf[0:3]))
            return cp.Problem(objective, con)
        elif pmark == 4:
            expression = 0
            for i in range(self.N):
                expression += cp.norm(x[4:6, i]) * (i / self.N)
            expression *= straight_fac
            expression += -z[self.N - 1] * self.N
            objective = cp.Minimize(expression)
            # objective = cp.Maximize(z[self.N - 1])
            return cp.Problem(objective, con)
    
    def solve(self, v_data, solver='ECOS', **kwargs):
        self._params['alpha'].value = 1 / 9.80665 / v_data['isp']
        self._params['G_max'].value = v_data['G_max']
        self._params['V_max'].value = v_data['V_max']
        self._params['y_gs_cot'].value = 1 / np.tan(v_data['y_gs'])
        self._params['p_cs_cos'].value = np.cos(v_data['p_cs'])
        self._params['m_wet'].value = v_data['m_wet']
        self._params['r1'].value = v_data['T_max'] * v_data['throt'][0]
        self._params['r2'].value = v_data['T_max'] * v_data['throt'][1]
        self._params['tf'].value = v_data['tf']
        self._params['straight_fac'].value = v_data['straight_fac']
        self._params['x0'].value = v_data['x0']
        self._params['xf'].value = v_data['xf']
        self._params['g'].value = v_data['g']
        t = np.linspace(0, (self.N - 1) * (v_data['tf'] / self.N), self.N)
        z0 = self._params['m_wet'].value - self._params['alpha'].value * self._params['r2'].value * t
        self._params['z0_inv'].value = 1 / z0
        self._params['z0_log'].value = np.log(z0)

        opt = self._p3.solve(solver, kwargs)
        if self._vars['x'].value is None:
            raise RuntimeError('Problem 3 failed.')
        x = self._vars['x'].value
        xf = self._params['xf'].value
        tf = self._params['tf'].value
        for i in range(x.shape[1]):
            if (np.linalg.norm(x[0:3, i] - xf[0:3]) + np.linalg.norm(x[3:6, i] - xf[3:6])) < 0.1:
                tf = i / x.shape[1] * tf
                break
        self._params['tf'].value = tf

        opt = self._p4.solve(solver, kwargs)
        if self._vars['x'].value is None:
            raise RuntimeError('Problem 4 failed.')
        x = self._vars['x'].value
        u = self._vars['u'].value
        s = self._vars['s'].value
        z = self._vars['z'].value
        m = np.exp(z)
        return tf, x, u, m, s, z
