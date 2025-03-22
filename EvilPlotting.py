import math
import numpy as np
import numpy.linalg as npl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def plot_run3D(tf, x, u, m, s, z, v_data):
    t = np.linspace(0, tf, num=len(m.T))

    r = np.array(x[0:3, :])
    v = np.array(x[3:6, :])
    z = np.array(z)
    s = np.array(s)
    u = np.array(u)
    m = np.array(m)

    r1 = v_data['T_max'] * v_data['throt'][0]
    r2 = v_data['T_max'] * v_data['throt'][1]

    if t.shape == () or r.shape == () or v.shape == () or u.shape == ():
        print('data actually empty')
        return

    Th = [np.linalg.norm(u[:, i]) * m[i] for i in range(len(v.T))]
    Th_ = [t_ for t_ in Th]
    for i in range(len(Th_) - 1):
        Th_[i] = (Th[i] + Th[i + 1]) / 2
    Th_[-1] *= 0
    vnorm = [np.linalg.norm(vel) for vel in v.T]

    print(f"Th0: {Th[0]} v0: {v[0, 0]}")
    print(f"Th1: {Th[1]} v0: {v[0, 1]}")
    print(f"Final mass: {m[-1]}")
    # u_dirs_1 = [90 - np.degrees(np.atan2(u[0,n], u[1,n])) for n in range(p.N)]
    # u_dirs_2 = [90 - np.degrees(np.atan2(u[0,n], u[2,n])) for n in range(p.N)]

    traj = plt.figure()
    ax = traj.add_axes(Axes3D(traj))
    ax.set_aspect('equal')

    r_ = np.linspace(0, max(max(r[1, :]), max(r[2, :])), 7)
    a_ = np.linspace(0, 2 * np.pi, 20)
    R, P = np.meshgrid(r_, a_)
    X, Y, Z = R * np.cos(P), R * np.sin(P), R * (np.tan(v_data['y_gs']))
    # X,Y,Z=R*np.cos(P), R*np.sin(P),((R**2 - 1)**2)

    # ax.plot(x(t),y(t),z(t),label='Flight Path')
    ax.plot(r[1, :], r[2, :], r[0, :], label='Flight Path')
    ax.plot(r[1, ::5], r[2, ::5], r[0, ::5], linestyle='None', marker='.')
    ax.plot_surface(X, Y, Z, cmap=plt.cm.YlGnBu_r, alpha=0.5)
    
    for n in range(0, u.shape[1]-1, 5):
        pos_x = r[1, n]
        pos_y = r[2, n]
        pos_z = r[0, n]
        
        thrust_x = u[1, n]
        thrust_y = u[2, n]
        thrust_z = u[0, n]
        
        ax.quiver(pos_x, pos_y, pos_z, 
                thrust_x, thrust_y, thrust_z, 
                length=80, normalize=True, color='orange')

    # Tweak the limits and add latex math labels.

    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
    ax.set_zlabel(r'$z$')

    ax.legend()

    f = plt.figure()
    f.add_subplot(321)

    plt.plot(t, vnorm)
    y = str(v_data['V_max'])
    x = np.array(range(0, int(max(t))))
    plt.plot(x, eval('0*x+' + y))
    plt.title('Velocity Magnitude (m/s)')

    plt.subplot(3, 2, 2)
    plt.plot(t, r[0, :])
    plt.title('Altitude (m)')

    plt.subplot(3, 2, 3)
    plt.plot(t, m)
    plt.title('Mass (kg)')

    plt.subplot(3, 2, 4)
    plt.plot(t, Th)
    # y = str(v_data['T_max'])
    x = np.array(range(0, int(max(t))))
    # print(eval('0*x+'+y))
    # plt.plot(x,eval('0*x+'+y))
    plt.plot(x, 0 * x + r1)
    plt.plot(x, 0 * x + r2)
    plt.plot(t, Th_)
    plt.plot(x, )
    plt.title('Thrust (N)')

    plt.subplot(3, 2, 5)
    u_angle = [np.degrees(math.acos(min(1, ui[0] / npl.norm(ui)))) for ui in u.T]
    plt.plot(x, 0 * x + np.degrees(v_data['p_cs']))
    plt.plot(t, u_angle)
    plt.title('Thrust angle')

    alpha = 1 / 9.80665 / v_data['isp']
    z0_term = (v_data['m_wet'] - alpha * r2)  # see ref [2], eq 34,35,36
    z1_term = (v_data['m_wet'] - alpha * r1)
    lim = []
    lim2 = []
    n = 0
    z = z.flatten()
    for t_ in t:
        if t_ > 0:
            try:
                v = r2 / (z0_term * t_) * (1 - (z[n] - np.log(z0_term * t_)))
            except ZeroDivisionError:
                v = 0
            lim.append(v)
            try:
                v = r1 / (z1_term * t_) * (1 - (z[n] - np.log(z0_term * t_)) + (1 / 2) * (z[n] - np.log(z0_term * t_)) ** 2)
            except ZeroDivisionError:
                v = 0
            lim2.append(v)
        else:
            lim.append(0)
            lim2.append(0)
        n += 1
    # lim = np.array(lim).flatten()
    plt.subplot(3, 2, 6)
    # plt.plot(t,lim)
    # plt.plot(t,lim2)
    s = s.flatten()
    if s.shape == (1, 65):
        s.reshape((65,))
        print('reshape', s)
    # print('s',s)
    plt.plot(t, s)
    plt.title('Sigma Slack')

    plt.show()
