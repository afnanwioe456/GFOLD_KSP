import krpc
import time
import simple_pid
import numpy as np
import numpy.linalg as npl
import EvilPlotting as plot

from GFOLD_solver import GFOLDSolver

DEBUG_LINE_THICKNESS = 0.7

def lerp(vec1, vec2, t):
    return t * vec2 + (1 - t) * vec1


def form_v3(f1, f2, f3):
    return np.array((f1, f2, f3))


def form_vec(vec):
    return np.array(vec)


def normalize(vec):
    return vec / npl.norm(vec)


def rotation_mat(src):
    (x, y, z, w) = src
    return np.array([
        [1 - 2 * y ** 2 - 2 * z ** 2, 2 * x * y + 2 * w * z, 2 * x * z - 2 * w * y],
        [2 * x * y - 2 * w * z, 1 - 2 * x ** 2 - 2 * z ** 2, 2 * y * z + 2 * w * x],
        [2 * x * z + 2 * w * y, 2 * y * z - 2 * w * x, 1 - 2 * x ** 2 - 2 * y ** 2]
    ])


def transform(vec, mat):
    res = np.array(vec) * mat
    return form_v3(res[0, 0], res[0, 1], res[0, 2])


def angle_around_axis(v1, v2, axis):
    axis = normalize(axis)
    v1 = normalize(np.cross(v1, axis))
    v2 = normalize(np.cross(v2, axis))
    direction = np.sign(np.dot(np.cross(v1, v2), axis))
    return direction * np.arccos(np.dot(v1, v2))


# 生成求解器所需的输入
def vessel_profile1(vessel_info):
    est_time = 0  # TODO: 处理点火延迟
    vel_est = vessel_info['vel'] + vessel_info['acceleration'] * est_time
    pos_est = vessel_info['error'] + vessel_info['vel'] * est_time + 0.5 * vessel_info['acceleration'] * est_time ** 2
    return {
        'isp': vessel_info['specific_impulse'],
        'G_max': params['G_max'],
        'V_max': params['V_max'],
        'y_gs': params['y_gs'] * deg2rad,
        'p_cs': max_tilt * 0.3,
        'm_wet': vessel_info['mass'],
        'T_max': vessel_info['max_thrust'],
        'throt': throttle_limit,
        'x0': np.array([pos_est[0], pos_est[1], pos_est[2], vel_est[0], vel_est[1], vel_est[2]]),
        'xf': np.zeros(6),
        'g': np.array([-g0, 0, 0]),
        'tf': params['tf'],
        'straight_fac': params['straight_fac'],
    }


def update_lines(x, u):
    m_u = vessel.max_thrust / vessel.mass
    for i in range(N - 1):
        lines[i].start = x[0:3, i]
        lines[i].end = x[0:3, i + 1]
        lines[i].thickness = DEBUG_LINE_THICKNESS
    for i in range(N):
        directions[i].start = x[0:3, i]
        directions[i].end = x[0:3, i] + u[:, i] * 20 / m_u
        directions[i].thickness = DEBUG_LINE_THICKNESS


# 找到路径上最近的点的索引值
def find_nearest_index(x, r):
    nearest_mag = np.inf
    nearest_i = 0
    for i in range(x.shape[1]):
        mag = npl.norm(x[0:3, i] - r)  # + npl.norm(x[3:6, i] - v) * 0.2
        if mag < nearest_mag:
            nearest_mag = mag
            nearest_i = i
    # 添加一个速度引起的偏移项, 用于采样插值
    v = x[3:6, nearest_i]
    v_norm = npl.norm(v)
    v_dir = v / v_norm
    frac = np.clip(np.dot(r - x[0:3, nearest_i], v_dir) / (gtf / N * v_norm), -0.5, 0.5)
    return nearest_i + frac


# 路径上采样
def sample_index(index):
    # if index >= N-1:
    if index >= N - 1:
        return form_v3(0, 0, 0), form_v3(0, 0, 0), form_v3(9.807, 0, 0)
    elif index <= 0:
        i = 0
        frac = index
    else:
        i = int(np.floor(index))
        frac = index - i
    x_i_s = lerp(gx[:, i], gx[:, i + 1], frac)
    u_i_s = lerp(gu[:, i], gu[:, i + 1], frac)
    if index < 0:
        u_i_s = gu[:, 1].copy()
    return x_i_s[0:3].copy(), x_i_s[3:6].copy(), u_i_s.copy()


def conic_clamp(target, min_mag, max_mag, max_t):
    # a_mag = npl.norm(target)
    hor_dir = form_v3(0., target[1], target[2])
    hor_dir /= npl.norm(hor_dir)
    # target_direction = target_a / a_mag
    a_hor = npl.norm(target[1:3])
    a_ver = target[0]

    if a_hor < min_mag * np.sin(max_t):
        a_ver_min = np.sqrt(min_mag ** 2 - a_hor ** 2)
    else:
        a_ver_min = np.cos(max_t) * min_mag

    if a_hor < max_mag * np.sin(max_t):
        a_ver_max = np.sqrt(max_mag ** 2 - a_hor ** 2)
    else:
        a_ver_max = np.cos(max_t) * max_mag
    a_ver = np.clip(a_ver, a_ver_min, a_ver_max)

    a_hor = min(a_hor, a_ver * np.tan(max_t))

    return hor_dir * a_hor + form_v3(a_ver, 0, 0)


# 求解最优路径
def solve_gfold(v_data):
    global gfold_path, n_i, nav_mode
    gfold_path = solver.solve(v_data)
    n_i = -100
    if gfold_path is not None:
        tf, x, u, m, s, z = gfold_path
        print("Slack0:", s[0])
        print("Velocity1:", x[3:6, 1])
        print("Final mass:", m[-1])
        # plot.plot_run3D(*gfold_path, v_data)
        if debug_lines:
            update_lines(x, u)
        print('gfold')
        nav_mode = 'gfold'
    conn.krpc.paused = False


if __name__ == '__main__':
    params = {}
    with open('params.txt', 'r', encoding='utf-8') as f:
        for line in f:
            pair = line.split('#')[0].split('=')
            if len(pair) == 2:
                key = pair[0].strip()
                value = eval(pair[1])
                params[key] = value

    deg2rad = np.pi / 180
    rad2deg = 180 / np.pi
    g0 = params['g0']

    conn = krpc.connect(name='PID')
    space_center = conn.space_center
    vessel = space_center.active_vessel
    body = vessel.orbit.body

    delta_time = 0.01

    # target
    target_lat = params['target_lat'] * deg2rad
    target_lon = params['target_lon'] * deg2rad
    target_height = params['target_height']
    target_axis = target_height + body.surface_height(target_lat * rad2deg, target_lon * rad2deg) + body.equatorial_radius
    target_body_pos = np.array((np.cos(target_lon) * np.cos(target_lat), np.sin(target_lat), np.sin(target_lon) * np.cos(target_lat))) * target_axis

    # limit
    max_tilt = params['max_tilt'] * deg2rad
    throttle_limit = params['throttle_limit']
    throttle_limit_ctrl = params['throttle_limit_ctrl']

    # rotation
    ctrl_x_rot = simple_pid.PID(Kp=params['ctrl_xz_rot.kp'], Kd=params['ctrl_xz_rot.kd'], differential_on_measurement=False)
    ctrl_y_avel_kp = params['ctrl_y_avel_kp']
    ctrl_z_rot = simple_pid.PID(Kp=params['ctrl_xz_rot.kp'], Kd=params['ctrl_xz_rot.kd'], differential_on_measurement=False)
    k_x = params['k_x']
    k_v = params['k_v']

    # final
    final_throttle = params['final_throttle']
    final_kp = params['final_kp']

    # time
    game_delta_time = 0.02
    game_prev_time = space_center.ut

    # references
    ref_local = vessel.reference_frame
    ref_surface = vessel.surface_reference_frame  # 地面参考系
    ref_body = body.reference_frame
    ref_target_temp = space_center.ReferenceFrame.create_relative(ref_body, position=target_body_pos)
    ref_target = space_center.ReferenceFrame.create_hybrid(ref_target_temp, rotation=ref_surface, velocity=ref_target_temp)

    prev_vel = form_vec(vessel.velocity(ref_target))
    N = 160
    solver = GFOLDSolver(N)
    gfold_path: None | tuple = None
    n_i = -1
    error = form_vec(vessel.position(ref_target))
    debug_lines = params['debug_lines']
    head_line = None
    target_line = None
    target2_line = None
    if debug_lines:
        lines = [conn.drawing.add_line((0, 0, 0), (0, 0, 0), ref_target) for i in range(N - 1)]
        directions = [conn.drawing.add_line((0, 0, 0), (1, 0, 0), ref_target) for i in range(N)]
        target_line = conn.drawing.add_line((0, 0, 0), (1, 0, 0), ref_target)
        target_line.color = (0, 0, 1)
        target2_line = conn.drawing.add_line((0, 0, 0), (1, 0, 0), ref_target)
        target2_line.color = (0, 0, 1)
        head_line = conn.drawing.add_line((0, 0, 0), (1, 0, 0), ref_target)
        head_line.color = (0, 1, 1)
        for line in directions:
            line.color = (0, 1, 0)

    nav_mode = 'none'
    onceFlag = False

    # 添加数据流
    ut_stream = conn.add_stream(getattr, space_center, 'ut')
    error_stream = conn.add_stream(vessel.position, ref_target)
    avel_stream = conn.add_stream(vessel.angular_velocity, ref_surface)
    vel_stream = conn.add_stream(vessel.velocity, ref_target)
    rotation_stream = conn.add_stream(vessel.rotation, ref_surface)
    mass_stream = conn.add_stream(getattr, vessel, 'mass')
    max_thrust_stream = conn.add_stream(getattr, vessel, 'max_thrust')
    flight = vessel.flight(ref_target)
    aerodynamic_force_stream = conn.add_stream(getattr, flight, 'aerodynamic_force')

    while True:
        time.sleep(delta_time)
        ut = ut_stream()
        game_delta_time = ut - game_prev_time
        if game_delta_time < 0.01:  # 意味着游戏中还没有经过一个物理帧，所以不进行计算
            continue

        # 取得载具数据
        error = form_vec(error_stream())  # 目标系里的偏差
        avel = form_vec(avel_stream())  # 地面系下角速度（等于目标系角速度
        vel = form_vec(vel_stream())  # 地面速度
        # rotation_local2srf = rotation_mat(form_vec(rotation_stream()))  # 机体系到地面系旋转矩阵
        # rotation_srf2local = npl.inv(rotation_local2srf)  # 地面系到机体系旋转矩阵
        mass = mass_stream()
        max_thrust = max_thrust_stream()

        if error[0] < params['start_altitude']:
            if not onceFlag:
                vessel_d = {}
                onceFlag = True
                vessel.control.throttle = 0.01
                time.sleep(3)
                # TODO: 推力修正
                corrected_thrust = vessel.max_thrust_at(1)
                corrected_isp = vessel.specific_impulse_at(1)
                vessel_d['max_thrust'] = corrected_thrust  # 修正推力以减缓引擎效率下降问题
                vessel_d['specific_impulse'] = corrected_isp
                vessel_d['error'] = form_vec(error_stream())  # 目标系里的偏差
                vessel_d['avel'] = form_vec(avel_stream())  # 地面系下角速度（等于目标系角速度
                vessel_d['vel'] = form_vec(vel_stream())  # 地面速度
                vessel_d['mass'] = mass_stream()
                conn.krpc.paused = True  # pause game
                vessel_d['acceleration'] = (vel - prev_vel) / game_delta_time
                v_data = vessel_profile1(vessel_d)
                print(v_data)
                solve_gfold(v_data)
                while conn.krpc.paused:
                    time.sleep(0.1)
                vessel.auto_pilot.engage()
                vessel.auto_pilot.reference_frame = ref_target

        if nav_mode == 'gfold':  # 跟随gfold路径
            gtf, gx, gu, gm, gs, gz = gfold_path  # g~: global naming issue
            # n_i = max(n_i - game_delta_time * 0.2 * N/tf, find_nearest_index(x, error))
            n_i = find_nearest_index(gx, error)
            (x_i, v_i, u_i) = sample_index(n_i)
            (x_i_, v_i_, u_i_) = sample_index(n_i + N / gtf)

            aerodynamic_force = np.array(aerodynamic_force_stream())
            # print(round(npl.norm(u_i_) * mass, 2), '\t', round(vessel.thrust), '\t', round(npl.norm(aerodynamic_force)))

            target_a = u_i + (v_i - vel) * k_v + (x_i - error) * k_x - aerodynamic_force / mass
            target_a_ = u_i_ + (v_i_ - vel) * k_v + (x_i - error) * k_x - aerodynamic_force / mass
            # print('Th:', u_i.round(2), target_a.round(2))

            if debug_lines:
                target_line.start = error
                target_line.end = (x_i[0], x_i[1], x_i[2])
                target_line.thickness = DEBUG_LINE_THICKNESS
                target2_line.start = error
                target2_line.end = (x_i_[0], x_i_[1], x_i_[2])
                target2_line.thickness = DEBUG_LINE_THICKNESS

            max_throttle_ctrl = throttle_limit_ctrl[1] * (max_thrust / mass)
            min_throttle_ctrl = throttle_limit_ctrl[0] * (max_thrust / mass)
            target_a = conic_clamp(target_a, min_throttle_ctrl, max_throttle_ctrl, max_tilt)
            target_a_ = conic_clamp(target_a_, min_throttle_ctrl, max_throttle_ctrl, max_tilt)
            if n_i < 0:
                target_a = np.array([g0, 0, 0]) + u_i
            target_direction = target_a_ / npl.norm(target_a_)
            target_throttle = npl.norm(target_a) / (max_thrust / mass)

            # target_u = conic_clamp(u_i, min_throttle_ctrl, max_throttle_ctrl, max_tilt)
            # throttle_u = npl.norm(target_u) / (max_thrust / mass)
            # print(v_i - vel)
            # print(target_throttle, '\t', throttle_u)
            # print('vx:', (v_i - vel).round(1), (x_i - error).round(1))

            if debug_lines:
                head_line.start = error
                head_line.end = error + target_direction * target_throttle * 20
                head_line.thickness = DEBUG_LINE_THICKNESS

            if n_i > 0:
                vessel.control.throttle = target_throttle

            if (N - n_i) * gtf / N < 10:
                vessel.control.gear = True
            if npl.norm(error[1:3]) < params['final_radius'] and npl.norm(error[0]) < params['final_height']:
                vessel.control.gear = True
                print('final')
                nav_mode = 'final'

        elif nav_mode == 'final':  # 最终降落阶段，直接pid
            max_acc = throttle_limit_ctrl[1] * (max_thrust / mass) - g0
            max_acc_low = throttle_limit_ctrl[1] * final_throttle * (max_thrust / mass) - g0
            est_h = error[0] - vel[0] ** 2 / (2 * max_acc)
            est_h_low = error[0] - vel[0] ** 2 / (2 * max_acc_low)
            est_h_center = (est_h + est_h_low) / 2
            vessel.control.throttle = np.clip(lerp(throttle_limit_ctrl[1] * final_throttle, throttle_limit_ctrl[1], -est_h_low / (est_h - est_h_low) * (1 + final_kp)), throttle_limit_ctrl[0],
                                              throttle_limit_ctrl[1])

            error_hor = form_v3(0, error[1], error[2])
            vel_hor = form_v3(0, vel[1], vel[2])
            ctrl_hor = -error_hor * 0.03 - vel_hor * 0.06
            target_direction = ctrl_hor + form_v3(1, 0, 0)
            target_direction /= npl.norm(target_direction)
            target_direction = conic_clamp(target_direction, 1, 1, max_tilt)
            
        else:
            target_direction = -vel

        # 变换到机体坐标系计算姿态控制，以下xyz均指机体系
        # target_direction_local = transform(target_direction, rotation_srf2local)  # 机体系的目标姿态的机体y轴指向
        # avel_local = transform(avel, rotation_srf2local)  # 机体系角速度
        # pid控制，roll直接消除角速度
        # control_pitch = -np.clip(ctrl_x_rot(angle_around_axis(target_direction_local, form_v3(0, 1, 0), form_v3(1, 0, 0)), game_delta_time), -1, 1)
        # control_yaw = -np.clip(ctrl_z_rot(angle_around_axis(target_direction_local, form_v3(0, 1, 0), form_v3(0, 0, 1)), game_delta_time), -1, 1)
        # control_roll = np.clip(avel_local[1] * ctrl_y_avel_kp, -1, 1)

        # 改用krpc内置PID控制器
        vessel.auto_pilot.target_direction = target_direction

        # vessel.control.pitch = control_pitch
        # vessel.control.yaw = control_yaw
        # vessel.control.roll = control_roll

        # 终止条件
        if (npl.norm(error[1:3]) < 3 and npl.norm(error[0]) < 1 and npl.norm(vel[1:3]) < 0.3 and npl.norm(vel[0]) < 0.5 and npl.norm(avel) < 0.2) or (vel[0] > 0 and npl.norm(error[0]) < 1):
            vessel.control.throttle = 0
            break

        prev_vel = vel
        game_prev_time = ut
