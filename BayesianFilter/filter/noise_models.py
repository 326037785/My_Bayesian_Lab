"""
过程噪声模型模块

提供共享的过程噪声协方差矩阵（Q矩阵）生成函数，
消除各滤波器中重复的 _get_process_noise_matrix 实现。

参考 filterpy 中 Q_discrete_white_noise() 的设计思路。
"""
import numpy as np


def get_process_noise_matrix(state_dim: int,
                             process_noise_std: float,
                             dt: float) -> np.ndarray:
    """获取过程噪声协方差矩阵 Q

    使用离散白噪声加速度模型（CV/CA模型），根据状态维度生成对应的Q矩阵。
    所有滤波器共享此函数，确保Q矩阵计算一致。

    Args:
        state_dim: 状态维度
            支持 4 (2D匀速: x,vx,y,vy)
            支持 6 (3D匀加速: x,vx,ax,y,vy,ay)
            其他维度使用默认白噪声模型
        process_noise_std: 过程噪声标准差
        dt: 时间步长

    Returns:
        过程噪声协方差矩阵 Q (state_dim x state_dim)
    """
    q = process_noise_std ** 2

    if state_dim == 4:
        # 2D 匀速运动模型 (x, vx, y, vy)
        return q * np.array([
            [dt**3/3, dt**2/2, 0, 0],
            [dt**2/2, dt, 0, 0],
            [0, 0, dt**3/3, dt**2/2],
            [0, 0, dt**2/2, dt]
        ])

    elif state_dim == 6:
        # 3D 匀加速运动模型 (x, vx, ax, y, vy, ay)
        dt2 = dt ** 2
        dt3 = dt ** 3
        dt4 = dt ** 4
        dt5 = dt ** 5
        return q * np.array([
            [dt5/20, dt4/8, dt3/6, 0, 0, 0],
            [dt4/8, dt3/3, dt2/2, 0, 0, 0],
            [dt3/6, dt2/2, dt, 0, 0, 0],
            [0, 0, 0, dt5/20, dt4/8, dt3/6],
            [0, 0, 0, dt4/8, dt3/3, dt2/2],
            [0, 0, 0, dt3/6, dt2/2, dt]
        ])

    else:
        # 默认：离散白噪声模型
        return q * dt * np.eye(state_dim)
