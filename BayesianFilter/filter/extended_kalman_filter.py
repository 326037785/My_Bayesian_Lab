"""
扩展卡尔曼滤波器 (EKF)
"""
import numpy as np
from typing import Optional, Callable
from .base_filter import BaseFilter
from .noise_models import get_process_noise_matrix


class ExtendedKalmanFilter(BaseFilter):
    """扩展卡尔曼滤波器
    
    适用于非线性系统
    
    状态转移: x(k) = f(x(k-1), u(k)) + w(k)
    观测方程: z(k) = h(x(k)) + v(k)
    
    其中：
    - f: 非线性状态转移函数
    - h: 非线性观测函数
    - w: 过程噪声，w ~ N(0, Q)
    - v: 观测噪声，v ~ N(0, R)
    
    EKF通过对非线性函数进行一阶泰勒展开（雅可比矩阵）来线性化
    """
    
    def __init__(self,
                 state_dim: int = 4,
                 measurement_dim: int = 2,
                 process_noise_std: float = 0.1,
                 measurement_noise_std: float = 1.0,
                 state_transition_func: Optional[Callable] = None,
                 measurement_func: Optional[Callable] = None,
                 state_transition_jacobian: Optional[Callable] = None,
                 measurement_jacobian: Optional[Callable] = None,
                 measurement_noise_matrix: Optional[np.ndarray] = None):
        """
        初始化扩展卡尔曼滤波器

        Args:
            state_dim: 状态维度
            measurement_dim: 观测维度
            process_noise_std: 过程噪声标准差
            measurement_noise_std: 观测噪声标准差（线性笛卡尔场景的简化设置）
            state_transition_func: 状态转移函数 f(x, dt)
            measurement_func: 观测函数 h(x)
            state_transition_jacobian: 状态转移雅可比矩阵 F(x, dt)
            measurement_jacobian: 观测雅可比矩阵 H(x)
            measurement_noise_matrix: 完整的观测噪声协方差矩阵 R
                当提供时，直接使用该矩阵（适用于极坐标等各向异性噪声）
                当为 None 时，回退到 diag(measurement_noise_std²)
        """
        super().__init__(state_dim, measurement_dim, process_noise_std)

        self.measurement_noise_std = measurement_noise_std

        # 设置函数
        self.f = state_transition_func if state_transition_func is not None else self._default_state_transition
        self.h = measurement_func if measurement_func is not None else self._default_measurement
        self.F_func = state_transition_jacobian if state_transition_jacobian is not None else self._default_F_jacobian
        self.H_func = measurement_jacobian if measurement_jacobian is not None else self._default_H_jacobian

        # 观测噪声协方差：优先使用完整矩阵（极坐标等各向异性场景），否则回退到各向同性默认值
        if measurement_noise_matrix is not None:
            self.R = np.asarray(measurement_noise_matrix, dtype=np.float64)
        else:
            self.R = np.eye(measurement_dim) * measurement_noise_std ** 2
    
    def _default_state_transition(self, state: np.ndarray, dt: float) -> np.ndarray:
        """默认状态转移函数（匀速模型）
        
        Args:
            state: 状态向量 [x, vx, y, vy]
            dt: 时间步长
            
        Returns:
            下一时刻状态
        """
        if self.state_dim == 4:
            x, vx, y, vy = state[0], state[1], state[2], state[3]
            return np.array([
                x + vx * dt,
                vx,
                y + vy * dt,
                vy
            ])
        elif self.state_dim == 6:
            x, vx, ax, y, vy, ay = state
            return np.array([
                x + vx * dt + 0.5 * ax * dt**2,
                vx + ax * dt,
                ax,
                y + vy * dt + 0.5 * ay * dt**2,
                vy + ay * dt,
                ay
            ])
        return state
    
    def _default_measurement(self, state: np.ndarray) -> np.ndarray:
        """默认观测函数（直接观测位置）
        
        Args:
            state: 状态向量
            
        Returns:
            观测向量 [x, y]
        """
        if self.state_dim >= 4:
            return np.array([state[0], state[2]])
        elif self.state_dim >= 2:
            return np.array([state[0], state[1]])
        return state[:self.measurement_dim]
    
    def _default_F_jacobian(self, state: np.ndarray, dt: float) -> np.ndarray:
        """默认状态转移雅可比矩阵
        
        Args:
            state: 状态向量
            dt: 时间步长
            
        Returns:
            雅可比矩阵 F
        """
        if self.state_dim == 4:
            F = np.array([
                [1, dt, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, dt],
                [0, 0, 0, 1]
            ])
        elif self.state_dim == 6:
            F = np.array([
                [1, dt, dt**2/2, 0, 0, 0],
                [0, 1, dt, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, 1, dt, dt**2/2],
                [0, 0, 0, 0, 1, dt],
                [0, 0, 0, 0, 0, 1]
            ])
        else:
            F = np.eye(self.state_dim)
        
        return F
    
    def _default_H_jacobian(self, state: np.ndarray) -> np.ndarray:
        """默认观测雅可比矩阵
        
        Args:
            state: 状态向量
            
        Returns:
            雅可比矩阵 H
        """
        if self.state_dim == 4:
            H = np.array([
                [1, 0, 0, 0],
                [0, 0, 1, 0]
            ])
        elif self.state_dim == 6:
            H = np.array([
                [1, 0, 0, 0, 0, 0],
                [0, 0, 0, 1, 0, 0]
            ])
        else:
            H = np.zeros((self.measurement_dim, self.state_dim))
            H[0, 0] = 1
            if self.measurement_dim > 1 and self.state_dim > 2:
                H[1, 2] = 1
        
        return H
    
    def _get_state_transition_matrix(self, dt: float) -> np.ndarray:
        """获取状态转移矩阵（雅可比矩阵）"""
        return self.F_func(self.state, dt)
    
    def _get_process_noise_matrix(self, dt: float) -> np.ndarray:
        """获取过程噪声协方差矩阵"""
        return get_process_noise_matrix(self.state_dim, self.process_noise_std, dt)
    
    def _get_measurement_matrix(self) -> np.ndarray:
        """获取观测矩阵（雅可比矩阵）"""
        return self.H_func(self.state)
    
    def predict(self, dt: float) -> None:
        """状态预测
        
        使用非线性状态转移函数进行预测
        
        Args:
            dt: 时间步长
        """
        if not self.initialized:
            raise RuntimeError("Filter not initialized")
        
        # 获取雅可比矩阵
        F = self.F_func(self.state, dt)
        Q = self._get_process_noise_matrix(dt)
        
        # 使用非线性函数进行状态预测
        self.state = self.f(self.state, dt)
        
        # 协方差预测
        self.covariance = F @ self.covariance @ F.T + Q
        
        # 确保协方差正定
        self.covariance = self._ensure_positive_definite(self.covariance)
    
    def update(self, measurement: np.ndarray,
               measurement_covariance: Optional[np.ndarray] = None) -> None:
        """量测更新
        
        使用非线性观测函数进行更新
        
        Args:
            measurement: 观测向量
            measurement_covariance: 观测噪声协方差（可选）
        """
        if not self.initialized:
            raise RuntimeError("Filter not initialized")
        
        # 使用提供的观测噪声协方差或默认值
        R = measurement_covariance if measurement_covariance is not None else self.R
        
        # 获取观测雅可比矩阵
        H = self.H_func(self.state)
        
        # 计算预测观测
        z_pred = self.h(self.state)
        
        # 计算新息
        z = np.asarray(measurement, dtype=np.float64)
        y = z - z_pred
        
        # 计算新息协方差
        S = H @ self.covariance @ H.T + R
        
        # 计算卡尔曼增益
        K = self.covariance @ H.T @ np.linalg.inv(S)
        
        # 状态更新
        self.state = self.state + K @ y
        
        # 协方差更新（使用Joseph形式）
        I_KH = np.eye(self.state_dim) - K @ H
        self.covariance = I_KH @ self.covariance @ I_KH.T + K @ R @ K.T
        
        # 确保协方差正定
        self.covariance = self._ensure_positive_definite(self.covariance)
        
        # 验证状态
        if not self._validate_state():
            print("Warning: Invalid state after update")
    
    def predict_and_update(self, dt: float, measurement: np.ndarray,
                           measurement_covariance: Optional[np.ndarray] = None) -> None:
        """预测并更新
        
        Args:
            dt: 时间步长
            measurement: 观测向量
            measurement_covariance: 观测噪声协方差（可选）
        """
        self.predict(dt)
        self.update(measurement, measurement_covariance)
    
    def set_state_transition_function(self, func: Callable, 
                                       jacobian: Optional[Callable] = None) -> None:
        """设置状态转移函数
        
        Args:
            func: 状态转移函数 f(x, dt)
            jacobian: 雅可比矩阵函数 F(x, dt)
        """
        self.f = func
        if jacobian is not None:
            self.F_func = jacobian
    
    def set_measurement_function(self, func: Callable,
                                  jacobian: Optional[Callable] = None) -> None:
        """设置观测函数

        Args:
            func: 观测函数 h(x)
            jacobian: 雅可比矩阵函数 H(x)
        """
        self.h = func
        if jacobian is not None:
            self.H_func = jacobian

    def set_measurement_noise(self, R: np.ndarray) -> None:
        """设置观测噪声协方差矩阵

        Args:
            R: 观测噪声协方差矩阵（必须与观测维度匹配）
        """
        self.R = np.asarray(R, dtype=np.float64)
    
    def get_predicted_measurement(self) -> np.ndarray:
        """获取预测观测
        
        Returns:
            预测观测向量
        """
        return self.h(self.state)
    
    def get_innovation(self, measurement: np.ndarray) -> np.ndarray:
        """计算新息
        
        Args:
            measurement: 观测向量
            
        Returns:
            新息向量
        """
        z = np.asarray(measurement, dtype=np.float64)
        z_pred = self.h(self.state)
        return z - z_pred
    
    def get_innovation_covariance(self, measurement_covariance: Optional[np.ndarray] = None) -> np.ndarray:
        """获取新息协方差矩阵
        
        Args:
            measurement_covariance: 观测噪声协方差（可选）
            
        Returns:
            新息协方差矩阵 S
        """
        H = self.H_func(self.state)
        R = measurement_covariance if measurement_covariance is not None else self.R
        return H @ self.covariance @ H.T + R
    
    def get_likelihood(self, measurement: np.ndarray,
                       measurement_covariance: Optional[np.ndarray] = None) -> float:
        """计算观测似然
        
        Args:
            measurement: 观测向量
            measurement_covariance: 观测噪声协方差（可选）
            
        Returns:
            似然值
        """
        innovation = self.get_innovation(measurement)
        S = self.get_innovation_covariance(measurement_covariance)
        
        n = len(innovation)
        S_det = np.linalg.det(S)
        S_inv = np.linalg.inv(S)
        
        exponent = -0.5 * innovation.T @ S_inv @ innovation
        coefficient = 1.0 / np.sqrt((2 * np.pi) ** n * S_det)
        
        return coefficient * np.exp(exponent)


# 一些常用的非线性观测模型

def polar_measurement_function(state: np.ndarray) -> np.ndarray:
    """极坐标观测函数
    
    状态: [x, vx, y, vy]
    观测: [range, bearing]
    
    Args:
        state: 状态向量
        
    Returns:
        极坐标观测 [range, bearing]
    """
    x, y = state[0], state[2]
    range_meas = np.sqrt(x**2 + y**2)
    bearing_meas = np.arctan2(y, x)
    return np.array([range_meas, bearing_meas])


def polar_measurement_jacobian(state: np.ndarray) -> np.ndarray:
    """极坐标观测雅可比矩阵
    
    Args:
        state: 状态向量
        
    Returns:
        雅可比矩阵 H
    """
    x, y = state[0], state[2]
    r = np.sqrt(x**2 + y**2)
    
    if r < 1e-6:
        return np.zeros((2, len(state)))
    
    H = np.zeros((2, len(state)))
    H[0, 0] = x / r  # d(range)/dx
    H[0, 2] = y / r  # d(range)/dy
    H[1, 0] = -y / r**2  # d(bearing)/dx
    H[1, 2] = x / r**2  # d(bearing)/dy
    
    return H


def range_only_measurement_function(state: np.ndarray) -> np.ndarray:
    """仅距离观测函数
    
    Args:
        state: 状态向量
        
    Returns:
        距离观测 [range]
    """
    x, y = state[0], state[2]
    return np.array([np.sqrt(x**2 + y**2)])


def range_only_measurement_jacobian(state: np.ndarray) -> np.ndarray:
    """仅距离观测雅可比矩阵
    
    Args:
        state: 状态向量
        
    Returns:
        雅可比矩阵 H
    """
    x, y = state[0], state[2]
    r = np.sqrt(x**2 + y**2)
    
    if r < 1e-6:
        return np.zeros((1, len(state)))
    
    H = np.zeros((1, len(state)))
    H[0, 0] = x / r
    H[0, 2] = y / r
    
    return H
