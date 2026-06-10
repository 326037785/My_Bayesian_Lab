"""
卡尔曼滤波器 (KF)
"""
import numpy as np
from typing import Optional
from .base_filter import BaseFilter
from .noise_models import get_process_noise_matrix


class KalmanFilter(BaseFilter):
    """卡尔曼滤波器
    
    适用于线性高斯系统
    
    状态转移: x(k) = F * x(k-1) + w(k)
    观测方程: z(k) = H * x(k) + v(k)
    
    其中：
    - F: 状态转移矩阵
    - H: 观测矩阵
    - w: 过程噪声，w ~ N(0, Q)
    - v: 观测噪声，v ~ N(0, R)
    """
    
    def __init__(self,
                 state_dim: int = 4,
                 measurement_dim: int = 2,
                 process_noise_std: float = 0.1,
                 measurement_noise_std: float = 1.0):
        """
        初始化卡尔曼滤波器
        
        Args:
            state_dim: 状态维度，默认4 [x, vx, y, vy]
            measurement_dim: 观测维度，默认2 [x, y]
            process_noise_std: 过程噪声标准差
            measurement_noise_std: 观测噪声标准差
        """
        super().__init__(state_dim, measurement_dim, process_noise_std)
        
        self.measurement_noise_std = measurement_noise_std
        
        # 默认观测噪声协方差
        self.R = np.eye(measurement_dim) * measurement_noise_std ** 2
        
        # 状态转移矩阵（匀速模型）
        self.F = self._get_state_transition_matrix(1.0)
        
        # 观测矩阵
        self.H = self._get_measurement_matrix()
    
    def _get_state_transition_matrix(self, dt: float) -> np.ndarray:
        """获取状态转移矩阵
        
        匀速运动模型:
        F = [[1, dt, 0,  0],
             [0,  1, 0,  0],
             [0,  0, 1, dt],
             [0,  0, 0,  1]]
        """
        if self.state_dim == 4:
            F = np.array([
                [1, dt, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, dt],
                [0, 0, 0, 1]
            ])
        elif self.state_dim == 6:
            # 匀加速模型 [x, vx, ax, y, vy, ay]
            F = np.array([
                [1, dt, dt**2/2, 0, 0, 0],
                [0, 1, dt, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, 1, dt, dt**2/2],
                [0, 0, 0, 0, 1, dt],
                [0, 0, 0, 0, 0, 1]
            ])
        else:
            # 默认使用单位矩阵
            F = np.eye(self.state_dim)
        
        return F
    
    def _get_process_noise_matrix(self, dt: float) -> np.ndarray:
        """获取过程噪声协方差矩阵
        
        使用共享的 get_process_noise_matrix 函数
        """
        return get_process_noise_matrix(self.state_dim, self.process_noise_std, dt)
    
    def _get_measurement_matrix(self) -> np.ndarray:
        """获取观测矩阵
        
        直接观测位置 [x, y]
        H = [[1, 0, 0, 0],
             [0, 0, 1, 0]]
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
            # 默认只观测前两个状态
            H = np.zeros((self.measurement_dim, self.state_dim))
            H[0, 0] = 1
            if self.measurement_dim > 1:
                H[1, 2] = 1 if self.state_dim > 2 else 1
        
        return H
    
    def predict(self, dt: float) -> None:
        """状态预测
        
        预测步骤：
        1. 状态预测: x_pred = F * x
        2. 协方差预测: P_pred = F * P * F' + Q
        
        Args:
            dt: 时间步长
        """
        if not self.initialized:
            raise RuntimeError("Filter not initialized")
        
        # 更新状态转移矩阵和过程噪声
        self.F = self._get_state_transition_matrix(dt)
        Q = self._get_process_noise_matrix(dt)
        
        # 状态预测
        self.state = self.F @ self.state
        
        # 协方差预测
        self.covariance = self.F @ self.covariance @ self.F.T + Q
        
        # 确保协方差正定
        self.covariance = self._ensure_positive_definite(self.covariance)
    
    def update(self, measurement: np.ndarray,
               measurement_covariance: Optional[np.ndarray] = None) -> None:
        """量测更新
        
        更新步骤：
       1. 计算卡尔曼增益: K = P * H' * (H * P * H' + R)^(-1)
        2. 状态更新: x = x + K * (z - H * x)
        3. 协方差更新: P = (I - K * H) * P
        
        Args:
            measurement: 观测向量
            measurement_covariance: 观测噪声协方差（可选）
        """
        if not self.initialized:
            raise RuntimeError("Filter not initialized")
        
        # 使用提供的观测噪声协方差或默认值
        R = measurement_covariance if measurement_covariance is not None else self.R
        
        # 计算新息（残差）
        z = np.asarray(measurement, dtype=np.float64)
        y = z - self.H @ self.state
        
        # 计算新息协方差
        S = self.H @ self.covariance @ self.H.T + R
        
        # 计算卡尔曼增益
        K = self.covariance @ self.H.T @ np.linalg.inv(S)
        
        # 状态更新
        self.state = self.state + K @ y
        
        # 协方差更新（使用Joseph形式保证正定性）
        I_KH = np.eye(self.state_dim) - K @ self.H
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
    
    def get_innovation(self, measurement: np.ndarray) -> np.ndarray:
        """计算新息（残差）
        
        Args:
            measurement: 观测向量
            
        Returns:
            新息向量
        """
        z = np.asarray(measurement, dtype=np.float64)
        return z - self.H @ self.state
    
    def get_innovation_covariance(self) -> np.ndarray:
        """获取新息协方差矩阵
        
        Returns:
            新息协方差矩阵 S
        """
        return self.H @ self.covariance @ self.H.T + self.R
    
    def get_likelihood(self, measurement: np.ndarray) -> float:
        """计算观测似然
        
        Args:
            measurement: 观测向量
            
        Returns:
            似然值
        """
        innovation = self.get_innovation(measurement)
        S = self.get_innovation_covariance()
        
        # 计算多元正态分布概率密度
        n = len(innovation)
        S_det = np.linalg.det(S)
        S_inv = np.linalg.inv(S)
        
        exponent = -0.5 * innovation.T @ S_inv @ innovation
        coefficient = 1.0 / np.sqrt((2 * np.pi) ** n * S_det)
        
        return coefficient * np.exp(exponent)
    
    def get_mahalanobis_distance(self, measurement: np.ndarray) -> float:
        """计算马氏距离
        
        Args:
            measurement: 观测向量
            
        Returns:
            马氏距离
        """
        innovation = self.get_innovation(measurement)
        S = self.get_innovation_covariance()
        S_inv = np.linalg.inv(S)
        
        return np.sqrt(innovation.T @ S_inv @ innovation)
    
    def set_measurement_noise(self, R: np.ndarray) -> None:
        """设置观测噪声协方差
        
        Args:
            R: 观测噪声协方差矩阵
        """
        self.R = np.asarray(R, dtype=np.float64)
    
    def set_process_noise_std(self, std: float) -> None:
        """设置过程噪声标准差
        
        Args:
            std: 过程噪声标准差
        """
        self.process_noise_std = std
