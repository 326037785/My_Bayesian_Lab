"""
容积卡尔曼滤波器 (CKF)
"""
import numpy as np
from typing import Optional, Callable, Tuple
from .base_filter import BaseFilter
from .noise_models import get_process_noise_matrix


class CubatureKalmanFilter(BaseFilter):
    """容积卡尔曼滤波器
    
    使用三阶球面-径向容积规则来近似非线性变换
    
    与UKF类似，但使用不同的Sigma点生成方法：
    - UKF: 使用缩放的Sigma点
    - CKF: 使用容积点（更规则的分布）
    
    优点：
    - 数值稳定性更好
    - 计算效率更高
    - 理论基础更坚实
    """
    
    def __init__(self,
                 state_dim: int = 4,
                 measurement_dim: int = 2,
                 process_noise_std: float = 0.1,
                 measurement_noise_std: float = 1.0,
                 state_transition_func: Optional[Callable] = None,
                 measurement_func: Optional[Callable] = None,
                 measurement_noise_matrix: Optional[np.ndarray] = None):
        """
        初始化容积卡尔曼滤波器

        Args:
            state_dim: 状态维度
            measurement_dim: 观测维度
            process_noise_std: 过程噪声标准差
            measurement_noise_std: 观测噪声标准差（线性笛卡尔场景的简化设置）
            state_transition_func: 状态转移函数 f(x, dt)
            measurement_func: 观测函数 h(x)
            measurement_noise_matrix: 完整的观测噪声协方差矩阵 R
                当提供时，直接使用该矩阵（适用于极坐标等各向异性噪声）
                当为 None 时，回退到 diag(measurement_noise_std²)
        """
        super().__init__(state_dim, measurement_dim, process_noise_std)

        self.measurement_noise_std = measurement_noise_std

        # 设置函数
        self.f = state_transition_func if state_transition_func is not None else self._default_state_transition
        self.h = measurement_func if measurement_func is not None else self._default_measurement

        # 观测噪声协方差：优先使用完整矩阵（极坐标等各向异性场景），否则回退到各向同性默认值
        if measurement_noise_matrix is not None:
            self.R = np.asarray(measurement_noise_matrix, dtype=np.float64)
        else:
            self.R = np.eye(measurement_dim) * measurement_noise_std ** 2
        
        # 计算容积点权重
        self.n_cubature = 2 * state_dim
        self.weight = 1.0 / self.n_cubature
        
        # 生成基础容积点
        self._generate_base_cubature_points()
    
    def _generate_base_cubature_points(self) -> None:
        """生成基础容积点
        
        容积点是单位向量的正负对
        """
        n = self.state_dim
        base_points = np.zeros((2 * n, n))
        
        for i in range(n):
            base_points[i, i] = 1.0
            base_points[n + i, i] = -1.0
        
        self.base_cubature_points = base_points
    
    def _default_state_transition(self, state: np.ndarray, dt: float) -> np.ndarray:
        """默认状态转移函数（匀速模型）"""
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
        """默认观测函数（直接观测位置）"""
        if self.state_dim >= 4:
            return np.array([state[0], state[2]])
        elif self.state_dim >= 2:
            return np.array([state[0], state[1]])
        return state[:self.measurement_dim]
    
    def _get_state_transition_matrix(self, dt: float) -> np.ndarray:
        """获取状态转移矩阵（用于基类接口）"""
        return np.eye(self.state_dim)
    
    def _get_process_noise_matrix(self, dt: float) -> np.ndarray:
        """获取过程噪声协方差矩阵"""
        return get_process_noise_matrix(self.state_dim, self.process_noise_std, dt)
    
    def _get_measurement_matrix(self) -> np.ndarray:
        """获取观测矩阵（用于基类接口）"""
        return np.eye(self.measurement_dim, self.state_dim)
    
    def _generate_cubature_points(self, mean: np.ndarray, 
                                   covariance: np.ndarray) -> np.ndarray:
        """生成容积点
        
        Args:
            mean: 均值向量
            covariance: 协方差矩阵
            
        Returns:
            容积点矩阵，形状为 (2*n, n)
        """
        n = len(mean)
        
        # 确保协方差正定
        covariance = self._ensure_positive_definite(covariance)
        
        # 计算矩阵平方根
        try:
            # Cholesky分解
            L = np.linalg.cholesky(covariance)
        except np.linalg.LinAlgError:
            # 如果Cholesky分解失败，使用特征值分解
            eigenvalues, eigenvectors = np.linalg.eigh(covariance)
            eigenvalues = np.maximum(eigenvalues, 1e-6)
            L = eigenvectors @ np.diag(np.sqrt(eigenvalues))
        
        # 生成容积点（向量化）
        scaled_points = np.sqrt(n) * L @ self.base_cubature_points.T
        cubature_points = mean + scaled_points.T
        
        return cubature_points
    
    def _cubature_transform(self, cubature_points: np.ndarray,
                            func: Callable, *args) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """容积变换
        
        Args:
            cubature_points: 容积点矩阵
            func: 非线性函数
            *args: 函数参数
            
        Returns:
            变换后的均值、协方差和变换后的容积点
        """
        n_cubature = cubature_points.shape[0]
        n_out = func(cubature_points[0], *args).shape[0]
        
        # 变换容积点（需要逐点应用非线性函数）
        transformed_points = np.zeros((n_cubature, n_out))
        for i in range(n_cubature):
            transformed_points[i] = func(cubature_points[i], *args)
        
        # 向量化计算均值
        mean = self.weight * np.sum(transformed_points, axis=0)
        
        # 向量化计算协方差
        diff = transformed_points - mean
        covariance = self.weight * (diff.T @ diff)
        
        return mean, covariance, transformed_points
    
    def predict(self, dt: float) -> None:
        """状态预测
        
        Args:
            dt: 时间步长
        """
        if not self.initialized:
            raise RuntimeError("Filter not initialized")
        
        # 生成容积点
        cubature_points = self._generate_cubature_points(self.state, self.covariance)
        
        # 通过状态转移函数传播容积点
        n_cubature = cubature_points.shape[0]
        propagated_points = np.zeros_like(cubature_points)
        for i in range(n_cubature):
            propagated_points[i] = self.f(cubature_points[i], dt)
        
        # 计算预测均值和协方差
        predicted_mean = np.zeros(self.state_dim)
        for i in range(n_cubature):
            predicted_mean += self.weight * propagated_points[i]
        
        predicted_cov = np.zeros((self.state_dim, self.state_dim))
        for i in range(n_cubature):
            diff = propagated_points[i] - predicted_mean
            predicted_cov += self.weight * np.outer(diff, diff)
        
        # 添加过程噪声
        Q = self._get_process_noise_matrix(dt)
        predicted_cov += Q
        
        # 更新状态
        self.state = predicted_mean
        self.covariance = self._ensure_positive_definite(predicted_cov)
    
    def update(self, measurement: np.ndarray,
               measurement_covariance: Optional[np.ndarray] = None) -> None:
        """量测更新
        
        Args:
            measurement: 观测向量
            measurement_covariance: 观测噪声协方差（可选）
        """
        if not self.initialized:
            raise RuntimeError("Filter not initialized")
        
        # 使用提供的观测噪声协方差或默认值
        R = measurement_covariance if measurement_covariance is not None else self.R
        
        # 生成容积点
        cubature_points = self._generate_cubature_points(self.state, self.covariance)
        
        # 通过观测函数传播容积点
        n_cubature = cubature_points.shape[0]
        measurement_points = np.zeros((n_cubature, self.measurement_dim))
        for i in range(n_cubature):
            measurement_points[i] = self.h(cubature_points[i])
        
        # 计算预测观测均值
        z_pred = np.zeros(self.measurement_dim)
        for i in range(n_cubature):
            z_pred += self.weight * measurement_points[i]
        
        # 计算协方差 Pxz 和 Pzz
        Pxz = np.zeros((self.state_dim, self.measurement_dim))
        Pzz = np.zeros((self.measurement_dim, self.measurement_dim))
        
        for i in range(n_cubature):
            dx = cubature_points[i] - self.state
            dz = measurement_points[i] - z_pred
            Pxz += self.weight * np.outer(dx, dz)
            Pzz += self.weight * np.outer(dz, dz)
        
        # 添加观测噪声
        Pzz += R
        
        # 计算卡尔曼增益
        K = Pxz @ np.linalg.inv(Pzz)
        
        # 计算新息
        z = np.asarray(measurement, dtype=np.float64)
        innovation = z - z_pred
        
        # 状态更新
        self.state = self.state + K @ innovation
        
        # 协方差更新
        self.covariance = self.covariance - K @ Pzz @ K.T
        
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
    
    def set_state_transition_function(self, func: Callable) -> None:
        """设置状态转移函数
        
        Args:
            func: 状态转移函数 f(x, dt)
        """
        self.f = func
    
    def set_measurement_function(self, func: Callable) -> None:
        """设置观测函数

        Args:
            func: 观测函数 h(x)
        """
        self.h = func

    def set_measurement_noise(self, R: np.ndarray) -> None:
        """设置观测噪声协方差矩阵

        Args:
            R: 观测噪声协方差矩阵（必须与观测维度匹配）
        """
        self.R = np.asarray(R, dtype=np.float64)

    def get_cubature_points(self) -> np.ndarray:
        """获取当前容积点
        
        Returns:
            容积点矩阵
        """
        return self._generate_cubature_points(self.state, self.covariance)
    
    def get_predicted_measurement(self) -> Tuple[np.ndarray, np.ndarray]:
        """获取预测观测和协方差
        
        Returns:
            (预测观测均值, 预测观测协方差)
        """
        cubature_points = self._generate_cubature_points(self.state, self.covariance)
        n_cubature = cubature_points.shape[0]
        
        measurement_points = np.zeros((n_cubature, self.measurement_dim))
        for i in range(n_cubature):
            measurement_points[i] = self.h(cubature_points[i])
        
        z_pred = np.zeros(self.measurement_dim)
        for i in range(n_cubature):
            z_pred += self.weight * measurement_points[i]
        
        Pzz = np.zeros((self.measurement_dim, self.measurement_dim))
        for i in range(n_cubature):
            dz = measurement_points[i] - z_pred
            Pzz += self.weight * np.outer(dz, dz)
        
        Pzz += self.R
        
        return z_pred, Pzz
    
    def get_innovation(self, measurement: np.ndarray) -> np.ndarray:
        """计算新息
        
        Args:
            measurement: 观测向量
            
        Returns:
            新息向量
        """
        z_pred, _ = self.get_predicted_measurement()
        z = np.asarray(measurement, dtype=np.float64)
        return z - z_pred
    
    def get_innovation_covariance(self, measurement_covariance: Optional[np.ndarray] = None) -> np.ndarray:
        """获取新息协方差矩阵
        
        Args:
            measurement_covariance: 观测噪声协方差（可选）
            
        Returns:
            新息协方差矩阵
        """
        _, Pzz = self.get_predicted_measurement()
        
        if measurement_covariance is not None:
            # 重新计算，使用提供的协方差
            cubature_points = self._generate_cubature_points(self.state, self.covariance)
            n_cubature = cubature_points.shape[0]
            
            measurement_points = np.zeros((n_cubature, self.measurement_dim))
            for i in range(n_cubature):
                measurement_points[i] = self.h(cubature_points[i])
            
            z_pred = np.zeros(self.measurement_dim)
            for i in range(n_cubature):
                z_pred += self.weight * measurement_points[i]
            
            Pzz = np.zeros((self.measurement_dim, self.measurement_dim))
            for i in range(n_cubature):
                dz = measurement_points[i] - z_pred
                Pzz += self.weight * np.outer(dz, dz)
            
            Pzz += measurement_covariance
        
        return Pzz
    
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


# 一些常用的非线性函数（可以与UKF共享）

def polar_measurement_function_ckf(state: np.ndarray) -> np.ndarray:
    """极坐标观测函数（用于CKF）
    
    Args:
        state: 状态向量 [x, vx, y, vy]
        
    Returns:
        极坐标观测 [range, bearing]
    """
    x, y = state[0], state[2]
    range_meas = np.sqrt(x**2 + y**2)
    bearing_meas = np.arctan2(y, x)
    return np.array([range_meas, bearing_meas])


def bearing_only_measurement_function(state: np.ndarray) -> np.ndarray:
    """仅角度观测函数
    
    Args:
        state: 状态向量 [x, vx, y, vy]
        
    Returns:
        角度观测 [bearing]
    """
    x, y = state[0], state[2]
    bearing = np.arctan2(y, x)
    return np.array([bearing])
