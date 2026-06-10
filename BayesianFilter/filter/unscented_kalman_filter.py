"""
无迹卡尔曼滤波器 (UKF)
"""
import numpy as np
from typing import Optional, Callable, Tuple
from .base_filter import BaseFilter
from .noise_models import get_process_noise_matrix


class UnscentedKalmanFilter(BaseFilter):
    """无迹卡尔曼滤波器
    
    使用Sigma点来捕获非线性变换后的均值和协方差
    
    优点：
    - 无需计算雅可比矩阵
    - 能够捕获高阶统计信息
    - 比EKF更准确
    
    参数：
    - alpha: 控制Sigma点的分布
    - beta: 包含先验知识（高斯分布时beta=2）
    - kappa: 二级缩放参数
    """
    
    def __init__(self,
                 state_dim: int = 4,
                 measurement_dim: int = 2,
                 process_noise_std: float = 0.1,
                 measurement_noise_std: float = 1.0,
                 state_transition_func: Optional[Callable] = None,
                 measurement_func: Optional[Callable] = None,
                 alpha: float = 1e-3,
                 beta: float = 2.0,
                 kappa: float = 0.0,
                 measurement_noise_matrix: Optional[np.ndarray] = None):
        """
        初始化无迹卡尔曼滤波器

        Args:
            state_dim: 状态维度
            measurement_dim: 观测维度
            process_noise_std: 过程噪声标准差
            measurement_noise_std: 观测噪声标准差（线性笛卡尔场景的简化设置）
            state_transition_func: 状态转移函数 f(x, dt)
            measurement_func: 观测函数 h(x)
            alpha: 控制Sigma点分布的参数
            beta: 包含先验知识的参数（高斯分布时beta=2）
            kappa: 二级缩放参数
            measurement_noise_matrix: 完整的观测噪声协方差矩阵 R
                当提供时，直接使用该矩阵（适用于极坐标等各向异性噪声）
                当为 None 时，回退到 diag(measurement_noise_std²)
        """
        super().__init__(state_dim, measurement_dim, process_noise_std)

        self.measurement_noise_std = measurement_noise_std

        # 设置函数
        self.f = state_transition_func if state_transition_func is not None else self._default_state_transition
        self.h = measurement_func if measurement_func is not None else self._default_measurement

        # UKF参数
        self.alpha = alpha
        self.beta = beta
        self.kappa = kappa

        # 计算缩放参数
        self.lambda_ = alpha ** 2 * (state_dim + kappa) - state_dim

        # 计算权重
        self._calculate_weights()

        # 观测噪声协方差：优先使用完整矩阵（极坐标等各向异性场景），否则回退到各向同性默认值
        if measurement_noise_matrix is not None:
            self.R = np.asarray(measurement_noise_matrix, dtype=np.float64)
        else:
            self.R = np.eye(measurement_dim) * measurement_noise_std ** 2
    
    def _calculate_weights(self) -> None:
        """计算Sigma点权重"""
        n = self.state_dim
        lambda_ = self.lambda_
        
        # 均值权重
        self.Wm = np.zeros(2 * n + 1)
        self.Wm[0] = lambda_ / (n + lambda_)
        self.Wm[1:] = 1.0 / (2 * (n + lambda_))
        
        # 协方差权重
        self.Wc = np.zeros(2 * n + 1)
        self.Wc[0] = lambda_ / (n + lambda_) + (1 - self.alpha ** 2 + self.beta)
        self.Wc[1:] = 1.0 / (2 * (n + lambda_))
    
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
        # UKF不使用显式的状态转移矩阵
        # 返回单位矩阵作为占位符
        return np.eye(self.state_dim)
    
    def _get_process_noise_matrix(self, dt: float) -> np.ndarray:
        """获取过程噪声协方差矩阵"""
        return get_process_noise_matrix(self.state_dim, self.process_noise_std, dt)
    
    def _get_measurement_matrix(self) -> np.ndarray:
        """获取观测矩阵（用于基类接口）"""
        # UKF不使用显式的观测矩阵
        # 返回单位矩阵作为占位符
        return np.eye(self.measurement_dim, self.state_dim)
    
    def _generate_sigma_points(self, mean: np.ndarray, 
                                covariance: np.ndarray) -> np.ndarray:
        """生成Sigma点
        
        参考MATLAB实现：
        Psqrtm = chol((n_x+lambda)*P)';  % 下三角矩阵
        X = repmat(m,[1 2*n_x+1]) + [zeros(n_x,1) -Psqrtm Psqrtm];
        
        Args:
            mean: 均值向量
            covariance: 协方差矩阵
            
        Returns:
            Sigma点矩阵，形状为 (2*n+1, n)
        """
        n = len(mean)
        lambda_ = self.lambda_
        
        # 确保协方差正定
        covariance = self._ensure_positive_definite(covariance)
        
        # 计算矩阵平方根（参考MATLAB: Psqrtm = chol((n_x+lambda)*P)'）
        # np.linalg.cholesky 返回下三角矩阵 L，使得 P = L @ L.T
        try:
            L = np.linalg.cholesky((n + lambda_) * covariance)
        except np.linalg.LinAlgError:
            # 如果Cholesky分解失败，使用特征值分解
            eigenvalues, eigenvectors = np.linalg.eigh(covariance)
            eigenvalues = np.maximum(eigenvalues, 1e-6)
            L = eigenvectors @ np.diag(np.sqrt(eigenvalues * (n + lambda_)))
        
        # 生成Sigma点（向量化）
        # 参考MATLAB: X = [zeros(n_x,1) -Psqrtm Psqrtm]
        # L 是下三角矩阵，每列对应一个偏移方向
        sigma_points = np.zeros((2 * n + 1, n))
        sigma_points[0] = mean
        sigma_points[1:n+1] = mean + L.T  # L.T 的每行是一个偏移方向
        sigma_points[n+1:] = mean - L.T
        
        return sigma_points
    
    def _unscented_transform(self, sigma_points: np.ndarray, 
                              func: Callable, *args) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """无迹变换
        
        Args:
            sigma_points: Sigma点矩阵
            func: 非线性函数
            *args: 函数参数
            
        Returns:
            变换后的均值、协方差和变换后的Sigma点
        """
        n_sigma = sigma_points.shape[0]
        n_out = func(sigma_points[0], *args).shape[0]
        
        # 变换Sigma点（需要逐点应用非线性函数）
        transformed_points = np.zeros((n_sigma, n_out))
        for i in range(n_sigma):
            transformed_points[i] = func(sigma_points[i], *args)
        
        # 向量化计算均值
        mean = np.dot(self.Wm, transformed_points)
        
        # 向量化计算协方差
        diff = transformed_points - mean
        covariance = diff.T @ (self.Wc[:, np.newaxis] * diff)
        
        return mean, covariance, transformed_points
    
    def predict(self, dt: float) -> None:
        """状态预测
        
        Args:
            dt: 时间步长
        """
        if not self.initialized:
            raise RuntimeError("Filter not initialized")
        
        # 生成Sigma点
        sigma_points = self._generate_sigma_points(self.state, self.covariance)
        
        # 通过状态转移函数传播Sigma点
        n_sigma = sigma_points.shape[0]
        propagated_points = np.zeros_like(sigma_points)
        for i in range(n_sigma):
            propagated_points[i] = self.f(sigma_points[i], dt)
        
        # 计算预测均值和协方差
        predicted_mean = np.zeros(self.state_dim)
        for i in range(n_sigma):
            predicted_mean += self.Wm[i] * propagated_points[i]
        
        predicted_cov = np.zeros((self.state_dim, self.state_dim))
        for i in range(n_sigma):
            diff = propagated_points[i] - predicted_mean
            predicted_cov += self.Wc[i] * np.outer(diff, diff)
        
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
        
        # 生成Sigma点
        sigma_points = self._generate_sigma_points(self.state, self.covariance)
        
        # 通过观测函数传播Sigma点
        n_sigma = sigma_points.shape[0]
        measurement_points = np.zeros((n_sigma, self.measurement_dim))
        for i in range(n_sigma):
            measurement_points[i] = self.h(sigma_points[i])
        
        # 计算预测观测均值
        z_pred = np.zeros(self.measurement_dim)
        for i in range(n_sigma):
            z_pred += self.Wm[i] * measurement_points[i]
        
        # 计算协方差 Pxz 和 Pzz
        Pxz = np.zeros((self.state_dim, self.measurement_dim))
        Pzz = np.zeros((self.measurement_dim, self.measurement_dim))
        
        for i in range(n_sigma):
            dx = sigma_points[i] - self.state
            dz = measurement_points[i] - z_pred
            Pxz += self.Wc[i] * np.outer(dx, dz)
            Pzz += self.Wc[i] * np.outer(dz, dz)
        
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
    
    def get_sigma_points(self) -> np.ndarray:
        """获取当前Sigma点
        
        Returns:
            Sigma点矩阵
        """
        return self._generate_sigma_points(self.state, self.covariance)
    
    def get_predicted_measurement(self) -> Tuple[np.ndarray, np.ndarray]:
        """获取预测观测和协方差
        
        Returns:
            (预测观测均值, 预测观测协方差)
        """
        sigma_points = self._generate_sigma_points(self.state, self.covariance)
        n_sigma = sigma_points.shape[0]
        
        measurement_points = np.zeros((n_sigma, self.measurement_dim))
        for i in range(n_sigma):
            measurement_points[i] = self.h(sigma_points[i])
        
        z_pred = np.zeros(self.measurement_dim)
        for i in range(n_sigma):
            z_pred += self.Wm[i] * measurement_points[i]
        
        Pzz = np.zeros((self.measurement_dim, self.measurement_dim))
        for i in range(n_sigma):
            dz = measurement_points[i] - z_pred
            Pzz += self.Wc[i] * np.outer(dz, dz)
        
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
            sigma_points = self._generate_sigma_points(self.state, self.covariance)
            n_sigma = sigma_points.shape[0]
            
            measurement_points = np.zeros((n_sigma, self.measurement_dim))
            for i in range(n_sigma):
                measurement_points[i] = self.h(sigma_points[i])
            
            z_pred = np.zeros(self.measurement_dim)
            for i in range(n_sigma):
                z_pred += self.Wm[i] * measurement_points[i]
            
            Pzz = np.zeros((self.measurement_dim, self.measurement_dim))
            for i in range(n_sigma):
                dz = measurement_points[i] - z_pred
                Pzz += self.Wc[i] * np.outer(dz, dz)
            
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


# 一些常用的非线性函数

def polar_measurement_function_ukf(state: np.ndarray) -> np.ndarray:
    """极坐标观测函数（用于UKF）
    
    Args:
        state: 状态向量 [x, vx, y, vy]
        
    Returns:
        极坐标观测 [range, bearing]
    """
    x, y = state[0], state[2]
    range_meas = np.sqrt(x**2 + y**2)
    bearing_meas = np.arctan2(y, x)
    return np.array([range_meas, bearing_meas])


def coordinated_turn_state_transition(state: np.ndarray, dt: float, 
                                        omega: float = 0.0) -> np.ndarray:
    """协调转弯状态转移函数（用于UKF）
    
    Args:
        state: 状态向量 [x, vx, y, vy]
        dt: 时间步长
        omega: 转弯角速度（如果为0，则使用状态中的角速度）
        
    Returns:
        下一时刻状态
    """
    x, vx, y, vy = state[0], state[1], state[2], state[3]
    
    if abs(omega) < 1e-6:
        # 匀速直线
        return np.array([
            x + vx * dt,
            vx,
            y + vy * dt,
            vy
        ])
    else:
        # 协调转弯
        sin_omega_dt = np.sin(omega * dt)
        cos_omega_dt = np.cos(omega * dt)
        
        new_x = x + (vx * sin_omega_dt - vy * (1 - cos_omega_dt)) / omega
        new_vx = vx * cos_omega_dt - vy * sin_omega_dt
        new_y = y + (vx * (1 - cos_omega_dt) + vy * sin_omega_dt) / omega
        new_vy = vx * sin_omega_dt + vy * cos_omega_dt
        
        return np.array([new_x, new_vx, new_y, new_vy])
