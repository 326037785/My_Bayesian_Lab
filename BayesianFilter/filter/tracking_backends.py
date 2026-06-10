"""
滤波器后端 - 正确的设计

核心思想：
1. 场景决定测量模型（线性 vs 非线性）
2. 滤波器封装预测和更新的数学细节
3. 数据关联只关心 predicted_meas 和 innovation_cov

场景类型：
- 线性跟踪：z = H @ x（直接矩阵提取）
- 雷达跟踪：z = [range, bearing] = h(x)（极坐标转换）
"""
import numpy as np
from typing import Optional, Callable
from .filter_backend import FilterBackend, PredictedState, cholesky_inv
from .kalman_filter import KalmanFilter
from .extended_kalman_filter import ExtendedKalmanFilter
from .unscented_kalman_filter import UnscentedKalmanFilter


def _ensure_symmetric(M: np.ndarray) -> np.ndarray:
    """确保矩阵对称"""
    return (M + M.T) / 2


class LinearKFBackend(FilterBackend):
    """线性卡尔曼滤波器后端
    
    适用于：线性跟踪问题
    测量模型：z = H @ x + v, v ~ N(0, R)
    
    典型场景：
    - 2D/3D笛卡尔坐标跟踪
    - 位置直接观测
    """
    
    def __init__(self, state_dim=4, meas_dim=2,
                 process_noise_std=0.1, measurement_noise_std=1.0,
                 H: Optional[np.ndarray] = None):
        self._state_dim = state_dim
        self._meas_dim = meas_dim
        self._kf = KalmanFilter(state_dim, meas_dim, 
                                process_noise_std, measurement_noise_std)
        # 测量矩阵：默认提取位置 [x, y]
        self._H = H if H is not None else self._kf.H
        self._R = np.eye(meas_dim) * measurement_noise_std ** 2
        self._log_2pi = meas_dim * np.log(2 * np.pi)
    
    @property
    def state_dim(self): return self._state_dim
    
    @property
    def meas_dim(self): return self._meas_dim
    
    @property
    def is_initialized(self): return self._kf.initialized
    
    def initialize(self, state, cov=None):
        self._kf.initialize(state, cov)
    
    def predict(self, dt):
        self._kf.predict(dt)
    
    def update(self, z):
        self._kf.update(z)
    
    def get_state(self):
        return self._kf.state.copy()
    
    def get_covariance(self):
        return self._kf.covariance.copy()
    
    def get_predicted_measurement(self):
        """线性测量：z_pred = H @ x"""
        return self._H @ self._kf.state
    
    def get_innovation_covariance(self):
        """S = H @ P @ H^T + R"""
        H = self._H
        P = self._kf.covariance
        return _ensure_symmetric(H @ P @ H.T + self._R)
    
    def compute_likelihood(self, z):
        z_pred = self.get_predicted_measurement()
        S = self.get_innovation_covariance()
        S_inv, det_S = cholesky_inv(S)
        innovation = z - z_pred
        log_lik = -0.5 * (self._log_2pi + np.log(det_S) + 
                          innovation.T @ S_inv @ innovation)
        return np.exp(log_lik)
    
    def compute_mahalanobis(self, z):
        z_pred = self.get_predicted_measurement()
        S = self.get_innovation_covariance()
        innovation = z - z_pred
        S_inv, _ = cholesky_inv(S)
        return np.sqrt(innovation.T @ S_inv @ innovation)


class EKFBackend(FilterBackend):
    """扩展卡尔曼滤波器后端
    
    适用于：非线性跟踪问题
    测量模型：z = h(x) + v, v ~ N(0, R)
    
    典型场景：
    - 雷达极坐标观测：z = [range, bearing] = [sqrt(x^2+y^2), atan2(y,x)]
    - 非线性传感器
    """
    
    def __init__(self, state_dim=4, meas_dim=2,
                 process_noise_std=0.1, measurement_noise_std=1.0,
                 h: Callable = None, H_jacobian: Callable = None,
                 R: np.ndarray = None):
        """
        Args:
            h: 非线性测量函数 h(x) -> z
            H_jacobian: 测量雅可比矩阵 H(x) -> dh/dx
            R: 观测噪声协方差矩阵。当提供时直接使用（极坐标场景）。
               当为 None 时回退到 diag(measurement_noise_std²)。
        """
        self._state_dim = state_dim
        self._meas_dim = meas_dim

        if h is None or H_jacobian is None:
            raise ValueError("EKF需要提供测量函数h和雅可比H_jacobian")

        self._h = h
        self._H_jacobian = H_jacobian
        self._R = np.asarray(R, dtype=np.float64) if R is not None \
            else np.eye(meas_dim) * measurement_noise_std ** 2
        self._log_2pi = meas_dim * np.log(2 * np.pi)
        
        # 内部状态
        self._state = None
        self._cov = None
        self._initialized = False
        
        # 过程噪声
        self._Q_func = lambda dt: self._get_process_noise(dt, state_dim, process_noise_std)
    
    def _get_process_noise(self, dt, dim, std):
        q = std ** 2
        if dim == 4:
            return q * np.array([
                [dt**3/3, dt**2/2, 0, 0],
                [dt**2/2, dt, 0, 0],
                [0, 0, dt**3/3, dt**2/2],
                [0, 0, dt**2/2, dt]
            ])
        return q * dt * np.eye(dim)
    
    @property
    def state_dim(self): return self._state_dim
    
    @property
    def meas_dim(self): return self._meas_dim
    
    @property
    def is_initialized(self): return self._initialized
    
    def initialize(self, state, cov=None):
        self._state = state.copy()
        self._cov = cov.copy() if cov is not None else np.eye(self._state_dim) * 100
        self._initialized = True
    
    def predict(self, dt):
        """EKF预测：线性化状态转移（匀速模型）"""
        F = np.array([
            [1, dt, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, dt],
            [0, 0, 0, 1]
        ])[:self._state_dim, :self._state_dim]
        
        Q = self._Q_func(dt)
        self._state = F @ self._state
        self._cov = _ensure_symmetric(F @ self._cov @ F.T + Q)
    
    def update(self, z):
        """EKF更新：使用雅可比线性化测量函数"""
        x = self._state
        P = self._cov
        
        # 预测观测
        z_pred = self._h(x)
        
        # 雅可比矩阵
        H = self._H_jacobian(x)
        
        # 新息协方差
        S = _ensure_symmetric(H @ P @ H.T + self._R)
        
        # 卡尔曼增益
        S_inv, _ = cholesky_inv(S)
        K = P @ H.T @ S_inv
        
        # 更新
        innovation = z - z_pred
        self._state = x + K @ innovation
        self._cov = _ensure_symmetric((np.eye(self._state_dim) - K @ H) @ P)
    
    def get_state(self):
        return self._state.copy()
    
    def get_covariance(self):
        return self._cov.copy()
    
    def get_predicted_measurement(self):
        """非线性测量：z_pred = h(x)"""
        return self._h(self._state)
    
    def get_innovation_covariance(self):
        """S = H(x) @ P @ H(x)^T + R"""
        H = self._H_jacobian(self._state)
        return _ensure_symmetric(H @ self._cov @ H.T + self._R)
    
    def compute_likelihood(self, z):
        z_pred = self.get_predicted_measurement()
        S = self.get_innovation_covariance()
        S_inv, det_S = cholesky_inv(S)
        innovation = z - z_pred
        log_lik = -0.5 * (self._log_2pi + np.log(det_S) + 
                          innovation.T @ S_inv @ innovation)
        return np.exp(log_lik)
    
    def compute_mahalanobis(self, z):
        z_pred = self.get_predicted_measurement()
        S = self.get_innovation_covariance()
        innovation = z - z_pred
        S_inv, _ = cholesky_inv(S)
        return np.sqrt(innovation.T @ S_inv @ innovation)


class UKFBackend(FilterBackend):
    """无迹卡尔曼滤波器后端
    
    适用于：非线性跟踪问题（无需计算雅可比）
    测量模型：z = h(x) + v, v ~ N(0, R)
    
    与EKF的区别：
    - EKF：需要显式计算雅可比矩阵 H(x)
    - UKF：通过sigma点变换隐式处理非线性
    """
    
    def __init__(self, state_dim=4, meas_dim=2,
                 process_noise_std=0.1, measurement_noise_std=1.0,
                 h: Callable = None,
                 alpha: float = 1e-3, beta: float = 2.0, kappa: float = 0.0,
                 R: np.ndarray = None):
        """
        Args:
            h: 非线性测量函数 h(x) -> z
            alpha, beta, kappa: UKF参数
            R: 观测噪声协方差矩阵。当提供时直接使用（极坐标场景）。
               当为 None 时回退到 diag(measurement_noise_std²)。
        """
        self._state_dim = state_dim
        self._meas_dim = meas_dim

        if h is None:
            raise ValueError("UKF需要提供测量函数h")

        self._h = h
        self._R = np.asarray(R, dtype=np.float64) if R is not None \
            else np.eye(meas_dim) * measurement_noise_std ** 2
        self._log_2pi = meas_dim * np.log(2 * np.pi)
        
        # UKF参数
        self._alpha = alpha
        self._beta = beta
        self._kappa = kappa
        
        # 内部状态
        self._state = None
        self._cov = None
        self._initialized = False
        
        # 过程噪声
        self._Q_func = lambda dt: self._get_process_noise(dt, state_dim, process_noise_std)
    
    def _get_process_noise(self, dt, dim, std):
        q = std ** 2
        if dim == 4:
            return q * np.array([
                [dt**3/3, dt**2/2, 0, 0],
                [dt**2/2, dt, 0, 0],
                [0, 0, dt**3/3, dt**2/2],
                [0, 0, dt**2/2, dt]
            ])
        return q * dt * np.eye(dim)
    
    def _generate_sigma_points(self, x, P):
        """生成sigma点"""
        n = len(x)
        lam = self._alpha**2 * (n + self._kappa) - n
        
        # Cholesky分解
        try:
            L = np.linalg.cholesky((n + lam) * P)
        except np.linalg.LinAlgError:
            # 如果分解失败，添加小扰动
            L = np.linalg.cholesky((n + lam) * P + np.eye(n) * 1e-6)
        
        sigma = np.zeros((2*n + 1, n))
        sigma[0] = x
        for i in range(n):
            sigma[i+1] = x + L[:, i]
            sigma[n+i+1] = x - L[:, i]
        
        # 权重
        Wm = np.zeros(2*n + 1)
        Wc = np.zeros(2*n + 1)
        Wm[0] = lam / (n + lam)
        Wc[0] = lam / (n + lam) + (1 - self._alpha**2 + self._beta)
        for i in range(1, 2*n + 1):
            Wm[i] = 0.5 / (n + lam)
            Wc[i] = 0.5 / (n + lam)
        
        return sigma, Wm, Wc
    
    @property
    def state_dim(self): return self._state_dim
    
    @property
    def meas_dim(self): return self._meas_dim
    
    @property
    def is_initialized(self): return self._initialized
    
    def initialize(self, state, cov=None):
        self._state = state.copy()
        self._cov = cov.copy() if cov is not None else np.eye(self._state_dim) * 100
        self._initialized = True
    
    def predict(self, dt):
        """UKF预测：sigma点变换"""
        F = np.array([
            [1, dt, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, dt],
            [0, 0, 0, 1]
        ])[:self._state_dim, :self._state_dim]
        
        Q = self._Q_func(dt)
        
        # 生成sigma点
        sigma, Wm, Wc = self._generate_sigma_points(self._state, self._cov)
        
        # 传播sigma点
        sigma_pred = np.array([F @ s for s in sigma])
        
        # 预测均值
        x_pred = np.sum(Wm[:, np.newaxis] * sigma_pred, axis=0)
        
        # 预测协方差
        P_pred = np.zeros_like(self._cov)
        for i in range(len(sigma)):
            diff = sigma_pred[i] - x_pred
            P_pred += Wc[i] * np.outer(diff, diff)
        
        self._state = x_pred
        self._cov = _ensure_symmetric(P_pred + Q)
    
    def update(self, z):
        """UKF更新：sigma点变换"""
        sigma, Wm, Wc = self._generate_sigma_points(self._state, self._cov)
        
        # 观测sigma点
        z_sigma = np.array([self._h(s) for s in sigma])
        
        # 预测观测
        z_pred = np.sum(Wm[:, np.newaxis] * z_sigma, axis=0)
        
        # 新息协方差
        S = np.zeros((self._meas_dim, self._meas_dim))
        for i in range(len(sigma)):
            diff = z_sigma[i] - z_pred
            S += Wc[i] * np.outer(diff, diff)
        S = _ensure_symmetric(S + self._R)
        
        # 互协方差
        Pxz = np.zeros((self._state_dim, self._meas_dim))
        for i in range(len(sigma)):
            dx = sigma[i] - self._state
            dz = z_sigma[i] - z_pred
            Pxz += Wc[i] * np.outer(dx, dz)
        
        # 卡尔曼增益
        S_inv, _ = cholesky_inv(S)
        K = Pxz @ S_inv
        
        # 更新
        innovation = z - z_pred
        self._state = self._state + K @ innovation
        self._cov = _ensure_symmetric(self._cov - K @ S @ K.T)
    
    def get_state(self):
        return self._state.copy()
    
    def get_covariance(self):
        return self._cov.copy()
    
    def get_predicted_measurement(self):
        """通过sigma点变换计算预测观测"""
        sigma, Wm, _ = self._generate_sigma_points(self._state, self._cov)
        z_sigma = np.array([self._h(s) for s in sigma])
        return np.sum(Wm[:, np.newaxis] * z_sigma, axis=0)
    
    def get_innovation_covariance(self):
        """通过sigma点变换计算新息协方差"""
        sigma, Wm, Wc = self._generate_sigma_points(self._state, self._cov)
        z_sigma = np.array([self._h(s) for s in sigma])
        z_pred = np.sum(Wm[:, np.newaxis] * z_sigma, axis=0)
        
        S = np.zeros((self._meas_dim, self._meas_dim))
        for i in range(len(sigma)):
            diff = z_sigma[i] - z_pred
            S += Wc[i] * np.outer(diff, diff)
        return _ensure_symmetric(S + self._R)
    
    def compute_likelihood(self, z):
        z_pred = self.get_predicted_measurement()
        S = self.get_innovation_covariance()
        S_inv, det_S = cholesky_inv(S)
        innovation = z - z_pred
        log_lik = -0.5 * (self._log_2pi + np.log(det_S) + 
                          innovation.T @ S_inv @ innovation)
        return np.exp(log_lik)
    
    def compute_mahalanobis(self, z):
        z_pred = self.get_predicted_measurement()
        S = self.get_innovation_covariance()
        innovation = z - z_pred
        S_inv, _ = cholesky_inv(S)
        return np.sqrt(innovation.T @ S_inv @ innovation)
