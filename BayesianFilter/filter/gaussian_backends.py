"""
高斯滤波器后端实现

支持 KF/EKF/UKF/CKF，通过 FilterBackend 统一接口。

数值稳定技巧（参考 MATLAB RFS Toolbox）：
1. Cholesky分解计算新息协方差逆和行列式
2. 对数空间计算似然
3. 协方差矩阵对称化：S = (S + S^T) / 2
"""
import numpy as np
from typing import Optional, List, Dict
from scipy.stats import multivariate_normal
from .filter_backend import (
    FilterBackend, MultiTargetBackend, PredictedState,
    cholesky_inv, logsumexp, mahalanobis_distance
)
from .kalman_filter import KalmanFilter
from .extended_kalman_filter import ExtendedKalmanFilter
from .unscented_kalman_filter import UnscentedKalmanFilter
from .cubature_kalman_filter import CubatureKalmanFilter


def _ensure_symmetric(M: np.ndarray) -> np.ndarray:
    """确保矩阵对称，避免数值漂移"""
    return (M + M.T) / 2


class KFBackend(FilterBackend):
    """卡尔曼滤波器后端
    
    线性高斯模型：z = H @ x + v, v ~ N(0, R)
    """
    
    def __init__(self, 
                 state_dim: int = 4,
                 meas_dim: int = 2,
                 process_noise_std: float = 0.1,
                 measurement_noise_std: float = 1.0):
        self._state_dim = state_dim
        self._meas_dim = meas_dim
        self._kf = KalmanFilter(
            state_dim=state_dim,
            measurement_dim=meas_dim,
            process_noise_std=process_noise_std,
            measurement_noise_std=measurement_noise_std
        )
        self._R = np.eye(meas_dim) * measurement_noise_std ** 2
        self._R_inv = np.linalg.inv(self._R)
        self._log_2pi = meas_dim * np.log(2 * np.pi)
    
    @property
    def state_dim(self) -> int:
        return self._state_dim
    
    @property
    def meas_dim(self) -> int:
        return self._meas_dim
    
    @property
    def is_initialized(self) -> bool:
        return self._kf.initialized
    
    def initialize(self, initial_state: np.ndarray, 
                   initial_covariance: Optional[np.ndarray] = None) -> None:
        self._kf.initialize(initial_state, initial_covariance)
    
    def predict(self, dt: float) -> None:
        self._kf.predict(dt)
    
    def update(self, measurement: np.ndarray) -> None:
        self._kf.update(measurement, measurement_covariance=self._R)

    def set_measurement_noise(self, R: np.ndarray) -> None:
        """设置观测噪声协方差矩阵"""
        self._R = np.asarray(R, dtype=np.float64)
        self._kf.set_measurement_noise(self._R)

    def get_state(self) -> np.ndarray:
        return self._kf.state.copy()

    def get_covariance(self) -> np.ndarray:
        return self._kf.covariance.copy()

    def get_predicted_measurement(self) -> np.ndarray:
        """z_pred = H @ x"""
        return self._kf.H @ self._kf.state

    def get_innovation_covariance(self) -> np.ndarray:
        """S = H @ P @ H^T + R，确保对称"""
        H = self._kf.H
        P = self._kf.covariance
        S = H @ P @ H.T + self._R
        return _ensure_symmetric(S)
    
    def compute_likelihood(self, measurement: np.ndarray) -> float:
        """计算似然，使用Cholesky分解保证数值稳定"""
        z_pred = self.get_predicted_measurement()
        S = self.get_innovation_covariance()
        
        # 使用Cholesky分解
        S_inv, det_S = cholesky_inv(S)
        innovation = measurement - z_pred
        
        # 对数似然：-0.5 * (d*log(2π) + log(det(S)) + innovation^T @ S^{-1} @ innovation)
        log_likelihood = -0.5 * (
            self._log_2pi + np.log(det_S) + innovation.T @ S_inv @ innovation
        )
        
        return np.exp(log_likelihood)
    
    def compute_mahalanobis(self, measurement: np.ndarray) -> float:
        """计算马氏距离"""
        z_pred = self.get_predicted_measurement()
        S = self.get_innovation_covariance()
        innovation = measurement - z_pred
        S_inv, _ = cholesky_inv(S)
        return mahalanobis_distance(innovation, S_inv)


class EKFBackend(FilterBackend):
    """扩展卡尔曼滤波器后端
    
    非线性模型：x' = f(x) + w, z = h(x) + v
    
    注意事项（参考MATLAB RFS Toolbox）：
    1. 协方差矩阵对称化：S = (S + S^T) / 2
    2. 使用Cholesky分解计算新息协方差逆
    """
    
    def __init__(self,
                 state_dim: int = 4,
                 meas_dim: int = 2,
                 process_noise_std: float = 0.1,
                 measurement_noise_std: float = 1.0,
                 state_transition_func=None,
                 measurement_func=None,
                 state_transition_jacobian=None,
                 measurement_jacobian=None,
                 R: Optional[np.ndarray] = None):
        self._state_dim = state_dim
        self._meas_dim = meas_dim
        self._ekf = ExtendedKalmanFilter(
            state_dim=state_dim,
            measurement_dim=meas_dim,
            process_noise_std=process_noise_std,
            measurement_noise_std=measurement_noise_std
        )
        # R矩阵：优先使用传入的完整矩阵（极坐标等），否则各向同性默认值
        self._R = np.asarray(R, dtype=np.float64) if R is not None \
            else np.eye(meas_dim) * measurement_noise_std ** 2
        self._ekf.set_measurement_noise(self._R)
        self._log_2pi = meas_dim * np.log(2 * np.pi)

        # 设置非线性函数（正确映射到EKF内部属性名）
        if state_transition_func is not None:
            self._ekf.f = state_transition_func
        if measurement_func is not None:
            self._ekf.h = measurement_func
        if state_transition_jacobian is not None:
            self._ekf.F_func = state_transition_jacobian
        if measurement_jacobian is not None:
            self._ekf.H_func = measurement_jacobian
    
    @property
    def state_dim(self) -> int:
        return self._state_dim
    
    @property
    def meas_dim(self) -> int:
        return self._meas_dim
    
    @property
    def is_initialized(self) -> bool:
        return self._ekf.initialized
    
    def initialize(self, initial_state: np.ndarray, 
                   initial_covariance: Optional[np.ndarray] = None) -> None:
        self._ekf.initialize(initial_state, initial_covariance)
    
    def predict(self, dt: float) -> None:
        self._ekf.predict(dt)
    
    def update(self, measurement: np.ndarray) -> None:
        self._ekf.update(measurement, measurement_covariance=self._R)

    def set_measurement_noise(self, R: np.ndarray) -> None:
        """设置观测噪声协方差矩阵"""
        self._R = np.asarray(R, dtype=np.float64)
        self._ekf.set_measurement_noise(self._R)

    def get_state(self) -> np.ndarray:
        return self._ekf.state.copy()

    def get_covariance(self) -> np.ndarray:
        return self._ekf.covariance.copy()

    def get_predicted_measurement(self) -> np.ndarray:
        """z_pred = h(x)"""
        return self._ekf.h(self._ekf.state)

    def get_innovation_covariance(self) -> np.ndarray:
        """S = H(x) @ P @ H(x)^T + R，确保对称"""
        H = self._ekf.H_func(self._ekf.state)
        P = self._ekf.covariance
        S = H @ P @ H.T + self._R
        return _ensure_symmetric(S)

    def compute_likelihood(self, measurement: np.ndarray) -> float:
        """计算似然，使用Cholesky分解"""
        z_pred = self.get_predicted_measurement()
        S = self.get_innovation_covariance()

        S_inv, det_S = cholesky_inv(S)
        innovation = measurement - z_pred

        log_likelihood = -0.5 * (
            self._log_2pi + np.log(det_S) + innovation.T @ S_inv @ innovation
        )

        return np.exp(log_likelihood)

    def compute_mahalanobis(self, measurement: np.ndarray) -> float:
        """计算马氏距离"""
        z_pred = self.get_predicted_measurement()
        S = self.get_innovation_covariance()
        innovation = measurement - z_pred
        S_inv, _ = cholesky_inv(S)
        return mahalanobis_distance(innovation, S_inv)


class UKFBackend(FilterBackend):
    """无迹卡尔曼滤波器后端
    
    非线性模型，使用sigma点变换
    
    注意事项（参考MATLAB RFS Toolbox ut.m）：
    1. 使用Cholesky分解计算sigma点：Psqrtm = chol((n+λ)P)'
    2. 权重调整：Wc[0] += (1 - alpha^2 + beta)
    3. 协方差矩阵对称化
    """
    
    def __init__(self,
                 state_dim: int = 4,
                 meas_dim: int = 2,
                 process_noise_std: float = 0.1,
                 measurement_noise_std: float = 1.0,
                 state_transition_func=None,
                 measurement_func=None,
                 alpha: float = 1e-3,
                 beta: float = 2.0,
                 kappa: float = 0.0,
                 R: Optional[np.ndarray] = None):
        self._state_dim = state_dim
        self._meas_dim = meas_dim
        self._ukf = UnscentedKalmanFilter(
            state_dim=state_dim,
            measurement_dim=meas_dim,
            process_noise_std=process_noise_std,
            measurement_noise_std=measurement_noise_std,
            alpha=alpha,
            beta=beta,
            kappa=kappa
        )
        # R矩阵：优先使用传入的完整矩阵（极坐标等），否则各向同性默认值
        self._R = np.asarray(R, dtype=np.float64) if R is not None \
            else np.eye(meas_dim) * measurement_noise_std ** 2
        self._ukf.set_measurement_noise(self._R)
        self._log_2pi = meas_dim * np.log(2 * np.pi)

        if state_transition_func is not None:
            self._ukf.f = state_transition_func
        if measurement_func is not None:
            self._ukf.h = measurement_func

        # 保存UKF参数用于手动计算
        self._alpha = alpha
        self._beta = beta
        self._kappa = kappa
    
    @property
    def state_dim(self) -> int:
        return self._state_dim
    
    @property
    def meas_dim(self) -> int:
        return self._meas_dim
    
    @property
    def is_initialized(self) -> bool:
        return self._ukf.initialized
    
    def initialize(self, initial_state: np.ndarray, 
                   initial_covariance: Optional[np.ndarray] = None) -> None:
        self._ukf.initialize(initial_state, initial_covariance)
    
    def predict(self, dt: float) -> None:
        self._ukf.predict(dt)
    
    def update(self, measurement: np.ndarray) -> None:
        self._ukf.update(measurement, measurement_covariance=self._R)

    def set_measurement_noise(self, R: np.ndarray) -> None:
        """设置观测噪声协方差矩阵"""
        self._R = np.asarray(R, dtype=np.float64)
        self._ukf.set_measurement_noise(self._R)

    def get_state(self) -> np.ndarray:
        return self._ukf.state.copy()
    
    def get_covariance(self) -> np.ndarray:
        return self._ukf.covariance.copy()
    
    def get_predicted_measurement(self) -> np.ndarray:
        """UKF的预测观测通过sigma点变换获得"""
        sigma_points = self._ukf._generate_sigma_points(
            self._ukf.state, self._ukf.covariance
        )
        transformed = np.array([self._ukf.h(sp) for sp in sigma_points])
        return np.dot(self._ukf.Wm, transformed)
    
    def get_innovation_covariance(self) -> np.ndarray:
        """UKF的新息协方差通过sigma点变换获得
        
        注意：Wc[0] 已包含 (1 - alpha^2 + beta) 调整
        """
        sigma_points = self._ukf._generate_sigma_points(
            self._ukf.state, self._ukf.covariance
        )
        transformed = np.array([self._ukf.h(sp) for sp in sigma_points])
        z_pred = np.dot(self._ukf.Wm, transformed)
        
        n_sigma = sigma_points.shape[0]
        S = np.zeros((self._meas_dim, self._meas_dim))
        for i in range(n_sigma):
            diff = transformed[i] - z_pred
            S += self._ukf.Wc[i] * np.outer(diff, diff)
        
        # 确保对称
        return _ensure_symmetric(S + self._R)
    
    def compute_likelihood(self, measurement: np.ndarray) -> float:
        """计算似然，使用Cholesky分解"""
        z_pred = self.get_predicted_measurement()
        S = self.get_innovation_covariance()
        
        S_inv, det_S = cholesky_inv(S)
        innovation = measurement - z_pred
        
        log_likelihood = -0.5 * (
            self._log_2pi + np.log(det_S) + innovation.T @ S_inv @ innovation
        )
        
        return np.exp(log_likelihood)
    
    def compute_mahalanobis(self, measurement: np.ndarray) -> float:
        """计算马氏距离"""
        z_pred = self.get_predicted_measurement()
        S = self.get_innovation_covariance()
        innovation = measurement - z_pred
        S_inv, _ = cholesky_inv(S)
        return mahalanobis_distance(innovation, S_inv)


class MultiTargetFilterManager(MultiTargetBackend):
    """多目标滤波器管理器
    
    管理多个目标的滤波器后端，提供统一的数据关联接口。
    """
    
    def __init__(self, backend_factory):
        """
        Args:
            backend_factory: 创建 FilterBackend 的工厂函数
                           例如: lambda: KFBackend(state_dim=4, meas_dim=2)
        """
        self._backend_factory = backend_factory
        self._targets: Dict[int, FilterBackend] = {}
    
    def add_target(self, target_id: int, 
                   initial_state: np.ndarray,
                   initial_covariance: Optional[np.ndarray] = None) -> None:
        backend = self._backend_factory()
        backend.initialize(initial_state, initial_covariance)
        self._targets[target_id] = backend
    
    def remove_target(self, target_id: int) -> None:
        if target_id in self._targets:
            del self._targets[target_id]
    
    def predict_all(self, dt: float) -> None:
        for backend in self._targets.values():
            backend.predict(dt)
    
    def get_predicted_states(self) -> List[PredictedState]:
        states = []
        for target_id, backend in self._targets.items():
            if backend.is_initialized:
                states.append(PredictedState(
                    state=backend.get_state(),
                    covariance=backend.get_covariance(),
                    predicted_meas=backend.get_predicted_measurement(),
                    innovation_cov=backend.get_innovation_covariance(),
                    target_id=target_id
                ))
        return states
    
    def update_target(self, target_id: int, 
                      measurement: np.ndarray) -> None:
        if target_id in self._targets:
            self._targets[target_id].update(measurement)
    
    def get_target_ids(self) -> List[int]:
        return list(self._targets.keys())
