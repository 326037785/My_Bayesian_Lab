"""
概率假设密度 (PHD) 滤波器

支持任意状态转移函数和滤波器类型（KF/EKF/UKF/CKF）
"""
import numpy as np
from typing import Optional, List, Dict, Tuple, Callable, Any
from dataclasses import dataclass
from scipy.stats import chi2


@dataclass
class GaussianComponent:
    """高斯分量"""
    weight: float
    mean: np.ndarray
    covariance: np.ndarray


class PHDFilter:
    """概率假设密度滤波器

    GM-PHD滤波器的Python实现，参考：
    B.-N. Vo and W.-K. Ma, "The Gaussian mixture Probability Hypothesis Density Filter,"
    IEEE Trans Signal Processing, Vol. 54, No. 11, pp. 4091-4104, 2006.

    支持任意状态转移函数和滤波器类型（KF/EKF/UKF/CKF）

    出生模型（Birth Model）说明：
    ---------------------------
    出生分量是GM-PHD的关键组成部分——它们代表监控区域中可能出现新目标的位置。
    由于出生分量高度依赖于实际场景（监控区域大小、目标运动特性等），应始终由
    外部指定，不要依赖默认值。

    三种指定方式（按优先级）：
    1. `birth_components`: 直接传入高斯分量列表（最精确，推荐生产使用）
    2. `surveillance_bounds` + `birth_covariance`: 根据边界自动生成网格分布的分量
    3. 都不传：回退到基于 `surveillance_area` 的粗略默认值（仅用于快速demo）
    """

    def __init__(self,
                 state_dim: int = 4,
                 measurement_dim: int = 2,
                 survival_probability: float = 0.99,
                 detection_probability: float = 0.98,
                 clutter_rate: float = 60.0,
                 surveillance_area: float = 1000000.0,
                 birth_weight: float = 0.03,
                 pruning_threshold: float = 1e-5,
                 merging_threshold: float = 4.0,
                 max_components: int = 100,
                 gating_threshold: float = 13.816,  # chi2inv(0.999, 2) (MATLAB default)
                 # 状态转移函数
                 state_transition_func: Optional[Callable] = None,
                 state_transition_jacobian: Optional[Callable] = None,
                 # 观测函数
                 measurement_func: Optional[Callable] = None,
                 measurement_jacobian: Optional[Callable] = None,
                 # 噪声矩阵
                 process_noise_matrix: Optional[np.ndarray] = None,
                 measurement_noise_matrix: Optional[np.ndarray] = None,
                 # 滤波器类型
                 filter_type: str = 'KF',
                 # === 出生模型（Birth Model）===
                 # 直接指定出生分量列表（优先级最高）
                 birth_components: Optional[List[GaussianComponent]] = None,
                 # 监控区域边界：((x_min, x_max), (y_min, y_max))，用于自动生成出生分量
                 surveillance_bounds: Optional[Tuple[Tuple[float, float],
                                                      Tuple[float, float]]] = None,
                 # 出生分量协方差矩阵（自动生成时使用）
                 birth_covariance: Optional[np.ndarray] = None,
                 # 自动生成的出生分量数量（默认4：上下左右各一，MATLAB风格）
                 n_birth_components: int = 4):
        """
        初始化PHD滤波器

        Args:
            state_dim: 状态维度
            measurement_dim: 观测维度
            survival_probability: 存活概率 P_S
            detection_probability: 检测概率 P_D
            clutter_rate: 杂波率（泊松分布的平均杂波数 lambda_c）
            surveillance_area: 监视区域面积（用于计算杂波密度 kappa = lambda_c / V）
            birth_weight: 出生分量权重（自动生成时使用）
            pruning_threshold: 剪枝阈值 T
            merging_threshold: 合并阈值 U（马氏距离平方）
            max_components: 最大分量数 J_max
            gating_threshold: 门限阈值（卡方分布分位数）
            state_transition_func: 状态转移函数 f(x, dt) -> x_next
            state_transition_jacobian: 状态转移雅可比矩阵 F(x, dt) -> F
            measurement_func: 观测函数 h(x) -> z
            measurement_jacobian: 观测雅可比矩阵 H(x) -> H
            process_noise_matrix: 过程噪声协方差矩阵 Q
            measurement_noise_matrix: 观测噪声协方差矩阵 R
            filter_type: 滤波器类型 ('KF', 'EKF', 'UKF', 'CKF')

            birth_components: 出生分量列表。若提供，直接使用（推荐方式）。
            surveillance_bounds: 监控区域边界，格式 ((x_min, x_max), (y_min, y_max))。
                用于自动生成出生分量时的网格分布。
            birth_covariance: 出生分量协方差矩阵。
                用于自动生成。应与场景大小匹配（位置方差 ~= (区域半径/3)^2）。
            n_birth_components: 自动生成的出生分量数量。默认9（3x3网格）。
        """
        self.state_dim = state_dim
        self.measurement_dim = measurement_dim
        self.P_S = survival_probability
        self.P_D = detection_probability
        self.Q_D = 1 - detection_probability
        self.clutter_rate = clutter_rate
        # 杂波密度 = 杂波率 / 监视区域面积
        self.clutter_density = clutter_rate / surveillance_area
        self.birth_weight = birth_weight
        self.pruning_threshold = pruning_threshold
        self.merging_threshold = merging_threshold
        self.max_components = max_components
        self.gating_threshold = gating_threshold
        self.filter_type = filter_type
        self._surveillance_bounds = surveillance_bounds

        # 设置状态转移函数
        self.f = state_transition_func if state_transition_func is not None else self._default_state_transition
        self.F_func = state_transition_jacobian if state_transition_jacobian is not None else self._default_F_jacobian

        # 设置观测函数
        self.h = measurement_func if measurement_func is not None else self._default_measurement
        self.H_func = measurement_jacobian if measurement_jacobian is not None else self._default_H_jacobian

        # 设置噪声矩阵
        self.Q = process_noise_matrix if process_noise_matrix is not None else self._default_process_noise(1.0)
        self.R = measurement_noise_matrix if measurement_noise_matrix is not None else self._default_measurement_noise()

        # === 出生分量初始化 ===
        if birth_components is not None:
            # 方式1：用户直接指定（推荐生产使用）
            self.birth_components = list(birth_components)
        else:
            # 方式2/3：自动生成
            self.birth_components = self._auto_generate_birth_components(
                surveillance_bounds=surveillance_bounds,
                birth_covariance=birth_covariance,
                n_components=n_birth_components,
                birth_weight=birth_weight
            )

        # PHD分量（初始化为空）
        self.components: List[GaussianComponent] = []
        self._initialize_components()
    
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
    
    def _default_F_jacobian(self, state: np.ndarray, dt: float) -> np.ndarray:
        """默认状态转移雅可比矩阵"""
        if self.state_dim == 4:
            return np.array([
                [1, dt, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, dt],
                [0, 0, 0, 1]
            ])
        elif self.state_dim == 6:
            return np.array([
                [1, dt, dt**2/2, 0, 0, 0],
                [0, 1, dt, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, 1, dt, dt**2/2],
                [0, 0, 0, 0, 1, dt],
                [0, 0, 0, 0, 0, 1]
            ])
        return np.eye(self.state_dim)
    
    def _default_measurement(self, state: np.ndarray) -> np.ndarray:
        """默认观测函数（直接观测位置）"""
        if self.state_dim >= 4:
            return np.array([state[0], state[2]])
        elif self.state_dim >= 2:
            return np.array([state[0], state[1]])
        return state[:self.measurement_dim]
    
    def _default_H_jacobian(self, state: np.ndarray) -> np.ndarray:
        """默认观测雅可比矩阵"""
        if self.state_dim == 4:
            return np.array([
                [1, 0, 0, 0],
                [0, 0, 1, 0]
            ])
        elif self.state_dim == 6:
            return np.array([
                [1, 0, 0, 0, 0, 0],
                [0, 0, 0, 1, 0, 0]
            ])
        H = np.zeros((self.measurement_dim, self.state_dim))
        H[0, 0] = 1
        if self.measurement_dim > 1 and self.state_dim > 2:
            H[1, 2] = 1
        return H
    
    def _default_process_noise(self, dt: float) -> np.ndarray:
        """默认过程噪声协方差矩阵"""
        q = 5.0 ** 2
        
        if self.state_dim == 4:
            return q * np.array([
                [dt**3/3, dt**2/2, 0, 0],
                [dt**2/2, dt, 0, 0],
                [0, 0, dt**3/3, dt**2/2],
                [0, 0, dt**2/2, dt]
            ])
        elif self.state_dim == 6:
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
        return q * dt * np.eye(self.state_dim)
    
    def _default_measurement_noise(self) -> np.ndarray:
        """默认观测噪声协方差矩阵"""
        return np.diag([10.0, 10.0]) ** 2
    
    @staticmethod
    def _auto_generate_birth_components(
            surveillance_bounds: Optional[Tuple[Tuple[float, float],
                                                  Tuple[float, float]]] = None,
            birth_covariance: Optional[np.ndarray] = None,
            n_components: int = 4,
            birth_weight: float = 0.03) -> List[GaussianComponent]:
        """根据场景参数自动生成出生分量

        策略（参考 MATLAB demo2 resize_cov）：
        1. 出生分量放在监控区域外围（上/右/下/左四边中点），而非网格
           —— 外围位置不与中心跟踪分量竞争测量，但大协方差覆盖整个场景。
        2. 协方差与场景大小成比例：pos_std = 0.5 × range，vel_std = 0.2 × range。
        3. 若 n_components > 4，额外分量均匀分布在外围环上。

        Args:
            surveillance_bounds: 监控区域边界 ((x_min, x_max), (y_min, y_max))。
            birth_covariance: 出生分量协方差。若 None，按场景动态计算。
            n_components: 分量数量（默认 4：上下左右各一）。
            birth_weight: 每个分量的总权重。

        Returns:
            出生分量列表
        """
        if surveillance_bounds is None:
            x_lim, y_lim = (-600.0, 600.0), (-600.0, 600.0)
        else:
            x_lim, y_lim = surveillance_bounds

        x_span = x_lim[1] - x_lim[0]
        y_span = y_lim[1] - y_lim[0]
        x_mid = (x_lim[0] + x_lim[1]) / 2.0
        y_mid = (y_lim[0] + y_lim[1]) / 2.0
        margin = 0.15  # 15% 外扩（MATLAB demo2）

        x_left = x_lim[0] - margin * x_span
        x_right = x_lim[1] + margin * x_span
        y_bottom = y_lim[0] - margin * y_span
        y_top = y_lim[1] + margin * y_span

        # 协方差：与场景大小成比例（MATLAB demo2 resize_cov 公式）
        if birth_covariance is None:
            pos_std = 0.5 * max(x_span, y_span)   # 大协方差覆盖整个场景
            vel_std = 0.2 * max(x_span, y_span)   # 速度不确定性
            birth_covariance = np.diag([pos_std, vel_std, pos_std, vel_std]) ** 2

        # 出生分量位置：外围环形分布
        if n_components <= 4:
            # 4 分量：上下左右四边中点
            positions = [
                (x_mid, y_top),       # 上
                (x_right, y_mid),     # 右
                (x_mid, y_bottom),    # 下
                (x_left, y_mid),      # 左
            ]
        else:
            # 多于 4：在外围环上均匀分布
            positions = []
            n_each = max(2, n_components // 4)
            # 上边
            for x in np.linspace(x_left, x_right, n_each + 2)[1:-1]:
                positions.append((x, y_top))
            # 右边
            for y in np.linspace(y_bottom, y_top, n_each + 2)[1:-1]:
                positions.append((x_right, y))
            # 下边
            for x in np.linspace(x_left, x_right, n_each + 2)[1:-1]:
                positions.append((x, y_bottom))
            # 左边
            for y in np.linspace(y_bottom, y_top, n_each + 2)[1:-1]:
                positions.append((x_left, y))
            positions = positions[:n_components]

        n_birth = len(positions)
        weight_per = birth_weight / n_birth

        components = []
        for x0, y0 in positions:
            components.append(GaussianComponent(
                weight=weight_per,
                mean=np.array([x0, 0.0, y0, 0.0]),
                covariance=birth_covariance.copy()
            ))

        return components

    def set_birth_components(self, components: List[GaussianComponent]) -> None:
        """手动设置出生分量（替换当前所有出生分量）

        典型用法：
            bounds = ((-500, 500), (-500, 500))
            pos_std = 200.0  # 根据场景调整
            vel_std = 20.0
            cov = np.diag([pos_std, vel_std, pos_std, vel_std]) ** 2
            comps = [GaussianComponent(weight=0.005, mean=np.array([x,0,y,0]), covariance=cov)
                     for x in [-400, 0, 400] for y in [-400, 0, 400]]
            phd.set_birth_components(comps)

        Args:
            components: 新的出生分量列表
        """
        self.birth_components = list(components)

    def _initialize_components(self) -> None:
        """初始化PHD分量"""
        self.components.append(GaussianComponent(
            weight=np.finfo(float).eps,
            mean=np.array([0.1, 0, 0.1, 0]),
            covariance=np.diag([1, 1, 1, 1]) ** 2
        ))
    
    def _predict_single(self, comp: GaussianComponent, dt: float) -> GaussianComponent:
        """预测单个分量
        
        Args:
            comp: 高斯分量
            dt: 时间步长
            
        Returns:
            预测后的分量
        """
        if self.filter_type == 'KF':
            # 线性卡尔曼滤波
            F = self.F_func(comp.mean, dt)
            predicted_mean = F @ comp.mean
            predicted_cov = F @ comp.covariance @ F.T + self.Q
        elif self.filter_type == 'EKF':
            # 扩展卡尔曼滤波
            predicted_mean = self.f(comp.mean, dt)
            F = self.F_func(comp.mean, dt)
            predicted_cov = F @ comp.covariance @ F.T + self.Q
        elif self.filter_type == 'UKF':
            # 无迹卡尔曼滤波
            predicted_mean, predicted_cov = self._ukf_predict(comp.mean, comp.covariance, dt)
        elif self.filter_type == 'CKF':
            # 容积卡尔曼滤波
            predicted_mean, predicted_cov = self._ckf_predict(comp.mean, comp.covariance, dt)
        else:
            raise ValueError(f"Unknown filter type: {self.filter_type}")
        
        # 确保协方差对称
        predicted_cov = (predicted_cov + predicted_cov.T) / 2
        
        return GaussianComponent(
            weight=self.P_S * comp.weight,
            mean=predicted_mean,
            covariance=predicted_cov
        )
    
    def _ukf_predict(self, mean: np.ndarray, covariance: np.ndarray, dt: float) -> Tuple[np.ndarray, np.ndarray]:
        """UKF预测"""
        n = len(mean)
        alpha = 1e-3
        beta = 2.0
        kappa = 0.0
        lambda_ = alpha ** 2 * (n + kappa) - n
        
        # 生成Sigma点
        try:
            L = np.linalg.cholesky((n + lambda_) * covariance)
        except np.linalg.LinAlgError:
            eigenvalues, eigenvectors = np.linalg.eigh(covariance)
            eigenvalues = np.maximum(eigenvalues, 1e-6)
            L = eigenvectors @ np.diag(np.sqrt(eigenvalues * (n + lambda_)))
        
        sigma_points = np.zeros((2 * n + 1, n))
        sigma_points[0] = mean
        sigma_points[1:n+1] = mean + L.T
        sigma_points[n+1:] = mean - L.T
        
        # 计算权重
        Wm = np.zeros(2 * n + 1)
        Wm[0] = lambda_ / (n + lambda_)
        Wm[1:] = 1.0 / (2 * (n + lambda_))
        
        Wc = np.zeros(2 * n + 1)
        Wc[0] = lambda_ / (n + lambda_) + (1 - alpha ** 2 + beta)
        Wc[1:] = 1.0 / (2 * (n + lambda_))
        
        # 传播Sigma点
        propagated_points = np.zeros_like(sigma_points)
        for i in range(2 * n + 1):
            propagated_points[i] = self.f(sigma_points[i], dt)
        
        # 计算预测均值和协方差
        predicted_mean = np.dot(Wm, propagated_points)
        diff = propagated_points - predicted_mean
        predicted_cov = diff.T @ (Wc[:, np.newaxis] * diff) + self.Q
        
        return predicted_mean, predicted_cov
    
    def _ckf_predict(self, mean: np.ndarray, covariance: np.ndarray, dt: float) -> Tuple[np.ndarray, np.ndarray]:
        """CKF预测"""
        n = len(mean)
        
        # 生成容积点
        try:
            L = np.linalg.cholesky(covariance)
        except np.linalg.LinAlgError:
            eigenvalues, eigenvectors = np.linalg.eigh(covariance)
            eigenvalues = np.maximum(eigenvalues, 1e-6)
            L = eigenvectors @ np.diag(np.sqrt(eigenvalues))
        
        # 2n个容积点
        xi = np.sqrt(n) * np.hstack([np.eye(n), -np.eye(n)])
        cubature_points = mean[:, np.newaxis] + L @ xi
        cubature_points = cubature_points.T
        
        # 传播容积点
        propagated_points = np.zeros_like(cubature_points)
        for i in range(2 * n):
            propagated_points[i] = self.f(cubature_points[i], dt)
        
        # 计算预测均值和协方差
        predicted_mean = np.mean(propagated_points, axis=0)
        diff = propagated_points - predicted_mean
        predicted_cov = diff.T @ diff / (2 * n) + self.Q
        
        return predicted_mean, predicted_cov
    
    def _update_single(self, z: np.ndarray, comp: GaussianComponent) -> Tuple[float, np.ndarray, np.ndarray]:
        """更新单个分量
        
        Args:
            z: 观测
            comp: 预测分量
            
        Returns:
            (似然, 更新均值, 更新协方差)
        """
        if self.filter_type == 'KF':
            # 线性卡尔曼滤波
            H = self.H_func(comp.mean)
            z_pred = H @ comp.mean
            S = self._ensure_symmetric_and_positive_definite(
                H @ comp.covariance @ H.T + self.R)
            K = comp.covariance @ H.T @ np.linalg.inv(S)
            innovation = z - z_pred
            updated_mean = comp.mean + K @ innovation
            updated_cov = (np.eye(self.state_dim) - K @ H) @ comp.covariance
        elif self.filter_type == 'EKF':
            # 扩展卡尔曼滤波
            z_pred = self.h(comp.mean)
            H = self.H_func(comp.mean)
            S = self._ensure_symmetric_and_positive_definite(
                H @ comp.covariance @ H.T + self.R)
            K = comp.covariance @ H.T @ np.linalg.inv(S)
            innovation = z - z_pred
            updated_mean = comp.mean + K @ innovation
            updated_cov = (np.eye(self.state_dim) - K @ H) @ comp.covariance
        elif self.filter_type == 'UKF':
            # 无迹卡尔曼滤波
            likelihood, updated_mean, updated_cov = self._ukf_update(z, comp.mean, comp.covariance)
            return likelihood, updated_mean, updated_cov
        elif self.filter_type == 'CKF':
            # 容积卡尔曼滤波
            likelihood, updated_mean, updated_cov = self._ckf_update(z, comp.mean, comp.covariance)
            return likelihood, updated_mean, updated_cov
        else:
            raise ValueError(f"Unknown filter type: {self.filter_type}")
        
        # 计算似然（对数空间）
        S_det = np.linalg.det(S)
        S_inv = np.linalg.inv(S)
        d = self.measurement_dim
        innovation = z - z_pred
        log_likelihood = -0.5 * d * np.log(2 * np.pi) - 0.5 * np.log(S_det) - 0.5 * innovation.T @ S_inv @ innovation
        likelihood = np.exp(log_likelihood)
        
        # 确保协方差对称
        updated_cov = (updated_cov + updated_cov.T) / 2
        
        return likelihood, updated_mean, updated_cov
    
    def _ukf_update(self, z: np.ndarray, mean: np.ndarray, covariance: np.ndarray) -> Tuple[float, np.ndarray, np.ndarray]:
        """UKF更新"""
        n = len(mean)
        alpha = 1e-3
        beta = 2.0
        kappa = 0.0
        lambda_ = alpha ** 2 * (n + kappa) - n
        
        # 生成Sigma点
        try:
            L = np.linalg.cholesky((n + lambda_) * covariance)
        except np.linalg.LinAlgError:
            eigenvalues, eigenvectors = np.linalg.eigh(covariance)
            eigenvalues = np.maximum(eigenvalues, 1e-6)
            L = eigenvectors @ np.diag(np.sqrt(eigenvalues * (n + lambda_)))
        
        sigma_points = np.zeros((2 * n + 1, n))
        sigma_points[0] = mean
        sigma_points[1:n+1] = mean + L.T
        sigma_points[n+1:] = mean - L.T
        
        # 计算权重
        Wm = np.zeros(2 * n + 1)
        Wm[0] = lambda_ / (n + lambda_)
        Wm[1:] = 1.0 / (2 * (n + lambda_))
        
        Wc = np.zeros(2 * n + 1)
        Wc[0] = lambda_ / (n + lambda_) + (1 - alpha ** 2 + beta)
        Wc[1:] = 1.0 / (2 * (n + lambda_))
        
        # 传播观测
        measurement_points = np.zeros((2 * n + 1, self.measurement_dim))
        for i in range(2 * n + 1):
            measurement_points[i] = self.h(sigma_points[i])
        
        # 计算预测观测
        z_pred = np.dot(Wm, measurement_points)
        
        # 计算协方差
        dz = measurement_points - z_pred
        dx = sigma_points - mean
        
        Pzz = dz.T @ (Wc[:, np.newaxis] * dz) + self.R
        Pxz = dx.T @ (Wc[:, np.newaxis] * dz)
        
        # 卡尔曼增益
        K = Pxz @ np.linalg.inv(Pzz)
        
        # 更新
        innovation = z - z_pred
        updated_mean = mean + K @ innovation
        updated_cov = covariance - K @ Pzz @ K.T
        
        # 计算似然
        Pzz_det = np.linalg.det(Pzz)
        Pzz_inv = np.linalg.inv(Pzz)
        d = self.measurement_dim
        log_likelihood = -0.5 * d * np.log(2 * np.pi) - 0.5 * np.log(Pzz_det) - 0.5 * innovation.T @ Pzz_inv @ innovation
        likelihood = np.exp(log_likelihood)
        
        return likelihood, updated_mean, updated_cov
    
    def _ckf_update(self, z: np.ndarray, mean: np.ndarray, covariance: np.ndarray) -> Tuple[float, np.ndarray, np.ndarray]:
        """CKF更新"""
        n = len(mean)
        
        # 生成容积点
        try:
            L = np.linalg.cholesky(covariance)
        except np.linalg.LinAlgError:
            eigenvalues, eigenvectors = np.linalg.eigh(covariance)
            eigenvalues = np.maximum(eigenvalues, 1e-6)
            L = eigenvectors @ np.diag(np.sqrt(eigenvalues))
        
        # 2n个容积点
        xi = np.sqrt(n) * np.hstack([np.eye(n), -np.eye(n)])
        cubature_points = mean[:, np.newaxis] + L @ xi
        cubature_points = cubature_points.T
        
        # 传播观测
        measurement_points = np.zeros((2 * n, self.measurement_dim))
        for i in range(2 * n):
            measurement_points[i] = self.h(cubature_points[i])
        
        # 计算预测观测
        z_pred = np.mean(measurement_points, axis=0)
        
        # 计算协方差
        dz = measurement_points - z_pred
        dx = cubature_points - mean
        
        Pzz = dz.T @ dz / (2 * n) + self.R
        Pxz = dx.T @ dz / (2 * n)
        
        # 卡尔曼增益
        K = Pxz @ np.linalg.inv(Pzz)
        
        # 更新
        innovation = z - z_pred
        updated_mean = mean + K @ innovation
        updated_cov = covariance - K @ Pzz @ K.T
        
        # 计算似然
        Pzz_det = np.linalg.det(Pzz)
        Pzz_inv = np.linalg.inv(Pzz)
        d = self.measurement_dim
        log_likelihood = -0.5 * d * np.log(2 * np.pi) - 0.5 * np.log(Pzz_det) - 0.5 * innovation.T @ Pzz_inv @ innovation
        likelihood = np.exp(log_likelihood)
        
        return likelihood, updated_mean, updated_cov
    
    def predict(self, dt: float) -> None:
        """PHD预测步骤
        
        1. 预测存活分量
        2. 添加固定出生分量
        
        Args:
            dt: 时间步长
        """
        # 更新过程噪声
        if self.filter_type == 'KF':
            self.Q = self._default_process_noise(dt)
        
        # 预测存活分量
        predicted_components = []
        for comp in self.components:
            predicted_comp = self._predict_single(comp, dt)
            predicted_components.append(predicted_comp)
        
        # 添加固定出生分量
        predicted_components.extend(self.birth_components)
        
        self.components = predicted_components
    
    def update(self, measurements: np.ndarray) -> None:
        """PHD更新步骤（向量化版本）
        
        1. 门限过滤观测
        2. 漏检分量
        3. 检测分量（向量化计算）
        4. 归一化
        
        Args:
            measurements: 观测矩阵，形状为 (n_measurements, measurement_dim)
        """
        if len(self.components) == 0:
            return
        
        n_meas = measurements.shape[0] if measurements.ndim > 1 else 0
        
        # 门限过滤观测
        if n_meas > 0:
            measurements = self._gate_measurements(measurements)
            n_meas = len(measurements)
        
        # 漏检分量
        updated_components = []
        for comp in self.components:
            updated_components.append(GaussianComponent(
                weight=self.Q_D * comp.weight,
                mean=comp.mean,
                covariance=comp.covariance
            ))
        
        # 检测分量
        if n_meas > 0:
            if self.filter_type in ('KF', 'EKF'):
                # KF/EKF: 向量化路径（使用雅可比预计算，高效）
                self._update_detection_kf_ekf(measurements, updated_components)
            elif self.filter_type in ('UKF', 'CKF'):
                # UKF/CKF: 逐分量路径（使用sigma/容积点，无雅可比）
                self._update_detection_ukf_ckf(measurements, updated_components)
        
        self.components = updated_components
        
        # 混合管理
        self._manage_components()
    
    @staticmethod
    def _ensure_symmetric_and_positive_definite(M: np.ndarray,
                                                  min_eig: float = 1e-3) -> np.ndarray:
        """确保矩阵对称正定，避免EKF/UKF线性化导致的数值不稳定

        当使用极坐标观测函数对远离传感器的分量做线性化时，
        H P H^T 可能产生不定矩阵——一阶泰勒展开在大协方差情况下
        无法准确捕获非线性观测函数的曲率。

        修复：对特征值做地板处理，保证最小特征值 >= min_eig。

        Args:
            M: 待修复矩阵
            min_eig: 最小允许特征值（默认 1e-3）

        Returns:
            对称正定的矩阵
        """
        M = (M + M.T) / 2
        try:
            eigenvalues = np.linalg.eigvalsh(M)
            min_actual = eigenvalues[0]
            if min_actual < min_eig:
                M = M + np.eye(len(M)) * (min_eig - min_actual)
        except np.linalg.LinAlgError:
            try:
                np.linalg.cholesky(M)
            except np.linalg.LinAlgError:
                M = M + np.eye(len(M)) * min_eig
        return M

    def _update_detection_kf_ekf(self, measurements: np.ndarray,
                                   updated_components: List[GaussianComponent]) -> None:
        """KF/EKF检测分量更新（向量化雅可比路径）

        Args:
            measurements: 门限过滤后的观测矩阵
            updated_components: 结果列表（原地修改）
        """
        n_comp = len(self.components)
        n_meas = len(measurements)

        # 预计算所有分量的预测观测和新息协方差
        z_preds = np.zeros((n_comp, self.measurement_dim))
        S_invs = np.zeros((n_comp, self.measurement_dim, self.measurement_dim))
        S_dets = np.zeros(n_comp)
        Ks = np.zeros((n_comp, self.state_dim, self.measurement_dim))
        Hs = np.zeros((n_comp, self.measurement_dim, self.state_dim))

        for i, comp in enumerate(self.components):
            # 预测观测
            z_preds[i] = self.h(comp.mean)

            # 新息协方差
            H = self.H_func(comp.mean)
            Hs[i] = H
            S = H @ comp.covariance @ H.T + self.R
            S = self._ensure_symmetric_and_positive_definite(S)

            # 使用Cholesky分解计算逆和行列式
            try:
                L = np.linalg.cholesky(S)
                S_inv = np.linalg.inv(L.T) @ np.linalg.inv(L)
                S_det = np.prod(np.diag(L)) ** 2
            except np.linalg.LinAlgError:
                S_inv = np.linalg.inv(S)
                S_det = np.linalg.det(S)

            S_invs[i] = S_inv
            S_dets[i] = S_det

            # 卡尔曼增益
            Ks[i] = comp.covariance @ H.T @ S_inv

        # 对每个观测更新所有分量
        for ell in range(n_meas):
            z = measurements[ell]
            detection_components = []

            for i, comp in enumerate(self.components):
                # 新息
                innovation = z - z_preds[i]

                # 似然（对数空间）
                d = self.measurement_dim
                log_likelihood = -0.5 * d * np.log(2 * np.pi) - 0.5 * np.log(S_dets[i]) \
                    - 0.5 * innovation.T @ S_invs[i] @ innovation
                likelihood = np.exp(log_likelihood)

                # 更新状态
                updated_mean = comp.mean + Ks[i] @ innovation
                updated_cov = (np.eye(self.state_dim) - Ks[i] @ Hs[i]) @ comp.covariance
                updated_cov = (updated_cov + updated_cov.T) / 2

                # 更新权重
                updated_weight = self.P_D * comp.weight * likelihood

                detection_components.append(GaussianComponent(
                    weight=updated_weight,
                    mean=updated_mean,
                    covariance=updated_cov
                ))

            # 归一化权重
            total_weight = sum(c.weight for c in detection_components)
            if total_weight > 0:
                for c in detection_components:
                    c.weight /= (self.clutter_density + total_weight)

            updated_components.extend(detection_components)

    def _update_detection_ukf_ckf(self, measurements: np.ndarray,
                                    updated_components: List[GaussianComponent]) -> None:
        """UKF/CKF检测分量更新（sigma点/容积点路径，无雅可比）

        Args:
            measurements: 门限过滤后的观测矩阵
            updated_components: 结果列表（原地修改）
        """
        for z in measurements:
            detection_components = []

            for comp in self.components:
                # 调用 _update_single 分派到 UKF/CKF 更新
                likelihood, updated_mean, updated_cov = self._update_single(z, comp)

                updated_weight = self.P_D * comp.weight * likelihood

                detection_components.append(GaussianComponent(
                    weight=updated_weight,
                    mean=updated_mean,
                    covariance=updated_cov
                ))

            # 归一化权重
            total_weight = sum(c.weight for c in detection_components)
            if total_weight > 0:
                for c in detection_components:
                    c.weight /= (self.clutter_density + total_weight)

            updated_components.extend(detection_components)

    def _gate_measurements(self, measurements: np.ndarray) -> np.ndarray:
        """门限过滤观测

        对于每个分量，计算观测与预测观测之间的马氏距离，保留距离小于阈值的观测。

        根据 filter_type 使用对应的新息协方差计算方法：
        - KF/EKF: 雅可比线性化 S = H P H^T + R
        - UKF: sigma点变换
        - CKF: 容积点变换

        Args:
            measurements: 观测矩阵

        Returns:
            过滤后的观测矩阵
        """
        if len(measurements) == 0:
            return measurements

        valid_indices = set()

        for comp in self.components:
            # 预测观测
            z_pred = self.h(comp.mean)

            # 新息协方差（根据滤波器类型）
            if self.filter_type in ('KF', 'EKF'):
                H = self.H_func(comp.mean)
                S = self._ensure_symmetric_and_positive_definite(
                    H @ comp.covariance @ H.T + self.R)
            elif self.filter_type == 'UKF':
                _, S = self._ukf_predict_measurement(comp.mean, comp.covariance)
            elif self.filter_type == 'CKF':
                _, S = self._ckf_predict_measurement(comp.mean, comp.covariance)
            else:
                H = self.H_func(comp.mean)
                S = self._ensure_symmetric_and_positive_definite(
                    H @ comp.covariance @ H.T + self.R)

            # 计算马氏距离
            try:
                S_inv = np.linalg.inv(S)
                for i, z in enumerate(measurements):
                    innovation = z - z_pred
                    distance = innovation.T @ S_inv @ innovation
                    if distance < self.gating_threshold:
                        valid_indices.add(i)
            except np.linalg.LinAlgError:
                continue

        if valid_indices:
            return measurements[list(valid_indices)]
        else:
            return np.array([]).reshape(0, self.measurement_dim)

    def _ukf_predict_measurement(self, mean: np.ndarray, covariance: np.ndarray) \
            -> Tuple[np.ndarray, np.ndarray]:
        """UKF预测测量均值和协方差（用于门限）"""
        n = len(mean)
        alpha = 1e-3
        beta = 2.0
        kappa = 0.0
        lambda_ = alpha ** 2 * (n + kappa) - n

        try:
            L = np.linalg.cholesky((n + lambda_) * covariance)
        except np.linalg.LinAlgError:
            eigenvalues, eigenvectors = np.linalg.eigh(covariance)
            eigenvalues = np.maximum(eigenvalues, 1e-6)
            L = eigenvectors @ np.diag(np.sqrt(eigenvalues * (n + lambda_)))

        sigma_points = np.zeros((2 * n + 1, n))
        sigma_points[0] = mean
        sigma_points[1:n+1] = mean + L.T
        sigma_points[n+1:] = mean - L.T

        Wm = np.zeros(2 * n + 1)
        Wc = np.zeros(2 * n + 1)
        Wm[0] = lambda_ / (n + lambda_)
        Wc[0] = lambda_ / (n + lambda_) + (1 - alpha ** 2 + beta)
        Wm[1:] = 1.0 / (2 * (n + lambda_))
        Wc[1:] = 1.0 / (2 * (n + lambda_))

        measurement_points = np.zeros((2 * n + 1, self.measurement_dim))
        for i in range(2 * n + 1):
            measurement_points[i] = self.h(sigma_points[i])

        z_pred = np.dot(Wm, measurement_points)
        dz = measurement_points - z_pred
        S = dz.T @ (Wc[:, np.newaxis] * dz) + self.R

        return z_pred, S

    def _ckf_predict_measurement(self, mean: np.ndarray, covariance: np.ndarray) \
            -> Tuple[np.ndarray, np.ndarray]:
        """CKF预测测量均值和协方差（用于门限）"""
        n = len(mean)

        try:
            L = np.linalg.cholesky(covariance)
        except np.linalg.LinAlgError:
            eigenvalues, eigenvectors = np.linalg.eigh(covariance)
            eigenvalues = np.maximum(eigenvalues, 1e-6)
            L = eigenvectors @ np.diag(np.sqrt(eigenvalues))

        xi = np.sqrt(n) * np.hstack([np.eye(n), -np.eye(n)])
        cubature_points = mean[:, np.newaxis] + L @ xi
        cubature_points = cubature_points.T

        measurement_points = np.zeros((2 * n, self.measurement_dim))
        for i in range(2 * n):
            measurement_points[i] = self.h(cubature_points[i])

        z_pred = np.mean(measurement_points, axis=0)
        dz = measurement_points - z_pred
        S = dz.T @ dz / (2 * n) + self.R

        return z_pred, S
    
    def _manage_components(self) -> None:
        """管理高斯分量：剪枝、合并、限制数量"""
        self.components = self._prune(self.components)
        self.components = self._merge(self.components)
        self.components = self._cap(self.components)
    
    def _prune(self, components: List[GaussianComponent]) -> List[GaussianComponent]:
        """剪枝：移除权重小于阈值的分量"""
        return [comp for comp in components if comp.weight > self.pruning_threshold]
    
    def _merge(self, components: List[GaussianComponent]) -> List[GaussianComponent]:
        """合并：合并相近的分量"""
        if not components:
            return []
        
        sorted_components = sorted(components, key=lambda c: c.weight, reverse=True)
        
        merged_components = []
        used = set()
        
        for i, comp_i in enumerate(sorted_components):
            if i in used:
                continue
            
            merge_set = [i]
            
            for j, comp_j in enumerate(sorted_components):
                if j <= i or j in used:
                    continue
                
                diff = comp_j.mean - comp_i.mean
                try:
                    inv_cov = np.linalg.inv(comp_i.covariance)
                    distance = diff.T @ inv_cov @ diff
                except np.linalg.LinAlgError:
                    distance = np.inf
                
                if distance <= self.merging_threshold:
                    merge_set.append(j)
                    used.add(j)
            
            # 合并分量
            merged_weight = sum(sorted_components[k].weight for k in merge_set)
            
            if merged_weight > 0:
                merged_mean = np.zeros(self.state_dim)
                for k in merge_set:
                    merged_mean += sorted_components[k].weight * sorted_components[k].mean
                merged_mean /= merged_weight
                
                merged_cov = np.zeros((self.state_dim, self.state_dim))
                for k in merge_set:
                    diff = sorted_components[k].mean - merged_mean
                    merged_cov += sorted_components[k].weight * (
                        sorted_components[k].covariance + np.outer(diff, diff)
                    )
                merged_cov /= merged_weight
                
                merged_components.append(GaussianComponent(
                    weight=merged_weight,
                    mean=merged_mean,
                    covariance=merged_cov
                ))
        
        return merged_components
    
    def _cap(self, components: List[GaussianComponent]) -> List[GaussianComponent]:
        """限制分量数量：保留权重最大的前 J_max 个分量

        参考 MATLAB gaus_cap.m：截断后重归一化权重以保持总 PHD 质量。
        这防止了因简单丢弃低权重分量而导致的基数漂移。
        """
        if len(components) <= self.max_components:
            return components

        sorted_components = sorted(components, key=lambda c: c.weight, reverse=True)
        kept = sorted_components[:self.max_components]

        # 重归一化：保留总 PHD 质量（参考 MATLAB gaus_cap）
        total_old = sum(c.weight for c in sorted_components)
        total_kept = sum(c.weight for c in kept)
        if total_kept > 0:
            scale = total_old / total_kept
            for c in kept:
                c.weight *= scale

        return kept

    def extract_states(self) -> Tuple[np.ndarray, int]:
        """提取目标状态

        标准方法（Vo & Ma 2006）：
        1. 总权重求和取整得到估计目标数 N = round(sum(w))
        2. 按权重比例将 N 个目标分配到各分量
        3. 每个分量分配到 round(w_i / sum(w) * N) 个目标

        Returns:
            (状态矩阵, 目标数量)
        """
        if not self.components:
            return np.array([]), 0

        total_weight = sum(comp.weight for comp in self.components)
        n_targets = round(total_weight)

        if n_targets <= 0:
            return np.array([]), 0

        states = []
        for comp in self.components:
            # 按权重比例分配
            n_comp = round(comp.weight / total_weight * n_targets)
            for _ in range(n_comp):
                states.append(comp.mean)

        if states:
            return np.array(states), len(states)
        return np.array([]), 0
    
    def get_expected_cardinality(self) -> float:
        """获取期望目标数"""
        return sum(comp.weight for comp in self.components)
    
    def get_components(self) -> List[GaussianComponent]:
        """获取所有分量"""
        return self.components.copy()
    
    def get_n_components(self) -> int:
        """获取分量数量"""
        return len(self.components)


class GMPHDFilter(PHDFilter):
    """高斯混合PHD滤波器（标准GM-PHD的别名）"""
    pass


class AdaptiveBirthPHDFilter(PHDFilter):
    """自适应出生PHD滤波器 (Vo 2012)

    参考：
    B. Ristic, D. Clark, B.-N. Vo, and B.-T. Vo, "Adaptive Target Birth
    Intensity for PHD and CPHD Filters," IEEE Trans. Aerospace and
    Electronic Systems, Vol. 48, No. 2, pp. 1656-1668, 2012.

    核心思想（Vo 2012 公式 15-17）：
    1. 对每个通过门限的量测 z，计算所有已知后验分量对其的解释程度
       surv_likelihood = Σ w_j * N(z; h(m_j), H_j P_j H_j^T + R)
    2. 若 surv_likelihood 小（无已知目标解释该量测），则可能是新生目标
       新生软权重: w_b = birth_rate / (clutter_density + surv_likelihood)
    3. 将新生分量预测到下一时刻后加入预测 PHD

    注意：此滤波器覆盖父类的 predict()，不再使用固定出生分量。
    """

    def __init__(self,
                 state_dim: int = 4,
                 measurement_dim: int = 2,
                 survival_probability: float = 0.99,
                 detection_probability: float = 0.98,
                 clutter_rate: float = 60.0,
                 surveillance_area: float = 1000000.0,
                 birth_weight: float = 0.03,
                 pruning_threshold: float = 1e-5,
                 merging_threshold: float = 4.0,
                 max_components: int = 100,
                 gating_threshold: float = 13.816,
                 # 状态转移函数
                 state_transition_func: Optional[Callable] = None,
                 state_transition_jacobian: Optional[Callable] = None,
                 # 观测函数
                 measurement_func: Optional[Callable] = None,
                 measurement_jacobian: Optional[Callable] = None,
                 # 噪声矩阵
                 process_noise_matrix: Optional[np.ndarray] = None,
                 measurement_noise_matrix: Optional[np.ndarray] = None,
                 # 滤波器类型
                 filter_type: str = 'KF',
                 # === 出生模型 ===
                 birth_components: Optional[List[GaussianComponent]] = None,
                 surveillance_bounds: Optional[Tuple[Tuple[float, float],
                                                      Tuple[float, float]]] = None,
                 birth_covariance: Optional[np.ndarray] = None,
                 n_birth_components: int = 4,
                 # 自适应新生速率 (Vo 2012: 期望每步新生目标数)
                 birth_rate: float = 0.2):
        """
        初始化自适应出生PHD滤波器

        Args:
            birth_rate: 期望每步新生目标数（Vo 2012 公式中的 birth_rate）。
                用于计算软权重 w_b = birth_rate / (clutter_density + surv_likelihood)。
                典型值：0.1-0.3。
            birth_covariance: 新生分量的协方差矩阵。
                默认 diag([10, 3.16, 10, 3.16])^2（与 MATLAB 一致）。
            其余参数同 PHDFilter。
        """
        self._birth_covariance = birth_covariance
        self.birth_rate = birth_rate

        super().__init__(
            state_dim=state_dim,
            measurement_dim=measurement_dim,
            survival_probability=survival_probability,
            detection_probability=detection_probability,
            clutter_rate=clutter_rate,
            surveillance_area=surveillance_area,
            birth_weight=birth_weight,
            pruning_threshold=pruning_threshold,
            merging_threshold=merging_threshold,
            max_components=max_components,
            gating_threshold=gating_threshold,
            state_transition_func=state_transition_func,
            state_transition_jacobian=state_transition_jacobian,
            measurement_func=measurement_func,
            measurement_jacobian=measurement_jacobian,
            process_noise_matrix=process_noise_matrix,
            measurement_noise_matrix=measurement_noise_matrix,
            filter_type=filter_type,
            birth_components=birth_components,
            surveillance_bounds=surveillance_bounds,
            birth_covariance=birth_covariance,
            n_birth_components=n_birth_components
        )

    def predict(self, dt: float, measurements: Optional[np.ndarray] = None) -> None:
        """PHD预测步骤（Vo 2012 自适应出生）

        使用当前时刻量测和上一时刻后验分量计算自适应出生强度，
        然后将新生分量预测到当前时刻后加入预测 PHD。

        Args:
            dt: 时间步长
            measurements: 当前时刻的观测矩阵
        """
        # === Step 1: 预测存活分量 ===
        predicted_components = []
        for comp in self.components:
            predicted_comp = self._predict_single(comp, dt)
            predicted_components.append(predicted_comp)

        # === Step 2: Vo 2012 自适应出生 ===
        if measurements is not None and len(measurements) > 0:
            adaptive_birth = self._compute_vo_adaptive_birth(
                measurements, self.components)
            # 将新生分量预测到当前时刻
            for comp in adaptive_birth:
                predicted_components.append(self._predict_single(comp, dt))
        else:
            # 无量测时回退到固定出生分量
            predicted_components.extend(self.birth_components)

        self.components = predicted_components

    def _compute_vo_adaptive_birth(
            self, measurements: np.ndarray,
            posterior_components: List[GaussianComponent]) -> List[GaussianComponent]:
        """自适应出生强度计算（门限筛选 + 均匀权重）

        策略（B.-N. Vo 早期实现）：
        1. 对每个分量计算预测量测和新息协方差
        2. 用量测与所有分量做马氏距离门限检测
        3. 未被任何分量门限"解释"的量测 → 新生候选
        4. 给候选量测分配均匀权重 → 出生分量

        这比 Vo 2012 的软权重公式更鲁棒，因为它不依赖 birth_rate/clutter_density
        的精确配比，而是直接通过门限判断量测是否被现有目标解释。

        Args:
            measurements: 当前时刻量测矩阵
            posterior_components: 上一时刻后验分量

        Returns:
            新生分量列表（时间对齐到当前时刻，需进一步预测）
        """
        n_meas = len(measurements)
        if n_meas == 0:
            return []

        # 对每个量测，检查是否落入任一分量的门限内
        is_explained = np.zeros(n_meas, dtype=bool)

        for comp in posterior_components:
            if comp.weight < self.pruning_threshold:
                continue
            try:
                z_pred = self.h(comp.mean)
                if self.filter_type in ('KF', 'EKF'):
                    H = self.H_func(comp.mean)
                    S = H @ comp.covariance @ H.T + self.R
                elif self.filter_type == 'UKF':
                    _, S = self._ukf_predict_measurement(
                        comp.mean, comp.covariance)
                elif self.filter_type == 'CKF':
                    _, S = self._ckf_predict_measurement(
                        comp.mean, comp.covariance)
                else:
                    continue
                S = self._ensure_symmetric_and_positive_definite(S)
                S_inv = np.linalg.inv(S)

                for i, z in enumerate(measurements):
                    if is_explained[i]:
                        continue
                    innovation = z - z_pred
                    mahal = innovation.T @ S_inv @ innovation
                    if mahal < self.gating_threshold:
                        is_explained[i] = True
            except (np.linalg.LinAlgError, ValueError):
                continue

        # 未被解释的量测 → 出生候选
        unexplained_idx = np.where(~is_explained)[0]
        n_birth = len(unexplained_idx)
        if n_birth == 0:
            return []

        # 限制每步最大出生分量数，防止早期步数爆炸
        max_birth = 10
        if n_birth > max_birth:
            # 随机选取或均匀采样
            step = n_birth / max_birth
            selected = [unexplained_idx[int(i * step)] for i in range(max_birth)]
            unexplained_idx = np.array(selected)
            n_birth = max_birth

        # 均匀权重：总出生权重平分（不超过总 weight 的合理比例）
        weight_per = min(self.birth_rate, 0.5) / n_birth

        # 出生协方差
        if self._birth_covariance is not None:
            birth_cov = self._birth_covariance.copy()
        else:
            birth_cov = np.diag([10.0, np.sqrt(10.0), 10.0, np.sqrt(10.0)]) ** 2

        birth_components = []
        for i in unexplained_idx:
            z = measurements[i]
            # 量测转状态：需要根据量测类型做坐标转换
            if self.state_dim == 4:
                if self.h is not self._default_measurement:
                    # 非线性量测（如极坐标）：将量测转换为笛卡尔位置
                    # 假设量测格式为 [range, bearing]
                    if self.measurement_dim == 2:
                        r, theta = z[0], z[1]
                        px = r * np.cos(theta)
                        py = r * np.sin(theta)
                        mean = np.array([px, 0.0, py, 0.0])
                    else:
                        mean = np.array([z[0], 0.0, z[1], 0.0])
                else:
                    # 线性量测：直接使用 z = [x, y]
                    mean = np.array([z[0], 0.0, z[1], 0.0])
            elif self.state_dim == 6:
                if self.h is not self._default_measurement:
                    r, theta = z[0], z[1]
                    px, py = r * np.cos(theta), r * np.sin(theta)
                    mean = np.array([px, 0.0, 0.0, py, 0.0, 0.0])
                else:
                    mean = np.array([z[0], 0.0, 0.0, z[1], 0.0, 0.0])
            else:
                mean = np.zeros(self.state_dim)
                mean[0] = z[0]
                if self.measurement_dim > 1 and self.state_dim > 2:
                    mean[2] = z[1]

            birth_components.append(GaussianComponent(
                weight=weight_per,
                mean=mean,
                covariance=birth_cov.copy()
            ))

        return birth_components
