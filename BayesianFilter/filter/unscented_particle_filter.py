"""
无迹粒子滤波器 (UPF)
"""
import numpy as np
from typing import Optional, Callable, List, Tuple
from .particle_filter import Particle, ParticleFilter
from .unscented_kalman_filter import UnscentedKalmanFilter
from .noise_models import get_process_noise_matrix


class UnscentedParticleFilter(ParticleFilter):
    """无迹粒子滤波器
    
    结合UKF和PF的优点：
    - 使用UKF生成更好的建议分布
    - 使用PF处理非高斯分布
    
    主要思想：
    1. 每个粒子运行一个UKF
    2. 使用UKF的预测均值和协方差生成建议分布
    3. 从建议分布中采样
    4. 计算重要性权重
    
    优点：
    - 比标准PF更准确
    - 能够处理强非线性
    - 减少粒子退化
    """
    
    def __init__(self,
                 state_dim: int = 4,
                 measurement_dim: int = 2,
                 n_particles: int = 500,
                 process_noise_std: float = 0.1,
                 measurement_noise_std: float = 1.0,
                 state_transition_func: Optional[Callable] = None,
                 measurement_func: Optional[Callable] = None,
                 likelihood_func: Optional[Callable] = None,
                 resampling_method: str = "systematic",
                 alpha: float = 1e-3,
                 beta: float = 2.0,
                 kappa: float = 0.0,
                 random_seed: Optional[int] = None):
        """
        初始化无迹粒子滤波器
        
        Args:
            state_dim: 状态维度
            measurement_dim: 观测维度
            n_particles: 粒子数量
            process_noise_std: 过程噪声标准差
            measurement_noise_std: 观测噪声标准差
            state_transition_func: 状态转移函数
            measurement_func: 观测函数
            likelihood_func: 似然函数
            resampling_method: 重采样方法
            alpha: UKF参数
            beta: UKF参数
            kappa: UKF参数
            random_seed: 随机种子
        """
        super().__init__(
            state_dim=state_dim,
            measurement_dim=measurement_dim,
            n_particles=n_particles,
            process_noise_std=process_noise_std,
            measurement_noise_std=measurement_noise_std,
            resampling_method=resampling_method,
            random_seed=random_seed
        )
        
        # 设置函数（使用默认值作为回退）
        if state_transition_func is not None:
            self.f = state_transition_func
        self.h = measurement_func if measurement_func is not None else self._default_measurement
        if likelihood_func is not None:
            self.likelihood = likelihood_func
        
        # UKF参数
        self.alpha = alpha
        self.beta = beta
        self.kappa = kappa
        
        # 为每个粒子创建UKF
        self.ukfs: List[UnscentedKalmanFilter] = []
        
        # 观测噪声协方差
        self.R = np.eye(measurement_dim) * measurement_noise_std ** 2
        
        # 数值稳定性常数
        self._EPS = 1e-12
        
        # 存储建议分布参数和先验信息（用于重要性权重计算）
        self._prev_states: List[Optional[np.ndarray]] = []   # x_{k-1} 上一时刻粒子状态
        self._ukf_pred_means: List[Optional[np.ndarray]] = []  # μ_pred UKF预测均值
        self._ukf_pred_covs: List[Optional[np.ndarray]] = []   # Σ_pred UKF预测协方差
        self._prev_dt: float = 0.0  # 上一预测步的时间步长

    def _create_ukf(self) -> UnscentedKalmanFilter:
        """创建UKF实例"""
        ukf = UnscentedKalmanFilter(
            state_dim=self.state_dim,
            measurement_dim=self.measurement_dim,
            process_noise_std=self.process_noise_std,
            measurement_noise_std=self.measurement_noise_std,
            alpha=self.alpha,
            beta=self.beta,
            kappa=self.kappa
        )
        
        # 设置自定义函数
        if hasattr(self, 'f') and self.f != self._default_state_transition:
            ukf.set_state_transition_function(self.f)
        if hasattr(self, 'h') and self.h != self._default_measurement:
            ukf.set_measurement_function(self.h)
        
        return ukf
    
    def _log_gaussian_pdf(self, x: np.ndarray, mean: np.ndarray,
                          cov: np.ndarray) -> float:
        """计算多元高斯对数概率密度函数值 log p(x | mean, cov)
        
        返回log域的值以避免高维或大马氏距离时的数值下溢。
        使用Cholesky分解求解马氏距离，确保数值稳定性。
        
        Args:
            x: 待评估的点
            mean: 高斯分布的均值向量
            cov: 高斯分布的协方差矩阵
            
        Returns:
            对数概率密度值 log p(x | mean, cov)
        """
        n = len(x)
        diff = x - mean
        
        # 确保协方差正定
        cov = self._ensure_positive_definite(cov)
        
        try:
            # Cholesky分解: cov = L @ L^T
            L = np.linalg.cholesky(cov)
            # 求解 L @ y = diff → y = L^{-1} @ diff
            y = np.linalg.solve(L, diff)
            # 马氏距离平方: diff^T @ cov^{-1} @ diff = y^T @ y
            mahalanobis = np.dot(y, y)
            # log|cov| = 2 * Σ log(L_{ii})
            log_det = 2 * np.sum(np.log(np.diag(L)))
        except np.linalg.LinAlgError:
            # Cholesky分解失败时回退到特征值分解
            eigenvalues, eigenvectors = np.linalg.eigh(cov)
            eigenvalues = np.maximum(eigenvalues, 1e-10)
            log_det = np.sum(np.log(eigenvalues))
            inv_cov = eigenvectors @ np.diag(1.0 / eigenvalues) @ eigenvectors.T
            mahalanobis = diff.T @ inv_cov @ diff
        
        # log p(x|mean,cov) = -0.5 * [n*log(2π) + log|cov| + (x-mean)^T cov^{-1} (x-mean)]
        return -0.5 * (n * np.log(2 * np.pi) + log_det + mahalanobis)
    
    def _default_state_transition(self, state: np.ndarray, dt: float) -> np.ndarray:
        """默认状态转移函数"""
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
        """默认观测函数"""
        if self.state_dim >= 4:
            return np.array([state[0], state[2]])
        elif self.state_dim >= 2:
            return np.array([state[0], state[1]])
        return state[:self.measurement_dim]
    
    def _default_likelihood(self, measurement: np.ndarray, 
                            state: np.ndarray) -> float:
        """默认似然函数"""
        # 预测观测
        z_pred = self.h(state)
        
        # 计算高斯似然
        diff = measurement - z_pred
        R = self.R
        
        n = self.measurement_dim
        R_det = np.linalg.det(R)
        R_inv = np.linalg.inv(R)
        
        exponent = -0.5 * diff.T @ R_inv @ diff
        coefficient = 1.0 / np.sqrt((2 * np.pi) ** n * R_det)
        
        return coefficient * np.exp(exponent)
    
    def _get_state_transition_matrix(self, dt: float) -> np.ndarray:
        """获取状态转移矩阵（用于基类接口）"""
        return np.eye(self.state_dim)
    
    def _get_process_noise_matrix(self, dt: float) -> np.ndarray:
        """获取过程噪声协方差矩阵"""
        return get_process_noise_matrix(self.state_dim, self.process_noise_std, dt)
    
    def _get_measurement_matrix(self) -> np.ndarray:
        """获取观测矩阵（用于基类接口）"""
        return np.eye(self.measurement_dim, self.state_dim)
    
    def initialize(self, 
                   initial_state: np.ndarray,
                   initial_covariance: Optional[np.ndarray] = None,
                   initial_particles: Optional[List[Particle]] = None) -> None:
        """初始化UPF
        
        Args:
            initial_state: 初始状态向量
            initial_covariance: 初始协方差矩阵
            initial_particles: 初始粒子集合
        """
        # 初始化UKF
        self.ukfs = []
        for _ in range(self.n_particles):
            ukf = self._create_ukf()
            self.ukfs.append(ukf)
        
        # 初始化粒子
        if initial_particles is not None:
            self.particles = initial_particles
        else:
            if initial_covariance is None:
                initial_covariance = np.eye(self.state_dim) * 100.0
            
            self.particles = []
            for i in range(self.n_particles):
                state = np.random.multivariate_normal(initial_state, initial_covariance)
                weight = 1.0 / self.n_particles
                self.particles.append(Particle(state, weight))
                
                # 初始化对应的UKF
                self.ukfs[i].initialize(state, initial_covariance)
        
        # 初始化存储列表（用于重要性权重计算）
        self._prev_states = [None] * self.n_particles
        self._ukf_pred_means = [None] * self.n_particles
        self._ukf_pred_covs = [None] * self.n_particles
        self._prev_dt = 0.0
        
        # 更新状态估计
        self._update_state_estimate()
        self.initialized = True
    
    def predict(self, dt: float) -> None:
        """状态预测
        
        使用UKF生成建议分布 q(x_k | x_{k-1}, z_k)。
        
        对于每个粒子：
        1. 运行UKF预测，得到建议分布 N(μ_pred, Σ_pred)
        2. 从建议分布中采样新粒子状态
        3. 保存建议分布参数和上一时刻状态，用于后续重要性权重计算
        
        Args:
            dt: 时间步长
        """
        if not self.initialized:
            raise RuntimeError("Filter not initialized")
        
        Q = self._get_process_noise_matrix(dt)
        self._prev_dt = dt
        
        # 确保存储列表长度与粒子数一致
        n = self.n_particles
        if len(self._prev_states) != n:
            self._prev_states = [None] * n
            self._ukf_pred_means = [None] * n
            self._ukf_pred_covs = [None] * n
        
        for i, particle in enumerate(self.particles):
            # 保存UKF当前状态 x_{k-1}（用于先验计算，与建议分布保持一致）
            # 注：UKF内部状态与粒子状态在update后可能不同，先验必须基于UKF的状态
            self._prev_states[i] = self.ukfs[i].get_state().copy() if self.ukfs[i].get_state() is not None else particle.state.copy()
            
            # 使用UKF预测
            self.ukfs[i].predict(dt)
            
            # 获取UKF的预测均值和协方差（建议分布参数）
            ukf_mean = self.ukfs[i].get_state()
            ukf_cov = self.ukfs[i].get_covariance()
            
            if ukf_mean is not None and ukf_cov is not None:
                # 保存建议分布参数 q(x_k | x_{k-1}, z_k) = N(μ_pred, Σ_pred)
                self._ukf_pred_means[i] = ukf_mean.copy()
                self._ukf_pred_covs[i] = ukf_cov.copy()
                
                # 从UKF的建议分布中采样新粒子状态 x_k ~ N(μ_pred, Σ_pred)
                particle.state = np.random.multivariate_normal(ukf_mean, ukf_cov)
            else:
                # UKF失败时，先验和建议分布相同（均使用状态转移模型）
                self._ukf_pred_means[i] = None
                self._ukf_pred_covs[i] = None
                noise = np.random.multivariate_normal(np.zeros(self.state_dim), Q)
                particle.state = self.f(particle.state, dt) + noise
    
    def update(self, measurement: np.ndarray,
               measurement_covariance: Optional[np.ndarray] = None) -> None:
        """量测更新
        
        重要性权重计算（核心UPF逻辑）:
        w_k ∝ p(z_k | x_k) × p(x_k | x_{k-1}) / q(x_k | x_{k-1}, z_k)
        
        其中:
          p(z_k | x_k)          = 似然 (likelihood) — 观测z_k在状态x_k下的概率
          p(x_k | x_{k-1})      = 先验 (prior) — 由状态转移模型 N(f(x_{k-1},dt), Q) 定义
          q(x_k | x_{k-1}, z_k) = 建议分布 (proposal) — 由UKF预测 N(μ_pred, Σ_pred) 定义
        
        注意:
        - 粒子状态x_k由predict()方法从建议分布中采样得到
        - UKF update仅用于更新UKF内部状态（供下一轮预测使用），不覆盖粒子状态
        - 重要性权重在log空间计算，使用log-sum-exp技巧归一化，避免数值下溢
        
        Args:
            measurement: 观测向量
            measurement_covariance: 观测噪声协方差
        """
        if not self.initialized:
            raise RuntimeError("Filter not initialized")
        
        z = np.asarray(measurement, dtype=np.float64)
        R = measurement_covariance if measurement_covariance is not None else self.R
        
        # log w_k = log p(z_k|x_k) + log p(x_k|x_{k-1}) - log q(x_k|x_{k-1},z_k)
        log_weights = np.zeros(self.n_particles)
        
        for i, particle in enumerate(self.particles):
            # 运行UKF update（更新UKF内部状态，用于下一轮predict）
            # 注意：不覆盖particle.state，保持predict()中从建议分布采样的结果
            self.ukfs[i].update(z, R)
            
            # ---------- 重要性权重计算 ----------
            # 1. log 似然: log p(z_k | x_k)
            likelihood = self.likelihood(z, particle.state)
            log_likelihood = np.log(likelihood + self._EPS)
            
            if (self._ukf_pred_means[i] is not None 
                    and self._ukf_pred_covs[i] is not None
                    and self._prev_states[i] is not None):
                
                # 2. log 先验: log p(x_k | x_{k-1}) = log N(x_k; f(x_{k-1}, dt), Q)
                prior_mean = self.f(self._prev_states[i], self._prev_dt)
                prior_cov = self._get_process_noise_matrix(self._prev_dt)
                log_prior = self._log_gaussian_pdf(
                    particle.state, prior_mean, prior_cov
                )
                
                # 3. log 建议分布: log q(x_k | x_{k-1}, z_k) = log N(x_k; μ_pred, Σ_pred)
                log_proposal = self._log_gaussian_pdf(
                    particle.state,
                    self._ukf_pred_means[i],
                    self._ukf_pred_covs[i]
                )
                
                # 组合: log w_k = log_likelihood + log_prior - log_proposal
                log_weights[i] = log_likelihood + log_prior - log_proposal
            else:
                # UKF失败时，先验和建议分布近似相等，权重 ≈ 似然
                log_weights[i] = log_likelihood
        
        # log-sum-exp 归一化: log w_i -= log Σ exp(log w_i)
        log_max = np.max(log_weights)
        log_sum_exp = log_max + np.log(
            np.sum(np.exp(log_weights - log_max))
        )
        log_weights -= log_sum_exp
        
        # 存储归一化后的权重（exp域）
        for i, particle in enumerate(self.particles):
            particle.weight = np.exp(log_weights[i])
        
        # 归一化权重（确保数值精度）
        self._normalize_weights()
        
        # 更新状态估计
        self._update_state_estimate()
        
        # 检查是否需要重采样
        if self._effective_sample_size() < self.n_particles / 2:
            self._resample()
    
    def predict_and_update(self, dt: float, measurement: np.ndarray,
                           measurement_covariance: Optional[np.ndarray] = None) -> None:
        """预测并更新
        
        Args:
            dt: 时间步长
            measurement: 观测向量
            measurement_covariance: 观测噪声协方差
        """
        self.predict(dt)
        self.update(measurement, measurement_covariance)
    
    def _resample(self) -> None:
        """重采样（重写父类方法以同步UKF）"""
        weights = np.array([p.weight for p in self.particles])
        weights = weights / np.sum(weights)
        
        # 系统重采样
        cumulative = np.cumsum(weights)
        u = np.random.uniform(0, 1.0 / self.n_particles)
        positions = u + np.arange(self.n_particles) / self.n_particles
        
        new_particles = []
        new_ukfs = []
        i, j = 0, 0
        
        while i < self.n_particles:
            if positions[i] < cumulative[j]:
                # 复制粒子
                new_particle = self.particles[j].copy()
                new_particle.weight = 1.0 / self.n_particles
                new_particles.append(new_particle)
                
                # 复制UKF
                new_ukf = self._create_ukf()
                new_ukf.initialize(
                    self.ukfs[j].get_state(),
                    self.ukfs[j].get_covariance()
                )
                new_ukfs.append(new_ukf)
                
                i += 1
            else:
                j += 1
        
        self.particles = new_particles
        self.ukfs = new_ukfs
    
    def get_ukf_estimates(self) -> List[Tuple[np.ndarray, np.ndarray]]:
        """获取所有UKF的估计
        
        Returns:
            列表，每个元素为 (状态估计, 协方差矩阵)
        """
        estimates = []
        for ukf in self.ukfs:
            state = ukf.get_state()
            cov = ukf.get_covariance()
            if state is not None and cov is not None:
                estimates.append((state, cov))
        return estimates
    
    def set_state_transition_function(self, func: Callable) -> None:
        """设置状态转移函数"""
        self.f = func
        for ukf in self.ukfs:
            ukf.set_state_transition_function(func)
    
    def set_measurement_function(self, func: Callable) -> None:
        """设置观测函数"""
        self.h = func
        for ukf in self.ukfs:
            ukf.set_measurement_function(func)


class SquareRootUnscentedParticleFilter(UnscentedParticleFilter):
    """平方根无迹粒子滤波器
    
    使用平方根UKF代替标准UKF，提高数值稳定性
    """
    
    def __init__(self,
                 state_dim: int = 4,
                 measurement_dim: int = 2,
                 n_particles: int = 500,
                 process_noise_std: float = 0.1,
                 measurement_noise_std: float = 1.0,
                 state_transition_func: Optional[Callable] = None,
                 measurement_func: Optional[Callable] = None,
                 likelihood_func: Optional[Callable] = None,
                 resampling_method: str = "systematic",
                 alpha: float = 1e-3,
                 beta: float = 2.0,
                 kappa: float = 0.0,
                 random_seed: Optional[int] = None):
        """
        初始化平方根无迹粒子滤波器
        """
        super().__init__(
            state_dim=state_dim,
            measurement_dim=measurement_dim,
            n_particles=n_particles,
            process_noise_std=process_noise_std,
            measurement_noise_std=measurement_noise_std,
            state_transition_func=state_transition_func,
            measurement_func=measurement_func,
            likelihood_func=likelihood_func,
            resampling_method=resampling_method,
            alpha=alpha,
            beta=beta,
            kappa=kappa,
            random_seed=random_seed
        )
    
    def _create_ukf(self) -> UnscentedKalmanFilter:
        """创建平方根UKF实例
        
        注：这里仍然使用标准UKF，但可以通过修改来实现真正的平方根UKF
        """
        # 目前使用标准UKF
        return super()._create_ukf()
