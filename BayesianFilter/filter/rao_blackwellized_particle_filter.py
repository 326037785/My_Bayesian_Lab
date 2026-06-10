"""
Rao-Blackwellized粒子滤波器 (RBPF)
"""
import numpy as np
from typing import Optional, Callable, List, Tuple
from .particle_filter import Particle, ParticleFilter
from .kalman_filter import KalmanFilter
from .noise_models import get_process_noise_matrix


class RaoBlackwellizedParticle(Particle):
    """Rao-Blackwellized粒子
    
    每个粒子包含：
    - 状态的线性部分（用卡尔曼滤波器维护）
    - 状态的非线性部分（用粒子表示）
    """
    
    def __init__(self, 
                 nonlinear_state: np.ndarray,
                 linear_state: np.ndarray,
                 linear_covariance: np.ndarray,
                 weight: float = 1.0):
        """
        初始化Rao-Blackwellized粒子
        
        Args:
            nonlinear_state: 非线性状态部分
            linear_state: 线性状态部分
            linear_covariance: 线性状态协方差
            weight: 粒子权重
        """
        super().__init__(nonlinear_state, weight)
        self.linear_state = linear_state.copy()
        self.linear_covariance = linear_covariance.copy()
    
    def copy(self) -> 'RaoBlackwellizedParticle':
        """复制粒子"""
        return RaoBlackwellizedParticle(
            self.state.copy(),
            self.linear_state.copy(),
            self.linear_covariance.copy(),
            self.weight
        )


class RaoBlackwellizedParticleFilter(ParticleFilter):
    """Rao-Blackwellized粒子滤波器
    
    将状态分解为线性和非线性部分：
    - 非线性部分：用粒子滤波器估计
    - 线性部分：用卡尔曼滤波器估计
    
    优点：
    - 降低粒子滤波器的维度
    - 提高估计精度
    - 减少粒子数量需求
    
    适用场景：
    - 状态可以分解为线性和非线性部分
    - 线性部分的维度较高
    """
    
    def __init__(self,
                 nonlinear_state_dim: int = 2,
                 linear_state_dim: int = 2,
                 measurement_dim: int = 2,
                 n_particles: int = 500,
                 process_noise_std: float = 0.1,
                 measurement_noise_std: float = 1.0,
                 nonlinear_transition_func: Optional[Callable] = None,
                 linear_transition_func: Optional[Callable] = None,
                 measurement_func: Optional[Callable] = None,
                 resampling_method: str = "systematic",
                 random_seed: Optional[int] = None):
        """
        初始化Rao-Blackwellized粒子滤波器
        
        Args:
            nonlinear_state_dim: 非线性状态维度
            linear_state_dim: 线性状态维度
            measurement_dim: 观测维度
            n_particles: 粒子数量
            process_noise_std: 过程噪声标准差
            measurement_noise_std: 观测噪声标准差
            nonlinear_transition_func: 非线性状态转移函数
            linear_transition_func: 线性状态转移函数
            measurement_func: 观测函数
            resampling_method: 重采样方法
            random_seed: 随机种子
        """
        # 总状态维度
        total_state_dim = nonlinear_state_dim + linear_state_dim
        
        super().__init__(
            state_dim=total_state_dim,
            measurement_dim=measurement_dim,
            n_particles=n_particles,
            process_noise_std=process_noise_std,
            measurement_noise_std=measurement_noise_std,
            resampling_method=resampling_method,
            random_seed=random_seed
        )
        
        self.nonlinear_state_dim = nonlinear_state_dim
        self.linear_state_dim = linear_state_dim
        
        # 设置函数
        self.f_nonlinear = nonlinear_transition_func if nonlinear_transition_func is not None else self._default_nonlinear_transition
        self.f_linear = linear_transition_func if linear_transition_func is not None else self._default_linear_transition
        self.h = measurement_func if measurement_func is not None else self._default_measurement
        
        # 线性部分的卡尔曼滤波器
        self.kf = KalmanFilter(
            state_dim=linear_state_dim,
            measurement_dim=measurement_dim,
            process_noise_std=process_noise_std,
            measurement_noise_std=measurement_noise_std
        )
    
    def _default_nonlinear_transition(self, nonlinear_state: np.ndarray, 
                                       dt: float) -> np.ndarray:
        """默认非线性状态转移函数"""
        # 简单的非线性转移（例如：转弯速率）
        return nonlinear_state
    
    def _default_linear_transition(self, linear_state: np.ndarray,
                                    nonlinear_state: np.ndarray,
                                    dt: float) -> Tuple[np.ndarray, np.ndarray]:
        """默认线性状态转移函数
        
        Returns:
            (状态转移矩阵F, 过程噪声矩阵Q)
        """
        # 假设线性部分是位置和速度
        F = np.array([
            [1, dt],
            [0, 1]
        ])
        
        q = self.process_noise_std ** 2
        Q = q * np.array([
            [dt**3/3, dt**2/2],
            [dt**2/2, dt]
        ])
        
        return F, Q
    
    def _default_measurement(self, state: np.ndarray) -> np.ndarray:
        """默认观测函数"""
        # 假设只观测线性部分
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
    
    def initialize(self, 
                   initial_state: np.ndarray,
                   initial_covariance: Optional[np.ndarray] = None,
                   initial_particles: Optional[List[RaoBlackwellizedParticle]] = None) -> None:
        """初始化RBPF
        
        Args:
            initial_state: 初始状态向量 [nonlinear_state, linear_state]
            initial_covariance: 初始协方差矩阵
            initial_particles: 初始粒子集合
        """
        if initial_particles is not None:
            self.particles = initial_particles
        else:
            # 分解初始状态
            nonlinear_state = initial_state[:self.nonlinear_state_dim]
            linear_state = initial_state[self.nonlinear_state_dim:]
            
            # 初始协方差
            if initial_covariance is not None:
                linear_cov = initial_covariance[self.nonlinear_state_dim:, 
                                               self.nonlinear_state_dim:]
            else:
                linear_cov = np.eye(self.linear_state_dim) * 100.0
            
            # 生成粒子
            self.particles = []
            for _ in range(self.n_particles):
                # 采样非线性状态
                nl_state = np.random.multivariate_normal(
                    nonlinear_state,
                    np.eye(self.nonlinear_state_dim) * 10.0
                )
                
                # 创建粒子
                particle = RaoBlackwellizedParticle(
                    nonlinear_state=nl_state,
                    linear_state=linear_state,
                    linear_covariance=linear_cov,
                    weight=1.0 / self.n_particles
                )
                self.particles.append(particle)
        
        # 更新状态估计
        self._update_state_estimate()
        self.initialized = True
    
    def _update_state_estimate(self) -> None:
        """从粒子集合更新状态估计"""
        if not self.particles:
            return
        
        # 加权平均
        nonlinear_states = []
        linear_states = []
        weights = []
        
        for particle in self.particles:
            if isinstance(particle, RaoBlackwellizedParticle):
                nonlinear_states.append(particle.state)
                linear_states.append(particle.linear_state)
                weights.append(particle.weight)
        
        weights = np.array(weights)
        weights = weights / np.sum(weights)
        
        # 计算加权均值
        nonlinear_mean = np.zeros(self.nonlinear_state_dim)
        linear_mean = np.zeros(self.linear_state_dim)
        
        for i in range(len(self.particles)):
            nonlinear_mean += weights[i] * nonlinear_states[i]
            linear_mean += weights[i] * linear_states[i]
        
        # 合并状态
        self.state = np.concatenate([nonlinear_mean, linear_mean])
        
        # 计算协方差（简化处理）
        self.covariance = np.eye(self.state_dim) * 10.0
    
    def predict(self, dt: float) -> None:
        """状态预测
        
        Args:
            dt: 时间步长
        """
        if not self.initialized:
            raise RuntimeError("Filter not initialized")
        
        Q = self._get_process_noise_matrix(dt)
        
        for particle in self.particles:
            if isinstance(particle, RaoBlackwellizedParticle):
                # 非线性部分：粒子传播
                noise_nl = np.random.multivariate_normal(
                    np.zeros(self.nonlinear_state_dim),
                    Q[:self.nonlinear_state_dim, :self.nonlinear_state_dim]
                )
                particle.state = self.f_nonlinear(particle.state, dt) + noise_nl
                
                # 线性部分：卡尔曼滤波器预测
                F, Q_linear = self.f_linear(particle.linear_state, particle.state, dt)
                
                # 状态预测
                particle.linear_state = F @ particle.linear_state
                
                # 协方差预测
                particle.linear_covariance = F @ particle.linear_covariance @ F.T + Q_linear
                
                # 确保协方差正定
                particle.linear_covariance = self._ensure_positive_definite(particle.linear_covariance)
    
    def update(self, measurement: np.ndarray,
               measurement_covariance: Optional[np.ndarray] = None) -> None:
        """量测更新
        
        Args:
            measurement: 观测向量
            measurement_covariance: 观测噪声协方差
        """
        if not self.initialized:
            raise RuntimeError("Filter not initialized")
        
        z = np.asarray(measurement, dtype=np.float64)
        R = measurement_covariance if measurement_covariance is not None else self.R
        
        for particle in self.particles:
            if isinstance(particle, RaoBlackwellizedParticle):
                # 线性部分：卡尔曼滤波器更新
                # 预测观测
                H = np.eye(self.measurement_dim, self.linear_state_dim)
                z_pred = H @ particle.linear_state
                
                # 新息
                innovation = z - z_pred
                
                # 新息协方差
                S = H @ particle.linear_covariance @ H.T + R
                
                # 卡尔曼增益
                K = particle.linear_covariance @ H.T @ np.linalg.inv(S)
                
                # 状态更新
                particle.linear_state = particle.linear_state + K @ innovation
                
                # 协方差更新
                I_KH = np.eye(self.linear_state_dim) - K @ H
                particle.linear_covariance = I_KH @ particle.linear_covariance @ I_KH.T + K @ R @ K.T
                
                # 确保协方差正定
                particle.linear_covariance = self._ensure_positive_definite(particle.linear_covariance)
                
                # 计算似然（用于权重更新）
                S_det = np.linalg.det(S)
                S_inv = np.linalg.inv(S)
                
                exponent = -0.5 * innovation.T @ S_inv @ innovation
                coefficient = 1.0 / np.sqrt((2 * np.pi) ** self.measurement_dim * S_det)
                
                likelihood = coefficient * np.exp(exponent)
                
                # 更新权重
                particle.weight *= likelihood
        
        # 归一化权重
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
        """重采样（重写父类方法以处理RBPF粒子）"""
        weights = np.array([p.weight for p in self.particles])
        weights = weights / np.sum(weights)
        
        # 系统重采样
        cumulative = np.cumsum(weights)
        u = np.random.uniform(0, 1.0 / self.n_particles)
        positions = u + np.arange(self.n_particles) / self.n_particles
        
        new_particles = []
        i, j = 0, 0
        while i < self.n_particles:
            if positions[i] < cumulative[j]:
                new_particle = self.particles[j].copy()
                new_particle.weight = 1.0 / self.n_particles
                new_particles.append(new_particle)
                i += 1
            else:
                j += 1
        
        self.particles = new_particles
    
    def get_linear_state_estimates(self) -> Tuple[np.ndarray, np.ndarray]:
        """获取线性状态估计
        
        Returns:
            (线性状态均值, 线性状态协方差)
        """
        if not self.particles:
            return np.zeros(self.linear_state_dim), np.eye(self.linear_state_dim)
        
        linear_states = []
        weights = []
        
        for particle in self.particles:
            if isinstance(particle, RaoBlackwellizedParticle):
                linear_states.append(particle.linear_state)
                weights.append(particle.weight)
        
        weights = np.array(weights)
        weights = weights / np.sum(weights)
        
        # 加权均值
        linear_mean = np.zeros(self.linear_state_dim)
        for i in range(len(linear_states)):
            linear_mean += weights[i] * linear_states[i]
        
        # 加权协方差
        linear_cov = np.zeros((self.linear_state_dim, self.linear_state_dim))
        for i in range(len(linear_states)):
            diff = linear_states[i] - linear_mean
            linear_cov += weights[i] * np.outer(diff, diff)
        
        return linear_mean, linear_cov
    
    def get_nonlinear_state_estimates(self) -> Tuple[np.ndarray, np.ndarray]:
        """获取非线性状态估计
        
        Returns:
            (非线性状态均值, 非线性状态协方差)
        """
        if not self.particles:
            return np.zeros(self.nonlinear_state_dim), np.eye(self.nonlinear_state_dim)
        
        nonlinear_states = []
        weights = []
        
        for particle in self.particles:
            if isinstance(particle, RaoBlackwellizedParticle):
                nonlinear_states.append(particle.state)
                weights.append(particle.weight)
        
        weights = np.array(weights)
        weights = weights / np.sum(weights)
        
        # 加权均值
        nonlinear_mean = np.zeros(self.nonlinear_state_dim)
        for i in range(len(nonlinear_states)):
            nonlinear_mean += weights[i] * nonlinear_states[i]
        
        # 加权协方差
        nonlinear_cov = np.zeros((self.nonlinear_state_dim, self.nonlinear_state_dim))
        for i in range(len(nonlinear_states)):
            diff = nonlinear_states[i] - nonlinear_mean
            nonlinear_cov += weights[i] * np.outer(diff, diff)
        
        return nonlinear_mean, nonlinear_cov
