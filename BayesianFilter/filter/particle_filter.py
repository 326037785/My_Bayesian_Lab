"""
粒子滤波器 (PF)
"""
import numpy as np
from typing import Optional, Callable, List, Tuple
from .base_filter import BaseFilter
from .noise_models import get_process_noise_matrix


class Particle:
    """粒子类"""
    
    def __init__(self, state: np.ndarray, weight: float = 1.0):
        """
        初始化粒子
        
        Args:
            state: 粒子状态
            weight: 粒子权重
        """
        self.state = np.asarray(state, dtype=np.float64)
        self.weight = weight
    
    def copy(self) -> 'Particle':
        """复制粒子"""
        return Particle(self.state.copy(), self.weight)


class ParticleFilter(BaseFilter):
    """标准粒子滤波器
    
    使用序贯蒙特卡洛方法进行状态估计
    
    主要步骤：
    1. 预测：根据状态转移模型传播粒子
    2. 更新：根据观测计算粒子权重
    3. 重采样：根据权重重采样粒子
    """
    
    def __init__(self,
                 state_dim: int = 4,
                 measurement_dim: int = 2,
                 n_particles: int = 1000,
                 process_noise_std: float = 0.1,
                 measurement_noise_std: float = 1.0,
                 state_transition_func: Optional[Callable] = None,
                 likelihood_func: Optional[Callable] = None,
                 resampling_method: str = "systematic",
                 random_seed: Optional[int] = None):
        """
        初始化粒子滤波器
        
        Args:
            state_dim: 状态维度
            measurement_dim: 观测维度
            n_particles: 粒子数量
            process_noise_std: 过程噪声标准差
            measurement_noise_std: 观测噪声标准差
            state_transition_func: 状态转移函数 f(x, dt, noise)
            likelihood_func: 似然函数 p(z|x)
            resampling_method: 重采样方法 ("multinomial", "systematic", "residual")
            random_seed: 随机种子
        """
        super().__init__(state_dim, measurement_dim, process_noise_std)
        
        self.n_particles = n_particles
        self.measurement_noise_std = measurement_noise_std
        self.resampling_method = resampling_method
        
        # 设置函数
        self.f = state_transition_func if state_transition_func is not None else self._default_state_transition
        self.likelihood = likelihood_func if likelihood_func is not None else self._default_likelihood
        
        # 粒子集合
        self.particles: List[Particle] = []
        
        # 预计算似然函数的缓存（避免重复矩阵求逆）
        self._R = np.eye(self.measurement_dim) * self.measurement_noise_std ** 2
        self._R_inv = np.linalg.inv(self._R)
        self._R_det = np.linalg.det(self._R)
        self._likelihood_coefficient = 1.0 / np.sqrt((2 * np.pi) ** self.measurement_dim * self._R_det)
        
        # 随机种子
        if random_seed is not None:
            np.random.seed(random_seed)
    
    def _default_state_transition(self, state: np.ndarray, dt: float, 
                                   noise: np.ndarray) -> np.ndarray:
        """默认状态转移函数（匀速模型）
        
        Args:
            state: 当前状态
            dt: 时间步长
            noise: 过程噪声
            
        Returns:
            下一时刻状态
        """
        if self.state_dim == 4:
            x, vx, y, vy = state[0], state[1], state[2], state[3]
            new_state = np.array([
                x + vx * dt,
                vx,
                y + vy * dt,
                vy
            ])
        elif self.state_dim == 6:
            x, vx, ax, y, vy, ay = state
            new_state = np.array([
                x + vx * dt + 0.5 * ax * dt**2,
                vx + ax * dt,
                ax,
                y + vy * dt + 0.5 * ay * dt**2,
                vy + ay * dt,
                ay
            ])
        else:
            new_state = state.copy()
        
        return new_state + noise
    
    def _default_likelihood(self, measurement: np.ndarray, 
                            state: np.ndarray) -> float:
        """默认似然函数（高斯似然）
        
        Args:
            measurement: 观测向量
            state: 状态向量
            
        Returns:
            似然值
        """
        # 提取位置
        if self.state_dim >= 4:
            predicted_meas = np.array([state[0], state[2]])
        else:
            predicted_meas = state[:self.measurement_dim]
        
        # 计算高斯似然（使用预计算的缓存值）
        diff = measurement - predicted_meas
        exponent = -0.5 * diff.T @ self._R_inv @ diff
        
        return self._likelihood_coefficient * np.exp(exponent)
    
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
        """初始化粒子滤波器
        
        Args:
            initial_state: 初始状态向量
            initial_covariance: 初始协方差矩阵
            initial_particles: 初始粒子集合（可选）
        """
        if initial_particles is not None:
            self.particles = initial_particles
        else:
            # 从高斯分布采样初始粒子
            if initial_covariance is None:
                initial_covariance = np.eye(self.state_dim) * 100.0
            
            self.particles = []
            for _ in range(self.n_particles):
                state = np.random.multivariate_normal(initial_state, initial_covariance)
                weight = 1.0 / self.n_particles
                self.particles.append(Particle(state, weight))
        
        # 计算初始状态估计
        self._update_state_estimate()
        self.initialized = True
    
    def _update_state_estimate(self) -> None:
        """从粒子集合更新状态估计"""
        if not self.particles:
            return
        
        # 加权平均
        states = np.array([p.state for p in self.particles])
        weights = np.array([p.weight for p in self.particles])
        
        # 归一化权重
        weights = weights / np.sum(weights)
        
        # 计算加权均值
        self.state = np.zeros(self.state_dim)
        for i in range(len(self.particles)):
            self.state += weights[i] * states[i]
        
        # 计算加权协方差
        self.covariance = np.zeros((self.state_dim, self.state_dim))
        for i in range(len(self.particles)):
            diff = states[i] - self.state
            self.covariance += weights[i] * np.outer(diff, diff)
    
    def predict(self, dt: float) -> None:
        """状态预测
        
        根据状态转移模型传播粒子
        
        Args:
            dt: 时间步长
        """
        if not self.initialized:
            raise RuntimeError("Filter not initialized")
        
        # 获取过程噪声协方差
        Q = self._get_process_noise_matrix(dt)
        
        # 传播每个粒子
        for particle in self.particles:
            # 生成过程噪声
            noise = np.random.multivariate_normal(np.zeros(self.state_dim), Q)
            
            # 状态转移
            particle.state = self.f(particle.state, dt, noise)
    
    def update(self, measurement: np.ndarray,
               measurement_covariance: Optional[np.ndarray] = None) -> None:
        """量测更新
        
        根据观测更新粒子权重
        
        Args:
            measurement: 观测向量
            measurement_covariance: 观测噪声协方差（可选）
        """
        if not self.initialized:
            raise RuntimeError("Filter not initialized")
        
        z = np.asarray(measurement, dtype=np.float64)
        
        # 更新每个粒子的权重
        for particle in self.particles:
            # 计算似然
            likelihood = self.likelihood(z, particle.state)
            particle.weight *= likelihood
        
        # 归一化权重
        self._normalize_weights()
        
        # 更新状态估计
        self._update_state_estimate()
        
        # 检查是否需要重采样
        if self._effective_sample_size() < self.n_particles / 2:
            self._resample()
        
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
    
    def _normalize_weights(self) -> None:
        """归一化粒子权重"""
        total_weight = sum(p.weight for p in self.particles)
        if total_weight > 0:
            for particle in self.particles:
                particle.weight /= total_weight
        else:
            # 如果所有权重为0，重新均匀分配
            for particle in self.particles:
                particle.weight = 1.0 / self.n_particles
    
    def _effective_sample_size(self) -> float:
        """计算有效样本数
        
        Returns:
            有效样本数 N_eff
        """
        weights = np.array([p.weight for p in self.particles])
        return 1.0 / np.sum(weights ** 2)
    
    def _resample(self) -> None:
        """重采样粒子"""
        if self.resampling_method == "multinomial":
            self._multinomial_resample()
        elif self.resampling_method == "systematic":
            self._systematic_resample()
        elif self.resampling_method == "residual":
            self._residual_resample()
        else:
            self._systematic_resample()
    
    def _multinomial_resample(self) -> None:
        """多项式重采样"""
        weights = np.array([p.weight for p in self.particles])
        weights = weights / np.sum(weights)
        
        # 采样索引
        indices = np.random.choice(
            self.n_particles, 
            size=self.n_particles, 
            p=weights
        )
        
        # 创建新粒子集合
        new_particles = []
        for idx in indices:
            new_particle = self.particles[idx].copy()
            new_particle.weight = 1.0 / self.n_particles
            new_particles.append(new_particle)
        
        self.particles = new_particles
    
    def _systematic_resample(self) -> None:
        """系统重采样"""
        weights = np.array([p.weight for p in self.particles])
        weights = weights / np.sum(weights)
        
        # 计算累积分布
        cumulative = np.cumsum(weights)
        
        # 生成均匀分布样本
        u = np.random.uniform(0, 1.0 / self.n_particles)
        positions = u + np.arange(self.n_particles) / self.n_particles
        
        # 重采样
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
    
    def _residual_resample(self) -> None:
        """残差重采样"""
        weights = np.array([p.weight for p in self.particles])
        weights = weights / np.sum(weights)
        
        # 计算每个粒子的复制次数
        n_copies = np.floor(self.n_particles * weights).astype(int)
        
        # 确定残差
        residual = weights - n_copies / self.n_particles
        residual = residual / np.sum(residual)
        
        # 第一阶段：确定性复制
        new_particles = []
        for i in range(self.n_particles):
            for _ in range(n_copies[i]):
                new_particle = self.particles[i].copy()
                new_particle.weight = 1.0 / self.n_particles
                new_particles.append(new_particle)
        
        # 第二阶段：多项式采样残差
        n_remaining = self.n_particles - len(new_particles)
        if n_remaining > 0:
            indices = np.random.choice(
                self.n_particles,
                size=n_remaining,
                p=residual
            )
            for idx in indices:
                new_particle = self.particles[idx].copy()
                new_particle.weight = 1.0 / self.n_particles
                new_particles.append(new_particle)
        
        self.particles = new_particles
    
    def get_particles(self) -> List[Particle]:
        """获取粒子集合"""
        return self.particles.copy()
    
    def get_particle_states(self) -> np.ndarray:
        """获取所有粒子的状态
        
        Returns:
            粒子状态数组，形状为 (n_particles, state_dim)
        """
        return np.array([p.state for p in self.particles])
    
    def get_particle_weights(self) -> np.ndarray:
        """获取所有粒子的权重
        
        Returns:
            粒子权重数组
        """
        return np.array([p.weight for p in self.particles])
    
    def get_position(self) -> np.ndarray:
        """获取位置估计（重写基类方法）"""
        if self.state is None:
            return np.zeros(2)
        
        if len(self.state) >= 4:
            return np.array([self.state[0], self.state[2]])
        elif len(self.state) >= 2:
            return np.array([self.state[0], self.state[1]])
        return self.state[:2]
    
    def get_position_uncertainty(self) -> float:
        """获取位置不确定性"""
        if self.covariance is None:
            return 0.0
        
        if self.state_dim >= 4:
            pos_cov = self.covariance[[0, 2], :][:, [0, 2]]
        else:
            pos_cov = self.covariance[:2, :2]
        
        return np.sqrt(np.trace(pos_cov))
    
    def set_state_transition_function(self, func: Callable) -> None:
        """设置状态转移函数"""
        self.f = func
    
    def set_likelihood_function(self, func: Callable) -> None:
        """设置似然函数"""
        self.likelihood = func
    
    def set_n_particles(self, n: int) -> None:
        """设置粒子数量"""
        self.n_particles = n
