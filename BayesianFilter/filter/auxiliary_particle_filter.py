"""
辅助粒子滤波器 (APF)
"""
import numpy as np
from typing import Optional, Callable, List
from .particle_filter import Particle, ParticleFilter


class AuxiliaryParticleFilter(ParticleFilter):
    """辅助粒子滤波器
    
    改进的粒子滤波器，通过引入辅助变量来减少权重方差
    
    主要思想：
    1. 预采样阶段：使用辅助变量选择"有前途"的粒子
    2. 重采样阶段：根据辅助权重重采样
    3. 传播阶段：传播重采样后的粒子
    4. 更新阶段：计算最终权重
    
    优点：
    - 减少权重方差
    - 提高采样效率
    - 更好的鲁棒性
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
        初始化辅助粒子滤波器
        
        Args:
            state_dim: 状态维度
            measurement_dim: 观测维度
            n_particles: 粒子数量
            process_noise_std: 过程噪声标准差
            measurement_noise_std: 观测噪声标准差
            state_transition_func: 状态转移函数
            likelihood_func: 似然函数
            resampling_method: 重采样方法
            random_seed: 随机种子
        """
        super().__init__(
            state_dim=state_dim,
            measurement_dim=measurement_dim,
            n_particles=n_particles,
            process_noise_std=process_noise_std,
            measurement_noise_std=measurement_noise_std,
            state_transition_func=state_transition_func,
            likelihood_func=likelihood_func,
            resampling_method=resampling_method,
            random_seed=random_seed
        )
        
        # 辅助变量
        self.auxiliary_indices: Optional[np.ndarray] = None
    
    def predict(self, dt: float) -> None:
        """状态预测（APF版本）
        
        APF的预测步骤：
        1. 计算辅助权重（基于预测观测）
        2. 根据辅助权重重采样
        3. 传播重采样后的粒子
        
        Args:
            dt: 时间步长
        """
        if not self.initialized:
            raise RuntimeError("Filter not initialized")
        
        # 获取过程噪声协方差
        Q = self._get_process_noise_matrix(dt)
        
        # 步骤1：计算每个粒子的预测均值（确定性传播）
        predicted_means = []
        for particle in self.particles:
            # 确定性传播（无噪声）
            mean = self.f(particle.state, dt, np.zeros(self.state_dim))
            predicted_means.append(mean)
        
        # 保存预测均值，供 update() 中的辅助似然计算使用
        self._predicted_means = predicted_means

        # 步骤2：计算辅助权重
        # 辅助权重 w_aux^(i) = w_{k-1}^{(i)} × p(z_k | μ_k^{(i)})
        # 其中 μ_k^{(i)} = f(x_{k-1}^{(i)}) 是粒子i的确定性预测均值
        # 观测似然 p(z_k | μ_k^{(i)}) 衡量预测粒子与当前观测的匹配程度
        auxiliary_weights = np.zeros(self.n_particles)
        measurement = getattr(self, '_last_measurement', None)
        for i, particle in enumerate(self.particles):
            if measurement is not None:
                # 辅助权重 = 先验权重 × 观测似然(在预测均值处)
                # 使用当前观测信息评估每个粒子的"前途"，优先保留与观测一致的粒子
                likelihood_at_pred_mean = self.likelihood(measurement, predicted_means[i])
                auxiliary_weights[i] = particle.weight * likelihood_at_pred_mean
            else:
                # 无观测信息时降级为先验权重（等价于标准粒子滤波）
                auxiliary_weights[i] = particle.weight
        
        # 归一化辅助权重
        auxiliary_weights = auxiliary_weights / np.sum(auxiliary_weights)
        
        # 步骤3：根据辅助权重重采样
        indices = np.random.choice(
            self.n_particles,
            size=self.n_particles,
            p=auxiliary_weights
        )
        
        # 保存辅助索引（用于后续权重计算）
        self.auxiliary_indices = indices
        
        # 步骤4：传播重采样后的粒子
        new_particles = []
        for i, idx in enumerate(indices):
            # 从原始粒子复制
            particle = self.particles[idx].copy()
            
            # 生成过程噪声
            noise = np.random.multivariate_normal(np.zeros(self.state_dim), Q)
            
            # 传播粒子
            particle.state = self.f(particle.state, dt, noise)
            
            # 重置权重
            particle.weight = 1.0 / self.n_particles
            
            new_particles.append(particle)
        
        self.particles = new_particles
    
    def update(self, measurement: np.ndarray,
               measurement_covariance: Optional[np.ndarray] = None) -> None:
        """量测更新（APF版本）
        
        APF的更新步骤：
        1. 计算最终权重 = 似然 / 辅助似然
        2. 归一化权重
        3. 更新状态估计
        
        Args:
            measurement: 观测向量
            measurement_covariance: 观测噪声协方差（可选）
        """
        if not self.initialized:
            raise RuntimeError("Filter not initialized")
        
        z = np.asarray(measurement, dtype=np.float64)
        
        # 保存当前观测，供下一次 predict() 中的辅助权重计算使用
        self._last_measurement = z
        
        # 计算最终权重
        # 最终权重 w_k^{(i)} = p(z_k | x_k^{(i)}) / p(z_k | μ_k^{(i)})
        # 其中 x_k^{(i)} 是传播后的粒子状态，μ_k^{(i)} 是原始粒子的预测均值
        # 分子：当前粒子状态的观测似然
        # 分母：原始粒子预测均值的观测似然（辅助似然）
        for i, particle in enumerate(self.particles):
            # 计算当前似然 p(z_k | x_k^{(i)})
            current_likelihood = self.likelihood(z, particle.state)
            
            # 计算辅助似然 p(z_k | μ_k^{(original_idx)})
            # 使用 predict() 中保存的预测均值，而非当前粒子状态
            if self.auxiliary_indices is not None and hasattr(self, '_predicted_means'):
                original_idx = self.auxiliary_indices[i]
                # 原始粒子的预测均值已经在 predict() 中计算并保存
                pred_mean = self._predicted_means[original_idx]
                auxiliary_likelihood = self.likelihood(z, pred_mean)
            else:
                # 降级：没有辅助信息时使用当前似然（权重为1）
                auxiliary_likelihood = current_likelihood
            
            # 最终权重 = 当前似然 / 辅助似然
            if auxiliary_likelihood > 0:
                particle.weight = current_likelihood / auxiliary_likelihood
            else:
                particle.weight = 0.0
        
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
        """预测并更新（APF版本）
        
        APF的预测-更新联合步骤：
        1. 在predict()中使用观测信息计算辅助权重
        2. 根据辅助权重重采样并传播粒子
        3. 在update()中计算最终权重 = p(z_k|x_k) / p(z_k|μ_k)
        
        Args:
            dt: 时间步长
            measurement: 观测向量
            measurement_covariance: 观测噪声协方差（可选）
        """
        # 先保存观测，供 predict() 中的辅助权重计算使用
        self._last_measurement = np.asarray(measurement, dtype=np.float64)
        self.predict(dt)
        self.update(measurement, measurement_covariance)


class AdaptiveAuxiliaryParticleFilter(AuxiliaryParticleFilter):
    """自适应辅助粒子滤波器
    
    根据粒子分布自适应调整粒子数量
    """
    
    def __init__(self,
                 state_dim: int = 4,
                 measurement_dim: int = 2,
                 n_particles: int = 1000,
                 min_particles: int = 100,
                 max_particles: int = 5000,
                 process_noise_std: float = 0.1,
                 measurement_noise_std: float = 1.0,
                 state_transition_func: Optional[Callable] = None,
                 likelihood_func: Optional[Callable] = None,
                 resampling_method: str = "systematic",
                 random_seed: Optional[int] = None):
        """
        初始化自适应辅助粒子滤波器
        
        Args:
            state_dim: 状态维度
            measurement_dim: 观测维度
            n_particles: 初始粒子数量
            min_particles: 最小粒子数量
            max_particles: 最大粒子数量
            process_noise_std: 过程噪声标准差
            measurement_noise_std: 观测噪声标准差
            state_transition_func: 状态转移函数
            likelihood_func: 似然函数
            resampling_method: 重采样方法
            random_seed: 随机种子
        """
        super().__init__(
            state_dim=state_dim,
            measurement_dim=measurement_dim,
            n_particles=n_particles,
            process_noise_std=process_noise_std,
            measurement_noise_std=measurement_noise_std,
            state_transition_func=state_transition_func,
            likelihood_func=likelihood_func,
            resampling_method=resampling_method,
            random_seed=random_seed
        )
        
        self.min_particles = min_particles
        self.max_particles = max_particles
        
        # 自适应参数
        self.ess_threshold = 0.5  # 有效样本数阈值
        self.adaptation_rate = 0.1  # 自适应速率
    
    def _adapt_particle_count(self) -> None:
        """自适应调整粒子数量"""
        ess = self._effective_sample_size()
        ess_ratio = ess / self.n_particles
        
        if ess_ratio < self.ess_threshold:
            # 有效样本数太低，增加粒子
            new_n = min(int(self.n_particles * (1 + self.adaptation_rate)), 
                       self.max_particles)
        elif ess_ratio > 0.8:
            # 有效样本数很高，可以减少粒子
            new_n = max(int(self.n_particles * (1 - self.adaptation_rate)),
                       self.min_particles)
        else:
            new_n = self.n_particles
        
        if new_n != self.n_particles:
            self._resize_particles(new_n)
    
    def _resize_particles(self, new_n: int) -> None:
        """调整粒子数量
        
        Args:
            new_n: 新的粒子数量
        """
        if new_n > self.n_particles:
            # 增加粒子
            n_additional = new_n - self.n_particles
            
            # 从现有粒子中采样新粒子
            weights = np.array([p.weight for p in self.particles])
            weights = weights / np.sum(weights)
            
            indices = np.random.choice(
                self.n_particles,
                size=n_additional,
                p=weights
            )
            
            for idx in indices:
                new_particle = self.particles[idx].copy()
                new_particle.weight = 1.0 / new_n
                self.particles.append(new_particle)
        
        elif new_n < self.n_particles:
            # 减少粒子
            weights = np.array([p.weight for p in self.particles])
            weights = weights / np.sum(weights)
            
            indices = np.random.choice(
                self.n_particles,
                size=new_n,
                p=weights,
                replace=False
            )
            
            self.particles = [self.particles[i] for i in indices]
        
        self.n_particles = new_n
        
        # 重新归一化权重
        self._normalize_weights()
    
    def update(self, measurement: np.ndarray,
               measurement_covariance: Optional[np.ndarray] = None) -> None:
        """量测更新（自适应版本）"""
        super().update(measurement, measurement_covariance)
        
        # 自适应调整粒子数量
        self._adapt_particle_count()
