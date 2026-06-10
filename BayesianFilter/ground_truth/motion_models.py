"""
运动模型实现
"""
import numpy as np
from typing import Optional, Tuple
from abc import ABC, abstractmethod


class MotionModelBase(ABC):
    """运动模型基类"""
    
    @abstractmethod
    def state_transition(self, state: np.ndarray, dt: float, 
                         process_noise: Optional[np.ndarray] = None) -> np.ndarray:
        """状态转移
        
        Args:
            state: 当前状态
            dt: 时间步长
            process_noise: 过程噪声（可选）
            
        Returns:
            下一时刻状态
        """
        pass
    
    @abstractmethod
    def get_state_transition_matrix(self, dt: float) -> np.ndarray:
        """获取状态转移矩阵
        
        Args:
            dt: 时间步长
            
        Returns:
            状态转移矩阵 F
        """
        pass
    
    @abstractmethod
    def get_process_noise_matrix(self, dt: float, 
                                  noise_std: float) -> np.ndarray:
        """获取过程噪声协方差矩阵
        
        Args:
            dt: 时间步长
            noise_std: 噪声标准差
            
        Returns:
            过程噪声协方差矩阵 Q
        """
        pass
    
    @abstractmethod
    def get_state_dimension(self) -> int:
        """获取状态维度"""
        pass


class ConstantVelocityModel(MotionModelBase):
    """匀速直线运动模型 (CV)
    
    状态向量: [x, vx, y, vy]
    """
    
    def state_transition(self, state: np.ndarray, dt: float,
                         process_noise: Optional[np.ndarray] = None) -> np.ndarray:
        F = self.get_state_transition_matrix(dt)
        new_state = F @ state
        
        if process_noise is not None:
            new_state += process_noise
            
        return new_state
    
    def get_state_transition_matrix(self, dt: float) -> np.ndarray:
        """状态转移矩阵
        
        F = [[1, dt, 0,  0],
             [0,  1, 0,  0],
             [0,  0, 1, dt],
             [0,  0, 0,  1]]
        """
        F = np.array([
            [1, dt, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, dt],
            [0, 0, 0, 1]
        ])
        return F
    
    def get_process_noise_matrix(self, dt: float, 
                                  noise_std: float) -> np.ndarray:
        """过程噪声协方差矩阵
        
        使用离散白噪声加速度模型
        """
        q = noise_std ** 2
        Q = q * np.array([
            [dt**3/3, dt**2/2, 0, 0],
            [dt**2/2, dt, 0, 0],
            [0, 0, dt**3/3, dt**2/2],
            [0, 0, dt**2/2, dt]
        ])
        return Q
    
    def get_state_dimension(self) -> int:
        return 4


class ConstantAccelerationModel(MotionModelBase):
    """匀加速运动模型 (CA)
    
    状态向量: [x, vx, ax, y, vy, ay]
    """
    
    def state_transition(self, state: np.ndarray, dt: float,
                         process_noise: Optional[np.ndarray] = None) -> np.ndarray:
        F = self.get_state_transition_matrix(dt)
        new_state = F @ state
        
        if process_noise is not None:
            new_state += process_noise
            
        return new_state
    
    def get_state_transition_matrix(self, dt: float) -> np.ndarray:
        """状态转移矩阵"""
        F = np.array([
            [1, dt, dt**2/2, 0, 0, 0],
            [0, 1, dt, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, 1, dt, dt**2/2],
            [0, 0, 0, 0, 1, dt],
            [0, 0, 0, 0, 0, 1]
        ])
        return F
    
    def get_process_noise_matrix(self, dt: float, 
                                  noise_std: float) -> np.ndarray:
        """过程噪声协方差矩阵"""
        q = noise_std ** 2
        dt2 = dt ** 2
        dt3 = dt ** 3
        dt4 = dt ** 4
        dt5 = dt ** 5
        
        Q = q * np.array([
            [dt5/20, dt4/8, dt3/6, 0, 0, 0],
            [dt4/8, dt3/3, dt2/2, 0, 0, 0],
            [dt3/6, dt2/2, dt, 0, 0, 0],
            [0, 0, 0, dt5/20, dt4/8, dt3/6],
            [0, 0, 0, dt4/8, dt3/3, dt2/2],
            [0, 0, 0, dt3/6, dt2/2, dt]
        ])
        return Q
    
    def get_state_dimension(self) -> int:
        return 6


class CoordinatedTurnModel(MotionModelBase):
    """协调转弯模型 (CT)
    
    状态向量: [x, vx, y, vy, omega]
    其中omega为转弯角速度
    """
    
    def __init__(self, known_turn_rate: Optional[float] = None):
        """
        初始化协调转弯模型
        
        Args:
            known_turn_rate: 已知的转弯角速度，如果为None则作为状态估计
        """
        self.known_turn_rate = known_turn_rate
    
    def state_transition(self, state: np.ndarray, dt: float,
                         process_noise: Optional[np.ndarray] = None) -> np.ndarray:
        x, vx, y, vy = state[0], state[1], state[2], state[3]
        
        if self.known_turn_rate is not None:
            omega = self.known_turn_rate
        else:
            omega = state[4] if len(state) > 0 else 0.0
        
        # 处理小角度情况
        if abs(omega) < 1e-6:
            new_x = x + vx * dt
            new_vx = vx
            new_y = y + vy * dt
            new_vy = vy
        else:
            sin_omega_dt = np.sin(omega * dt)
            cos_omega_dt = np.cos(omega * dt)
            
            new_x = x + (vx * sin_omega_dt - vy * (1 - cos_omega_dt)) / omega
            new_vx = vx * cos_omega_dt - vy * sin_omega_dt
            new_y = y + (vx * (1 - cos_omega_dt) + vy * sin_omega_dt) / omega
            new_vy = vx * sin_omega_dt + vy * cos_omega_dt
        
        if self.known_turn_rate is not None:
            new_state = np.array([new_x, new_vx, new_y, new_vy])
        else:
            new_omega = omega
            new_state = np.array([new_x, new_vx, new_y, new_vy, new_omega])
        
        if process_noise is not None:
            new_state += process_noise
            
        return new_state
    
    def get_state_transition_matrix(self, dt: float) -> np.ndarray:
        """状态转移矩阵（雅可比矩阵）"""
        # 简化版本，使用单位矩阵加上时间项
        if self.known_turn_rate is not None:
            F = np.array([
                [1, dt, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, dt],
                [0, 0, 0, 1]
            ])
        else:
            F = np.array([
                [1, dt, 0, 0, 0],
                [0, 1, 0, 0, 0],
                [0, 0, 1, dt, 0],
                [0, 0, 0, 1, 0],
                [0, 0, 0, 0, 1]
            ])
        return F
    
    def get_process_noise_matrix(self, dt: float, 
                                  noise_std: float) -> np.ndarray:
        """过程噪声协方差矩阵"""
        q = noise_std ** 2
        
        if self.known_turn_rate is not None:
            Q = q * np.array([
                [dt**3/3, dt**2/2, 0, 0],
                [dt**2/2, dt, 0, 0],
                [0, 0, dt**3/3, dt**2/2],
                [0, 0, dt**2/2, dt]
            ])
        else:
            Q = q * np.array([
                [dt**3/3, dt**2/2, 0, 0, 0],
                [dt**2/2, dt, 0, 0, 0],
                [0, 0, dt**3/3, dt**2/2, 0],
                [0, 0, dt**2/2, dt, 0],
                [0, 0, 0, 0, dt]
            ])
        return Q
    
    def get_state_dimension(self) -> int:
        if self.known_turn_rate is not None:
            return 4
        return 5


class RandomWalkModel(MotionModelBase):
    """随机游走模型 (RW)
    
    状态向量: [x, y]
    """
    
    def state_transition(self, state: np.ndarray, dt: float,
                         process_noise: Optional[np.ndarray] = None) -> np.ndarray:
        new_state = state.copy()
        
        if process_noise is not None:
            new_state += process_noise
            
        return new_state
    
    def get_state_transition_matrix(self, dt: float) -> np.ndarray:
        """状态转移矩阵（单位矩阵）"""
        return np.eye(2)
    
    def get_process_noise_matrix(self, dt: float, 
                                  noise_std: float) -> np.ndarray:
        """过程噪声协方差矩阵"""
        q = noise_std ** 2
        Q = q * dt * np.eye(2)
        return Q
    
    def get_state_dimension(self) -> int:
        return 2


def create_motion_model(model_type: str, **kwargs) -> MotionModelBase:
    """创建运动模型的工厂函数
    
    Args:
        model_type: 模型类型 ("CV", "CA", "CT", "RW")
        **kwargs: 模型参数
        
    Returns:
        运动模型实例
    """
    models = {
        "CV": ConstantVelocityModel,
        "CA": ConstantAccelerationModel,
        "CT": CoordinatedTurnModel,
        "RW": RandomWalkModel
    }
    
    if model_type not in models:
        raise ValueError(f"Unknown model type: {model_type}")
    
    return models[model_type](**kwargs)
