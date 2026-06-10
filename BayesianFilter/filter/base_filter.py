"""
滤波器基类
"""
import numpy as np
from typing import Optional, List, Dict, Any
from abc import ABC, abstractmethod
from dataclasses import dataclass


@dataclass
class FilterState:
    """滤波器状态"""
    timestamp: float
    state: np.ndarray  # 状态估计
    covariance: np.ndarray  # 协方差矩阵
    target_id: Optional[int] = None


class BaseFilter(ABC):
    """滤波器基类
    
    所有滤波器必须实现以下接口：
    - predict: 状态预测
    - update: 量测更新
    - initialize: 初始化
    """
    
    def __init__(self, 
                 state_dim: int,
                 measurement_dim: int,
                 process_noise_std: float = 0.1):
        """
        初始化滤波器
        
        Args:
            state_dim: 状态维度
            measurement_dim: 观测维度
            process_noise_std: 过程噪声标准差
        """
        self.state_dim = state_dim
        self.measurement_dim = measurement_dim
        self.process_noise_std = process_noise_std
        
        # 滤波器状态
        self.state: Optional[np.ndarray] = None
        self.covariance: Optional[np.ndarray] = None
        self.initialized = False
        
        # 历史记录
        self.state_history: List[FilterState] = []
    
    @abstractmethod
    def predict(self, dt: float) -> None:
        """状态预测
        
        Args:
            dt: 时间步长
        """
        pass
    
    @abstractmethod
    def update(self, measurement: np.ndarray, 
               measurement_covariance: Optional[np.ndarray] = None) -> None:
        """量测更新
        
        Args:
            measurement: 观测向量
            measurement_covariance: 观测噪声协方差（可选，使用默认值）
        """
        pass
    
    @abstractmethod
    def _get_state_transition_matrix(self, dt: float) -> np.ndarray:
        """获取状态转移矩阵
        
        Args:
            dt: 时间步长
            
        Returns:
            状态转移矩阵 F
        """
        pass
    
    @abstractmethod
    def _get_process_noise_matrix(self, dt: float) -> np.ndarray:
        """获取过程噪声协方差矩阵
        
        Args:
            dt: 时间步长
            
        Returns:
            过程噪声协方差矩阵 Q
        """
        pass
    
    @abstractmethod
    def _get_measurement_matrix(self) -> np.ndarray:
        """获取观测矩阵
        
        Returns:
            观测矩阵 H
        """
        pass
    
    def initialize(self, 
                   initial_state: np.ndarray,
                   initial_covariance: Optional[np.ndarray] = None) -> None:
        """初始化滤波器
        
        Args:
            initial_state: 初始状态向量
            initial_covariance: 初始协方差矩阵
        """
        self.state = np.asarray(initial_state, dtype=np.float64)
        
        if initial_covariance is not None:
            self.covariance = np.asarray(initial_covariance, dtype=np.float64)
        else:
            # 默认使用较大的初始协方差
            self.covariance = np.eye(self.state_dim) * 100.0
        
        self.initialized = True
        self._record_state(0.0)
    
    def _record_state(self, timestamp: float) -> None:
        """记录当前状态到历史"""
        if self.state is not None and self.covariance is not None:
            self.state_history.append(FilterState(
                timestamp=timestamp,
                state=self.state.copy(),
                covariance=self.covariance.copy()
            ))
    
    def get_state(self) -> Optional[np.ndarray]:
        """获取当前状态估计"""
        return self.state.copy() if self.state is not None else None
    
    def get_covariance(self) -> Optional[np.ndarray]:
        """获取当前协方差矩阵"""
        return self.covariance.copy() if self.covariance is not None else None
    
    def get_position(self) -> Optional[np.ndarray]:
        """获取当前位置估计
        
        Returns:
            位置向量 [x, y]
        """
        if self.state is None:
            return None
        
        # 假设状态格式为 [x, vx, y, vy] 或 [x, vx, ax, y, vy, ay]
        if len(self.state) >= 4:
            return np.array([self.state[0], self.state[2]])
        elif len(self.state) >= 2:
            return np.array([self.state[0], self.state[1]])
        return None
    
    def get_velocity(self) -> Optional[np.ndarray]:
        """获取当前速度估计
        
        Returns:
            速度向量 [vx, vy]
        """
        if self.state is None or len(self.state) < 4:
            return None
        
        return np.array([self.state[1], self.state[3]])
    
    def get_position_uncertainty(self) -> Optional[float]:
        """获取位置不确定性（标准差）
        
        Returns:
            位置不确定性
        """
        if self.covariance is None:
            return None
        
        # 提取位置部分的协方差
        if self.state_dim >= 4:
            pos_cov = self.covariance[[0, 2], :][:, [0, 2]]
        elif self.state_dim >= 2:
            pos_cov = self.covariance[:2, :2]
        else:
            return None
        
        return np.sqrt(np.trace(pos_cov))
    
    def get_state_history(self) -> List[FilterState]:
        """获取状态历史"""
        return self.state_history.copy()
    
    def get_estimated_positions(self) -> np.ndarray:
        """获取估计位置历史
        
        Returns:
            位置数组，形状为 (n_steps, 2)
        """
        positions = []
        for state in self.state_history:
            if len(state.state) >= 4:
                positions.append([state.state[0], state.state[2]])
            elif len(state.state) >= 2:
                positions.append([state.state[0], state.state[1]])
        return np.array(positions) if positions else np.array([])
    
    def reset(self) -> None:
        """重置滤波器"""
        self.state = None
        self.covariance = None
        self.initialized = False
        self.state_history = []
    
    def _validate_state(self) -> bool:
        """验证状态是否有效"""
        if self.state is None or self.covariance is None:
            return False
        
        # 检查NaN
        if np.any(np.isnan(self.state)) or np.any(np.isnan(self.covariance)):
            return False
        
        # 检查无穷大
        if np.any(np.isinf(self.state)) or np.any(np.isinf(self.covariance)):
            return False
        
        return True
    
    def _ensure_positive_definite(self, matrix: np.ndarray) -> np.ndarray:
        """确保矩阵正定
        
        Args:
            matrix: 输入矩阵
            
        Returns:
            正定矩阵
        """
        # 对称化
        matrix = (matrix + matrix.T) / 2
        
        # 检查特征值
        eigenvalues = np.linalg.eigvalsh(matrix)
        if np.all(eigenvalues > 0):
            return matrix
        
        # 如果不是正定，进行修正
        min_eigenvalue = np.min(eigenvalues)
        if min_eigenvalue <= 0:
            # 添加小的正数到对角线
            offset = abs(min_eigenvalue) + 1e-6
            matrix += offset * np.eye(matrix.shape[0])
        
        return matrix
