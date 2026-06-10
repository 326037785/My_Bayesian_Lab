"""
目标基类
"""
import numpy as np
from typing import Optional, List, Dict, Any
from dataclasses import dataclass
from enum import Enum


class MotionModel(Enum):
    """运动模型类型"""
    CONSTANT_VELOCITY = "CV"      # 匀速直线
    CONSTANT_ACCELERATION = "CA"  # 匀加速
    COORDINATED_TURN = "CT"       # 协调转弯
    RANDOM_WALK = "RW"            # 随机游走


@dataclass
class TargetState:
    """目标状态"""
    timestamp: float
    state: np.ndarray  # 状态向量
    target_id: int


class Target:
    """目标类
    
    表示一个运动目标，包含其轨迹和属性
    """
    
    def __init__(self, 
                 target_id: int,
                 initial_state: np.ndarray,
                 birth_time: float,
                 death_time: float,
                 motion_model: MotionModel = MotionModel.CONSTANT_VELOCITY):
        """
        初始化目标
        
        Args:
            target_id: 目标ID
            initial_state: 初始状态向量
            birth_time: 目标出现时间
            death_time: 目标消失时间
            motion_model: 运动模型类型
        """
        self.target_id = target_id
        self.initial_state = initial_state.copy()
        self.birth_time = birth_time
        self.death_time = death_time
        self.motion_model = motion_model
        
        # 存储轨迹历史
        self.trajectory: List[TargetState] = []
        
    def is_alive(self, time: float) -> bool:
        """检查目标在指定时间是否存活
        
        Args:
            time: 时间戳
            
        Returns:
            目标是否存活
        """
        return self.birth_time <= time <= self.death_time
    
    def add_state(self, timestamp: float, state: np.ndarray):
        """添加状态到轨迹
        
        Args:
            timestamp: 时间戳
            state: 状态向量
        """
        self.trajectory.append(TargetState(
            timestamp=timestamp,
            state=state.copy(),
            target_id=self.target_id
        ))
    
    def get_state_at_time(self, time: float) -> Optional[np.ndarray]:
        """获取指定时间的状态
        
        Args:
            time: 时间戳
            
        Returns:
            状态向量，如果不存在则返回None
        """
        for state in self.trajectory:
            if abs(state.timestamp - time) < 1e-6:
                return state.state
        return None
    
    def get_position_at_time(self, time: float) -> Optional[np.ndarray]:
        """获取指定时间的位置
        
        Args:
            time: 时间戳
            
        Returns:
            位置向量 [x, y] 或 [x, y, z]，如果不存在则返回None
        """
        state = self.get_state_at_time(time)
        if state is None:
            return None
        
        # 根据状态向量提取位置
        # 假设状态向量格式: [x, vx, y, vy] 或 [x, vx, ax, y, vy, ay]
        if len(state) >= 4:
            return np.array([state[0], state[2]])
        return None
    
    def get_trajectory_positions(self) -> np.ndarray:
        """获取整个轨迹的位置序列
        
        Returns:
            位置序列，形状为 (n_steps, 2)
        """
        positions = []
        for ts in self.trajectory:
            pos = self.get_position_at_time(ts.timestamp)
            if pos is not None:
                positions.append(pos)
        return np.array(positions) if positions else np.array([])
    
    def __repr__(self) -> str:
        return (f"Target(id={self.target_id}, "
                f"birth={self.birth_time}, "
                f"death={self.death_time}, "
                f"model={self.motion_model.value})")
