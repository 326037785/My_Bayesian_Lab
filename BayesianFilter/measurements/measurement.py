"""
观测数据结构
"""
import numpy as np
from typing import Optional, List
from dataclasses import dataclass
from enum import Enum


class MeasurementType(Enum):
    """观测类型"""
    LINEAR = "linear"          # 线性观测 [x, y]
    POLAR = "polar"            # 极坐标观测 [range, bearing]
    RANGE_ONLY = "range_only"  # 仅距离观测
    BEARING_ONLY = "bearing_only"  # 仅角度观测


@dataclass
class Measurement:
    """单个观测数据
    
    Attributes:
        timestamp: 观测时间戳
        measurement: 观测向量
        measurement_type: 观测类型
        target_id: 真实目标ID（用于评估，实际跟踪中未知）
        is_clutter: 是否为杂波
        covariance: 观测噪声协方差矩阵
    """
    timestamp: float
    measurement: np.ndarray
    measurement_type: MeasurementType
    target_id: Optional[int] = None
    is_clutter: bool = False
    covariance: Optional[np.ndarray] = None
    
    def __post_init__(self):
        """初始化后处理"""
        self.measurement = np.asarray(self.measurement, dtype=np.float64)
        if self.covariance is not None:
            self.covariance = np.asarray(self.covariance, dtype=np.float64)
    
    @property
    def dimension(self) -> int:
        """观测维度"""
        return len(self.measurement)
    
    def to_cartesian(self) -> np.ndarray:
        """将观测转换为笛卡尔坐标
        
        Returns:
            笛卡尔坐标 [x, y]
        """
        if self.measurement_type == MeasurementType.LINEAR:
            return self.measurement.copy()
        elif self.measurement_type == MeasurementType.POLAR:
            r, theta = self.measurement[0], self.measurement[1]
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            return np.array([x, y])
        else:
            raise NotImplementedError(
                f"Cannot convert {self.measurement_type} to cartesian"
            )
    
    def __repr__(self) -> str:
        return (f"Measurement(time={self.timestamp}, "
                f"type={self.measurement_type.value}, "
                f"target_id={self.target_id}, "
                f"is_clutter={self.is_clutter})")


class MeasurementSet:
    """观测集合
    
    某一时刻的所有观测
    """
    
    def __init__(self, timestamp: float):
        """
        初始化观测集合
        
        Args:
            timestamp: 时间戳
        """
        self.timestamp = timestamp
        self.measurements: List[Measurement] = []
    
    def add_measurement(self, measurement: Measurement):
        """添加观测"""
        self.measurements.append(measurement)
    
    @property
    def size(self) -> int:
        """观测数量"""
        return len(self.measurements)
    
    @property
    def n_clutter(self) -> int:
        """杂波数量"""
        return sum(1 for m in self.measurements if m.is_clutter)
    
    @property
    def n_target_measurements(self) -> int:
        """目标观测数量"""
        return sum(1 for m in self.measurements if not m.is_clutter)
    
    def get_measurements_array(self) -> np.ndarray:
        """获取观测数组
        
        Returns:
            观测数组，形状为 (n_measurements, measurement_dim)
        """
        if not self.measurements:
            return np.array([])
        return np.array([m.measurement for m in self.measurements])
    
    def get_target_measurements(self) -> List[Measurement]:
        """获取目标观测"""
        return [m for m in self.measurements if not m.is_clutter]
    
    def get_clutter_measurements(self) -> List[Measurement]:
        """获取杂波观测"""
        return [m for m in self.measurements if m.is_clutter]
    
    def __len__(self) -> int:
        return self.size
    
    def __iter__(self):
        return iter(self.measurements)
    
    def __getitem__(self, index: int) -> Measurement:
        return self.measurements[index]
