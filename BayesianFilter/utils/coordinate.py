"""
坐标转换工具
"""
import numpy as np
from typing import Tuple


def cartesian_to_polar(x: float, y: float) -> Tuple[float, float]:
    """笛卡尔坐标转极坐标
    
    Args:
        x: x坐标
        y: y坐标
        
    Returns:
        (距离, 角度) 角度单位为弧度
    """
    r = np.sqrt(x ** 2 + y ** 2)
    theta = np.arctan2(y, x)
    return r, theta


def polar_to_cartesian(r: float, theta: float) -> Tuple[float, float]:
    """极坐标转笛卡尔坐标
    
    Args:
        r: 距离
        theta: 角度（弧度）
        
    Returns:
        (x, y)
    """
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y


def cartesian_to_spherical(x: float, y: float, z: float) -> Tuple[float, float, float]:
    """笛卡尔坐标转球坐标
    
    Args:
        x, y, z: 笛卡尔坐标
        
    Returns:
        (距离, 方位角, 俯仰角)
    """
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    azimuth = np.arctan2(y, x)
    elevation = np.arccos(z / r) if r > 0 else 0.0
    return r, azimuth, elevation


def spherical_to_cartesian(r: float, azimuth: float, elevation: float) -> Tuple[float, float, float]:
    """球坐标转笛卡尔坐标
    
    Args:
        r: 距离
        azimuth: 方位角（弧度）
        elevation: 俯仰角（弧度）
        
    Returns:
        (x, y, z)
    """
    x = r * np.sin(elevation) * np.cos(azimuth)
    y = r * np.sin(elevation) * np.sin(azimuth)
    z = r * np.cos(elevation)
    return x, y, z


def state_to_measurement_polar(state: np.ndarray) -> np.ndarray:
    """将状态向量 [x, vx, y, vy] 转换为极坐标观测 [r, theta]
    
    Args:
        state: 状态向量 [x, vx, y, vy]
        
    Returns:
        极坐标观测 [r, theta]
    """
    x, y = state[0], state[2]
    r, theta = cartesian_to_polar(x, y)
    return np.array([r, theta])


def measurement_polar_to_cartesian(measurement: np.ndarray) -> np.ndarray:
    """将极坐标观测 [r, theta] 转换为笛卡尔坐标 [x, y]
    
    Args:
        measurement: 极坐标观测 [r, theta]
        
    Returns:
        笛卡尔坐标 [x, y]
    """
    r, theta = measurement[0], measurement[1]
    x, y = polar_to_cartesian(r, theta)
    return np.array([x, y])
