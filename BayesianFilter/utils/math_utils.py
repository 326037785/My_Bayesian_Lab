"""
数学工具函数
"""
import numpy as np
from typing import Tuple, Optional


def normalize_angle(angle: float) -> float:
    """将角度归一化到[-pi, pi]范围"""
    return (angle + np.pi) % (2 * np.pi) - np.pi


def wrap_angle(angle: float) -> float:
    """将角度包装到[0, 2*pi]范围"""
    return angle % (2 * np.pi)


def gaussian_pdf(x: np.ndarray, mean: np.ndarray, cov: np.ndarray) -> float:
    """计算多元高斯概率密度函数
    
    Args:
        x: 观测值
        mean: 均值
        cov: 协方差矩阵
        
    Returns:
        概率密度值
    """
    n = len(mean)
    diff = x - mean
    cov_inv = np.linalg.inv(cov)
    cov_det = np.linalg.det(cov)
    
    exponent = -0.5 * diff.T @ cov_inv @ diff
    coefficient = 1.0 / np.sqrt((2 * np.pi) ** n * cov_det)
    
    return coefficient * np.exp(exponent)


def multivariate_normal_sample(mean: np.ndarray, cov: np.ndarray, 
                               n_samples: int = 1) -> np.ndarray:
    """多元正态分布采样
    
    Args:
        mean: 均值向量
        cov: 协方差矩阵
        n_samples: 采样数量
        
    Returns:
        采样结果，形状为 (n_samples, len(mean))
    """
    return np.random.multivariate_normal(mean, cov, n_samples)


def cholesky_sample(mean: np.ndarray, cov: np.ndarray) -> np.ndarray:
    """使用Cholesky分解进行采样
    
    Args:
        mean: 均值向量
        cov: 协方差矩阵
        
    Returns:
        采样结果
    """
    L = np.linalg.cholesky(cov)
    z = np.random.randn(len(mean))
    return mean + L @ z


def mahalanobis_distance(x: np.ndarray, mean: np.ndarray, 
                         cov: np.ndarray) -> float:
    """计算马氏距离
    
    Args:
        x: 观测值
        mean: 均值
        cov: 协方差矩阵
        
    Returns:
        马氏距离
    """
    diff = x - mean
    cov_inv = np.linalg.inv(cov)
    return np.sqrt(diff.T @ cov_inv @ diff)


def nearest_positive_definite(A: np.ndarray) -> np.ndarray:
    """找到最近的正定矩阵
    
    Args:
        A: 输入矩阵
        
    Returns:
        最近的正定矩阵
    """
    B = (A + A.T) / 2
    _, s, V = np.linalg.svd(B)
    H = V.T @ np.diag(s) @ V
    A2 = (B + H) / 2
    A3 = (A2 + A2.T) / 2
    
    if is_positive_definite(A3):
        return A3
    
    spacing = np.spacing(np.linalg.norm(A))
    identity = np.eye(A.shape[0])
    k = 1
    while not is_positive_definite(A3):
        mineig = np.min(np.real(np.linalg.eigvals(A3)))
        A3 += identity * (-mineig * k ** 2 + spacing)
        k += 1
    
    return A3


def is_positive_definite(A: np.ndarray) -> bool:
    """检查矩阵是否正定"""
    try:
        np.linalg.cholesky(A)
        return True
    except np.linalg.LinAlgError:
        return False
