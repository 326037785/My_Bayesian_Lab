"""
滤波器后端抽象层

设计目标：
1. 统一预测观测、新息协方差、门限计算接口
2. 数据关联算法只需依赖 FilterBackend 接口，不关心具体滤波器实现
3. 支持高斯族（KF/EKF/UKF/CKF）和粒子滤波器族（PF/APF/UPF）无缝切换

数值稳定技巧（参考MATLAB RFS Toolbox）：
1. Cholesky分解计算新息协方差逆：S = LL^T, S^{-1} = L^{-T} L^{-1}
2. 对数空间计算似然，避免数值下溢
3. logsumexp技巧归一化权重
"""
import numpy as np
from typing import Optional, List, Dict, Tuple
from abc import ABC, abstractmethod
from dataclasses import dataclass


def cholesky_inv(M: np.ndarray) -> Tuple[np.ndarray, float]:
    """使用Cholesky分解计算矩阵逆和行列式
    
    比直接 np.linalg.inv 更数值稳定。
    
    Args:
        M: 正定矩阵
        
    Returns:
        (M_inv, det_M): 逆矩阵和行列式
    """
    try:
        L = np.linalg.cholesky(M)
        # det(M) = det(L)^2 = prod(diag(L))^2
        det_M = np.prod(np.diag(L)) ** 2
        # M^{-1} = L^{-T} L^{-1}
        L_inv = np.linalg.inv(L)
        M_inv = L_inv.T @ L_inv
        return M_inv, det_M
    except np.linalg.LinAlgError:
        # 如果Cholesky分解失败，回退到标准方法
        return np.linalg.inv(M), np.linalg.det(M)


def logsumexp(w: np.ndarray) -> float:
    """log-sum-exp技巧，避免数值下溢
    
    Args:
        w: 对数权重数组
        
    Returns:
        log(sum(exp(w)))
    """
    if np.all(w == -np.inf):
        return -np.inf
    val = np.max(w)
    return np.log(np.sum(np.exp(w - val))) + val


def mahalanobis_distance(innovation: np.ndarray, 
                         S_inv: np.ndarray) -> float:
    """计算马氏距离 d = sqrt(innovation^T @ S^{-1} @ innovation)
    
    Args:
        innovation: 新息向量
        S_inv: 新息协方差逆矩阵
        
    Returns:
        马氏距离
    """
    return np.sqrt(innovation.T @ S_inv @ innovation)


@dataclass
class PredictedState:
    """预测状态（用于数据关联）"""
    state: np.ndarray              # 状态向量
    covariance: np.ndarray         # 协方差矩阵
    predicted_meas: np.ndarray     # 预测观测 z_pred = h(x)
    innovation_cov: np.ndarray     # 新息协方差 S = H P H^T + R
    target_id: Optional[int] = None


class FilterBackend(ABC):
    """滤波器后端抽象接口
    
    所有滤波器（KF/EKF/UKF/CKF/PF）必须实现此接口。
    数据关联算法通过此接口获取预测信息，不直接访问滤波器内部。
    """
    
    @property
    @abstractmethod
    def state_dim(self) -> int:
        """状态维度"""
        pass
    
    @property
    @abstractmethod
    def meas_dim(self) -> int:
        """观测维度"""
        pass
    
    @property
    @abstractmethod
    def is_initialized(self) -> bool:
        """是否已初始化"""
        pass
    
    @abstractmethod
    def initialize(self, 
                   initial_state: np.ndarray,
                   initial_covariance: Optional[np.ndarray] = None) -> None:
        """初始化滤波器"""
        pass
    
    @abstractmethod
    def predict(self, dt: float) -> None:
        """状态预测"""
        pass
    
    @abstractmethod
    def update(self, measurement: np.ndarray) -> None:
        """量测更新（单个观测）"""
        pass
    
    @abstractmethod
    def get_state(self) -> np.ndarray:
        """获取当前状态估计"""
        pass
    
    @abstractmethod
    def get_covariance(self) -> np.ndarray:
        """获取当前协方差矩阵"""
        pass
    
    @abstractmethod
    def get_predicted_measurement(self) -> np.ndarray:
        """获取预测观测 z_pred = h(x_predict)
        
        对于线性模型: z_pred = H @ x
        对于非线性模型: z_pred = h(x)
        """
        pass
    
    @abstractmethod
    def get_innovation_covariance(self) -> np.ndarray:
        """获取新息协方差 S = H P H^T + R
        
        对于非线性模型: S = H_jacobian @ P @ H_jacobian^T + R
        """
        pass
    
    @abstractmethod
    def compute_likelihood(self, measurement: np.ndarray) -> float:
        """计算观测似然 p(z|x)
        
        Args:
            measurement: 观测向量
            
        Returns:
            似然值
        """
        pass
    
    @abstractmethod
    def compute_mahalanobis(self, measurement: np.ndarray) -> float:
        """计算马氏距离 d = sqrt(innovation^T @ S^{-1} @ innovation)
        
        Args:
            measurement: 观测向量
            
        Returns:
            马氏距离
        """
        pass
    
    def gating_test(self, measurement: np.ndarray, 
                    threshold: float = 9.21) -> bool:
        """门限测试
        
        Args:
            measurement: 观测向量
            threshold: 卡方阈值（自由度=2，99%置信度=9.21）
            
        Returns:
            True if 通过门限
        """
        d_sq = self.compute_mahalanobis(measurement) ** 2
        return d_sq < threshold


class MultiTargetBackend(ABC):
    """多目标滤波器后端接口
    
    管理多个目标的滤波器，提供统一的数据关联接口。
    """
    
    @abstractmethod
    def add_target(self, target_id: int, 
                   initial_state: np.ndarray,
                   initial_covariance: Optional[np.ndarray] = None) -> None:
        """添加新目标"""
        pass
    
    @abstractmethod
    def remove_target(self, target_id: int) -> None:
        """删除目标"""
        pass
    
    @abstractmethod
    def predict_all(self, dt: float) -> None:
        """预测所有目标"""
        pass
    
    @abstractmethod
    def get_predicted_states(self) -> List[PredictedState]:
        """获取所有目标的预测状态（用于数据关联）"""
        pass
    
    @abstractmethod
    def update_target(self, target_id: int, 
                      measurement: np.ndarray) -> None:
        """更新指定目标"""
        pass
    
    @abstractmethod
    def get_target_ids(self) -> List[int]:
        """获取所有目标ID"""
        pass
