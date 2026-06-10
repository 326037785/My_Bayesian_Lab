"""
数据关联基类
"""
import numpy as np
from typing import Optional, List, Dict, Tuple, Set
from abc import ABC, abstractmethod
from dataclasses import dataclass


@dataclass
class AssociationResult:
    """关联结果
    
    Attributes:
        associations: 关联字典，键为观测索引，值为目标索引
        unassociated_measurements: 未关联的观测索引集合
        unassociated_targets: 未关联的目标索引集合
        association_matrix: 关联矩阵（可选）
    """
    associations: Dict[int, int]
    unassociated_measurements: Set[int]
    unassociated_targets: Set[int]
    association_matrix: Optional[np.ndarray] = None


class BaseAssociation(ABC):
    """数据关联基类
    
    所有关联算法必须实现以下接口：
    - associate: 执行数据关联
    - compute_association_matrix: 计算关联矩阵
    """
    
    def __init__(self, 
                 gating_threshold: float = 9.21,
                 use_mahalanobis: bool = True):
        """
        初始化数据关联
        
        Args:
            gating_threshold: 门限阈值（卡方分布，自由度=2，95%置信度=5.99，99%置信度=9.21）
            use_mahalanobis: 是否使用马氏距离（否则使用欧氏距离）
        """
        self.gating_threshold = gating_threshold
        self.use_mahalanobis = use_mahalanobis
    
    @abstractmethod
    def associate(self, 
                  measurements: np.ndarray,
                  predicted_measurements: np.ndarray,
                  measurement_covariances: Optional[List[np.ndarray]] = None,
                  innovation_covariances: Optional[List[np.ndarray]] = None) -> AssociationResult:
        """执行数据关联
        
        Args:
            measurements: 观测矩阵，形状为 (n_measurements, measurement_dim)
            predicted_measurements: 预测观测矩阵，形状为 (n_targets, measurement_dim)
            measurement_covariances: 观测噪声协方差列表
            innovation_covariances: 新息协方差列表
            
        Returns:
            关联结果
        """
        pass
    
    def compute_association_matrix(self,
                                    measurements: np.ndarray,
                                    predicted_measurements: np.ndarray,
                                    measurement_covariances: Optional[List[np.ndarray]] = None,
                                    innovation_covariances: Optional[List[np.ndarray]] = None) -> np.ndarray:
        """计算关联矩阵
        
        Args:
            measurements: 观测矩阵
            predicted_measurements: 预测观测矩阵
            measurement_covariances: 观测噪声协方差列表
            innovation_covariances: 新息协方差列表
            
        Returns:
            关联矩阵，形状为 (n_measurements, n_targets)
        """
        n_meas = measurements.shape[0]
        n_targets = predicted_measurements.shape[0]
        
        # 初始化关联矩阵
        association_matrix = np.full((n_meas, n_targets), np.inf)
        
        for i in range(n_meas):
            for j in range(n_targets):
                # 计算距离
                if self.use_mahalanobis and innovation_covariances is not None:
                    distance = self._mahalanobis_distance(
                        measurements[i],
                        predicted_measurements[j],
                        innovation_covariances[j]
                    )
                else:
                    distance = self._euclidean_distance(
                        measurements[i],
                        predicted_measurements[j]
                    )
                
                # 应用门限
                if distance <= self.gating_threshold:
                    association_matrix[i, j] = distance
        
        return association_matrix
    
    def _euclidean_distance(self, meas: np.ndarray, pred: np.ndarray) -> float:
        """计算欧氏距离"""
        return np.sqrt(np.sum((meas - pred) ** 2))
    
    def _mahalanobis_distance(self, meas: np.ndarray, pred: np.ndarray,
                               covariance: np.ndarray) -> float:
        """计算马氏距离"""
        diff = meas - pred
        try:
            cov_inv = np.linalg.inv(covariance)
            return np.sqrt(diff.T @ cov_inv @ diff)
        except np.linalg.LinAlgError:
            # 如果矩阵奇异，使用欧氏距离
            return self._euclidean_distance(meas, pred)
    
    def _gating(self, distance: float) -> bool:
        """门限判断
        
        Args:
            distance: 距离值
            
        Returns:
            是否在门限内
        """
        return distance <= self.gating_threshold
    
    def _chi2_quantile(self, p: float, df: int) -> float:
        """卡方分布分位数
        
        Args:
            p: 置信水平
            df: 自由度
            
        Returns:
            分位数
        """
        from scipy.stats import chi2
        return chi2.ppf(p, df)
    
    def get_gating_region(self, 
                           predicted_measurement: np.ndarray,
                           innovation_covariance: np.ndarray) -> Tuple[np.ndarray, float]:
        """获取门限区域
        
        Args:
            predicted_measurement: 预测观测
            innovation_covariance: 新息协方差
            
        Returns:
            (中心, 半径)
        """
        # 计算门限区域的半径
        radius = np.sqrt(self.gating_threshold)
        
        return predicted_measurement, radius


class GlobalNearestNeighbor(BaseAssociation):
    """全局最近邻 (GNN)
    
    使用匈牙利算法求解最优关联
    """
    
    def __init__(self, 
                 gating_threshold: float = 9.21,
                 use_mahalanobis: bool = True,
                 max_cost: float = 1e6):
        """
        初始化GNN
        
        Args:
            gating_threshold: 门限阈值
            use_mahalanobis: 是否使用马氏距离
            max_cost: 最大代价（用于填充无效关联）
        """
        super().__init__(gating_threshold, use_mahalanobis)
        self.max_cost = max_cost
    
    def associate(self, 
                  measurements: np.ndarray,
                  predicted_measurements: np.ndarray,
                  measurement_covariances: Optional[List[np.ndarray]] = None,
                  innovation_covariances: Optional[List[np.ndarray]] = None) -> AssociationResult:
        """执行GNN关联
        
        使用匈牙利算法求解最优分配
        """
        n_meas = measurements.shape[0]
        n_targets = predicted_measurements.shape[0]
        
        if n_meas == 0 or n_targets == 0:
            return AssociationResult(
                associations={},
                unassociated_measurements=set(range(n_meas)),
                unassociated_targets=set(range(n_targets))
            )
        
        # 计算关联矩阵
        association_matrix = self.compute_association_matrix(
            measurements,
            predicted_measurements,
            measurement_covariances,
            innovation_covariances
        )
        
        # 创建代价矩阵（用于匈牙利算法）
        cost_matrix = np.copy(association_matrix)
        cost_matrix[cost_matrix == np.inf] = self.max_cost
        
        # 使用匈牙利算法求解最优分配
        from scipy.optimize import linear_sum_assignment
        row_indices, col_indices = linear_sum_assignment(cost_matrix)
        
        # 构建关联结果
        associations = {}
        unassociated_measurements = set(range(n_meas))
        unassociated_targets = set(range(n_targets))
        
        for row, col in zip(row_indices, col_indices):
            if association_matrix[row, col] < np.inf:
                associations[row] = col
                unassociated_measurements.discard(row)
                unassociated_targets.discard(col)
        
        return AssociationResult(
            associations=associations,
            unassociated_measurements=unassociated_measurements,
            unassociated_targets=unassociated_targets,
            association_matrix=association_matrix
        )
