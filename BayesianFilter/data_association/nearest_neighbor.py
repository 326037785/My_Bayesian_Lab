"""
最近邻数据关联算法
"""
import numpy as np
from typing import Optional, List, Dict, Set
from .base_association import BaseAssociation, AssociationResult


class NearestNeighborAssociation(BaseAssociation):
    """最近邻数据关联 (NN)
    
    将每个观测关联到最近的目标
    
    特点：
    - 简单高效
    - 贪心策略，可能不是全局最优
    - 适用于目标稀疏的场景
    """
    
    def __init__(self, 
                 gating_threshold: float = 9.21,
                 use_mahalanobis: bool = True):
        """
        初始化最近邻关联
        
        Args:
            gating_threshold: 门限阈值
            use_mahalanobis: 是否使用马氏距离
        """
        super().__init__(gating_threshold, use_mahalanobis)
    
    def associate(self, 
                  measurements: np.ndarray,
                  predicted_measurements: np.ndarray,
                  measurement_covariances: Optional[List[np.ndarray]] = None,
                  innovation_covariances: Optional[List[np.ndarray]] = None) -> AssociationResult:
        """执行最近邻关联
        
        Args:
            measurements: 观测矩阵，形状为 (n_measurements, measurement_dim)
            predicted_measurements: 预测观测矩阵，形状为 (n_targets, measurement_dim)
            measurement_covariances: 观测噪声协方差列表
            innovation_covariances: 新息协方差列表
            
        Returns:
            关联结果
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
        
        # 初始化结果
        associations = {}
        unassociated_measurements = set(range(n_meas))
        unassociated_targets = set(range(n_targets))
        
        # 对每个观测，找到最近的目标
        for i in range(n_meas):
            # 找到门限内的目标
            valid_targets = np.where(association_matrix[i] < np.inf)[0]
            
            if len(valid_targets) > 0:
                # 找到最近的目标
                best_target = valid_targets[np.argmin(association_matrix[i, valid_targets])]
                
                # 检查目标是否已被关联
                if best_target in unassociated_targets:
                    associations[i] = best_target
                    unassociated_measurements.discard(i)
                    unassociated_targets.discard(best_target)
        
        return AssociationResult(
            associations=associations,
            unassociated_measurements=unassociated_measurements,
            unassociated_targets=unassociated_targets,
            association_matrix=association_matrix
        )


class KNearestNeighborAssociation(BaseAssociation):
    """K最近邻数据关联 (KNN)
    
    将每个观测关联到K个最近的目标，使用概率加权
    
    特点：
    - 考虑多个可能的关联
    - 使用概率加权，更鲁棒
    - 适用于目标密集的场景
    """
    
    def __init__(self, 
                 k: int = 3,
                 gating_threshold: float = 9.21,
                 use_mahalanobis: bool = True):
        """
        初始化K最近邻关联
        
        Args:
            k: 考虑的最近邻数量
            gating_threshold: 门限阈值
            use_mahalanobis: 是否使用马氏距离
        """
        super().__init__(gating_threshold, use_mahalanobis)
        self.k = k
    
    def associate(self, 
                  measurements: np.ndarray,
                  predicted_measurements: np.ndarray,
                  measurement_covariances: Optional[List[np.ndarray]] = None,
                  innovation_covariances: Optional[List[np.ndarray]] = None) -> AssociationResult:
        """执行K最近邻关联
        
        Args:
            measurements: 观测矩阵
            predicted_measurements: 预测观测矩阵
            measurement_covariances: 观测噪声协方差列表
            innovation_covariances: 新息协方差列表
            
        Returns:
            关联结果
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
        
        # 初始化结果
        associations = {}
        unassociated_measurements = set(range(n_meas))
        unassociated_targets = set(range(n_targets))
        
        # 对每个观测，找到K个最近的目标
        for i in range(n_meas):
            # 找到门限内的目标
            valid_targets = np.where(association_matrix[i] < np.inf)[0]
            
            if len(valid_targets) > 0:
                # 按距离排序
                sorted_indices = np.argsort(association_matrix[i, valid_targets])
                
                # 取前K个
                k_nearest = valid_targets[sorted_indices[:self.k]]
                
                # 选择最近的且未被关联的目标
                for target_idx in k_nearest:
                    if target_idx in unassociated_targets:
                        associations[i] = target_idx
                        unassociated_measurements.discard(i)
                        unassociated_targets.discard(target_idx)
                        break
        
        return AssociationResult(
            associations=associations,
            unassociated_measurements=unassociated_measurements,
            unassociated_targets=unassociated_targets,
            association_matrix=association_matrix
        )
    
    def associate_with_probabilities(self,
                                      measurements: np.ndarray,
                                      predicted_measurements: np.ndarray,
                                      measurement_covariances: Optional[List[np.ndarray]] = None,
                                      innovation_covariances: Optional[List[np.ndarray]] = None) -> Dict[int, Dict[int, float]]:
        """执行K最近邻关联，返回概率
        
        Args:
            measurements: 观测矩阵
            predicted_measurements: 预测观测矩阵
            measurement_covariances: 观测噪声协方差列表
            innovation_covariances: 新息协方差列表
            
        Returns:
            字典，键为观测索引，值为目标概率字典
        """
        n_meas = measurements.shape[0]
        n_targets = predicted_measurements.shape[0]
        
        if n_meas == 0 or n_targets == 0:
            return {}
        
        # 计算关联矩阵
        association_matrix = self.compute_association_matrix(
            measurements,
            predicted_measurements,
            measurement_covariances,
            innovation_covariances
        )
        
        # 计算每个观测关联到每个目标的概率
        association_probs = {}
        
        for i in range(n_meas):
            # 找到门限内的目标
            valid_targets = np.where(association_matrix[i] < np.inf)[0]
            
            if len(valid_targets) > 0:
                # 计算距离（转换为相似度）
                distances = association_matrix[i, valid_targets]
                
                # 使用指数核转换为相似度
                similarities = np.exp(-0.5 * distances)
                
                # 归一化得到概率
                total_similarity = np.sum(similarities)
                if total_similarity > 0:
                    probabilities = similarities / total_similarity
                else:
                    probabilities = np.ones(len(valid_targets)) / len(valid_targets)
                
                # 存储概率
                association_probs[i] = {}
                for j, target_idx in enumerate(valid_targets):
                    association_probs[i][target_idx] = probabilities[j]
        
        return association_probs


class ProbabilisticNearestNeighbor(BaseAssociation):
    """概率最近邻关联 (PNN)
    
    使用概率方法选择关联，而不是简单的最近邻
    
    特点：
    - 考虑距离的概率分布
    - 更鲁棒，减少错误关联
    - 适用于噪声较大的场景
    """
    
    def __init__(self, 
                 gating_threshold: float = 9.21,
                 use_mahalanobis: bool = True,
                 temperature: float = 1.0):
        """
        初始化概率最近邻关联
        
        Args:
            gating_threshold: 门限阈值
            use_mahalanobis: 是否使用马氏距离
            temperature: 温度参数（控制选择的随机性）
        """
        super().__init__(gating_threshold, use_mahalanobis)
        self.temperature = temperature
    
    def associate(self, 
                  measurements: np.ndarray,
                  predicted_measurements: np.ndarray,
                  measurement_covariances: Optional[List[np.ndarray]] = None,
                  innovation_covariances: Optional[List[np.ndarray]] = None) -> AssociationResult:
        """执行概率最近邻关联
        
        Args:
            measurements: 观测矩阵
            predicted_measurements: 预测观测矩阵
            measurement_covariances: 观测噪声协方差列表
            innovation_covariances: 新息协方差列表
            
        Returns:
            关联结果
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
        
        # 初始化结果
        associations = {}
        unassociated_measurements = set(range(n_meas))
        unassociated_targets = set(range(n_targets))
        
        # 对每个观测，使用概率方法选择目标
        for i in range(n_meas):
            # 找到门限内的目标
            valid_targets = np.where(association_matrix[i] < np.inf)[0]
            
            if len(valid_targets) > 0:
                # 计算距离
                distances = association_matrix[i, valid_targets]
                
                # 使用softmax计算概率
                exp_distances = np.exp(-distances / self.temperature)
                probabilities = exp_distances / np.sum(exp_distances)
                
                # 根据概率选择目标
                selected_idx = np.random.choice(len(valid_targets), p=probabilities)
                target_idx = valid_targets[selected_idx]
                
                # 检查目标是否已被关联
                if target_idx in unassociated_targets:
                    associations[i] = target_idx
                    unassociated_measurements.discard(i)
                    unassociated_targets.discard(target_idx)
        
        return AssociationResult(
            associations=associations,
            unassociated_measurements=unassociated_measurements,
            unassociated_targets=unassociated_targets,
            association_matrix=association_matrix
        )
