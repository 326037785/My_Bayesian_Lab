"""
性能评估指标
"""
import numpy as np
from typing import Optional, List, Tuple
from scipy.optimize import linear_sum_assignment


def RMSE(estimates: np.ndarray, truths: np.ndarray) -> float:
    """计算均方根误差 (RMSE)
    
    Args:
        estimates: 估计值数组，形状为 (n_samples, dim)
        truths: 真实值数组，形状为 (n_samples, dim)
        
    Returns:
        RMSE值
    """
    if len(estimates) == 0 or len(truths) == 0:
        return 0.0
    
    # 确保形状匹配
    estimates = np.asarray(estimates)
    truths = np.asarray(truths)
    
    if estimates.shape != truths.shape:
        raise ValueError(f"Shape mismatch: {estimates.shape} vs {truths.shape}")
    
    # 计算RMSE
    errors = estimates - truths
    mse = np.mean(np.sum(errors ** 2, axis=1))
    return np.sqrt(mse)


def RMSE_position(estimates: np.ndarray, truths: np.ndarray) -> float:
    """计算位置RMSE
    
    Args:
        estimates: 估计状态数组，形状为 (n_samples, state_dim)
        truths: 真实状态数组，形状为 (n_samples, state_dim)
        
    Returns:
        位置RMSE值
    """
    if len(estimates) == 0 or len(truths) == 0:
        return 0.0
    
    estimates = np.asarray(estimates)
    truths = np.asarray(truths)
    
    # 提取位置分量（假设状态为 [x, vx, y, vy] 或 [x, vx, ax, y, vy, ay]）
    if estimates.shape[1] >= 4:
        est_pos = estimates[:, [0, 2]]
        true_pos = truths[:, [0, 2]]
    elif estimates.shape[1] >= 2:
        est_pos = estimates[:, :2]
        true_pos = truths[:, :2]
    else:
        est_pos = estimates
        true_pos = truths
    
    return RMSE(est_pos, true_pos)


def OSPA(true_sets: List[np.ndarray], 
         estimated_sets: List[np.ndarray],
         c: float = 100.0,
         p: float = 2.0) -> Tuple[float, List[float], List[float]]:
    """计算最优子模式分配 (OSPA) 指标
    
    OSPA是多目标跟踪性能的综合指标，考虑：
    - 定位误差（估计目标与真实目标的距离）
    - 基数误差（目标数量的差异）
    
    Args:
        true_sets: 真实目标集合列表，每个元素为 (n_targets, dim) 数组
        estimated_sets: 估计目标集合列表，每个元素为 (n_estimates, dim) 数组
        c: 截断参数（惩罚基数误差的上限）
        p: 距离范数阶数（通常为2）
        
    Returns:
        (平均OSPA值, 每时刻的OSPA值列表, 每时刻的基数误差列表)
    """
    n_times = len(true_sets)
    ospa_values = []
    cardinality_errors = []
    
    for t in range(n_times):
        true_set = true_sets[t]
        est_set = estimated_sets[t]
        
        n_true = len(true_set) if true_set.ndim > 1 else 0
        n_est = len(est_set) if est_set.ndim > 1 else 0
        
        if n_true == 0 and n_est == 0:
            ospa_values.append(0.0)
            cardinality_errors.append(0.0)
            continue
        
        if n_true == 0 or n_est == 0:
            # 纯基数误差
            card_error = c * abs(n_true - n_est) ** (1.0 / p)
            ospa_values.append(card_error)
            cardinality_errors.append(card_error)
            continue
        
        # 计算距离矩阵
        distance_matrix = np.zeros((n_true, n_est))
        for i in range(n_true):
            for j in range(n_est):
                distance_matrix[i, j] = np.linalg.norm(true_set[i] - est_set[j], p)
        
        # 截断距离
        distance_matrix = np.minimum(distance_matrix, c)
        
        # 使用匈牙利算法求解最优分配
        row_indices, col_indices = linear_sum_assignment(distance_matrix)
        
        # 计算定位误差
        localization_error = 0.0
        if len(row_indices) > 0:
            localization_error = np.sum(distance_matrix[row_indices, col_indices] ** p)
        
        # 计算基数误差
        cardinality_error = c ** p * abs(n_true - n_est)
        
        # 计算OSPA
        total_error = (localization_error + cardinality_error) / max(n_true, n_est)
        ospa_value = total_error ** (1.0 / p)
        
        ospa_values.append(ospa_value)
        cardinality_errors.append(cardinality_error ** (1.0 / p))
    
    avg_ospa = np.mean(ospa_values) if ospa_values else 0.0
    
    return avg_ospa, ospa_values, cardinality_errors


def GOSPA(true_sets: List[np.ndarray],
          estimated_sets: List[np.ndarray],
          alpha: float = 2.0,
          c: float = 100.0,
          p: float = 2.0) -> Tuple[float, List[float], List[float], List[float]]:
    """计算广义最优子模式分配 (GOSPA) 指标
    
    GOSPA是OSPA的改进版本，将误差分解为：
    - 定位误差
    - 漏检误差（真实目标未被检测到）
    - 虚警误差（估计目标不存在）
    
    Args:
        true_sets: 真实目标集合列表
        estimated_sets: 估计目标集合列表
        alpha: 权衡参数（控制基数误差的重要性）
        c: 截断参数
        p: 距离范数阶数
        
    Returns:
        (平均GOSPA值, 每时刻的GOSPA值, 每时刻的定位误差, 每时刻的基数误差)
    """
    n_times = len(true_sets)
    gospa_values = []
    localization_errors = []
    missed_errors = []
    false_alarm_errors = []
    
    for t in range(n_times):
        true_set = true_sets[t]
        est_set = estimated_sets[t]
        
        n_true = len(true_set) if true_set.ndim > 1 else 0
        n_est = len(est_set) if est_set.ndim > 1 else 0
        
        if n_true == 0 and n_est == 0:
            gospa_values.append(0.0)
            localization_errors.append(0.0)
            missed_errors.append(0.0)
            false_alarm_errors.append(0.0)
            continue
        
        # 计算距离矩阵
        if n_true > 0 and n_est > 0:
            distance_matrix = np.zeros((n_true, n_est))
            for i in range(n_true):
                for j in range(n_est):
                    distance_matrix[i, j] = np.linalg.norm(true_set[i] - est_set[j], p)
            
            # 使用匈牙利算法求解最优分配
            row_indices, col_indices = linear_sum_assignment(distance_matrix)
            
            # 计算分配的距离
            assigned_distances = distance_matrix[row_indices, col_indices]
            
            # 计算定位误差（截断）
            localization_error = np.sum(np.minimum(assigned_distances, c) ** p)
            
            # 计算漏检误差
            n_assigned = len(row_indices)
            n_missed = n_true - n_assigned
            missed_error = (c ** p) * n_missed
            
            # 计算虚警误差
            n_false_alarm = n_est - n_assigned
            false_alarm_error = (c ** p) * n_false_alarm
        else:
            localization_error = 0.0
            missed_error = (c ** p) * n_true
            false_alarm_error = (c ** p) * n_est
        
        # 计算GOSPA
        total_error = localization_error + alpha * (missed_error + false_alarm_error)
        gospa_value = (total_error / max(n_true, n_est)) ** (1.0 / p) if max(n_true, n_est) > 0 else 0.0
        
        gospa_values.append(gospa_value)
        localization_errors.append(localization_error ** (1.0 / p))
        missed_errors.append(missed_error ** (1.0 / p))
        false_alarm_errors.append(false_alarm_error ** (1.0 / p))
    
    avg_gospa = np.mean(gospa_values) if gospa_values else 0.0
    
    return avg_gospa, gospa_values, localization_errors, missed_errors, false_alarm_errors


def compute_ospa_for_single_pair(true_set: np.ndarray,
                                  estimated_set: np.ndarray,
                                  c: float = 100.0,
                                  p: float = 2.0) -> Tuple[float, float, float]:
    """计算单对目标集合的OSPA
    
    Args:
        true_set: 真实目标集合，形状为 (n_true, dim)
        estimated_set: 估计目标集合，形状为 (n_est, dim)
        c: 截断参数
        p: 距离范数阶数
        
    Returns:
        (OSPA值, 定位误差, 基数误差)
    """
    n_true = len(true_set) if true_set.ndim > 1 else 0
    n_est = len(estimated_set) if estimated_set.ndim > 1 else 0
    
    if n_true == 0 and n_est == 0:
        return 0.0, 0.0, 0.0
    
    if n_true == 0 or n_est == 0:
        card_error = c * abs(n_true - n_est) ** (1.0 / p)
        return card_error, 0.0, card_error
    
    # 计算距离矩阵
    distance_matrix = np.zeros((n_true, n_est))
    for i in range(n_true):
        for j in range(n_est):
            distance_matrix[i, j] = np.linalg.norm(true_set[i] - estimated_set[j], p)
    
    # 截断距离
    distance_matrix = np.minimum(distance_matrix, c)
    
    # 使用匈牙利算法求解最优分配
    row_indices, col_indices = linear_sum_assignment(distance_matrix)
    
    # 计算定位误差
    localization_error = 0.0
    if len(row_indices) > 0:
        localization_error = np.sum(distance_matrix[row_indices, col_indices] ** p)
    
    # 计算基数误差
    cardinality_error = c ** p * abs(n_true - n_est)
    
    # 计算OSPA
    total_error = (localization_error + cardinality_error) / max(n_true, n_est)
    ospa_value = total_error ** (1.0 / p)
    
    return ospa_value, localization_error ** (1.0 / p), cardinality_error ** (1.0 / p)


class MetricTracker:
    """指标跟踪器
    
    用于累积计算和跟踪性能指标
    """
    
    def __init__(self):
        """初始化指标跟踪器"""
        self.rmse_values = []
        self.ospa_values = []
        self.gospa_values = []
        self.localization_errors = []
        self.cardinality_errors = []
        self.missed_errors = []
        self.false_alarm_errors = []
    
    def update(self,
               true_set: np.ndarray,
               estimated_set: np.ndarray,
               ospa_c: float = 100.0,
               ospa_p: float = 2.0,
               gospa_alpha: float = 2.0,
               gospa_c: float = 100.0,
               gospa_p: float = 2.0) -> None:
        """更新指标
        
        Args:
            true_set: 真实目标集合
            estimated_set: 估计目标集合
            ospa_c: OSPA截断参数
            ospa_p: OSPA距离范数阶数
            gospa_alpha: GOSPA权衡参数
            gospa_c: GOSPA截断参数
            gospa_p: GOSPA距离范数阶数
        """
        # 计算RMSE
        n_true = len(true_set) if true_set.ndim > 1 else 0
        n_est = len(estimated_set) if estimated_set.ndim > 1 else 0
        
        if n_true > 0 and n_est > 0:
            # 使用最近邻匹配计算RMSE
            distance_matrix = np.zeros((n_true, n_est))
            for i in range(n_true):
                for j in range(n_est):
                    distance_matrix[i, j] = np.linalg.norm(true_set[i] - estimated_set[j])
            
            row_indices, col_indices = linear_sum_assignment(distance_matrix)
            
            if len(row_indices) > 0:
                rmse = np.sqrt(np.mean(distance_matrix[row_indices, col_indices] ** 2))
                self.rmse_values.append(rmse)
        
        # 计算OSPA
        ospa, loc_err, card_err = compute_ospa_for_single_pair(
            true_set, estimated_set, ospa_c, ospa_p
        )
        self.ospa_values.append(ospa)
        self.localization_errors.append(loc_err)
        self.cardinality_errors.append(card_err)
        
        # 计算GOSPA
        gospa, gospa_values, loc_errors, missed, false_alarm = GOSPA(
            [true_set], [estimated_set], gospa_alpha, gospa_c, gospa_p
        )
        self.gospa_values.append(gospa_values[0])
        self.missed_errors.append(missed[0])
        self.false_alarm_errors.append(false_alarm[0])
    
    def get_summary(self) -> dict:
        """获取指标摘要
        
        Returns:
            指摘摘要字典
        """
        summary = {}
        
        if self.rmse_values:
            summary['rmse_mean'] = np.mean(self.rmse_values)
            summary['rmse_std'] = np.std(self.rmse_values)
        
        if self.ospa_values:
            summary['ospa_mean'] = np.mean(self.ospa_values)
            summary['ospa_std'] = np.std(self.ospa_values)
        
        if self.gospa_values:
            summary['gospa_mean'] = np.mean(self.gospa_values)
            summary['gospa_std'] = np.std(self.gospa_values)
        
        if self.localization_errors:
            summary['localization_error_mean'] = np.mean(self.localization_errors)
        
        if self.cardinality_errors:
            summary['cardinality_error_mean'] = np.mean(self.cardinality_errors)
        
        if self.missed_errors:
            summary['missed_error_mean'] = np.mean(self.missed_errors)
        
        if self.false_alarm_errors:
            summary['false_alarm_error_mean'] = np.mean(self.false_alarm_errors)
        
        return summary
    
    def reset(self) -> None:
        """重置指标跟踪器"""
        self.rmse_values = []
        self.ospa_values = []
        self.gospa_values = []
        self.localization_errors = []
        self.cardinality_errors = []
        self.missed_errors = []
        self.false_alarm_errors = []
