"""
联合概率数据关联 (JPDA)

优化特性：
1. Murty's k-best 分配算法：高效找到k个最优分配，避免枚举所有假设
2. 对数空间计算：全程log域避免数值下溢
3. 假设管理：prune/cap/merge 控制假设数量
"""
import heapq
import numpy as np
from typing import Optional, List, Dict, Set, Tuple
from collections import defaultdict
from dataclasses import dataclass
from scipy.optimize import linear_sum_assignment
from .base_association import BaseAssociation, AssociationResult


# =========================================================================
# Murty's k-best 2D 分配算法
# =========================================================================


def k_best_2d_assignments(cost_matrix: np.ndarray, k: int) -> List[Tuple[float, np.ndarray, np.ndarray]]:
    """Murty's k-best 2D分配算法
    
    高效寻找k个最小代价的分配方案。核心思想：
    1. 用匈牙利算法找到最优分配
    2. 将解空间逐层划分（forbid一对+fix之前所有对）
    3. 用优先队列管理子问题，每次取最优解
    4. 重复直到找到k个解或空间耗尽
    
    对数空间：代价矩阵和返回的代价均在log域（越小越好）
    
    Args:
        cost_matrix: (n_rows, n_cols) 代价矩阵（log域，小=好）
        k: 返回的最优分配数量
        
    Returns:
        List of (total_cost, row_indices, column_indices) 按代价升序排列
        其中 row_indices, column_indices 是匈牙利算法的行/列分配索引
    """
    n_rows, n_cols = cost_matrix.shape
    
    def _solve_with_constraints(C: np.ndarray,
                                fixed_rows: Dict[int, int],
                                fixed_cols: Dict[int, int],
                                forbidden: Set[Tuple[int, int]]) -> Optional[Tuple[float, np.ndarray, np.ndarray]]:
        """在固定分配和禁止分配约束下求解2D分配问题"""
        C_mod = C.astype(np.float64, copy=True) if not np.issubdtype(C.dtype, np.floating) else C.copy()
        
        # 应用固定分配：固定行和列只能配对彼此
        for r, c in fixed_rows.items():
            C_mod[r, :] = np.inf
            C_mod[:, c] = np.inf
            C_mod[r, c] = C[r, c]  # 恢复固定对的原始代价
        
        # 应用禁止分配
        for r, c in forbidden:
            if r not in fixed_rows and c not in fixed_cols:
                C_mod[r, c] = np.inf
        
        # 可行性检查：每行每列至少有一个有限值
        if not np.all(np.any(np.isfinite(C_mod), axis=1)):
            return None
        if not np.all(np.any(np.isfinite(C_mod), axis=0)):
            return None
        
        try:
            row_ind, col_ind = linear_sum_assignment(C_mod)
            total_cost = C_mod[row_ind, col_ind].sum()
            if not np.isfinite(total_cost):
                return None
            return total_cost, row_ind, col_ind
        except (ValueError, RuntimeError):
            return None
    
    # ---- 步骤1: 求解原始问题 ----
    result = _solve_with_constraints(cost_matrix, {}, {}, set())
    if result is None:
        return []
    
    best_cost, best_row, best_col = result
    unique_id = 0
    
    # 结果列表: [(cost, row_ind, col_ind), ...]
    results = [(best_cost, best_row.copy(), best_col.copy())]
    
    # 最小堆: (cost, unique_id, fixed_rows, fixed_cols, forbidden, row_ind, col_ind)
    heap = []
    
    # ---- 划分原始最优解 ----
    fixed_rows: Dict[int, int] = {}
    fixed_cols: Dict[int, int] = {}
    
    for idx in range(len(best_row)):
        r = int(best_row[idx])
        c = int(best_col[idx])
        
        # 子问题: forbid (r,c), 保持之前所有 fixed
        child_forbidden = {(r, c)}
        child_result = _solve_with_constraints(
            cost_matrix, fixed_rows, fixed_cols, child_forbidden
        )
        if child_result is not None:
            child_cost, child_row, child_col = child_result
            heapq.heappush(heap, (
                child_cost, unique_id,
                dict(fixed_rows), dict(fixed_cols),
                child_forbidden, child_row, child_col
            ))
            unique_id += 1
        
        # 将 (r,c) 加入固定集，继续划分
        fixed_rows[r] = c
        fixed_cols[c] = r
    
    # ---- 步骤2: 从堆中依次提取k个最优 ----
    while len(results) < k and heap:
        item_cost, _uid, item_fr, item_fc, item_fb, item_row, item_col = heapq.heappop(heap)
        results.append((item_cost, item_row.copy(), item_col.copy()))
        
        # 进一步划分这个解
        new_fixed_rows = dict(item_fr)
        new_fixed_cols = dict(item_fc)
        
        for idx in range(len(item_row)):
            r = int(item_row[idx])
            c = int(item_col[idx])
            
            if r in new_fixed_rows:
                continue  # 已固定的分配不需要再划分
            
            # 子问题: 在已有禁止集基础上，额外 forbid (r,c)
            child_forbidden = set(item_fb) | {(r, c)}
            child_result = _solve_with_constraints(
                cost_matrix, new_fixed_rows, new_fixed_cols, child_forbidden
            )
            if child_result is not None:
                child_cost, child_row, child_col = child_result
                heapq.heappush(heap, (
                    child_cost, unique_id,
                    dict(new_fixed_rows), dict(new_fixed_cols),
                    child_forbidden, child_row, child_col
                ))
                unique_id += 1
            
            # 将 (r,c) 加入新固定集
            new_fixed_rows[r] = c
            new_fixed_cols[c] = r
    
    return results


class JPDAFilter(BaseAssociation):
    """联合概率数据关联滤波器
    
    JPDA考虑所有可能的关联假设，计算每个目标的关联概率
    
    特点：
    - 考虑所有可能的关联假设
    - 计算联合关联概率
    - 适用于目标密集的场景
    - 计算复杂度较高
    """
    
    def __init__(self, 
                 gating_threshold: float = 9.21,
                 use_mahalanobis: bool = True,
                 detection_probability: float = 0.9,
                 clutter_density: float = 1e-4):
        """
        初始化JPDA滤波器
        
        Args:
            gating_threshold: 门限阈值
            use_mahalanobis: 是否使用马氏距离
            detection_probability: 检测概率
            clutter_density: 杂波密度
        """
        super().__init__(gating_threshold, use_mahalanobis)
        self.detection_probability = detection_probability
        self.clutter_density = clutter_density
    
    def associate(self, 
                  measurements: np.ndarray,
                  predicted_measurements: np.ndarray,
                  measurement_covariances: Optional[List[np.ndarray]] = None,
                  innovation_covariances: Optional[List[np.ndarray]] = None) -> AssociationResult:
        """执行JPDA关联
        
        Args:
            measurements: 观测矩阵
            predicted_measurements: 预测观测矩阵
            measurement_covariances: 观测噪声协方差列表
            innovation_covariances: 新息协方差列表
            
        Returns:
            关联结果（包含关联概率）
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
        
        # 获取有效关联（门限内的关联）
        valid_associations = self._get_valid_associations(association_matrix)
        
        # 计算关联概率
        association_probs = self._compute_association_probabilities(
            n_meas, n_targets, valid_associations,
            association_matrix, innovation_covariances
        )
        
        # 构建关联结果
        associations = {}
        unassociated_measurements = set(range(n_meas))
        unassociated_targets = set(range(n_targets))
        
        # 对每个目标，选择概率最高的观测
        for j in range(n_targets):
            best_meas = None
            best_prob = 0.0
            
            for i in range(n_meas):
                if (i, j) in association_probs:
                    prob = association_probs[(i, j)]
                    if prob > best_prob:
                        best_prob = prob
                        best_meas = i
            
            if best_meas is not None and best_prob > 0.5:
                associations[best_meas] = j
                unassociated_measurements.discard(best_meas)
                unassociated_targets.discard(j)
        
        return AssociationResult(
            associations=associations,
            unassociated_measurements=unassociated_measurements,
            unassociated_targets=unassociated_targets,
            association_matrix=association_matrix
        )
    
    def _get_valid_associations(self, 
                                 association_matrix: np.ndarray) -> List[Tuple[int, int]]:
        """获取有效关联
        
        Args:
            association_matrix: 关联矩阵
            
        Returns:
            有效关联列表 [(观测索引, 目标索引), ...]
        """
        valid = []
        n_meas, n_targets = association_matrix.shape
        
        for i in range(n_meas):
            for j in range(n_targets):
                if association_matrix[i, j] < np.inf:
                    valid.append((i, j))
        
        return valid
    
    def _compute_log_likelihood_matrix(
        self,
        n_meas: int,
        n_targets: int,
        association_matrix: np.ndarray,
        innovation_covariances: Optional[List[np.ndarray]] = None,
    ) -> np.ndarray:
        """向量化计算对数似然矩阵
        
        对于每个观测-目标对 (i, j)，计算:
            log(N(z_i; Hx_j, S_j)) = -0.5*d² - 0.5*log(det(S_j))
        
        Args:
            n_meas: 观测数量
            n_targets: 目标数量
            association_matrix: 关联矩阵 (n_meas, n_targets)，元素为马氏距离或inf
            innovation_covariances: 新息协方差列表
            
        Returns:
            对数似然矩阵 (n_meas, n_targets)
        """
        # 预计算每个目标的新息协方差行列式的对数 log(det(S_j))
        log_det_s = np.zeros(n_targets)
        if innovation_covariances is not None:
            for j in range(min(n_targets, len(innovation_covariances))):
                det_val = np.linalg.det(innovation_covariances[j])
                if det_val > 0:
                    log_det_s[j] = np.log(det_val)
        
        # 向量化计算对数似然矩阵
        # 使用NumPy广播: association_matrix 是 (n_meas, n_targets)
        # 对每个元素计算 -0.5 * d² - 0.5 * log_det_s[j]
        finite_mask = np.isfinite(association_matrix)
        log_likelihood = np.full((n_meas, n_targets), -np.inf)
        log_likelihood[finite_mask] = (
            -0.5 * association_matrix[finite_mask] ** 2
            - 0.5 * log_det_s[np.where(finite_mask)[1]]
        )
        
        return log_likelihood
    
    def _cluster_targets(
        self,
        n_targets: int,
        valid_associations: List[Tuple[int, int]],
    ) -> List[List[int]]:
        """将目标分组为独立簇
        
        如果两个目标共享同一个观测（同一观测在两个目标的门限内），
        则它们属于同一个簇。不同簇之间完全独立，可以分别计算关联概率。
        
        Args:
            n_targets: 目标数量
            valid_associations: 有效关联列表 [(观测索引, 目标索引)]
            
        Returns:
            簇列表，每个簇是目标索引的列表
        """
        # 构建观测 → 目标列表 的映射
        meas_to_targets = defaultdict(set)
        for i, j in valid_associations:
            meas_to_targets[i].add(j)
        
        # 构建目标邻接图: 共享观测的目标相连
        target_adj = defaultdict(set)
        for targets in meas_to_targets.values():
            targets_list = list(targets)
            for idx1 in range(len(targets_list)):
                t1 = targets_list[idx1]
                for idx2 in range(idx1 + 1, len(targets_list)):
                    t2 = targets_list[idx2]
                    target_adj[t1].add(t2)
                    target_adj[t2].add(t1)
        
        # BFS 查找连通分量
        visited = set()
        clusters = []
        all_targets = set(range(n_targets))
        
        for j in range(n_targets):
            if j in visited:
                continue
            component = []
            stack = [j]
            visited.add(j)
            while stack:
                node = stack.pop()
                component.append(node)
                for neighbor in target_adj.get(node, set()):
                    if neighbor not in visited:
                        visited.add(neighbor)
                        stack.append(neighbor)
            clusters.append(component)
        
        # 确保所有目标都被覆盖（孤立目标没有有效关联）
        # visited 已经包含了所有可达目标，但 {孤立目标} ∩ visited 可能不全
        for j in range(n_targets):
            if j not in visited:
                clusters.append([j])
                visited.add(j)
        
        return clusters
    
    def _get_cluster_measurements(
        self,
        cluster_targets: List[int],
        valid_associations: List[Tuple[int, int]],
    ) -> List[int]:
        """获取与簇中任意目标有关联的所有观测
        
        Args:
            cluster_targets: 簇中的目标索引列表
            valid_associations: 有效关联列表
            
        Returns:
            观测索引列表（排序去重）
        """
        target_set = set(cluster_targets)
        measurements = sorted({i for i, j in valid_associations if j in target_set})
        return measurements
    
    def _exact_cluster_probabilities(
        self,
        cluster_targets: List[int],
        cluster_measurements: List[int],
        log_likelihood: np.ndarray,
        log_pd: float,
        log_1mpd: float,
        log_lambda: float,
    ) -> Dict[Tuple[int, int], float]:
        """对小型簇进行精确枚举计算边际关联概率
        
        使用迭代式动态规划构建所有有效假设，避免递归开销。
        适用于目标数 ≤ 6 的簇。
        
        Args:
            cluster_targets: 簇中的目标索引（全局索引）
            cluster_measurements: 簇中的观测索引（全局索引）
            log_likelihood: 全局对数似然矩阵
            log_pd: log(P_D)
            log_1mpd: log(1-P_D)  
            log_lambda: log(λ)
            
        Returns:
            边际关联概率字典 {(观测索引, 目标索引): 概率}
        """
        n_targets = len(cluster_targets)
        n_meas = len(cluster_measurements)
        
        if n_targets == 0 or n_meas == 0:
            return {}
        
        # 构建全局→局部索引映射
        target_set = set(cluster_targets)
        meas_set = set(cluster_measurements)
        
        # 预提取子矩阵: local_log_likelihood[local_i, local_j]
        local_ll = np.full((n_meas, n_targets), -np.inf)
        for local_j, global_j in enumerate(cluster_targets):
            for local_i, global_i in enumerate(cluster_measurements):
                local_ll[local_i, local_j] = log_likelihood[global_i, global_j]
        
        # 构建每个目标的可用观测列表（局部索引）
        gates = [[] for _ in range(n_targets)]
        for local_j in range(n_targets):
            for local_i in range(n_meas):
                if local_ll[local_i, local_j] > -np.inf:
                    gates[local_j].append(local_i)
        
        # 迭代式假设枚举（使用列表扩展，避免递归）
        # 每个假设表示为数组: assignment[local_j] = local_i (>=0) or -1 (未检测)
        # 使用 set 记录已用观测以强制一对一约束
        
        hypotheses_list = [np.full(n_targets, -1, dtype=np.int32)]
        used_sets = [set()]
        
        for local_j in range(n_targets):
            meas_options = gates[local_j]
            new_hypotheses = []
            new_used_sets = []
            
            for hyp_idx in range(len(hypotheses_list)):
                hyp = hypotheses_list[hyp_idx]
                used = used_sets[hyp_idx]
                
                # 分支1: 未检测
                new_hyp = hyp.copy()
                new_hyp[local_j] = -1
                new_hypotheses.append(new_hyp)
                new_used_sets.append(used.copy())
                
                # 分支2: 每个可用观测
                for local_i in meas_options:
                    if local_i not in used:
                        new_hyp = hyp.copy()
                        new_hyp[local_j] = local_i
                        new_hypotheses.append(new_hyp)
                        new_used_sets.append(used | {local_i})
            
            hypotheses_list = new_hypotheses
            used_sets = new_used_sets
        
        n_hypotheses = len(hypotheses_list)
        if n_hypotheses == 0:
            return {}
        
        # 向量化计算所有假设的对数概率
        # 构建分配矩阵 (n_hypotheses, n_targets)
        assign_matrix = np.stack(hypotheses_list, axis=0)  # (n_hyp, n_targets)
        
        # δ = 每个假设中已分配的目标数
        delta = np.sum(assign_matrix >= 0, axis=1)  # (n_hyp,)
        
        # 计算对数概率: λ^φ × P_D^δ × (1-P_D)^(T-δ) × ∏ L
        log_probs = np.zeros(n_hypotheses)
        
        # φ = n_meas - δ  (簇内虚警数)
        phi = n_meas - delta
        if log_lambda > -np.inf:
            log_probs += phi * log_lambda
        else:
            log_probs[phi > 0] = -np.inf  # λ=0且有虚警时概率为0
        
        if log_pd > -np.inf:
            log_probs += delta * log_pd
        
        undetected = n_targets - delta
        if log_1mpd > -np.inf:
            log_probs += undetected * log_1mpd
        
        # 高斯似然乘积: 对每个假设，累加已分配观测的对数似然
        # 使用向量化索引
        valid_mask = assign_matrix >= 0  # (n_hyp, n_targets)
        if np.any(valid_mask):
            for local_j in range(n_targets):
                col_assignments = assign_matrix[:, local_j]  # (n_hyp,)
                col_valid = col_assignments >= 0
                if np.any(col_valid):
                    log_probs[col_valid] += local_ll[col_assignments[col_valid], local_j]
        
        # 归一化（对数空间 → 概率空间，防止下溢）
        max_log = np.max(log_probs)
        if max_log == -np.inf or not np.isfinite(max_log):
            return {}
        
        probs = np.exp(log_probs - max_log)
        total_prob = np.sum(probs)
        
        if total_prob <= 0:
            return {}
        
        probs /= total_prob
        
        # 计算边际关联概率
        # β_{global_i, global_j} = Σ_{假设中 global_j 分配了 global_i} prob(假设)
        marginal_probs = {}
        
        for local_j, global_j in enumerate(cluster_targets):
            col_assignments = assign_matrix[:, local_j]  # (n_hyp,)
            for local_i, global_i in enumerate(cluster_measurements):
                if local_ll[local_i, local_j] <= -np.inf:
                    continue  # 无效关联
                match_mask = col_assignments == local_i
                if np.any(match_mask):
                    beta = np.sum(probs[match_mask])
                    if beta > 1e-10:
                        marginal_probs[(global_i, global_j)] = beta
        
        return marginal_probs
    
    def _gibbs_cluster_probabilities(
        self,
        cluster_targets: List[int],
        cluster_measurements: List[int],
        log_likelihood: np.ndarray,
        log_pd: float,
        log_1mpd: float,
        log_lambda: float,
    ) -> Dict[Tuple[int, int], float]:
        """使用Gibbs采样计算大型簇的边际关联概率（优化版）
        
        关键优化:
        1. 自适应迭代次数: 根据簇大小调整采样量
        2. 预计算权重矩阵: 避免在热循环中重复计算
        3. 布尔数组替代set: O(1)可用性检查
        4. 预生成随机排列: 减少每轮开销
        
        Args:
            cluster_targets: 簇中的目标索引（全局索引）
            cluster_measurements: 簇中的观测索引（全局索引）
            log_likelihood: 全局对数似然矩阵
            log_pd: log(P_D)
            log_1mpd: log(1-P_D)
            log_lambda: log(λ)
            
        Returns:
            边际关联概率字典 {(观测索引, 目标索引): 概率}
        """
        n_targets = len(cluster_targets)
        n_meas = len(cluster_measurements)
        
        if n_targets == 0 or n_meas == 0:
            return {}
        
        # 索引映射
        meas_to_global = {local_i: global_i for local_i, global_i in enumerate(cluster_measurements)}
        
        # --- 预计算局部矩阵 ---
        local_ll = np.full((n_meas, n_targets), -np.inf)
        for local_j, global_j in enumerate(cluster_targets):
            for local_i, global_i in enumerate(cluster_measurements):
                local_ll[local_i, local_j] = log_likelihood[global_i, global_j]
        
        # 有效关联布尔掩码
        valid_mask = local_ll > -np.inf  # (n_meas, n_targets)
        
        # --- 预计算权重矩阵 ---
        # weight_mat[j, i] = log(P(θ[j]=i | θ_{-j}) 的对数非归一化权重)
        # 形状: (n_targets, n_meas + 1), 列0 = 未检测, 列1..n_meas = 观测0..n_meas-1
        # 无效关联填 -np.inf
        
        if log_lambda > -np.inf:
            log_lambda_eff = log_lambda
        else:
            log_lambda_eff = -100.0  # λ=0: 测量权重极大
        
        weight_mat = np.full((n_targets, n_meas + 1), -np.inf, dtype=np.float64)
        weight_mat[:, 0] = log_1mpd  # 未检测
        
        # 向量化填充有效关联的权重
        for local_j in range(n_targets):
            valid_row = np.where(valid_mask[:, local_j])[0]
            if len(valid_row) > 0:
                weight_mat[local_j, valid_row + 1] = (
                    log_pd + local_ll[valid_row, local_j] - log_lambda_eff
                )
        
        # --- 自适应采样参数 ---
        n_samples = max(500, min(5000, n_targets * 120))
        n_burnin = max(200, n_samples // 4)
        n_total = n_samples + n_burnin
        
        # --- 预生成随机排列（大幅减少Python层开销） ---
        all_orders = [np.random.permutation(n_targets).astype(np.int32) 
                      for _ in range(n_total)]
        
        # --- 初始化状态 ---
        assignment = np.full(n_targets, -1, dtype=np.int32)
        used = np.zeros(n_meas, dtype=bool)  # 布尔数组: O(1)可用性检查
        counts = np.zeros((n_meas, n_targets))
        
        # hard_neg 用于屏蔽不可用观测: (n_meas+1,) 列0=未检测总是可用
        hard_neg = np.full(n_meas + 1, -np.inf)
        hard_neg[0] = 0.0  # 不屏蔽未检测
        
        # --- 主采样循环 ---
        for iteration in range(n_total):
            order = all_orders[iteration]
            
            for idx in range(n_targets):
                local_j = order[idx]
                
                # 移除旧分配
                old_assign = assignment[local_j]
                if old_assign >= 0:
                    used[old_assign] = False
                
                # 获取该目标的预计算权重行
                w_row = weight_mat[local_j]  # (n_meas+1,), 视图
                
                # 屏蔽已用的观测: 构建临时权重数组
                # 将已用观测的对应列设为 -inf
                if np.any(used):
                    w = w_row.copy()
                    w[1:][used] = -np.inf
                else:
                    w = w_row
                
                # 从条件分布采样
                w_max = np.max(w)
                if w_max <= -np.inf or not np.isfinite(w_max):
                    # 所有候选权重无效 → 均匀采样
                    valid_candidates = np.where(np.isfinite(w))[0]
                    if len(valid_candidates) == 0:
                        chosen = -1  # 退化为未检测
                    else:
                        chosen_idx = valid_candidates[
                            np.random.randint(len(valid_candidates))
                        ]
                        chosen = chosen_idx - 1 if chosen_idx > 0 else -1
                else:
                    # softmax 归一化
                    probs = np.exp(w - w_max)
                    probs /= np.sum(probs)
                    # 使用累积分布采样（比 np.random.choice 更快）
                    r = np.random.random()
                    cumsum = 0.0
                    chosen_idx = 0
                    for ci in range(n_meas + 1):
                        cumsum += probs[ci]
                        if r < cumsum:
                            chosen_idx = ci
                            break
                    chosen = chosen_idx - 1 if chosen_idx > 0 else -1
                
                # 更新状态
                assignment[local_j] = chosen
                if chosen >= 0:
                    used[chosen] = True
            
            # 烧入期后收集样本
            if iteration >= n_burnin:
                for local_j in range(n_targets):
                    aj = assignment[local_j]
                    if aj >= 0:
                        counts[aj, local_j] += 1
        
        # --- 计算边际概率 ---
        marginal_probs = {}
        for local_j, global_j in enumerate(cluster_targets):
            valid_local_is = np.where(valid_mask[:, local_j])[0]
            for local_i in valid_local_is:
                prob = counts[local_i, local_j] / n_samples
                if prob > 1e-10:
                    marginal_probs[(meas_to_global[local_i], global_j)] = prob
        
        return marginal_probs
    
    def _kbest_cluster_probabilities(
        self,
        cluster_targets: List[int],
        cluster_measurements: List[int],
        log_likelihood: np.ndarray,
        log_pd: float,
        log_1mpd: float,
        log_lambda: float,
    ) -> Dict[Tuple[int, int], float]:
        """使用 k-best 分配计算簇的边际关联概率（核心优化）
        
        用Murty's k-best代替Gibbs采样或精确枚举：
        1. 构建方形代价矩阵 (M+T × M+T)，其中：
           - 左上 M×T = -log_likelihood - log_pd  (检测分配代价)
           - 右下 T×T 对角 = -log_1mpd  (漏检代价)
           - 右上 M×M 对角 = -log_lambda (杂波代价)
           - 其余 = Inf (无效分配)
        2. 运行Murty求k个最小代价分配（=最大后验概率）
        3. 从k个最好分配中计算边际关联概率
        
        Args:
            cluster_targets: 簇中的目标索引（全局索引）
            cluster_measurements: 簇中的观测索引（全局索引）
            log_likelihood: 全局对数似然矩阵
            log_pd: log(P_D)
            log_1mpd: log(1-P_D)
            log_lambda: log(λ)
            
        Returns:
            边际关联概率字典 {(观测索引, 目标索引): 概率}
        """
        n_targets = len(cluster_targets)
        n_meas = len(cluster_measurements)
        
        if n_targets == 0 or n_meas == 0:
            return {}
        
        # 索引映射
        meas_to_global = {local_i: global_i for local_i, global_i in enumerate(cluster_measurements)}
        target_to_global = {local_j: global_j for local_j, global_j in enumerate(cluster_targets)}
        
        # --- 构建局部对数似然子矩阵 ---
        local_ll = np.full((n_meas, n_targets), -np.inf)
        for local_j, global_j in enumerate(cluster_targets):
            for local_i, global_i in enumerate(cluster_measurements):
                local_ll[local_i, local_j] = log_likelihood[global_i, global_j]
        
        # --- 构建方形代价矩阵 (n_meas + n_targets) × (n_meas + n_targets) ---
        # 上半行 = 测量, 下半行 = 漏检虚拟行
        # 左半列 = 目标, 右半列 = 杂波虚拟列
        n_total = n_meas + n_targets
        C = np.full((n_total, n_total), np.inf)
        
        # 1) 检测分配: 测量i→目标j, 代价 = -log_likelihood - log_pd
        finite_mask = np.isfinite(local_ll)
        for local_i in range(n_meas):
            for local_j in range(n_targets):
                if finite_mask[local_i, local_j]:
                    C[local_i, local_j] = -local_ll[local_i, local_j] - log_pd
        
        # 2) 杂波: 测量i→杂波虚拟列T+i, 代价 = -log_lambda
        if np.isfinite(log_lambda):
            for local_i in range(n_meas):
                C[local_i, n_targets + local_i] = -log_lambda
        else:
            # λ=0: 不允许杂波
            pass  # 保留Inf
        
        # 3) 漏检: 虚拟行M+j→目标j, 代价 = -log(1-P_D)
        if np.isfinite(log_1mpd):
            for local_j in range(n_targets):
                C[n_meas + local_j, local_j] = -log_1mpd
        
        # --- 运行 k-best ---
        k = min(100, max(10, n_targets * n_meas * 2))
        k_results = k_best_2d_assignments(C, k)
        
        if len(k_results) == 0:
            return {}
        
        # --- 将k个分配转换回假设并计算概率 ---
        # 每个假设: assignments[local_j] = local_i 或 -1(漏检)
        # 使用对数概率累加
        hypotheses_log_probs = np.zeros(len(k_results))
        hypotheses_assign = np.full((len(k_results), n_targets), -1, dtype=np.int32)
        
        for h_idx, (cost, row_ind, col_ind) in enumerate(k_results):
            # 构建从col(row分配)到target的映射
            # 匈牙利返回的是 row→col 分配
            assigned_cols = set()
            for r, c in zip(row_ind, col_ind):
                r_int, c_int = int(r), int(c)
                if r_int < n_meas and c_int < n_targets:
                    # 测量→目标
                    hypotheses_assign[h_idx, c_int] = r_int
                    assigned_cols.add(c_int)
            
            # 计算对数概率（验证用，也可以用代价直接换算）
            # 实际上代价已经是最小化 -log(P)，所以 log_prob = -cost + const
            delta = np.sum(hypotheses_assign[h_idx] >= 0)  # 已分配目标数
            phi = n_meas - delta  # 杂波数
            undetected = n_targets - delta
            
            log_prob = 0.0
            if phi > 0 and np.isfinite(log_lambda):
                log_prob += phi * log_lambda
            if delta > 0 and np.isfinite(log_pd):
                log_prob += delta * log_pd
            if undetected > 0 and np.isfinite(log_1mpd):
                log_prob += undetected * log_1mpd
            
            # 累加分配测量的对数似然
            for local_j in range(n_targets):
                local_i = hypotheses_assign[h_idx, local_j]
                if local_i >= 0:
                    log_prob += local_ll[local_i, local_j]
            
            hypotheses_log_probs[h_idx] = log_prob
        
        # --- 对数空间归一化 ---
        max_log = np.max(hypotheses_log_probs)
        if not np.isfinite(max_log):
            return {}
        
        probs = np.exp(hypotheses_log_probs - max_log)
        total_prob = np.sum(probs)
        if total_prob <= 0:
            return {}
        probs /= total_prob
        
        # --- 计算边际关联概率 ---
        marginal_probs = {}
        for local_j, global_j in enumerate(cluster_targets):
            for local_i, global_i in enumerate(cluster_measurements):
                if not finite_mask[local_i, local_j]:
                    continue
                match_mask = hypotheses_assign[:, local_j] == local_i
                if np.any(match_mask):
                    beta = np.sum(probs[match_mask])
                    if beta > 1e-10:
                        marginal_probs[(global_i, global_j)] = beta
        
        return marginal_probs

    def _compute_association_probabilities(self,
                                            n_meas: int,
                                            n_targets: int,
                                            valid_associations: List[Tuple[int, int]],
                                            association_matrix: np.ndarray,
                                            innovation_covariances: Optional[List[np.ndarray]] = None) -> Dict[Tuple[int, int], float]:
        """计算JPDA联合关联概率（优化版本 v2）
        
        使用以下优化策略:
        1. 向量化计算所有关联对的似然（NumPy广播）
        2. 目标聚类: 将不共享观测的目标分成独立簇
        3. 小簇精确枚举（≤6个目标）: 迭代式DP代替递归
        4. 大簇 k-best 分配: 用Murty算法找k个最优分配近似边际概率
           （替代原Gibbs采样，避免采样噪声且更高效）
        5. 超大簇 Gibbs 采样（>20个目标）: 保留作为后备
        6. 全程对数空间避免数值下溢
        
        关联假设 θ 的后验概率:
            P(θ|Z) ∝ λ^φ × P_D^δ × (1-P_D)^{(n-δ)} × ∏_{j:θ_j>0} N(z_{θ_j}; Hx_j, S_j)
        
        Args:
            n_meas: 观测数量
            n_targets: 目标数量
            valid_associations: 有效关联列表 [(观测索引, 目标索引)]
            association_matrix: 关联矩阵 (n_meas, n_targets)
            innovation_covariances: 新息协方差列表
                
        Returns:
            关联概率字典，键为 (观测索引, 目标索引)，值为概率
        """
        # --- Step 1: 向量化计算对数似然矩阵 ---
        log_likelihood = self._compute_log_likelihood_matrix(
            n_meas, n_targets, association_matrix, innovation_covariances
        )
        
        # --- Step 2: 目标聚类 ---
        clusters = self._cluster_targets(n_targets, valid_associations)
        
        # 预计算对数参数
        log_pd = np.log(self.detection_probability) if self.detection_probability > 0 else -np.inf
        log_1mpd = np.log(1.0 - self.detection_probability) if self.detection_probability < 1.0 else -np.inf
        log_lambda = np.log(self.clutter_density) if self.clutter_density > 0 else -np.inf
        
        # --- Step 3: 每个簇独立计算 ---
        marginal_probs = {}
        
        for cluster_targets in clusters:
            # 获取簇内的观测
            cluster_measurements = self._get_cluster_measurements(
                cluster_targets, valid_associations
            )
            
            n_ct = len(cluster_targets)
            
            if n_ct <= 6:
                # 小簇: 精确枚举
                probs = self._exact_cluster_probabilities(
                    cluster_targets, cluster_measurements,
                    log_likelihood, log_pd, log_1mpd, log_lambda
                )
            elif n_ct <= 20:
                # 中簇: k-best 分配（Murty算法），高效且确定性
                probs = self._kbest_cluster_probabilities(
                    cluster_targets, cluster_measurements,
                    log_likelihood, log_pd, log_1mpd, log_lambda
                )
            else:
                # 超大簇: Gibbs采样（后备）
                probs = self._gibbs_cluster_probabilities(
                    cluster_targets, cluster_measurements,
                    log_likelihood, log_pd, log_1mpd, log_lambda
                )
            
            marginal_probs.update(probs)
        
        return marginal_probs
    
    def compute_joint_probabilities(self,
                                     n_meas: int,
                                     n_targets: int,
                                     valid_associations: List[Tuple[int, int]],
                                     distances: Dict[Tuple[int, int], float]) -> Dict[Tuple[int, int], float]:
        """计算联合关联概率（完整版本）
        
        Args:
            n_meas: 观测数量
            n_targets: 目标数量
            valid_associations: 有效关联列表
            distances: 距离字典
            
        Returns:
            关联概率字典
        """
        # 枚举所有可能的关联假设
        hypotheses = self._enumerate_hypotheses(n_meas, n_targets, valid_associations)
        
        # 计算每个假设的概率
        hypothesis_probs = []
        for hypothesis in hypotheses:
            prob = self._compute_hypothesis_probability(hypothesis, distances)
            hypothesis_probs.append(prob)
        
        # 归一化
        total_prob = sum(hypothesis_probs)
        if total_prob > 0:
            hypothesis_probs = [p / total_prob for p in hypothesis_probs]
        
        # 计算边际概率
        association_probs = {}
        for i in range(n_meas):
            for j in range(n_targets):
                prob = 0.0
                for k, hypothesis in enumerate(hypotheses):
                    if (i, j) in hypothesis:
                        prob += hypothesis_probs[k]
                if prob > 0:
                    association_probs[(i, j)] = prob
        
        return association_probs
    
    def _enumerate_hypotheses(self,
                               n_meas: int,
                               n_targets: int,
                               valid_associations: List[Tuple[int, int]]) -> List[List[Tuple[int, int]]]:
        """枚举所有可能的关联假设
        
        Args:
            n_meas: 观测数量
            n_targets: 目标数量
            valid_associations: 有效关联列表
            
        Returns:
            假设列表
        """
        # 简化版本：只考虑一对一关联
        hypotheses = []
        
        # 使用回溯法枚举
        def backtrack(current_hypothesis, used_meas, used_targets, start_idx):
            hypotheses.append(current_hypothesis.copy())
            
            for idx in range(start_idx, len(valid_associations)):
                i, j = valid_associations[idx]
                
                if i not in used_meas and j not in used_targets:
                    current_hypothesis.append((i, j))
                    used_meas.add(i)
                    used_targets.add(j)
                    
                    backtrack(current_hypothesis, used_meas, used_targets, idx + 1)
                    
                    current_hypothesis.pop()
                    used_meas.remove(i)
                    used_targets.remove(j)
        
        backtrack([], set(), set(), 0)
        
        return hypotheses
    
    def _compute_hypothesis_probability(self,
                                         hypothesis: List[Tuple[int, int]],
                                         distances: Dict[Tuple[int, int], float]) -> float:
        """计算假设概率
        
        Args:
            hypothesis: 关联假设
            distances: 距离字典
            
        Returns:
            假设概率
        """
        prob = 1.0
        
        for i, j in hypothesis:
            if (i, j) in distances:
                # 使用高斯似然
                dist = distances[(i, j)]
                prob *= np.exp(-0.5 * dist)
        
        return prob


class FastJPDA(JPDAFilter):
    """快速JPDA（基于 k-best 分配优化）
    
    使用 Murty's k-best 分配算法替代原始枚举：
    - 避免枚举所有指数级假设
    - 只保留k个最可能的分配用于边际概率计算
    - k自适应（默认=min(max_hypotheses, 20*n_targets)）
    """
    
    def __init__(self, 
                 gating_threshold: float = 9.21,
                 use_mahalanobis: bool = True,
                 detection_probability: float = 0.9,
                 clutter_density: float = 1e-4,
                 max_hypotheses: int = 100,
                 k_best_factor: int = 20):
        """
        初始化快速JPDA（k-best版）
        
        Args:
            gating_threshold: 门限阈值
            use_mahalanobis: 是否使用马氏距离
            detection_probability: 检测概率
            clutter_density: 杂波密度
            max_hypotheses: k-best最大假设数
            k_best_factor: k = min(max_hypotheses, n_targets * k_best_factor)
        """
        super().__init__(gating_threshold, use_mahalanobis, detection_probability, clutter_density)
        self.max_hypotheses = max_hypotheses
        self.k_best_factor = k_best_factor
    
    def _compute_association_probabilities(self,
                                            n_meas: int,
                                            n_targets: int,
                                            valid_associations: List[Tuple[int, int]],
                                            association_matrix: np.ndarray,
                                            innovation_covariances: Optional[List[np.ndarray]] = None) -> Dict[Tuple[int, int], float]:
        """用 k-best 分配计算关联概率（全程对数空间）
        
        直接构建全局代价矩阵，用Murty求k个最优分配，
        从分配中计算边际概率。完全避免假设枚举。
        """
        log_likelihood = self._compute_log_likelihood_matrix(
            n_meas, n_targets, association_matrix, innovation_covariances
        )
        
        clusters = self._cluster_targets(n_targets, valid_associations)
        
        log_pd = np.log(self.detection_probability) if self.detection_probability > 0 else -np.inf
        log_1mpd = np.log(1.0 - self.detection_probability) if self.detection_probability < 1.0 else -np.inf
        log_lambda = np.log(self.clutter_density) if self.clutter_density > 0 else -np.inf
        
        marginal_probs = {}
        
        for cluster_targets in clusters:
            cluster_measurements = self._get_cluster_measurements(
                cluster_targets, valid_associations
            )
            
            n_ct = len(cluster_targets)
            n_cm = len(cluster_measurements)
            
            # k自适应: 小簇多k, 大簇少k
            k = min(self.max_hypotheses, max(10, n_ct * self.k_best_factor, n_cm * 5))
            
            if n_ct <= 6:
                probs = self._exact_cluster_probabilities(
                    cluster_targets, cluster_measurements,
                    log_likelihood, log_pd, log_1mpd, log_lambda
                )
            else:
                probs = self._kbest_cluster_probabilities(
                    cluster_targets, cluster_measurements,
                    log_likelihood, log_pd, log_1mpd, log_lambda
                )
            
            marginal_probs.update(probs)
        
        return marginal_probs
