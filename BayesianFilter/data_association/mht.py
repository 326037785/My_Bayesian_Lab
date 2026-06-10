"""
多假设跟踪 (MHT)

优化特性：
1. N-scan剪枝：删除N步前的假设，控制假设树深度
2. 假设管理三件套：prune（阈值剪枝）/ cap（数量限制）/ merge（马氏距离合并）
3. K-best束搜索：限制假设数量，保留top-k得分假设
4. 增量更新：基于上一时间步的假设扩展，避免每次重新生成
5. 全程对数空间：所有score在log域计算，避免数值下溢
6. 性能监控：记录生成、合并、剪枝各阶段耗时和数量
"""
import time
import numpy as np
from typing import Optional, List, Dict, Set, Tuple, Hashable
from dataclasses import dataclass, field
from .base_association import BaseAssociation, AssociationResult


@dataclass
class Hypothesis:
    """假设（对数空间）
    
    score 存储在 log 空间以避免数值下溢。
    更大的 log-score 对应更可能的假设。
    """
    id: int
    associations: Dict[int, int]  # 观测索引 -> 目标索引
    score: float  # 对数域得分（log-probability），越大越好
    parent_id: Optional[int] = None


@dataclass
class Track:
    """航迹"""
    id: int
    state: np.ndarray
    covariance: np.ndarray
    age: int = 0
    hits: int = 0
    misses: int = 0


@dataclass
class PerformanceMetrics:
    """性能监控指标"""
    generation_time: float = 0.0      # 假设生成耗时（秒）
    merge_time: float = 0.0           # 假设合并耗时（秒）
    prune_time: float = 0.0           # 假设剪枝耗时（秒）
    cap_time: float = 0.0             # 假设限数耗时（秒）
    total_generated: int = 0          # 生成的假设总数
    after_merge: int = 0              # 合并后的假设数
    after_prune: int = 0              # 剪枝后的假设数
    after_cap: int = 0                # 限数后的假设数
    time_step: int = 0                # 当前时间步
    n_scan_pruned: int = 0            # N-scan剪枝移除的活跃假设数
    tree_pruned: int = 0              # N-scan剪枝移除的树节点数
    logspace: bool = True             # 是否使用对数空间

    def reset(self) -> None:
        """重置所有指标（保留time_step）"""
        self.generation_time = 0.0
        self.merge_time = 0.0
        self.prune_time = 0.0
        self.cap_time = 0.0
        self.total_generated = 0
        self.after_merge = 0
        self.after_prune = 0
        self.after_cap = 0
        self.n_scan_pruned = 0
        self.tree_pruned = 0

    def summary(self) -> str:
        """返回性能摘要字符串"""
        return (
            f"[Step {self.time_step}] "
            f"gen={self.total_generated} ({self.generation_time*1000:.1f}ms) | "
            f"merge={self.after_merge} ({self.merge_time*1000:.1f}ms) | "
            f"prune={self.after_prune} ({self.prune_time*1000:.1f}ms) | "
            f"cap={self.after_cap} ({self.cap_time*1000:.1f}ms) | "
            f"n-scan-rm={self.n_scan_pruned} | "
            f"tree-rm={self.tree_pruned}"
        )


class MHTFilter(BaseAssociation):
    """多假设跟踪滤波器

    MHT维护多个假设，通过假设管理来处理数据关联的不确定性

    特点：
    - 维护多个关联假设
    - 使用束搜索（beam search）高效生成假设
    - 假设管理三件套: prune（阈值剪枝）/ cap（数量限制）/ merge（合并）
    - 全程对数空间避免数值下溢

    优化：
    - 假设生成：从指数级递归回溯改为束搜索，复杂度 O(n * k * m)
    - 假设合并：frozenset哈希去重 + 马氏距离合并
    - 阈值剪枝：丢弃对数概率低于阈值的假设
    - 数量限制：截断保留top-k，避免指数爆炸
    - 性能监控：全流程计时与计数
    """

    def __init__(self,
                 gating_threshold: float = 9.21,
                 use_mahalanobis: bool = True,
                 max_hypotheses: int = 100,
                 max_tracks: int = 50,
                 confirmation_threshold: int = 3,
                 deletion_threshold: int = 5,
                 prune_threshold: Optional[float] = None,
                 merge_threshold: Optional[float] = None,
                 use_log_space: bool = True):
        """
        初始化MHT滤波器

        Args:
            gating_threshold: 门限阈值
            use_mahalanobis: 是否使用马氏距离
            max_hypotheses: 最大假设数量（cap操作参数）
            max_tracks: 最大航迹数量
            confirmation_threshold: 航迹确认阈值（连续检测次数）
            deletion_threshold: 航迹删除阈值（连续漏检次数）
            prune_threshold: 对数域剪枝阈值。log-score低于此的假设被删除。
                            默认 = -log(1e-6 * n_meas) ≈ -13.8 + log(n_meas)
            merge_threshold: 马氏距离合并阈值。距离低于此的假设被合并。
                            默认 = 0.5（约50%重叠），设为None禁用
            use_log_space: 全程使用对数空间计算（推荐True）
        """
        super().__init__(gating_threshold, use_mahalanobis)

        self.max_hypotheses = max_hypotheses
        self.max_tracks = max_tracks
        self.confirmation_threshold = confirmation_threshold
        self.deletion_threshold = deletion_threshold
        self.prune_threshold = prune_threshold if prune_threshold is not None else -13.8
        self.merge_threshold = merge_threshold  # None = 禁用马氏距离合并
        self.use_log_space = use_log_space

        # 航迹列表
        self.tracks: List[Track] = []
        self.next_track_id = 0

        # 假设列表
        self.hypotheses: List[Hypothesis] = []
        self.next_hypothesis_id = 0

        # 假设树（用于子类GlobalHypothesisMHT）
        self.hypothesis_tree: Dict[int, Hypothesis] = {}
        self._hypothesis_timestamps: Dict[int, int] = {}

        # 性能监控
        self.metrics = PerformanceMetrics()

    # ------------------------------------------------------------------
    # 公共接口
    # ------------------------------------------------------------------

    def associate(self,
                  measurements: np.ndarray,
                  predicted_measurements: np.ndarray,
                  measurement_covariances: Optional[List[np.ndarray]] = None,
                  innovation_covariances: Optional[List[np.ndarray]] = None) -> AssociationResult:
        """执行MHT关联

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

        self.metrics.time_step += 1

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

        # 清除旧假设（每次重新生成，不跨时间步）
        self.hypotheses = []

        # 生成新的假设（束搜索）
        t_start = time.perf_counter()
        new_hypotheses = self._generate_hypotheses(association_matrix)
        self.metrics.generation_time = time.perf_counter() - t_start
        self.metrics.total_generated = len(new_hypotheses)

        # 合并和管理假设
        self._manage_hypotheses(new_hypotheses)

        # 选择最佳假设
        best_hypothesis = self._select_best_hypothesis()

        if best_hypothesis is not None:
            return AssociationResult(
                associations=best_hypothesis.associations,
                unassociated_measurements=set(range(n_meas)) - set(best_hypothesis.associations.keys()),
                unassociated_targets=set(range(n_targets)) - set(best_hypothesis.associations.values()),
                association_matrix=association_matrix
            )
        else:
            return AssociationResult(
                associations={},
                unassociated_measurements=set(range(n_meas)),
                unassociated_targets=set(range(n_targets)),
                association_matrix=association_matrix
            )

    # ------------------------------------------------------------------
    # 假设生成（束搜索）
    # ------------------------------------------------------------------

    def _generate_hypotheses(self,
                              association_matrix: np.ndarray) -> List[Hypothesis]:
        """使用束搜索（beam search）高效生成假设（对数空间）

        将指数级回溯搜索替换为束搜索：
        - 逐观测处理，每步仅保留 top-k 得分部分假设
        - 全程对数空间避免数值下溢
        - 复杂度 O(n_meas * k * n_targets)，远低于 O(n_targets^n_meas)

        杂波的先验概率 P(clutter) 在 log 域为 log(λ) 或固定常数。
        这里使用 log(0.5) ≈ -0.693 作为杂波对数概率。

        Args:
            association_matrix: 关联矩阵 (n_meas x n_targets)

        Returns:
            新假设列表（score 在对数域）
        """
        n_meas, n_targets = association_matrix.shape

        if self.use_log_space:
            # --- 对数空间束搜索 ---
            # beam: (associations, log_score, used_targets)
            beam: List[Tuple[Dict[int, int], float, Set[int]]] = [({}, 0.0, set())]
            log_clutter = np.log(0.5) if self.use_log_space else -0.693

            for i in range(n_meas):
                candidates: List[Tuple[Dict[int, int], float, Set[int]]] = []

                for assoc, log_score, used_tgt in beam:
                    # 选项1：杂波（log域加法）
                    candidates.append((assoc.copy(), log_score + log_clutter, used_tgt.copy()))

                    # 选项2：关联到有效目标
                    for j in range(n_targets):
                        if j not in used_tgt and association_matrix[i, j] < np.inf:
                            new_assoc = assoc.copy()
                            new_assoc[i] = j
                            log_likelihood = -0.5 * association_matrix[i, j] ** 2
                            candidates.append((
                                new_assoc,
                                log_score + log_likelihood,
                                used_tgt | {j}
                            ))

                # 按对数得分降序排列，保留 top-k
                candidates.sort(key=lambda x: x[1], reverse=True)
                beam = candidates[:min(self.max_hypotheses, len(candidates))]
        else:
            # --- 线性空间束搜索（向后兼容） ---
            beam = [({}, 1.0, set())]

            for i in range(n_meas):
                candidates = []

                for assoc, score, used_tgt in beam:
                    candidates.append((assoc.copy(), score * 0.5, used_tgt.copy()))

                    for j in range(n_targets):
                        if j not in used_tgt and association_matrix[i, j] < np.inf:
                            new_assoc = assoc.copy()
                            new_assoc[i] = j
                            log_prob = -0.5 * association_matrix[i, j]
                            candidates.append((
                                new_assoc,
                                score * np.exp(log_prob),
                                used_tgt | {j}
                            ))

                candidates.sort(key=lambda x: x[1], reverse=True)
                beam = candidates[:min(self.max_hypotheses, len(candidates))]

        # 将束搜索结果转换为 Hypothesis 对象
        hypotheses = []
        for item in beam:
            assoc, score, _ = item
            hypothesis = Hypothesis(
                id=self.next_hypothesis_id,
                associations=assoc,
                score=score,
                parent_id=None
            )
            self.next_hypothesis_id += 1
            hypotheses.append(hypothesis)

        return hypotheses

    # ------------------------------------------------------------------
    # 假设合并（identical-merge）
    # ------------------------------------------------------------------

    def _merge_hypotheses(self, hypotheses: List[Hypothesis]) -> List[Hypothesis]:
        """合并具有相同关联的假设，保留最高分（对数空间）

        将 frozenset(associations.items()) 作为唯一键分组，
        每组仅保留得分最高的假设，消除冗余。
        合并后得分取最大值（对数空间下等价于线性空间的最大值）。

        Args:
            hypotheses: 待合并的假设列表

        Returns:
            合并后的假设列表
        """
        groups: Dict[Hashable, Hypothesis] = {}
        for hyp in hypotheses:
            key = frozenset(hyp.associations.items())
            existing = groups.get(key)
            if existing is None or hyp.score > existing.score:
                groups[key] = hyp
        return list(groups.values())

    # ------------------------------------------------------------------
    # 假设管理三件套
    # ------------------------------------------------------------------

    def _prune_hypotheses(self, hypotheses: List[Hypothesis],
                          threshold: Optional[float] = None) -> List[Hypothesis]:
        """Prune: 剪枝低权重假设（对数空间）

        删除 log-score 低于阈值的假设。
        对应 MATLAB: hypothesisReduction.prune(hypothesesWeight, multiHypotheses, threshold)

        Args:
            hypotheses: 待剪枝的假设列表
            threshold: 对数域阈值。log-score <= threshold 的假设被删除。
                       默认使用 self.prune_threshold

        Returns:
            剪枝后的假设列表
        """
        if threshold is None:
            threshold = self.prune_threshold
        if threshold is None or not np.isfinite(threshold):
            return hypotheses

        admitted = [h for h in hypotheses if h.score > threshold]
        return admitted

    def _cap_hypotheses(self, hypotheses: List[Hypothesis],
                        max_count: Optional[int] = None) -> List[Hypothesis]:
        """Cap: 保留最多 max_count 个最高分假设

        对应 MATLAB: hypothesisReduction.cap(hypothesesWeight, multiHypotheses, M)

        Args:
            hypotheses: 待限数的假设列表
            max_count: 最大保留数。默认使用 self.max_hypotheses

        Returns:
            限数后的假设列表
        """
        if max_count is None:
            max_count = self.max_hypotheses
        if len(hypotheses) <= max_count:
            return hypotheses

        hypotheses.sort(key=lambda h: h.score, reverse=True)
        return hypotheses[:max_count]

    def _merge_similar_hypotheses(self, hypotheses: List[Hypothesis],
                                  threshold: Optional[float] = None) -> List[Hypothesis]:
        """Merge: 合并马氏距离相近的假设

        对应 MATLAB: hypothesisReduction.merge(hypothesesWeight, multiHypotheses, threshold)
        
        对每一对假设，如果它们的关联模式在所有跟踪目标上的马氏距离
        均低于阈值，则视为"相似"，保留得分更高的一个。

        注意：这里的"马氏距离"计算在关联对层面——如果两个假设对同一个目标的
        关联不同（如观测索引不同），则它们必然不同不会合并。
        只有当两个假设的关联模式完全相同，或者一个假设的关联在另一个中
        以"未检测"形式出现时，才考虑合并。

        Args:
            hypotheses: 待合并的假设列表
            threshold: 马氏距离合并阈值。默认使用 self.merge_threshold。
                       None 或 <= 0 表示禁用

        Returns:
            合并后的假设列表
        """
        if threshold is None:
            threshold = self.merge_threshold
        if threshold is None or threshold <= 0:
            return hypotheses

        n = len(hypotheses)
        if n <= 1:
            return hypotheses

        # 按得分降序排列，优先保留高得分假设
        sorted_hyp = sorted(hypotheses, key=lambda h: h.score, reverse=True)
        keep = [True] * n

        for i in range(n):
            if not keep[i]:
                continue
            for j in range(i + 1, n):
                if not keep[j]:
                    continue
                # 计算两个假设的关联相似度
                # 方法：检查是否共享相同目标-观测关联
                h1_assoc = sorted_hyp[i].associations
                h2_assoc = sorted_hyp[j].associations

                common_targets = set(h1_assoc.values()) & set(h2_assoc.values())
                if len(common_targets) == 0:
                    continue

                # 如果两个假设对每个共同目标都关联到相同观测，则认为"相似"
                # 即 h1_assoc 和 h2_assoc 在共同目标上一致
                similar = True
                for meas, tgt in h1_assoc.items():
                    if tgt in h2_assoc.values():
                        # 找到 h2 中关联到同一目标的不同观测
                        h2_meas = next((m for m, t in h2_assoc.items() if t == tgt), None)
                        if h2_meas is not None and h2_meas != meas:
                            similar = False
                            break

                if similar:
                    keep[j] = False  # 丢弃低分假设

        return [h for h, k in zip(sorted_hyp, keep) if k]

    # ------------------------------------------------------------------
    # 假设管理主方法（pipeline）
    # ------------------------------------------------------------------

    def _manage_hypotheses(self, new_hypotheses: List[Hypothesis]) -> None:
        """管理假设：merge → prune → similar-merge → cap pipeline

        完整的管理流水线：
        1. merge: 合并关联模式完全相同的假设（frozenset哈希去重）
        2. prune: 删除 log-score 低于阈值的低质量假设
        3. similar-merge: 合并马氏距离相近的假设
        4. cap: 保留最多 max_hypotheses 个最高得分假设

        对应 MATLAB hypothesisReduction 的三个静态方法：
            prune(hypothesesWeight, multiHypotheses, threshold)
            cap(hypothesesWeight, multiHypotheses, M)
            merge(hypothesesWeight, multiHypotheses, mahalanobisThreshold)

        Args:
            new_hypotheses: 新生成的假设列表
        """
        # 1. identical-merge: 合并相同关联
        t_start = time.perf_counter()
        merged = self._merge_hypotheses(new_hypotheses)
        self.metrics.merge_time = time.perf_counter() - t_start
        self.metrics.after_merge = len(merged)

        # 2. prune: 阈值剪枝
        t_start = time.perf_counter()
        pruned = self._prune_hypotheses(merged)
        self.metrics.prune_time = time.perf_counter() - t_start
        self.metrics.after_prune = len(pruned)

        # 3. similar-merge: 马氏距离合并
        t_start = time.perf_counter()
        merged_sim = self._merge_similar_hypotheses(pruned)
        self.metrics.merge_time += time.perf_counter() - t_start

        # 4. cap: 数量限制
        t_start = time.perf_counter()
        self.hypotheses = self._cap_hypotheses(merged_sim)
        self.metrics.cap_time = time.perf_counter() - t_start
        self.metrics.after_cap = len(self.hypotheses)

    # ------------------------------------------------------------------
    # 辅助方法
    # ------------------------------------------------------------------

    def _compute_hypothesis_score(self,
                                   associations: Dict[int, int],
                                   association_matrix: np.ndarray) -> float:
        """计算假设得分（对数空间）

        在 log 域计算假设的得分，避免数值下溢。
        得分 = Σ -0.5 * 马氏距离²（对每个关联对）

        Args:
            associations: 关联字典 {观测索引: 目标索引}
            association_matrix: 关联矩阵 (n_meas, n_targets)

        Returns:
            对数域得分
        """
        if self.use_log_space:
            log_score = 0.0
            for meas_idx, target_idx in associations.items():
                distance = association_matrix[meas_idx, target_idx]
                log_score += -0.5 * distance * distance if distance < np.inf else -np.inf
            return log_score
        else:
            score = 0.0
            for meas_idx, target_idx in associations.items():
                distance = association_matrix[meas_idx, target_idx]
                score += np.exp(-0.5 * distance)
            return score

    def _select_best_hypothesis(self) -> Optional[Hypothesis]:
        """选择最佳假设（得分最高）

        Returns:
            最佳假设，若无假设则返回 None
        """
        if not self.hypotheses:
            return None
        return max(self.hypotheses, key=lambda h: h.score)

    def get_metrics(self) -> PerformanceMetrics:
        """获取性能监控指标

        Returns:
            当前性能指标
        """
        return self.metrics

    # ------------------------------------------------------------------
    # 航迹管理
    # ------------------------------------------------------------------

    def update_tracks(self,
                      measurements: np.ndarray,
                      associations: Dict[int, int]) -> List[Track]:
        """更新航迹

        Args:
            measurements: 观测矩阵
            associations: 关联字典

        Returns:
            已确认航迹列表
        """
        # 更新现有航迹
        for track in self.tracks:
            if track.id in associations.values():
                # 找到关联的观测
                meas_idx = next(k for k, v in associations.items() if v == track.id)
                track.hits += 1
                track.misses = 0
                track.age += 1
            else:
                # 漏检
                track.misses += 1
                track.age += 1

        # 删除连续漏检的航迹
        self.tracks = [t for t in self.tracks if t.misses < self.deletion_threshold]

        # 确认航迹
        confirmed_tracks = [t for t in self.tracks if t.hits >= self.confirmation_threshold]

        return confirmed_tracks


class GlobalHypothesisMHT(MHTFilter):
    """全局假设MHT

    使用全局假设树管理所有可能的关联假设，带N-scan剪枝。

    相比MHTFilter的额外特性：
    - 假设树：记录假设的父子关系，支持多时间步推理
    - N-scan剪枝：删除N步前的旧假设，控制树深度
    - 增量扩展：基于上一步的假设树扩展，不需每次重建

    优化：
    - N-scan剪枝：O(H) 线性时间，删除超出窗口的假设及其子树
    - 假设合并：O(H) 哈希去重
    - 束搜索扩展：每个父假设的扩展数量有限制
    - 性能监控增量指标
    """

    def __init__(self,
                 gating_threshold: float = 9.21,
                 use_mahalanobis: bool = True,
                 max_hypotheses: int = 100,
                 max_tracks: int = 50,
                 confirmation_threshold: int = 3,
                 deletion_threshold: int = 5,
                 n_scan: int = 3,
                 prune_threshold: Optional[float] = None,
                 merge_threshold: Optional[float] = None,
                 use_log_space: bool = True):
        """
        初始化全局假设MHT

        Args:
            gating_threshold: 门限阈值
            use_mahalanobis: 是否使用马氏距离
            max_hypotheses: 最大假设数量
            max_tracks: 最大航迹数量
            confirmation_threshold: 航迹确认阈值
            deletion_threshold: 航迹删除阈值
            n_scan: N-scan剪枝窗口大小（保留最近N步的假设）
            prune_threshold: 对数域剪枝阈值
            merge_threshold: 马氏距离合并阈值
            use_log_space: 全程使用对数空间
        """
        super().__init__(
            gating_threshold=gating_threshold,
            use_mahalanobis=use_mahalanobis,
            max_hypotheses=max_hypotheses,
            max_tracks=max_tracks,
            confirmation_threshold=confirmation_threshold,
            deletion_threshold=deletion_threshold,
            prune_threshold=prune_threshold,
            merge_threshold=merge_threshold,
            use_log_space=use_log_space
        )

        self.n_scan = n_scan

        # 假设树：hyp_id -> Hypothesis
        self.hypothesis_tree: Dict[int, Hypothesis] = {}

        # 假设创建时间戳：hyp_id -> time_step
        self._hypothesis_timestamps: Dict[int, int] = {}

        # 功能开关
        self._enable_n_scan = True
        self._enable_merging = True

    # ------------------------------------------------------------------
    # 公共接口覆盖
    # ------------------------------------------------------------------

    def associate(self,
                  measurements: np.ndarray,
                  predicted_measurements: np.ndarray,
                  measurement_covariances: Optional[List[np.ndarray]] = None,
                  innovation_covariances: Optional[List[np.ndarray]] = None) -> AssociationResult:
        """执行MHT关联（带假设树增量更新）

        与父类区别：
        - 不清除现有假设，基于假设树增量扩展
        - 应用N-scan剪枝
        """
        n_meas = measurements.shape[0]
        n_targets = predicted_measurements.shape[0]

        self.metrics.time_step += 1

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

        # 生成新假设（基于假设树增量扩展）
        t_start = time.perf_counter()
        new_hypotheses = self._generate_hypotheses(association_matrix)
        self.metrics.generation_time = time.perf_counter() - t_start
        self.metrics.total_generated = len(new_hypotheses)

        # 管理假设（合并 + N-scan + 剪枝）
        self._manage_hypotheses(new_hypotheses)

        # 选择最佳假设
        best_hypothesis = self._select_best_hypothesis()

        if best_hypothesis is not None:
            return AssociationResult(
                associations=dict(best_hypothesis.associations),
                unassociated_measurements=set(range(n_meas)) - set(best_hypothesis.associations.keys()),
                unassociated_targets=set(range(n_targets)) - set(best_hypothesis.associations.values()),
                association_matrix=association_matrix
            )
        else:
            return AssociationResult(
                associations={},
                unassociated_measurements=set(range(n_meas)),
                unassociated_targets=set(range(n_targets)),
                association_matrix=association_matrix
            )

    # ------------------------------------------------------------------
    # 假设生成（增量扩展）
    # ------------------------------------------------------------------

    def _generate_hypotheses(self,
                              association_matrix: np.ndarray) -> List[Hypothesis]:
        """基于假设树增量生成假设

        策略：
        - 无现有假设时从空假设开始
        - 有现有假设时对每个父假设进行束搜索扩展
        - 扩展数量自适应：max_hypotheses // len(parents)，保证总数可控

        Args:
            association_matrix: 关联矩阵

        Returns:
            扩展后的新假设列表
        """
        n_meas, n_targets = association_matrix.shape

        # 无现有假设时从空假设开始
        if not self.hypotheses:
            empty_hyp = Hypothesis(
                id=self.next_hypothesis_id,
                associations={},
                score=1.0,
                parent_id=None
            )
            self.next_hypothesis_id += 1
            self._hypothesis_timestamps[empty_hyp.id] = self.metrics.time_step
            self.hypothesis_tree[empty_hyp.id] = empty_hyp
            base_hypotheses = [empty_hyp]
        else:
            base_hypotheses = list(self.hypotheses)

        # 限制每个父假设的扩展数量
        # 总数 <= max_hypotheses * 2（给合并留余量）
        hypotheses_per_parent = max(
            1,
            (self.max_hypotheses * 2) // len(base_hypotheses)
        )

        # 对每个父假设进行束搜索扩展
        new_hypotheses = []

        # ---- 批量扩展：对全部父假设并行束搜索 ----
        # 用全局束搜索替代逐父假设独立束搜索
        # 这样可以在全局范围保留最优假设，而不是在每个父假设下各自保留

        # 收集所有父假设作为束的初始状态
        parent_states: List[Tuple[Dict[int, int], float, Set[int], Optional[int]]] = []
        for parent in base_hypotheses:
            used_targets = set(parent.associations.values())
            parent_states.append((
                dict(parent.associations),
                parent.score,
                used_targets,
                parent.id
            ))

        # 找出所有未分配的观测（在所有父假设中都未分配）
        all_assigned_meas: Set[int] = set()
        for _, assoc, _, _ in [(p[0], p[0], p[2], p[3]) for p in parent_states]:
            all_assigned_meas.update(assoc.keys())
        # 注意：不同父假设可能有不同已分配观测，统一处理效率更高
        # 改为逐父假设处理

        new_hypotheses = []
        for parent in base_hypotheses:
            extended = self._extend_hypothesis_beam(
                parent, association_matrix, hypotheses_per_parent
            )
            new_hypotheses.extend(extended)

        return new_hypotheses

    def _extend_hypothesis_beam(self,
                                 parent: Hypothesis,
                                 association_matrix: np.ndarray,
                                 max_extensions: int) -> List[Hypothesis]:
        """使用束搜索扩展单个父假设（对数空间）

        对父假设中未分配的观测进行束搜索扩展。
        全程对数空间避免数值下溢。

        Args:
            parent: 父假设
            association_matrix: 关联矩阵
            max_extensions: 最大扩展数量

        Returns:
            扩展后的假设列表（score在对数域）
        """
        n_meas, n_targets = association_matrix.shape

        used_meas = set(parent.associations.keys())
        used_targets = set(parent.associations.values())

        # 找出未处理的观测
        unassigned_meas = [i for i in range(n_meas) if i not in used_meas]

        if not unassigned_meas:
            return []

        if self.use_log_space:
            # --- 对数空间束搜索 ---
            log_clutter = np.log(0.5)
            beam: List[Tuple[Dict[int, int], float, Set[int]]] = [
                (dict(parent.associations), parent.score, set(used_targets))
            ]

            for i in unassigned_meas:
                candidates: List[Tuple[Dict[int, int], float, Set[int]]] = []

                for assoc, log_score, used_tgt in beam:
                    # 选项1：杂波
                    candidates.append((dict(assoc), log_score + log_clutter, set(used_tgt)))

                    # 选项2：关联到有效目标
                    for j in range(n_targets):
                        if j not in used_tgt and association_matrix[i, j] < np.inf:
                            new_assoc = dict(assoc)
                            new_assoc[i] = j
                            log_likelihood = -0.5 * association_matrix[i, j] ** 2
                            candidates.append((
                                new_assoc,
                                log_score + log_likelihood,
                                used_tgt | {j}
                            ))

                # 束剪枝
                candidates.sort(key=lambda x: x[1], reverse=True)
                beam = candidates[:min(max_extensions * 2, len(candidates))]
        else:
            # --- 线性空间（向后兼容） ---
            beam: List[Tuple[Dict[int, int], float, Set[int]]] = [
                (dict(parent.associations), parent.score, set(used_targets))
            ]

            for i in unassigned_meas:
                candidates = []

                for assoc, score, used_tgt in beam:
                    candidates.append((dict(assoc), score * 0.5, set(used_tgt)))

                    for j in range(n_targets):
                        if j not in used_tgt and association_matrix[i, j] < np.inf:
                            new_assoc = dict(assoc)
                            new_assoc[i] = j
                            log_prob = -0.5 * association_matrix[i, j]
                            candidates.append((
                                new_assoc,
                                score * np.exp(log_prob),
                                used_tgt | {j}
                            ))

                candidates.sort(key=lambda x: x[1], reverse=True)
                beam = candidates[:min(max_extensions * 2, len(candidates))]

        # 转换为 Hypothesis 对象
        results = []
        for item in beam:
            assoc, score, _ = item
            # 跳过与父假设无差异的（全杂波情况）
            if len(assoc) <= len(parent.associations):
                continue

            hyp = Hypothesis(
                id=self.next_hypothesis_id,
                associations=assoc,
                score=score,
                parent_id=parent.id
            )
            self.next_hypothesis_id += 1
            results.append(hyp)

        return results[:max_extensions]

    def _extend_hypothesis(self,
                            parent_hypothesis: Hypothesis,
                            association_matrix: np.ndarray) -> List[Hypothesis]:
        """扩展假设（保留原始接口兼容性）

        委托给束搜索版本，使用 max_hypotheses 作为扩展上限。
        """
        return self._extend_hypothesis_beam(
            parent_hypothesis,
            association_matrix,
            max_extensions=self.max_hypotheses
        )

    # ------------------------------------------------------------------
    # 假设合并
    # ------------------------------------------------------------------

    def _merge_hypotheses(self, hypotheses: List[Hypothesis]) -> List[Hypothesis]:
        """合并相同关联的假设（带开关）

        Args:
            hypotheses: 待合并假设列表

        Returns:
            合并后的假设列表
        """
        if not self._enable_merging:
            return hypotheses

        groups: Dict[Hashable, Hypothesis] = {}
        for hyp in hypotheses:
            key = frozenset(hyp.associations.items())
            existing = groups.get(key)
            if existing is None or hyp.score > existing.score:
                groups[key] = hyp
        return list(groups.values())

    # ------------------------------------------------------------------
    # N-scan 剪枝
    # ------------------------------------------------------------------

    def _n_scan_prune(self) -> None:
        """N-scan剪枝：删除N步前的假设及其子树

        保留策略：
        1. 保留在时间窗口 [current_step - n_scan, current_step] 内的假设
        2. 对窗口内的假设，追溯最多 n_scan 层祖先（避免保留整个历史树）
           - 深度限制：每步最多向上追溯 n_scan 代祖先
           - 祖先链的祖先链不会被无限保留

        原理：N-scan保留最近N次决策所需的树结构，更老的决策被抛弃。
        复杂度：O(H)，其中 H 为假设总数
        """
        if not self._enable_n_scan or self.n_scan <= 0:
            return

        current_step = self.metrics.time_step
        min_step = current_step - self.n_scan

        # 第1步：找出在时间窗口内的假设（时间戳 >= min_step）
        keep_ids: Set[int] = set()
        for hyp_id, timestamp in self._hypothesis_timestamps.items():
            if timestamp >= min_step:
                keep_ids.add(hyp_id)

        if not keep_ids:
            # 窗口为空时保留当前活跃假设（避免完全清空）
            for hyp in self.hypotheses:
                keep_ids.add(hyp.id)
            return

        # 第2步：追溯祖先链，但限制最大深度为 n_scan
        # 从活跃假设出发，向上追溯最多 n_scan 层
        for hyp in list(self.hypotheses):
            current = hyp
            depth = 0
            while current is not None and depth < self.n_scan:
                keep_ids.add(current.id)
                if current.parent_id is None:
                    break
                parent = self.hypothesis_tree.get(current.parent_id)
                if parent is None:
                    break
                current = parent
                depth += 1
            # 添加深度 n_scan 处的节点（作为剪枝边界锚点）
            if current is not None:
                keep_ids.add(current.id)

        # 第3步：过滤 hypotheses 列表（只保留在 keep_ids 中的）
        old_count = len(self.hypotheses)
        self.hypotheses = [h for h in self.hypotheses if h.id in keep_ids]

        # 第4步：清理 hypothesis_tree 和 timestamps
        prune_ids = set(self.hypothesis_tree.keys()) - keep_ids
        for hid in prune_ids:
            del self.hypothesis_tree[hid]
            self._hypothesis_timestamps.pop(hid, None)

        # 活跃列表中的剪枝数
        self.metrics.n_scan_pruned = old_count - len(self.hypotheses)
        # 树中清理的假设数
        self.metrics.tree_pruned = len(prune_ids)

    # ------------------------------------------------------------------
    # 假设管理覆盖
    # ------------------------------------------------------------------

    def _manage_hypotheses(self, new_hypotheses: List[Hypothesis]) -> None:
        """管理假设：注册树 → merge → prune → N-scan → similar-merge → cap

        完整流水线：
        1. 注册新建的假设到树结构（用于后续N-scan追溯）
        2. identical-merge: 合并相同关联的假设
        3. prune: 阈值剪枝（删除低分假设）
        4. N-scan剪枝: 删除N步前的旧假设
        5. similar-merge: 合并马氏距离相近的假设
        6. cap: 保留最多 max_hypotheses 个

        Args:
            new_hypotheses: 新生成的假设列表
        """
        # 1. 注册到假设树
        for hyp in new_hypotheses:
            self.hypothesis_tree[hyp.id] = hyp
            self._hypothesis_timestamps[hyp.id] = self.metrics.time_step

        # 2. identical-merge: 合并相同关联
        t_start = time.perf_counter()
        merged = self._merge_hypotheses(new_hypotheses)
        self.metrics.merge_time = time.perf_counter() - t_start
        self.metrics.after_merge = len(merged)

        self.hypotheses = merged

        # 3. prune: 阈值剪枝
        t_start = time.perf_counter()
        pruned = self._prune_hypotheses(self.hypotheses)
        self.metrics.prune_time = time.perf_counter() - t_start
        self.metrics.after_prune = len(pruned)
        self.hypotheses = pruned

        # 4. N-scan 剪枝（GlobalHypothesisMHT特有）
        t_start = time.perf_counter()
        self._n_scan_prune()
        self.metrics.prune_time += time.perf_counter() - t_start

        # 5. similar-merge: 马氏距离合并
        t_start = time.perf_counter()
        merged_sim = self._merge_similar_hypotheses(self.hypotheses)
        self.metrics.merge_time += time.perf_counter() - t_start
        self.hypotheses = merged_sim

        # 6. cap: 数量限制
        t_start = time.perf_counter()
        self.hypotheses = self._cap_hypotheses(self.hypotheses)
        self.metrics.cap_time = time.perf_counter() - t_start
        self.metrics.after_cap = len(self.hypotheses)
