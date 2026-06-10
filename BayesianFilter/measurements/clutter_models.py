"""
杂波模型（严格分离概念）

杂波 = 虚假观测点，与目标状态无关
杂波模型 = 杂波数量模型 + 杂波位置分布

1. 杂波数量模型：
   - PoissonCount: N_c ~ Poisson(λ)
   - FixedCount: N_c = 固定值
   - BinomialCount: N_c ~ Binomial(n, p)

2. 杂波位置分布：
   - UniformDistribution: 均匀分布
   - GaussianDistribution: 高斯分布（围绕目标）
   - NonUniformDistribution: 非均匀强度分布

注意：
- "泊松"描述杂波数量，不描述位置分布
- "高斯杂波"描述虚假观测点的位置分布，不等价于高斯测量噪声
"""
import numpy as np
from typing import Optional, List, Tuple, Dict
from abc import ABC, abstractmethod


# ============================================================================
# 杂波数量模型
# ============================================================================

class ClutterCountModel(ABC):
    """杂波数量模型基类"""
    
    @abstractmethod
    def sample_count(self) -> int:
        """采样杂波数量"""
        pass
    
    @abstractmethod
    def get_expected_count(self) -> float:
        """获取期望杂波数量"""
        pass


class PoissonCount(ClutterCountModel):
    """泊松杂波数量模型
    
    N_c ~ Poisson(λ)
    最常用的杂波数量模型
    """
    
    def __init__(self, clutter_rate: float = 5.0, random_seed: Optional[int] = None):
        """
        Args:
            clutter_rate: λ，平均杂波数量
            random_seed: 随机种子
        """
        self.clutter_rate = clutter_rate
        if random_seed is not None:
            np.random.seed(random_seed)
    
    def sample_count(self) -> int:
        return np.random.poisson(self.clutter_rate)
    
    def get_expected_count(self) -> float:
        return self.clutter_rate


class FixedCount(ClutterCountModel):
    """固定杂波数量模型
    
    N_c = 固定值
    适用于：已知杂波数量的场景
    """
    
    def __init__(self, n_clutter: int = 5):
        self.n_clutter = n_clutter
    
    def sample_count(self) -> int:
        return self.n_clutter
    
    def get_expected_count(self) -> float:
        return float(self.n_clutter)


class BinomialCount(ClutterCountModel):
    """二项分布杂波数量模型
    
    N_c ~ Binomial(n, p)
    适用于：有限次试验的杂波
    """
    
    def __init__(self, n_trials: int = 10, p: float = 0.5, random_seed: Optional[int] = None):
        """
        Args:
            n_trials: 试验次数
            p: 每次试验产生杂波的概率
            random_seed: 随机种子
        """
        self.n_trials = n_trials
        self.p = p
        if random_seed is not None:
            np.random.seed(random_seed)
    
    def sample_count(self) -> int:
        return np.random.binomial(self.n_trials, self.p)
    
    def get_expected_count(self) -> float:
        return self.n_trials * self.p


# ============================================================================
# 杂波位置分布
# ============================================================================

class ClutterDistribution(ABC):
    """杂波位置分布基类"""
    
    @abstractmethod
    def sample_positions(self, 
                         n_samples: int,
                         target_positions: Dict[int, np.ndarray]) -> List[np.ndarray]:
        """采样杂波位置
        
        Args:
            n_samples: 采样数量
            target_positions: 目标位置字典（用于确定采样区域）
            
        Returns:
            杂波位置列表，每个为 [x, y]
        """
        pass


class UniformDistribution(ClutterDistribution):
    """均匀位置分布
    
    杂波在目标包围盒内均匀分布
    p_c(z) = 1/Area
    """
    
    def __init__(self, padding: float = 200.0):
        """
        Args:
            padding: 包围盒边距（米）
        """
        self.padding = padding
    
    def sample_positions(self, 
                         n_samples: int,
                         target_positions: Dict[int, np.ndarray]) -> List[np.ndarray]:
        if not target_positions or n_samples == 0:
            return []
        
        # 计算目标包围盒
        positions = np.array(list(target_positions.values()))
        x_min, x_max = positions[:, 0].min() - self.padding, positions[:, 0].max() + self.padding
        y_min, y_max = positions[:, 1].min() - self.padding, positions[:, 1].max() + self.padding
        
        # 在包围盒内均匀采样
        clutter_list = []
        for _ in range(n_samples):
            x = np.random.uniform(x_min, x_max)
            y = np.random.uniform(y_min, y_max)
            clutter_list.append(np.array([x, y]))
        
        return clutter_list


class GaussianDistribution(ClutterDistribution):
    """高斯位置分布
    
    杂波围绕每个目标呈高斯分布
    p_c(z) = Σ N(z; x_i, Σ_c)
    
    注意：这是杂波位置分布，不是测量噪声！
    """
    
    def __init__(self, clutter_std: float = 50.0):
        """
        Args:
            clutter_std: 杂波分布标准差（米）
        """
        self.clutter_std = clutter_std
    
    def sample_positions(self, 
                         n_samples: int,
                         target_positions: Dict[int, np.ndarray]) -> List[np.ndarray]:
        if not target_positions or n_samples == 0:
            return []
        
        clutter_list = []
        n_targets = len(target_positions)
        samples_per_target = n_samples // n_targets
        
        for target_id, pos in target_positions.items():
            # 每个目标周围的杂波数量
            n_target_clutter = samples_per_target if target_id < n_targets - 1 else n_samples - len(clutter_list)
            
            # 在目标周围高斯采样
            for _ in range(n_target_clutter):
                noise = np.random.randn(2) * self.clutter_std
                clutter_list.append(pos + noise)
        
        return clutter_list


class NonUniformDistribution(ClutterDistribution):
    """非均匀位置分布
    
    杂波密度随位置变化
    p_c(z) = λ(z) / ∫λ(z)dz
    
    适用于：地形影响、城市环境等
    """
    
    def __init__(self, 
                 hotspots: Optional[List[Tuple[np.ndarray, float, float]]] = None,
                 padding: float = 200.0):
        """
        Args:
            hotspots: 热点列表 [(位置, 半径, 相对密度), ...]
            padding: 包围盒边距
        """
        self.hotspots = hotspots or []
        self.padding = padding
    
    def sample_positions(self, 
                         n_samples: int,
                         target_positions: Dict[int, np.ndarray]) -> List[np.ndarray]:
        if not target_positions or n_samples == 0:
            return []
        
        # 计算包围盒
        positions = np.array(list(target_positions.values()))
        x_min, x_max = positions[:, 0].min() - self.padding, positions[:, 0].max() + self.padding
        y_min, y_max = positions[:, 1].min() - self.padding, positions[:, 1].max() + self.padding
        
        clutter_list = []
        
        # 基础杂波（均匀分布）
        n_base = n_samples // 2
        for _ in range(n_base):
            x = np.random.uniform(x_min, x_max)
            y = np.random.uniform(y_min, y_max)
            clutter_list.append(np.array([x, y]))
        
        # 热点杂波
        n_hotspot = n_samples - n_base
        if self.hotspots and n_hotspot > 0:
            # 按密度加权分配热点杂波
            total_density = sum(density for _, _, density in self.hotspots)
            for center, radius, density in self.hotspots:
                n_center = int(n_hotspot * density / total_density)
                for _ in range(n_center):
                    angle = np.random.uniform(0, 2 * np.pi)
                    r = np.random.uniform(0, radius)
                    x = center[0] + r * np.cos(angle)
                    y = center[1] + r * np.sin(angle)
                    clutter_list.append(np.array([x, y]))
        
        return clutter_list


# ============================================================================
# 组合杂波模型
# ============================================================================

class ClutterModel:
    """组合杂波模型
    
    杂波 = 杂波数量模型 + 杂波位置分布
    
    使用示例：
        # 泊松数量 + 均匀分布
        clutter = ClutterModel(PoissonCount(5.0), UniformDistribution(200.0))
        
        # 泊松数量 + 高斯分布（围绕目标）
        clutter = ClutterModel(PoissonCount(3.0), GaussianDistribution(50.0))
        
        # 固定数量 + 非均匀分布
        clutter = ClutterModel(FixedCount(10), NonUniformDistribution(hotspots))
    """
    
    def __init__(self, 
                 count_model: ClutterCountModel,
                 distribution: ClutterDistribution):
        """
        Args:
            count_model: 杂波数量模型
            distribution: 杂波位置分布
        """
        self.count_model = count_model
        self.distribution = distribution
    
    def generate_clutter(self, 
                         timestamp: float,
                         target_positions: Dict[int, np.ndarray]) -> List[np.ndarray]:
        """生成杂波
        
        Args:
            timestamp: 时间戳
            target_positions: 目标位置字典 {target_id: [x, y]}
            
        Returns:
            杂波观测列表，每个为 [x, y]
        """
        # 1. 采样杂波数量
        n_clutter = self.count_model.sample_count()
        
        # 2. 采样杂波位置
        clutter_positions = self.distribution.sample_positions(n_clutter, target_positions)
        
        return clutter_positions
    
    def get_expected_clutter_count(self) -> float:
        """获取期望杂波数量"""
        return self.count_model.get_expected_count()


# ============================================================================
# 便捷工厂函数
# ============================================================================

def create_uniform_clutter(clutter_rate: float = 5.0, 
                           padding: float = 200.0) -> ClutterModel:
    """创建均匀杂波模型（泊松数量 + 均匀分布）"""
    return ClutterModel(PoissonCount(clutter_rate), UniformDistribution(padding))


def create_gaussian_clutter(clutter_rate: float = 3.0,
                            clutter_std: float = 50.0) -> ClutterModel:
    """创建高斯杂波模型（泊松数量 + 高斯分布）"""
    return ClutterModel(PoissonCount(clutter_rate), GaussianDistribution(clutter_std))


def create_nonuniform_clutter(base_rate: float = 2.0,
                               hotspots: Optional[List[Tuple[np.ndarray, float, float]]] = None) -> ClutterModel:
    """创建非均匀杂波模型（泊松数量 + 非均匀分布）"""
    return ClutterModel(PoissonCount(base_rate), NonUniformDistribution(hotspots))
