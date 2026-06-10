"""
噪声模型
"""
import numpy as np
from typing import Optional, Tuple
from abc import ABC, abstractmethod


class NoiseModel(ABC):
    """噪声模型基类"""
    
    @abstractmethod
    def sample(self, size: int = 1) -> np.ndarray:
        """生成噪声样本
        
        Args:
            size: 样本数量
            
        Returns:
            噪声样本
        """
        pass
    
    @abstractmethod
    def get_covariance(self) -> np.ndarray:
        """获取噪声协方差矩阵"""
        pass


class GaussianNoise(NoiseModel):
    """高斯噪声模型"""
    
    def __init__(self, 
                 covariance: Optional[np.ndarray] = None,
                 std_dev: Optional[np.ndarray] = None):
        """
        初始化高斯噪声
        
        Args:
            covariance: 协方差矩阵
            std_dev: 标准差向量（会被转换为对角协方差矩阵）
        """
        if covariance is not None:
            self.covariance = np.asarray(covariance, dtype=np.float64)
            self.std_dev = np.sqrt(np.diag(self.covariance))
        elif std_dev is not None:
            self.std_dev = np.asarray(std_dev, dtype=np.float64)
            self.covariance = np.diag(self.std_dev ** 2)
        else:
            raise ValueError("Must provide either covariance or std_dev")
        
        self.dimension = len(self.std_dev)
    
    def sample(self, size: int = 1) -> np.ndarray:
        """生成高斯噪声样本
        
        Args:
            size: 样本数量
            
        Returns:
            噪声样本，形状为 (size, dimension)
        """
        if size == 1:
            return np.random.multivariate_normal(
                np.zeros(self.dimension), self.covariance
            )
        return np.random.multivariate_normal(
            np.zeros(self.dimension), self.covariance, size
        )
    
    def get_covariance(self) -> np.ndarray:
        return self.covariance.copy()
    
    def __repr__(self) -> str:
        return f"GaussianNoise(dimension={self.dimension}, std_dev={self.std_dev})"


class LinearMeasurementNoise(GaussianNoise):
    """线性观测噪声
    
    直接在笛卡尔坐标系中添加噪声
    """
    
    def __init__(self, 
                 x_std: float = 1.0, 
                 y_std: float = 1.0,
                 correlation: float = 0.0):
        """
        初始化线性观测噪声
        
        Args:
            x_std: x方向标准差
            y_std: y方向标准差
            correlation: x和y的相关系数
        """
        # 构建协方差矩阵
        cov = np.array([
            [x_std ** 2, correlation * x_std * y_std],
            [correlation * x_std * y_std, y_std ** 2]
        ])
        super().__init__(covariance=cov)
        self.x_std = x_std
        self.y_std = y_std
        self.correlation = correlation


class PolarMeasurementNoise(GaussianNoise):
    """极坐标观测噪声
    
    在极坐标系中添加噪声（距离和角度）
    """
    
    def __init__(self, 
                 range_std: float = 10.0,
                 bearing_std: float = 0.01,
                 correlation: float = 0.0):
        """
        初始化极坐标观测噪声
        
        Args:
            range_std: 距离标准差（米）
            bearing_std: 角度标准差（弧度）
            correlation: 距离和角度的相关系数
        """
        # 构建协方差矩阵
        cov = np.array([
            [range_std ** 2, correlation * range_std * bearing_std],
            [correlation * range_std * bearing_std, bearing_std ** 2]
        ])
        super().__init__(covariance=cov)
        self.range_std = range_std
        self.bearing_std = bearing_std
        self.correlation = correlation
    
    def add_noise_to_polar(self, 
                            range_true: float, 
                            bearing_true: float) -> Tuple[float, float]:
        """在极坐标中添加噪声
        
        Args:
            range_true: 真实距离
            bearing_true: 真实角度
            
        Returns:
            带噪声的距离和角度
        """
        noise = self.sample()
        range_meas = range_true + noise[0]
        bearing_meas = bearing_true + noise[1]
        
        # 确保距离非负
        range_meas = max(0.0, range_meas)
        
        # 归一化角度到[-pi, pi]
        bearing_meas = (bearing_meas + np.pi) % (2 * np.pi) - np.pi
        
        return range_meas, bearing_meas


class NonGaussianNoise(NoiseModel):
    """非高斯噪声模型（重尾噪声）"""
    
    def __init__(self, 
                 dimension: int = 2,
                 scale: float = 1.0,
                 degrees_of_freedom: float = 3.0):
        """
        初始化非高斯噪声（t分布噪声）
        
        Args:
            dimension: 噪声维度
            scale: 尺度参数
            degrees_of_freedom: 自由度（越小尾部越重）
        """
        self.dimension = dimension
        self.scale = scale
        self.dof = degrees_of_freedom
        self._covariance = (scale ** 2 * self.dof / (self.dof - 2)) * np.eye(dimension)
    
    def sample(self, size: int = 1) -> np.ndarray:
        """生成t分布噪声样本
        
        Args:
            size: 样本数量
            
        Returns:
            噪声样本
        """
        # 使用t分布生成重尾噪声
        samples = np.random.standard_t(self.dof, size=(size, self.dimension))
        return self.scale * samples
    
    def get_covariance(self) -> np.ndarray:
        return self._covariance.copy()


def create_noise_model(noise_type: str, **kwargs) -> NoiseModel:
    """创建噪声模型的工厂函数
    
    Args:
        noise_type: 噪声类型 ("linear", "polar", "nongaussian")
        **kwargs: 模型参数
        
    Returns:
        噪声模型实例
    """
    noise_models = {
        "linear": LinearMeasurementNoise,
        "polar": PolarMeasurementNoise,
        "nongaussian": NonGaussianNoise
    }
    
    if noise_type not in noise_models:
        raise ValueError(f"Unknown noise type: {noise_type}")
    
    return noise_models[noise_type](**kwargs)
