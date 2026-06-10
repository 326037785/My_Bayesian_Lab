"""
观测模拟器
"""
import numpy as np
from typing import Optional, List, Dict, Tuple
from .measurement import Measurement, MeasurementType, MeasurementSet
from .noise_models import NoiseModel, LinearMeasurementNoise, PolarMeasurementNoise
from .clutter_models import ClutterModel, create_uniform_clutter


class MeasurementSimulator:
    """观测模拟器
    
    整合目标观测生成和杂波生成
    """
    
    def __init__(self,
                 measurement_type: MeasurementType = MeasurementType.LINEAR,
                 noise_model: Optional[NoiseModel] = None,
                 clutter_model: Optional[ClutterModel] = None,
                 detection_probability: float = 0.9,
                 random_seed: Optional[int] = None):
        """
        初始化观测模拟器
        
        Args:
            measurement_type: 观测类型
            noise_model: 噪声模型
            clutter_model: 杂波模型
            detection_probability: 目标检测概率
            random_seed: 随机种子
        """
        self.measurement_type = measurement_type
        self.detection_probability = detection_probability
        
        # 设置默认噪声模型
        if noise_model is None:
            if measurement_type == MeasurementType.LINEAR:
                self.noise_model = LinearMeasurementNoise(x_std=1.0, y_std=1.0)
            elif measurement_type == MeasurementType.POLAR:
                self.noise_model = PolarMeasurementNoise(range_std=10.0, bearing_std=0.01)
            else:
                self.noise_model = LinearMeasurementNoise(x_std=1.0, y_std=1.0)
        else:
            self.noise_model = noise_model
        
        # 设置默认杂波模型
        if clutter_model is None:
            self.clutter_model = create_uniform_clutter(clutter_rate=5.0, padding=200.0)
        else:
            self.clutter_model = clutter_model
        
        # 默认监视区域
        self.surveillance_region = ((-1000, 1000), (-1000, 1000))
        
        if random_seed is not None:
            np.random.seed(random_seed)
    
    def set_surveillance_region(self, 
                                 x_range: Tuple[float, float],
                                 y_range: Tuple[float, float]):
        """设置监视区域
        
        Args:
            x_range: x范围 (min, max)
            y_range: y范围 (min, y_max)
        """
        self.surveillance_region = (x_range, y_range)
    
    def generate_measurements(self,
                               timestamp: float,
                               target_states: Dict[int, np.ndarray]) -> MeasurementSet:
        """生成一个时刻的观测
        
        Args:
            timestamp: 时间戳
            target_states: 目标状态字典，键为目标ID，值为状态向量
            
        Returns:
            观测集合
        """
        measurement_set = MeasurementSet(timestamp)
        
        # 生成目标观测
        for target_id, state in target_states.items():
            # 检测概率
            if np.random.random() > self.detection_probability:
                continue
            
            # 根据观测类型生成观测
            if self.measurement_type == MeasurementType.LINEAR:
                meas = self._generate_linear_measurement(state)
            elif self.measurement_type == MeasurementType.POLAR:
                meas = self._generate_polar_measurement(state)
            else:
                meas = self._generate_linear_measurement(state)
            
            measurement = Measurement(
                timestamp=timestamp,
                measurement=meas,
                measurement_type=self.measurement_type,
                target_id=target_id,
                is_clutter=False,
                covariance=self.noise_model.get_covariance()
            )
            measurement_set.add_measurement(measurement)
        
        # 生成杂波（动态包围盒：基于当前时刻的目标位置）
        target_positions = {}
        for target_id, state in target_states.items():
            if len(state) >= 4:
                target_positions[target_id] = np.array([state[0], state[2]])
            else:
                target_positions[target_id] = state[:2]
        
        clutter_measurements = self.clutter_model.generate_clutter(
            timestamp, target_positions
        )
        
        for clutter in clutter_measurements:
            # 杂波坐标与测量类型保持一致
            if self.measurement_type == MeasurementType.POLAR:
                # 将笛卡尔杂波 [x, y] 转换为极坐标 [range, bearing]
                r = np.sqrt(clutter[0]**2 + clutter[1]**2)
                theta = np.arctan2(clutter[1], clutter[0])
                clutter_meas = np.array([r, theta])
                clutter_type = MeasurementType.POLAR
                clutter_cov = self.noise_model.get_covariance()
            else:
                clutter_meas = clutter
                clutter_type = MeasurementType.LINEAR
                clutter_cov = None

            measurement = Measurement(
                timestamp=timestamp,
                measurement=clutter_meas,
                measurement_type=clutter_type,
                target_id=None,
                is_clutter=True,
                covariance=clutter_cov
            )
            measurement_set.add_measurement(measurement)
        
        return measurement_set
    
    def _generate_linear_measurement(self, state: np.ndarray) -> np.ndarray:
        """生成线性观测
        
        Args:
            state: 目标状态向量 [x, vx, y, vy] 或 [x, vx, ax, y, vy, ay]
            
        Returns:
            观测向量 [x, y]
        """
        # 提取真实位置
        if len(state) >= 4:
            x_true, y_true = state[0], state[2]
        else:
            x_true, y_true = state[0], state[1]
        
        true_position = np.array([x_true, y_true])
        
        # 添加噪声
        noise = self.noise_model.sample()
        
        return true_position + noise
    
    def _generate_polar_measurement(self, state: np.ndarray) -> np.ndarray:
        """生成极坐标观测
        
        Args:
            state: 目标状态向量
            
        Returns:
            观测向量 [range, bearing]
        """
        # 提取真实位置
        if len(state) >= 4:
            x_true, y_true = state[0], state[2]
        else:
            x_true, y_true = state[0], state[1]
        
        # 转换为极坐标
        range_true = np.sqrt(x_true ** 2 + y_true ** 2)
        bearing_true = np.arctan2(y_true, x_true)
        
        # 添加极坐标噪声
        if isinstance(self.noise_model, PolarMeasurementNoise):
            range_meas, bearing_meas = self.noise_model.add_noise_to_polar(
                range_true, bearing_true
            )
        else:
            noise = self.noise_model.sample()
            range_meas = range_true + noise[0]
            bearing_meas = bearing_true + noise[1]
        
        return np.array([range_meas, bearing_meas])
    
    def generate_scenario_measurements(self,
                                        scenario_data: Dict[int, List[Tuple[float, np.ndarray]]],
                                        time_step: float = 1.0) -> Dict[float, MeasurementSet]:
        """生成整个场景的观测
        
        Args:
            scenario_data: 场景数据，格式为 {target_id: [(timestamp, state), ...]}
            time_step: 时间步长
            
        Returns:
            字典，键为时间戳，值为观测集合
        """
        # 获取所有时间戳
        all_timestamps = set()
        for trajectory in scenario_data.values():
            for timestamp, _ in trajectory:
                all_timestamps.add(timestamp)
        
        all_timestamps = sorted(all_timestamps)
        
        # 生成每个时刻的观测
        measurements_by_time = {}
        
        for t in all_timestamps:
            # 获取当前时刻存活的目标状态
            target_states = {}
            for target_id, trajectory in scenario_data.items():
                for timestamp, state in trajectory:
                    if abs(timestamp - t) < 1e-6:
                        target_states[target_id] = state
                        break
            
            # 生成观测
            measurement_set = self.generate_measurements(t, target_states)
            measurements_by_time[t] = measurement_set
        
        return measurements_by_time
    
    def generate_measurements_for_tracking(self,
                                            scenario_data: Dict[int, List[Tuple[float, np.ndarray]]]) -> List[MeasurementSet]:
        """生成用于跟踪的观测序列
        
        Args:
            scenario_data: 场景数据
            
        Returns:
            观测集合列表（按时间排序）
        """
        measurements_by_time = self.generate_scenario_measurements(scenario_data)
        return [measurements_by_time[t] for t in sorted(measurements_by_time.keys())]
