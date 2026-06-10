"""
场景管理器
"""
import numpy as np
from typing import List, Dict, Optional, Tuple
from .target import Target, MotionModel
from .motion_models import create_motion_model, MotionModelBase


class ScenarioManager:
    """场景管理器
    
    管理多个目标的轨迹生成
    """
    
    def __init__(self, 
                 time_step: float = 1.0,
                 process_noise_std: float = 0.1,
                 random_seed: Optional[int] = None):
        """
        初始化场景管理器
        
        Args:
            time_step: 时间步长（秒）
            process_noise_std: 过程噪声标准差
            random_seed: 随机种子
        """
        self.time_step = time_step
        self.process_noise_std = process_noise_std
        self.targets: List[Target] = []
        self.current_time = 0.0
        
        if random_seed is not None:
            np.random.seed(random_seed)
    
    def add_target(self, 
                   target_id: int,
                   initial_state: np.ndarray,
                   birth_time: float,
                   death_time: float,
                   motion_model: str = "CV",
                   **model_kwargs) -> Target:
        """添加目标
        
        Args:
            target_id: 目标ID
            initial_state: 初始状态向量
            birth_time: 目标出现时间
            death_time: 目标消失时间
            motion_model: 运动模型类型 ("CV", "CA", "CT", "RW")
            **model_kwargs: 模型参数
            
        Returns:
            创建的目标对象
        """
        model_enum = MotionModel(motion_model)
        target = Target(
            target_id=target_id,
            initial_state=initial_state,
            birth_time=birth_time,
            death_time=death_time,
            motion_model=model_enum
        )
        # 存储模型参数，供 generate_scenario 创建运动模型时使用
        target._model_kwargs = model_kwargs
        self.targets.append(target)
        return target
    
    def generate_scenario(self, duration: float) -> Dict[int, List[Tuple[float, np.ndarray]]]:
        """生成整个场景
        
        Args:
            duration: 场景总时长（秒）
            
        Returns:
            字典，键为目标ID，值为(时间戳, 状态)列表
        """
        time_steps = np.arange(0, duration + self.time_step, self.time_step)
        scenario_data = {}
        
        for target in self.targets:
            # 检查是否为混合机动目标（使用 _maneuver_schedule 自定义轨迹生成）
            if hasattr(target, '_maneuver_schedule') and target._maneuver_schedule:
                trajectory = self._generate_mixed_trajectory(target, time_steps)
            else:
                trajectory = []
                model_kwargs = getattr(target, '_model_kwargs', {})
                motion_model = create_motion_model(target.motion_model.value, **model_kwargs)
                
                for t in time_steps:
                    if not target.is_alive(t):
                        continue
                    
                    if len(trajectory) == 0:
                        # 第一个时间步，使用初始状态
                        state = target.initial_state.copy()
                    else:
                        # 使用运动模型进行状态转移
                        prev_state = trajectory[-1][1]
                        process_noise = self._generate_process_noise(motion_model)
                        state = motion_model.state_transition(
                            prev_state, self.time_step, process_noise
                        )
                    
                    trajectory.append((t, state))
                    target.add_state(t, state)
            
            scenario_data[target.target_id] = trajectory
        
        return scenario_data
    
    def _generate_process_noise(self, motion_model: MotionModelBase) -> np.ndarray:
        """生成过程噪声
        
        Args:
            motion_model: 运动模型
            
        Returns:
            过程噪声向量
        """
        Q = motion_model.get_process_noise_matrix(
            self.time_step, self.process_noise_std
        )
        return np.random.multivariate_normal(
            np.zeros(Q.shape[0]), Q
        )
    
    def get_all_positions_at_time(self, time: float) -> Dict[int, np.ndarray]:
        """获取指定时间所有存活目标的位置
        
        Args:
            time: 时间戳
            
        Returns:
            字典，键为目标ID，值为位置向量
        """
        positions = {}
        for target in self.targets:
            if target.is_alive(time):
                pos = target.get_position_at_time(time)
                if pos is not None:
                    positions[target.target_id] = pos
        return positions
    
    def get_all_states_at_time(self, time: float) -> Dict[int, np.ndarray]:
        """获取指定时间所有存活目标的状态
        
        Args:
            time: 时间戳
            
        Returns:
            字典，键为目标ID，值为状态向量
        """
        states = {}
        for target in self.targets:
            if target.is_alive(time):
                state = target.get_state_at_time(time)
                if state is not None:
                    states[target.target_id] = state
        return states
    
    def get_target_by_id(self, target_id: int) -> Optional[Target]:
        """根据目标ID获取目标对象
        
        Args:
            target_id: 目标ID
            
        Returns:
            目标对象，如果不存在返回 None
        """
        for target in self.targets:
            if target.target_id == target_id:
                return target
        return None
    
    # ------------------------------------------------------------------
    # 场景创建方法
    # ------------------------------------------------------------------
    
    def create_linear_scenario(self, 
                               n_targets: int = 1,
                               duration: float = 100.0,
                               x_range: Tuple[float, float] = (-500, 500),
                               y_range: Tuple[float, float] = (-500, 500),
                               speed_range: Tuple[float, float] = (5, 20)) -> None:
        """创建线性运动场景（所有目标使用 CV 匀速直线模型）
        
        Args:
            n_targets: 目标数量
            duration: 场景时长
            x_range: x坐标范围
            y_range: y坐标范围
            speed_range: 速度范围
        """
        for i in range(n_targets):
            # 随机初始位置
            x0 = np.random.uniform(*x_range)
            y0 = np.random.uniform(*y_range)
            
            # 随机速度和方向
            speed = np.random.uniform(*speed_range)
            angle = np.random.uniform(0, 2 * np.pi)
            vx0 = speed * np.cos(angle)
            vy0 = speed * np.sin(angle)
            
            initial_state = np.array([x0, vx0, y0, vy0])
            
            # 随机出生和死亡时间
            birth_time = 0.0
            death_time = duration
            
            self.add_target(
                target_id=i,
                initial_state=initial_state,
                birth_time=birth_time,
                death_time=death_time,
                motion_model="CV"
            )
    
    def create_maneuvering_scenario(self,
                                     n_targets: int = 1,
                                     duration: float = 100.0,
                                     n_maneuvers: int = 3) -> None:
        """创建机动目标场景
        
        为每个目标随机选择运动模型（CV/CA/CT），不同目标可使用不同模型。
        相比于线性场景，目标运动更加多样化。
        
        Args:
            n_targets: 目标数量
            duration: 场景时长
            n_maneuvers: 机动次数（保留参数，当前为每个目标使用一种随机模型）
        """
        model_choices = ["CV", "CA", "CT"]
        for i in range(n_targets):
            # 随机选择运动模型
            model_type = np.random.choice(model_choices)
            
            x0 = np.random.uniform(-500, 500)
            y0 = np.random.uniform(-500, 500)
            vx0 = np.random.uniform(-10, 10)
            vy0 = np.random.uniform(-10, 10)
            
            # 根据模型类型创建相应维度的初始状态
            if model_type == "CV":
                initial_state = np.array([x0, vx0, y0, vy0])
                model_kwargs = {}
            elif model_type == "CA":
                initial_state = np.array([x0, vx0, 0.0, y0, vy0, 0.0])
                model_kwargs = {}
            elif model_type == "CT":
                turn_rate = np.random.uniform(-0.2, 0.2)
                initial_state = np.array([x0, vx0, y0, vy0])
                model_kwargs = {"known_turn_rate": turn_rate}
            else:
                initial_state = np.array([x0, vx0, y0, vy0])
                model_kwargs = {}
            
            self.add_target(
                target_id=i,
                initial_state=initial_state,
                birth_time=0.0,
                death_time=duration,
                motion_model=model_type,
                **model_kwargs
            )
    
    def create_coordinated_turn_scenario(self,
                                          n_targets: int = 1,
                                          duration: float = 100.0,
                                          x_range: Tuple[float, float] = (-500, 500),
                                          y_range: Tuple[float, float] = (-500, 500),
                                          speed_range: Tuple[float, float] = (5, 20),
                                          turn_rate: float = 0.1) -> None:
        """创建协调转弯场景
        
        所有目标使用 CT 模型进行匀速转弯运动。
        每个目标以随机初始航向开始，并随机选择转弯方向（左转/右转）。
        
        Args:
            n_targets: 目标数量
            duration: 场景时长
            x_range: x坐标范围
            y_range: y坐标范围
            speed_range: 速度范围
            turn_rate: 转弯速率 (rad/s)，正值左转，负值右转
        """
        for i in range(n_targets):
            x0 = np.random.uniform(*x_range)
            y0 = np.random.uniform(*y_range)
            
            speed = np.random.uniform(*speed_range)
            angle = np.random.uniform(0, 2 * np.pi)
            vx0 = speed * np.cos(angle)
            vy0 = speed * np.sin(angle)
            
            # CT 模型使用 4D 状态 [x, vx, y, vy]（通过 known_turn_rate 指定转弯速率）
            initial_state = np.array([x0, vx0, y0, vy0])
            
            # 随机选择转弯方向，使部分目标左转、部分右转
            target_turn_rate = turn_rate * np.random.choice([-1, 1])
            
            self.add_target(
                target_id=i,
                initial_state=initial_state,
                birth_time=0.0,
                death_time=duration,
                motion_model="CT",
                known_turn_rate=target_turn_rate
            )
    
    def create_random_walk_scenario(self,
                                     n_targets: int = 1,
                                     duration: float = 100.0,
                                     x_range: Tuple[float, float] = (-500, 500),
                                     y_range: Tuple[float, float] = (-500, 500),
                                     walk_std: float = 5.0) -> None:
        """创建随机游走场景
        
        所有目标使用 RW 模型进行随机游走运动。
        目标状态仅包含位置 [x, y]，每一步添加高斯随机扰动。
        
        Args:
            n_targets: 目标数量
            duration: 场景时长
            x_range: x坐标范围
            y_range: y坐标范围
            walk_std: 游走标准差（控制随机运动的幅度，值越大运动越剧烈）
        """
        for i in range(n_targets):
            x0 = np.random.uniform(*x_range)
            y0 = np.random.uniform(*y_range)
            
            # RW 模型使用 2D 状态 [x, y]
            initial_state = np.array([x0, y0])
            
            self.add_target(
                target_id=i,
                initial_state=initial_state,
                birth_time=0.0,
                death_time=duration,
                motion_model="RW"
            )
    
    def create_mixed_maneuver_scenario(self,
                                        n_targets: int = 1,
                                        duration: float = 100.0,
                                        n_maneuvers: int = 3,
                                        x_range: Tuple[float, float] = (-500, 500),
                                        y_range: Tuple[float, float] = (-500, 500),
                                        speed_range: Tuple[float, float] = (5, 20)) -> None:
        """创建混合机动场景
        
        每个目标在运动过程中切换不同运动模型（CV / CT / RW），
        实现直线-转弯-随机游走等复杂机动模式。
        目标按照机动段顺序自动切换模型并处理状态维度转换。
        
        Args:
            n_targets: 目标数量
            duration: 场景时长
            n_maneuvers: 机动段数量（每个目标包含的机动段数）
            x_range: x坐标范围
            y_range: y坐标范围
            speed_range: 速度范围
        """
        for i in range(n_targets):
            x0 = np.random.uniform(*x_range)
            y0 = np.random.uniform(*y_range)
            
            speed = np.random.uniform(*speed_range)
            angle = np.random.uniform(0, 2 * np.pi)
            vx0 = speed * np.cos(angle)
            vy0 = speed * np.sin(angle)
            
            # 生成机动段序列：第一段固定为 CV，后续随机选择 CT 或 RW
            segment_duration = duration / n_maneuvers
            segments = []
            
            # 第一段：匀速直线
            segments.append(("CV", {}, segment_duration))
            
            # 后续段：随机选择 CT 或 RW，确保不连续重复
            model_pool = ["CT", "RW", "CV"]
            for j in range(1, n_maneuvers):
                prev_model = segments[-1][0]
                available = [m for m in model_pool if m != prev_model]
                model_type = np.random.choice(available)
                
                if model_type == "CT":
                    turn_rate = np.random.uniform(-0.2, 0.2)
                    model_kwargs = {"known_turn_rate": turn_rate}
                else:
                    model_kwargs = {}
                
                segments.append((model_type, model_kwargs, segment_duration))
            
            # 初始状态使用 CV 模型的 4D 状态 [x, vx, y, vy]
            initial_state = np.array([x0, vx0, y0, vy0])
            
            # 创建目标并附加机动段信息
            target = self.add_target(
                target_id=i,
                initial_state=initial_state,
                birth_time=0.0,
                death_time=duration,
                motion_model="CV"  # 第一段为 CV
            )
            # 在目标对象上附加机动段调度（由 _generate_mixed_trajectory 使用）
            target._maneuver_schedule = segments
    
    def create_multi_model_scenario(self,
                                     n_targets: int = 3,
                                     duration: float = 100.0,
                                     x_range: Tuple[float, float] = (-500, 500),
                                     y_range: Tuple[float, float] = (-500, 500),
                                     speed_range: Tuple[float, float] = (5, 20)) -> None:
        """创建多模型场景
        
        多个目标各自使用不同的运动模型（CV/CA/CT/RW），
        用于测试滤波器对不同运动模型的适应能力。
        模型按顺序循环分配，确保每个模型至少有一个目标使用。
        
        Args:
            n_targets: 目标数量
            duration: 场景时长
            x_range: x坐标范围
            y_range: y坐标范围
            speed_range: 速度范围
        """
        model_pool = ["CV", "CA", "CT", "RW"]
        
        for i in range(n_targets):
            # 循环分配不同模型
            model_type = model_pool[i % len(model_pool)]
            
            x0 = np.random.uniform(*x_range)
            y0 = np.random.uniform(*y_range)
            speed = np.random.uniform(*speed_range)
            angle = np.random.uniform(0, 2 * np.pi)
            vx0 = speed * np.cos(angle)
            vy0 = speed * np.sin(angle)
            
            # 根据模型类型创建相应维度的初始状态
            if model_type == "CV":
                initial_state = np.array([x0, vx0, y0, vy0])
                model_kwargs = {}
            elif model_type == "CA":
                initial_state = np.array([x0, vx0, 0.0, y0, vy0, 0.0])
                model_kwargs = {}
            elif model_type == "CT":
                turn_rate = np.random.uniform(-0.2, 0.2)
                initial_state = np.array([x0, vx0, y0, vy0])
                model_kwargs = {"known_turn_rate": turn_rate}
            elif model_type == "RW":
                initial_state = np.array([x0, y0])
                model_kwargs = {}
            
            self.add_target(
                target_id=i,
                initial_state=initial_state,
                birth_time=0.0,
                death_time=duration,
                motion_model=model_type,
                **model_kwargs
            )
    
    # ------------------------------------------------------------------
    # 混合机动场景辅助方法
    # ------------------------------------------------------------------
    
    def _generate_mixed_trajectory(self, target: Target,
                                    time_steps: np.ndarray) -> List[Tuple[float, np.ndarray]]:
        """生成混合机动目标的轨迹
        
        根据目标附加的 _maneuver_schedule 属性，在不同时间段自动切换运动模型，
        并处理状态维度转换（如 CV 4D ↔ RW 2D）。
        
        Args:
            target: 目标对象（需包含 _maneuver_schedule 属性）
            time_steps: 时间步数组
            
        Returns:
            (时间戳, 状态) 列表
        """
        trajectory = []
        current_state = None
        current_model_type = None
        current_model = None
        segments = target._maneuver_schedule
        
        for t in time_steps:
            if not target.is_alive(t):
                continue
            
            # 确定当前时间所属的机动段
            elapsed = t - target.birth_time
            model_type, model_kwargs = self._get_model_at_time(segments, elapsed)
            
            if len(trajectory) == 0:
                # 第一个时间步，使用初始状态
                current_state = target.initial_state.copy()
                current_model_type = model_type
                current_model = create_motion_model(model_type, **model_kwargs)
            else:
                # 检查运动模型是否切换
                if model_type != current_model_type:
                    # 模型切换：转换状态维度
                    current_state = self._convert_state_for_model(
                        current_state, current_model_type, model_type
                    )
                    current_model_type = model_type
                    current_model = create_motion_model(model_type, **model_kwargs)
                
                # 使用当前模型进行状态转移
                process_noise = self._generate_process_noise(current_model)
                current_state = current_model.state_transition(
                    current_state, self.time_step, process_noise
                )
            
            trajectory.append((t, current_state.copy()))
            target.add_state(t, current_state.copy())
        
        return trajectory
    
    @staticmethod
    def _get_model_at_time(segments: List[Tuple],
                            elapsed: float) -> Tuple[str, Dict]:
        """获取指定时间对应的运动模型
        
        遍历机动段列表，找到 elapsed 时间所在的段。
        
        Args:
            segments: 机动段列表，每项为 (model_type, kwargs, duration)
            elapsed: 从目标出生开始的经过时间
            
        Returns:
            (模型类型, 参数字典)
        """
        cumulative = 0.0
        for model_type, model_kwargs, duration in segments:
            if elapsed < cumulative + duration:
                return model_type, model_kwargs
            cumulative += duration
        # 超出所有段范围，使用最后一段
        return segments[-1][0], segments[-1][1]
    
    @staticmethod
    def _convert_state_for_model(state: np.ndarray,
                                  from_model: str,
                                  to_model: str) -> np.ndarray:
        """在不同运动模型之间转换状态向量维度
        
        支持的转换：
        - CV/CT(4D)  ↔ CA(6D)    : 添加/删除加速度分量
        - CV/CT/CA   → RW(2D)    : 提取位置分量
        - RW(2D)     → CV/CT/CA  : 添加零速度/加速度分量
        
        Args:
            state: 当前状态向量
            from_model: 当前模型类型
            to_model: 目标模型类型
            
        Returns:
            转换后的状态向量
        """
        if from_model == to_model:
            return state.copy()
        
        # CV/CT(known) 4D ↔ CA 6D
        if from_model in ("CV", "CT") and to_model == "CA":
            return np.array([state[0], state[1], 0.0, state[2], state[3], 0.0])
        if from_model == "CA" and to_model in ("CV", "CT"):
            return np.array([state[0], state[1], state[3], state[4]])
        
        # CV/CT/CA → RW 2D（仅提取位置）
        if from_model in ("CV", "CT") and to_model == "RW":
            return np.array([state[0], state[2]])
        if from_model == "CA" and to_model == "RW":
            return np.array([state[0], state[3]])
        
        # RW 2D → CV/CT/CA（添加零速度/加速度）
        if from_model == "RW" and to_model in ("CV", "CT"):
            return np.array([state[0], 0.0, state[1], 0.0])
        if from_model == "RW" and to_model == "CA":
            return np.array([state[0], 0.0, 0.0, state[1], 0.0, 0.0])
        
        # 未覆盖的转换，保持原始状态
        return state.copy()
    
    def create_multi_target_appearance_scenario(self,
                                                 n_targets: int = 5,
                                                 duration: float = 200.0,
                                                 avg_lifetime: float = 50.0) -> None:
        """创建多目标出现消失场景
        
        目标初始位置与新生分量均值接近（均匀分布在场景中）
        
        Args:
            n_targets: 目标数量
            duration: 场景总时长
            avg_lifetime: 平均生存时间
        """
        # 新生分量均值位置（均匀分布在 [-600, 600] x [-600, 600]）
        birth_positions = [
            (-300, -300),
            (-300, 300),
            (300, -300),
            (300, 300)
        ]
        
        for i in range(n_targets):
            # 随机出生时间
            birth_time = np.random.uniform(0, duration * 0.3)
            
            # 随机生存时间（指数分布）
            lifetime = np.random.exponential(avg_lifetime)
            death_time = min(birth_time + lifetime, duration)
            
            # 选择一个新生分量位置作为基础，添加随机偏移
            base_pos = birth_positions[i % len(birth_positions)]
            x0 = base_pos[0] + np.random.uniform(-200, 200)
            y0 = base_pos[1] + np.random.uniform(-200, 200)
            vx0 = np.random.uniform(-10, 10)
            vy0 = np.random.uniform(-10, 10)
            
            initial_state = np.array([x0, vx0, y0, vy0])
            
            self.add_target(
                target_id=i,
                initial_state=initial_state,
                birth_time=birth_time,
                death_time=death_time,
                motion_model="CV"
            )
