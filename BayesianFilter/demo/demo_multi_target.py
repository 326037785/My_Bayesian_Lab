"""

多目标跟踪演示"""

import sys

import os

from pathlib import Path



# 添加项目根目录到Python路径

PROJECT_ROOT = Path(__file__).parent.parent

sys.path.insert(0, str(PROJECT_ROOT))



import numpy as np

import matplotlib.pyplot as plt

from typing import Optional, Dict, List



from ground_truth import ScenarioManager

from measurements import create_uniform_clutter, LinearMeasurementNoise, MeasurementSimulator, MeasurementType

from filter import KalmanFilter

from data_association import (BaseAssociation, NearestNeighborAssociation,

                              KNearestNeighborAssociation, JPDAFilter,

                              MHTFilter, PHDFilter, AdaptiveBirthPHDFilter)

from visualize import ScenarioVisualizer, MetricTracker





def run_multi_target_demo(duration: float = 100.0,

                           n_targets: int = 3,

                           time_step: float = 1.0,

                           process_noise_std: float = 0.1,

                           measurement_noise_std: float = 1.0,

                           clutter_rate: float = 5.0,

                           detection_probability: float = 0.9,

                           random_seed: int = 42,

                           show_plot: bool = True,

                           save_path: Optional[str] = None,

                           association_method: str = "knn",

                           association_params: Optional[Dict] = None,

                           association: Optional[BaseAssociation] = None) -> Dict:

    """运行多目标跟踪演示    

    演示使用卡尔曼滤波器和数据关联进行多目标跟踪

    

    Args:

        duration: 场景时长（秒，        n_targets: 目标数量

        time_step: 时间步长

        process_noise_std: 过程噪声标准巨        measurement_noise_std: 观测噪声标准巨        clutter_rate: 杂波率        detection_probability: 检测概率        random_seed: 随机种子

        show_plot: 是否显示图形

        save_path: 保存路径

        association_method: 关联算法名称 ("nn", "knn", "jpda", "mht")

        association_params: 关联算法参数字典（覆盖默认值）

        association: 自定义关联算法对象（优先级最高）

        

    Returns:

        结果字典

        

    Examples:

        默认 KNN 关联:

        >>> results = run_multi_target_demo()

        

        使用最近邻关联:

        >>> results = run_multi_target_demo(association_method="nn")

        

        使用 JPDA 关联并自定义参数:

        >>> results = run_multi_target_demo(

        ...     association_method="jpda",

        ...     association_params={"clutter_density": 1e-3}

        ... )

        

        使用自定义关联算法对�?

        >>> from data_association import MHTFilter

        >>> my_mht = MHTFilter(max_hypotheses=200)

        >>> results = run_multi_target_demo(association=my_mht)

    """

    # 创建数据关联器或 PHD 滤波器
    _use_phd = False
    _phd_filter = None

    

    if association is None:

        # 根据 association_method 创建关联器
        _base_params = {"gating_threshold": 9.21}
        if association_method == "nn":

            _params = {**_base_params, "use_mahalanobis": True}

            if association_params:

                _params.update(association_params)

            association = NearestNeighborAssociation(**_params)

            _method_label = "NN"

        elif association_method == "knn":

            _params = {**_base_params, "k": 3, "use_mahalanobis": True}

            if association_params:

                _params.update(association_params)

            association = KNearestNeighborAssociation(**_params)

            _method_label = "KNN"

        elif association_method == "jpda":

            _params = {**_base_params, "use_mahalanobis": True,

                       "detection_probability": detection_probability,

                       "clutter_density": 1e-4}

            if association_params:

                _params.update(association_params)

            association = JPDAFilter(**_params)

            _method_label = "JPDA"

        elif association_method == "mht":

            _params = {**_base_params, "use_mahalanobis": True}

            if association_params:

                _params.update(association_params)

            association = MHTFilter(**_params)

            _method_label = "MHT"

        elif association_method == "phd":

            # PHD 滤波器：使用默认参数初始化（surveillance_bounds 将在
            # 测量模拟器创建后通过 set_birth_components 更新）
            _phd_params = {
                "state_dim": 4,
                "measurement_dim": 2,
                "detection_probability": detection_probability,
                "clutter_rate": clutter_rate,
                "surveillance_area": 1200 * 1200,
            }

            if association_params:

                _phd_params.update(association_params)

            _phd_filter = PHDFilter(**_phd_params)

            _use_phd = True

            _method_label = "PHD"

        else:

            raise ValueError(

                f"Unknown association_method: {association_method}. "

                f"Choose from: 'nn', 'knn', 'jpda', 'mht', 'phd'"

            )

    else:

        _method_label = type(association).__name__



    print("=" * 60)

    print(f"多目标跟踪演示(KF + {_method_label} Association)")

    print("=" * 60)

    

    # 设置随机种子

    np.random.seed(random_seed)

    

    # 创建场景管理器
    scenario_manager = ScenarioManager(
        time_step=time_step,

        process_noise_std=process_noise_std,

        random_seed=random_seed

    )

    

    # 创建多目标场景
    scenario_manager.create_multi_target_appearance_scenario(
        n_targets=n_targets,

        duration=duration,

        avg_lifetime=50.0

    )

    

    # 生成场景数据

    print(f"生成{n_targets}个目标的真实轨迹...")

    scenario_data = scenario_manager.generate_scenario(duration)

    

    # 创建观测模拟器
    noise_model = LinearMeasurementNoise(
        x_std=measurement_noise_std,

        y_std=measurement_noise_std

    )

    clutter_model = create_uniform_clutter(clutter_rate=clutter_rate)

    

    measurement_simulator = MeasurementSimulator(

        measurement_type=MeasurementType.LINEAR,

        noise_model=noise_model,

        clutter_model=clutter_model,

        detection_probability=detection_probability,

        random_seed=random_seed

    )

    

    # 设置监视区域

    measurement_simulator.set_surveillance_region(

        x_range=(-600, 600),

        y_range=(-600, 600)

    )

    # PHD: 根据实际监视区域配置出生分量

    if _use_phd:

        bounds = ((-600, 600), (-600, 600))

        # 动态生成出生分量：16个分量覆盖 1200x1200 区域

        _phd_filter.birth_components = PHDFilter._auto_generate_birth_components(

            surveillance_bounds=bounds,

            n_components=4,

            birth_weight=0.03

        )

        _phd_filter.clutter_density = clutter_rate / (1200 * 1200)

        _phd_filter.R = np.diag([measurement_noise_std**2, measurement_noise_std**2])

    # 生成观测

    print("生成观测数据...")

    measurements_by_time = measurement_simulator.generate_scenario_measurements(

        scenario_data, time_step

    )

    

    # 创建多个卡尔曼滤波器（每个目标一个）- 仅在�?PHD 模式使用

    filters: Dict[int, KalmanFilter] = {}

    

    # 运行滤波

    print("运行多目标跟踪..")

    estimated_trajectories: Dict[int, List[np.ndarray]] = {}

    

    time_steps = sorted(measurements_by_time.keys())

    

    for t_idx, t in enumerate(time_steps):

        measurement_set = measurements_by_time[t]

        

        if _use_phd:

            # ===== PHD 模式：使用PHD 滤波器直接处理=====
            _phd_filter.predict(time_step)

            if measurement_set.size > 0:
                meas_array = measurement_set.get_measurements_array()
                _phd_filter.update(meas_array)

            

            # 以PHD 获取估计的目标状态
            phd_states, n_phd_targets = _phd_filter.extract_states()
            

            # 射PHD 估计的状态转换为轨迹

            if n_phd_targets > 0:

                for i, state in enumerate(phd_states):

                    if i not in estimated_trajectories:

                        estimated_trajectories[i] = []

                    # 提取位置 [x, y]（状态为 [x, vx, y, vy]，
                    position = np.array([state[0], state[2]])
                    estimated_trajectories[i].append(position)

        else:

            # ===== 标准模式：KF + 数据关联 =====

            # 预测所有滤波器

            for target_id, kf in filters.items():

                kf.predict(time_step)

            

            # 数据关联

            if measurement_set.size > 0 and len(filters) > 0:

                meas_array = measurement_set.get_measurements_array()

                

                # 预测观测

                predicted_measurements = []

                innovation_covariances = []

                

                for target_id in sorted(filters.keys()):

                    kf = filters[target_id]

                    z_pred = kf.state[[0, 2]]  # KF: z_pred = [x, y]

                    predicted_measurements.append(z_pred)

                    

                    S = kf.get_innovation_covariance()

                    innovation_covariances.append(S)

                

                predicted_measurements = np.array(predicted_measurements)

                

                # 执行关联

                result = association.associate(

                    meas_array,

                    predicted_measurements,

                    innovation_covariances=innovation_covariances

                )

                

                # 根据关联结果更新滤波器
                for meas_idx, target_idx in result.associations.items():
                    if target_idx < len(sorted(filters.keys())):

                        target_id = sorted(filters.keys())[target_idx]

                        kf = filters[target_id]

                        meas_cov = measurement_set.measurements[meas_idx].covariance if meas_idx < len(measurement_set.measurements) else None
                        kf.update(meas_array[meas_idx], measurement_covariance=meas_cov)

            

            # 记录估计结果

            for target_id, kf in filters.items():

                if target_id not in estimated_trajectories:

                    estimated_trajectories[target_id] = []

                estimated_trajectories[target_id].append(kf.get_position())

            

            # 检查是否有新目标出率
            for target in scenario_manager.targets:
                if target.is_alive(t) and target.target_id not in filters:

                    # 创建新的滤波器
                    initial_state = target.get_state_at_time(t)
                    if initial_state is not None:

                        kf = KalmanFilter(

                            state_dim=4,

                            measurement_dim=2,

                            process_noise_std=process_noise_std,

                            measurement_noise_std=measurement_noise_std

                        )

                        kf.initialize(initial_state, np.eye(4) * 100.0)

                        filters[target.target_id] = kf

            

            # 检查是否有目标消失

            for target_id in list(filters.keys()):

                target = scenario_manager.get_target_by_id(target_id)

                if target is not None and not target.is_alive(t):

                    del filters[target_id]

    

    # 准备可视化数据
    true_trajectories = {}
    for target in scenario_manager.targets:

        positions = target.get_trajectory_positions()

        if len(positions) > 0:

            true_trajectories[target.target_id] = positions

    

    # 转换估计轨迹

    est_trajectories = {}

    for target_id, positions in estimated_trajectories.items():

        if positions:

            est_trajectories[target_id] = np.array([p for p in positions if p is not None])

    

    # 计算性能指标

    print("计算性能指标...")

    metric_tracker = MetricTracker()

    

    # 正确的OSPA计算：每个时刻独立计算（参考MATLAB RFS Toolbox，
    for t_idx, t in enumerate(time_steps):
        # 收集该时刻所有真实目标的位置

        true_positions_at_t = []

        for target in scenario_manager.targets:

            if target.is_alive(t):

                pos = target.get_position_at_time(t)

                if pos is not None:

                    true_positions_at_t.append(pos)

        

        # 收集该时刻所有估计目标的位置

        est_positions_at_t = []

        for target_id, positions in estimated_trajectories.items():

            if t_idx < len(positions):

                est_positions_at_t.append(positions[t_idx])

        

        # 转换为数统
        true_arr = np.array(true_positions_at_t) if true_positions_at_t else np.array([]).reshape(0, 2)
        est_arr = np.array(est_positions_at_t) if est_positions_at_t else np.array([]).reshape(0, 2)

        

        # 更新指标（MetricTracker内部会处理空集合的情况）

        metric_tracker.update(true_arr, est_arr)

    

    metrics_summary = metric_tracker.get_summary()

    

    # 打印结果

    print("\n性能指标:")

    for key, value in metrics_summary.items():

        print(f"  {key}: {value:.4f}")

    

    # 可视区
    if show_plot:
        print("\n生成可视区..")

        visualizer = ScenarioVisualizer()

        

        # 准备观测数据

        measurement_list = [measurements_by_time[t].get_measurements_array() 

                           for t in time_steps]

        

        # 1. 绘制跟踪场景

        fig = visualizer.plot_scenario(

            true_trajectories=true_trajectories,

            measurements=measurement_list,

            estimated_trajectories=est_trajectories,

            title=f"Multi-Target Tracking (KF + {_method_label})",

            show_measurements=True,

            show_estimates=True,

            show_true=True,

            save_path=save_path

        )

        

        # 2. 绘制OSPA误差曲线（参考MATLAB RFS Toolbox，
        ospa_values = np.array(metric_tracker.ospa_values)
        localization_errors = np.array(metric_tracker.localization_errors)

        cardinality_errors = np.array(metric_tracker.cardinality_errors)

        

        fig_ospa = visualizer.plot_ospa_over_time(

            ospa_values=ospa_values,

            localization_errors=localization_errors,

            cardinality_errors=cardinality_errors,

            time_steps=np.array(time_steps),

            title=f"OSPA Error (KF + {_method_label})",

            save_path=save_path.replace('.png', '_ospa.png') if save_path else None

        )

        

        # 3. 绘制基数估计

        true_cardinality = []

        estimated_cardinality = []

        for t_idx, t in enumerate(time_steps):

            # 真实基数

            n_true = sum(1 for target in scenario_manager.targets if target.is_alive(t))

            true_cardinality.append(n_true)

            

            # 估计基数

            n_est = sum(1 for positions in estimated_trajectories.values() 

                       if t_idx < len(positions))

            estimated_cardinality.append(n_est)

        

        fig_card = visualizer.plot_cardinality_over_time(

            true_cardinality=np.array(true_cardinality),

            estimated_cardinality=np.array(estimated_cardinality),

            time_steps=np.array(time_steps),

            title=f"Cardinality Estimation (KF + {_method_label})",

            save_path=save_path.replace('.png', '_cardinality.png') if save_path else None

        )

        

        plt.show()

    

    print("\n演示完成!")

    

    return {

        'scenario_manager': scenario_manager,

        'filters': filters,

        'estimated_trajectories': est_trajectories,

        'metrics': metrics_summary,

        'true_trajectories': true_trajectories,

        'measurements': measurements_by_time

    }





if __name__ == '__main__':

    # 运行演示

    results = run_multi_target_demo(

        duration=50.0,

        n_targets=3,

        show_plot=True

    )



