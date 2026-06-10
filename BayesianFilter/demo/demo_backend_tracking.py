"""

演示：使用FilterBackend 架构进行多目标跟踪

关键优势，1. 数据关联算法只依资FilterBackend 接口，不关心具体滤波器2. 可以轻松切换 KF/EKF/UKF/CKF，无需修改关联代码

3. 预测观测、新息协方差、门限计算统一在后端中

"""

import sys

from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))



import numpy as np

from typing import Dict, List



from filter import (

    KFBackend, EKFBackend, UKFBackend,

    MultiTargetFilterManager, PredictedState

)

from data_association import KNearestNeighborAssociation, JPDAFilter, MHTFilter

from ground_truth import ScenarioManager

from measurements import create_uniform_clutter, LinearMeasurementNoise, MeasurementSimulator, MeasurementType





def run_tracking_with_backend(

    backend_type: str = "kf",

    association_method: str = "knn",

    n_targets: int = 3,

    duration: float = 50.0,

    time_step: float = 1.0,

    process_noise_std: float = 0.1,

    measurement_noise_std: float = 1.0

) -> Dict:

    """使用指定后端运行多目标跟踪    

    Args:

        backend_type: 后端类型 ("kf", "ekf", "ukf")

        association_method: 关联方法 ("knn", "jpda", "mht")

        n_targets: 目标数量

        duration: 场景时长

        time_step: 时间步长

        process_noise_std: 过程噪声标准巨        measurement_noise_std: 观测噪声标准巨        

    Returns:

        跟踪结果字典

    """

    print(f"=" * 60)

    print(f"多目标跟踪 {backend_type.upper()} + {association_method.upper()}")

    print(f"=" * 60)

    

    # 1. 创建后端工厂

    if backend_type == "kf":

        backend_factory = lambda: KFBackend(

            state_dim=4,

            meas_dim=2,

            process_noise_std=process_noise_std,

            measurement_noise_std=measurement_noise_std

        )

    elif backend_type == "ekf":

        # EKF 使用非线性观测模型（示例：极坐标，
        from filter.extended_kalman_filter import polar_measurement_function, polar_measurement_jacobian
        backend_factory = lambda: EKFBackend(

            state_dim=4,

            meas_dim=2,

            process_noise_std=process_noise_std,

            measurement_noise_std=measurement_noise_std,

            measurement_func=polar_measurement_function,

            measurement_jacobian=polar_measurement_jacobian

        )

    elif backend_type == "ukf":

        from filter.unscented_kalman_filter import polar_measurement_function_ukf

        backend_factory = lambda: UKFBackend(

            state_dim=4,

            meas_dim=2,

            process_noise_std=process_noise_std,

            measurement_noise_std=measurement_noise_std,

            measurement_func=polar_measurement_function_ukf

        )

    else:

        raise ValueError(f"Unknown backend: {backend_type}")

    

    # 2. 创建多目标管理器

    manager = MultiTargetFilterManager(backend_factory)

    

    # 3. 创建数据关联器（只依资PredictedState，不依赖具体滤波器）

    if association_method == "knn":

        association = KNearestNeighborAssociation(k=3)

    elif association_method == "jpda":

        association = JPDAFilter()

    elif association_method == "mht":

        association = MHTFilter()

    else:

        raise ValueError(f"Unknown association: {association_method}")

    

    # 4. 创建场景

    np.random.seed(42)

    scenario_manager = ScenarioManager(

        time_step=time_step,

        process_noise_std=process_noise_std,

        random_seed=42

    )

    scenario_manager.create_linear_scenario(n_targets=n_targets, duration=duration)

    scenario_data = scenario_manager.generate_scenario(duration)

    

    # 5. 创建观测模拟器
    noise_model = LinearMeasurementNoise(
        x_std=measurement_noise_std,

        y_std=measurement_noise_std

    )

    clutter_model = create_uniform_clutter(clutter_rate=5.0)

    measurement_simulator = MeasurementSimulator(

        measurement_type=MeasurementType.LINEAR,

        noise_model=noise_model,

        clutter_model=clutter_model,

        detection_probability=0.9,

        random_seed=42

    )

    measurement_simulator.set_surveillance_region(

        x_range=(-600, 600),

        y_range=(-600, 600)

    )

    measurements_by_time = measurement_simulator.generate_scenario_measurements(

        scenario_data, time_step

    )

    

    # 6. 运行跟踪

    print("运行跟踪...")

    time_steps = sorted(measurements_by_time.keys())

    estimated_trajectories: Dict[int, List[np.ndarray]] = {}

    

    for t in time_steps:

        measurement_set = measurements_by_time[t]

        

        # 预测所有目标
        manager.predict_all(time_step)
        

        # 检查新目标出现

        for target in scenario_manager.targets:

            if target.is_alive(t) and target.target_id not in manager.get_target_ids():

                initial_state = target.get_state_at_time(t)

                if initial_state is not None:

                    manager.add_target(

                        target.target_id,

                        initial_state,

                        np.eye(4) * 100.0

                    )

        

        # 数据关联

        if measurement_set.size > 0 and len(manager.get_target_ids()) > 0:

            meas_array = measurement_set.get_measurements_array()

            

            # 获取预测状态（统一接口，不依赖具体滤波器）

            predicted_states = manager.get_predicted_states()

            

            # 构建关联所需的数据
            predicted_measurements = np.array([ps.predicted_meas for ps in predicted_states])
            innovation_covariances = [ps.innovation_cov for ps in predicted_states]

            

            # 执行关联

            result = association.associate(

                meas_array,

                predicted_measurements,

                innovation_covariances=innovation_covariances

            )

            

            # 更新关联的目标
            target_ids = manager.get_target_ids()
            for meas_idx, target_idx in result.associations.items():

                if target_idx < len(target_ids):

                    manager.update_target(target_ids[target_idx], meas_array[meas_idx])

        

        # 记录轨迹

        for ps in manager.get_predicted_states():

            if ps.target_id not in estimated_trajectories:

                estimated_trajectories[ps.target_id] = []

            estimated_trajectories[ps.target_id].append(

                np.array([ps.state[0], ps.state[2]])  # [x, y]

            )

        

        # 检查目标消�?
        for target in scenario_manager.targets:
            if not target.is_alive(t) and target.target_id in manager.get_target_ids():

                manager.remove_target(target.target_id)

    

    # 7. 计算性能指标

    from visualize import MetricTracker

    metric_tracker = MetricTracker()

    

    for target in scenario_manager.targets:

        if target.target_id in estimated_trajectories:

            true_pos = target.get_trajectory_positions()

            est_pos = np.array(estimated_trajectories[target.target_id])

            min_len = min(len(true_pos), len(est_pos))

            if min_len > 0:

                metric_tracker.update(true_pos[:min_len], est_pos[:min_len])

    

    metrics = metric_tracker.get_summary()

    

    print(f"\n性能指标:")

    for key, value in metrics.items():

        print(f"  {key}: {value:.4f}")

    

    return {

        'backend': backend_type,

        'association': association_method,

        'metrics': metrics,

        'estimated_trajectories': estimated_trajectories

    }





if __name__ == "__main__":

    # 演示：轻松切换不同后端
    backends = ["kf", "ukf"]
    associations = ["knn", "jpda"]

    

    for backend in backends:

        for assoc in associations:

            try:

                result = run_tracking_with_backend(

                    backend_type=backend,

                    association_method=assoc,

                    n_targets=3,

                    duration=30.0

                )

            except Exception as e:

                print(f"Error with {backend}+{assoc}: {e}")



