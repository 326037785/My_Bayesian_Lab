"""
非线性跟踪演示示例
"""
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
from measurements import create_uniform_clutter, PolarMeasurementNoise, MeasurementSimulator, MeasurementType
from filter import ExtendedKalmanFilter, UnscentedKalmanFilter, CubatureKalmanFilter
from filter.extended_kalman_filter import polar_measurement_function, polar_measurement_jacobian
from filter.unscented_kalman_filter import polar_measurement_function_ukf
from filter.cubature_kalman_filter import polar_measurement_function_ckf
from visualize import ScenarioVisualizer, MetricTracker


def run_nonlinear_tracking_demo(duration: float = 100.0,
                                 n_targets: int = 1,
                                 time_step: float = 1.0,
                                 process_noise_std: float = 0.1,
                                 range_std: float = 10.0,
                                 bearing_std: float = 0.01,
                                 clutter_rate: float = 5.0,
                                 detection_probability: float = 0.9,
                                 filter_type: str = 'UKF',
                                 random_seed: int = 42,
                                 show_plot: bool = True,
                                 save_path: Optional[str] = None) -> Dict:
    """运行非线性跟踪演    
    演示使用非线性滤波器进行极坐标观测下的目标跟    
    Args:
        duration: 场景时长（秒        n_targets: 目标数量
        time_step: 时间步长
        process_noise_std: 过程噪声标准        range_std: 距离观测噪声标准        bearing_std: 角度观测噪声标准        clutter_rate: 杂波        detection_probability: 检测概        filter_type: 滤波器类('EKF', 'UKF', 'CKF')
        random_seed: 随机种子
        show_plot: 是否显示图形
        save_path: 保存路径
        
    Returns:
        结果字典
    """
    print("=" * 60)
    print(f"非线性跟踪演示({filter_type})")
    print("=" * 60)
    
    # 设置随机种子
    np.random.seed(random_seed)
    
    # 创建场景管理
    scenario_manager = ScenarioManager(
        time_step=time_step,
        process_noise_std=process_noise_std,
        random_seed=random_seed
    )
    
    # 创建线性场
    scenario_manager.create_linear_scenario(
        n_targets=n_targets,
        duration=duration,
        x_range=(-500, 500),
        y_range=(-500, 500),
        speed_range=(5, 20)
    )
    
    # 生成场景数据
    print("生成真实轨迹...")
    scenario_data = scenario_manager.generate_scenario(duration)
    
    # 创建观测模拟器（极坐标）
    noise_model = PolarMeasurementNoise(
        range_std=range_std,
        bearing_std=bearing_std
    )
    clutter_model = create_uniform_clutter(clutter_rate=clutter_rate)
    
    measurement_simulator = MeasurementSimulator(
        measurement_type=MeasurementType.POLAR,
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
    
    # 生成观测
    print("生成极坐标观测数..")
    measurements_by_time = measurement_simulator.generate_scenario_measurements(
        scenario_data, time_step
    )
    
    # 极坐标观测噪声协方差（正确的各向异性矩阵）
    polar_R = np.diag([range_std**2, bearing_std**2])

    # 创建滤波
    if filter_type == 'EKF':
        kf = ExtendedKalmanFilter(
            state_dim=4,
            measurement_dim=2,
            process_noise_std=process_noise_std,
            measurement_noise_std=range_std,
            state_transition_func=None,  # 使用默认
            measurement_func=polar_measurement_function,
            state_transition_jacobian=None,  # 使用默认
            measurement_jacobian=polar_measurement_jacobian,
            measurement_noise_matrix=polar_R  # 正确的极坐标噪声协方差
        )
    elif filter_type == 'UKF':
        kf = UnscentedKalmanFilter(
            state_dim=4,
            measurement_dim=2,
            process_noise_std=process_noise_std,
            measurement_noise_std=range_std,
            state_transition_func=None,  # 使用默认
            measurement_func=polar_measurement_function_ukf,
            measurement_noise_matrix=polar_R  # 正确的极坐标噪声协方差
        )
    elif filter_type == 'CKF':
        kf = CubatureKalmanFilter(
            state_dim=4,
            measurement_dim=2,
            process_noise_std=process_noise_std,
            measurement_noise_std=range_std,
            state_transition_func=None,  # 使用默认
            measurement_func=polar_measurement_function_ckf,
            measurement_noise_matrix=polar_R  # 正确的极坐标噪声协方差
        )
    else:
        raise ValueError(f"Unknown filter type: {filter_type}")
    
    # 初始化滤波器
    initial_target = list(scenario_manager.targets)[0]
    initial_state = initial_target.initial_state
    initial_covariance = np.eye(4) * 100.0
    kf.initialize(initial_state, initial_covariance)
    
    # 运行滤波
    print(f"运行{filter_type}滤波...")
    estimated_states = []
    estimated_positions = []
    
    time_steps = sorted(measurements_by_time.keys())
    
    for t in time_steps:
        measurement_set = measurements_by_time[t]

        # 预测
        kf.predict(time_step)

        # 更新
        if measurement_set.size > 0:
            # 获取极坐标观测及其协方差
            meas_array = measurement_set.get_measurements_array()
            # 获取第一个测量（非杂波）的协方差作为R
            target_meas = [m for m in measurement_set.measurements if not m.is_clutter]
            meas_R = target_meas[0].covariance if target_meas else polar_R

            # 计算预测观测（极坐标）
            if filter_type == 'EKF':
                z_pred = polar_measurement_function(kf.state)
            elif filter_type == 'UKF':
                z_pred = polar_measurement_function_ukf(kf.state)
            else:
                z_pred = polar_measurement_function_ckf(kf.state)

            # 使用马氏距离进行门限检测
            best_distance = float('inf')
            best_idx = -1
            for i, z in enumerate(meas_array):
                # 计算新息
                innovation = z - z_pred
                # 使用预测的新息协方差
                S = kf.get_innovation_covariance(meas_R)
                try:
                    S_inv = np.linalg.inv(S)
                    mahal_dist = innovation.T @ S_inv @ innovation
                    if mahal_dist < best_distance:
                        best_distance = mahal_dist
                        best_idx = i
                except np.linalg.LinAlgError:
                    continue

            # 门限检测：卡方分布 2 DOF, 99% 置信度 -> 9.21
            if best_idx >= 0 and best_distance < 9.21:
                kf.update(meas_array[best_idx], measurement_covariance=meas_R)

        # 记录状态
        estimated_states.append(kf.get_state())
        estimated_positions.append(kf.get_position())
    
    # 准备可视化数据
    true_trajectories = {}
    for target in scenario_manager.targets:
        positions = target.get_trajectory_positions()
        if len(positions) > 0:
            true_trajectories[target.target_id] = positions
    
    estimated_trajectories = {}
    if estimated_positions:
        est_pos_array = np.array([p for p in estimated_positions if p is not None])
        if len(est_pos_array) > 0:
            estimated_trajectories[0] = est_pos_array
    
    # 计算性能指标
    print("计算性能指标...")
    metric_tracker = MetricTracker()
    
    for target in scenario_manager.targets:
        true_positions = target.get_trajectory_positions()
        if len(true_positions) > 0 and len(estimated_positions) > 0:
            min_len = min(len(true_positions), len(estimated_positions))
            true_positions = true_positions[:min_len]
            est_positions = np.array([p for p in estimated_positions[:min_len] if p is not None])
            
            if len(est_positions) > 0:
                min_len = min(len(true_positions), len(est_positions))
                metric_tracker.update(true_positions[:min_len], est_positions[:min_len])
    
    metrics_summary = metric_tracker.get_summary()
    
    # 打印结果
    print("\n性能指标:")
    for key, value in metrics_summary.items():
        print(f"  {key}: {value:.4f}")
    
    # 可视
    if show_plot:
        print("\n生成可视..")
        visualizer = ScenarioVisualizer()
        
        # 准备观测数据（转换为笛卡尔坐标用于可视化
        measurement_list = []
        for t in time_steps:
            meas_array = measurements_by_time[t].get_measurements_array()
            if len(meas_array) > 0:
                # 转换极坐标到笛卡尔坐
                cartesian_meas = []
                for m in meas_array:
                    r, theta = m[0], m[1]
                    x = r * np.cos(theta)
                    y = r * np.sin(theta)
                    cartesian_meas.append([x, y])
                measurement_list.append(np.array(cartesian_meas))
            else:
                measurement_list.append(np.array([]))
        
        fig = visualizer.plot_scenario(
            true_trajectories=true_trajectories,
            measurements=measurement_list,
            estimated_trajectories=estimated_trajectories,
            title=f"Nonlinear Tracking Demo ({filter_type})",
            show_measurements=True,
            show_estimates=True,
            show_true=True,
            save_path=save_path
        )
        
        plt.show()
    
    print("\n演示完成!")
    
    return {
        'scenario_manager': scenario_manager,
        'filter': kf,
        'estimated_states': estimated_states,
        'estimated_positions': estimated_positions,
        'metrics': metrics_summary,
        'true_trajectories': true_trajectories,
        'measurements': measurements_by_time
    }


if __name__ == '__main__':
    # 运行演示
    results = run_nonlinear_tracking_demo(
        duration=50.0,
        n_targets=32,
        filter_type='EKF',
        show_plot=True
    )

