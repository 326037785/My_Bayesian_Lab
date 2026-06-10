"""
性能对比演示
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
from measurements import create_uniform_clutter, LinearMeasurementNoise, PolarMeasurementNoise, MeasurementSimulator, MeasurementType
from filter import KalmanFilter, ExtendedKalmanFilter, UnscentedKalmanFilter, CubatureKalmanFilter
from visualize import PerformanceVisualizer, MetricTracker


def polar_to_cartesian(measurements: np.ndarray) -> np.ndarray:
    """将极坐标观测 [range, bearing] 转换为笛卡尔坐标 [x, y]
    
    Args:
        measurements: 极坐标观测数组，形状上(N, 2)，每行为 [range, bearing]
        
    Returns:
        笛卡尔坐标数组，形状上(N, 2)，每行为 [x, y]
    """
    ranges = measurements[:, 0]
    bearings = measurements[:, 1]
    x = ranges * np.cos(bearings)
    y = ranges * np.sin(bearings)
    return np.column_stack([x, y])


def run_performance_comparison_demo(duration: float = 100.0,
                                     n_targets: int = 1,
                                     time_step: float = 1.0,
                                     process_noise_std: float = 0.1,
                                     measurement_noise_std: float = 1.0,
                                     clutter_rate: float = 5.0,
                                     detection_probability: float = 0.9,
                                     random_seed: int = 42,
                                     show_plot: bool = True,
                                     save_path: Optional[str] = None,
                                     scenario_type: str = "linear",
                                     polar_range_std: float = 10.0,
                                     polar_bearing_std: float = 0.01) -> Dict:
    """运行性能对比演示
    
    比较不同滤波器的性能
    
    Args:
        duration: 场景时长（秒，        n_targets: 目标数量
        time_step: 时间步长
        process_noise_std: 过程噪声标准巨        measurement_noise_std: 线性观测噪声标准差
        clutter_rate: 杂波率        detection_probability: 检测概率        random_seed: 随机种子
        show_plot: 是否显示图形
        save_path: 保存路径
        scenario_type: 场景类型，linear"（线性观测）或"nonlinear"（极坐标观测，        polar_range_std: 极坐标观测距离标准差（仅 nonlinear 场景，        polar_bearing_std: 极坐标观测角度标准差（仅 nonlinear 场景，        
    Returns:
        结果字典
    """
    print("=" * 60)
    print("性能对比演示")
    print("=" * 60)
    
    # 设置随机种子
    np.random.seed(random_seed)
    
    # 创建场景管理器
    scenario_manager = ScenarioManager(
        time_step=time_step,
        process_noise_std=process_noise_std,
        random_seed=random_seed
    )
    
    # 创建场景（目标运动为线性CV模型，
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
    
    # 根据场景类型创建观测模拟器
    if scenario_type == 'linear':
        noise_model = LinearMeasurementNoise(
            x_std=measurement_noise_std,
            y_std=measurement_noise_std
        )
        meas_type = MeasurementType.LINEAR
    else:
        noise_model = PolarMeasurementNoise(
            range_std=polar_range_std,
            bearing_std=polar_bearing_std
        )
        meas_type = MeasurementType.POLAR
    
    clutter_model = create_uniform_clutter(clutter_rate=clutter_rate)
    
    measurement_simulator = MeasurementSimulator(
        measurement_type=meas_type,
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
    print("生成观测数据...")
    measurements_by_time = measurement_simulator.generate_scenario_measurements(
        scenario_data, time_step
    )
    
    # 获取真实轨迹
    true_trajectories = {}
    for target in scenario_manager.targets:
        positions = target.get_trajectory_positions()
        if len(positions) > 0:
            true_trajectories[target.target_id] = positions
    
    # 定义要比较的滤波器
    if scenario_type == 'linear':
        # 线性场景：KF最优，EKF/UKF/CKF等价于KF
        filters_to_compare = {
            'KF': KalmanFilter(
                state_dim=4,
                measurement_dim=2,
                process_noise_std=process_noise_std,
                measurement_noise_std=measurement_noise_std
            ),
            'EKF': ExtendedKalmanFilter(
                state_dim=4,
                measurement_dim=2,
                process_noise_std=process_noise_std,
                measurement_noise_std=measurement_noise_std
            ),
            'UKF': UnscentedKalmanFilter(
                state_dim=4,
                measurement_dim=2,
                process_noise_std=process_noise_std,
                measurement_noise_std=measurement_noise_std
            ),
            'CKF': CubatureKalmanFilter(
                state_dim=4,
                measurement_dim=2,
                process_noise_std=process_noise_std,
                measurement_noise_std=measurement_noise_std
            )
        }
    else:
        # 非线性场景（极坐标观测）：KF不适用（线性模型无法处理 z=[range,bearing]）
        from filter.extended_kalman_filter import polar_measurement_function, polar_measurement_jacobian
        from filter.unscented_kalman_filter import polar_measurement_function_ukf
        from filter.cubature_kalman_filter import polar_measurement_function_ckf

        polar_R = np.diag([polar_range_std**2, polar_bearing_std**2])

        filters_to_compare = {
            'EKF': ExtendedKalmanFilter(
                state_dim=4,
                measurement_dim=2,
                process_noise_std=process_noise_std,
                measurement_noise_std=measurement_noise_std,
                measurement_func=polar_measurement_function,
                measurement_jacobian=polar_measurement_jacobian,
                measurement_noise_matrix=polar_R
            ),
            'UKF': UnscentedKalmanFilter(
                state_dim=4,
                measurement_dim=2,
                process_noise_std=process_noise_std,
                measurement_noise_std=measurement_noise_std,
                measurement_func=polar_measurement_function_ukf,
                measurement_noise_matrix=polar_R
            ),
            'CKF': CubatureKalmanFilter(
                state_dim=4,
                measurement_dim=2,
                process_noise_std=process_noise_std,
                measurement_noise_std=measurement_noise_std,
                measurement_func=polar_measurement_function_ckf,
                measurement_noise_matrix=polar_R
            )
        }
    
    # 运行所有滤波器
    results = {}
    time_steps = sorted(measurements_by_time.keys())
    
    for filter_name, kf in filters_to_compare.items():
        print(f"\n运行{filter_name}滤波...")
        
        # 初始化滤波器
        initial_target = list(scenario_manager.targets)[0]
        initial_state = initial_target.initial_state
        initial_covariance = np.eye(4) * 100.0
        kf.initialize(initial_state, initial_covariance)
        
        # 运行滤波
        estimated_positions = []
        
        for t in time_steps:
            measurement_set = measurements_by_time[t]
            
            # 预测
            kf.predict(time_step)
            
            # 更新
            if measurement_set.size > 0:
                meas_array = measurement_set.get_measurements_array()

                # 获取测量噪声协方差
                target_meas = [m for m in measurement_set.measurements if not m.is_clutter]
                meas_R = target_meas[0].covariance if target_meas else None

                # 根据场景类型计算预测观测
                if scenario_type == 'linear':
                    # 线性场景：所有滤波器使用线性测量模型 z = [x, y]
                    z_pred = kf.state[[0, 2]]
                else:
                    # 非线性场景（极坐标观测）：使用各滤波器的观测函数
                    z_pred = kf.h(kf.state)  # 直接调用观测函数 h(x) -> z

                # 使用马氏距离进行门限检测
                best_distance = float('inf')
                best_idx = -1
                for i, z in enumerate(meas_array):
                    innovation = z - z_pred
                    try:
                        S = kf.get_innovation_covariance(meas_R)
                        S_inv = np.linalg.inv(S)
                        mahal_dist = innovation.T @ S_inv @ innovation
                        if mahal_dist < best_distance:
                            best_distance = mahal_dist
                            best_idx = i
                    except np.linalg.LinAlgError:
                        continue

                # 门限：卡方分布 2 DOF, 99% -> 9.21
                if best_idx >= 0 and best_distance < 9.21:
                    kf.update(meas_array[best_idx], measurement_covariance=meas_R)
            
            # 记录状态
            estimated_positions.append(kf.get_position())
        
        # 存储结果
        results[filter_name] = {
            'estimated_positions': estimated_positions,
            'filter': kf
        }
    
    # 计算性能指标
    print("\n计算性能指标...")
    metrics_dict = {}
    
    for filter_name, result in results.items():
        metric_tracker = MetricTracker()
        
        for target in scenario_manager.targets:
            true_positions = target.get_trajectory_positions()
            est_positions = [p for p in result['estimated_positions'] if p is not None]
            
            if len(true_positions) > 0 and len(est_positions) > 0:
                min_len = min(len(true_positions), len(est_positions))
                metric_tracker.update(true_positions[:min_len], np.array(est_positions[:min_len]))
        
        metrics_dict[filter_name] = metric_tracker.get_summary()
    
    # 打印结果
    print("\n性能指标对比:")
    print("-" * 60)
    print(f"{'算法':<10} {'RMSE':<15} {'OSPA':<15} {'GOSPA':<15}")
    print("-" * 60)
    
    for filter_name, metrics in metrics_dict.items():
        rmse = metrics.get('rmse_mean', 0)
        ospa = metrics.get('ospa_mean', 0)
        gospa = metrics.get('gospa_mean', 0)
        print(f"{filter_name:<10} {rmse:<15.4f} {ospa:<15.4f} {gospa:<15.4f}")
    
    # 可视区
    if show_plot:
        print("\n生成可视区..")
        
        # 准备RMSE数据
        rmse_dict = {}
        for filter_name, result in results.items():
            # 计算每个时间步的RMSE
            rmse_values = []
            for target in scenario_manager.targets:
                true_positions = target.get_trajectory_positions()
                est_positions = [p for p in result['estimated_positions'] if p is not None]
                
                if len(true_positions) > 0 and len(est_positions) > 0:
                    min_len = min(len(true_positions), len(est_positions))
                    for i in range(min_len):
                        error = np.linalg.norm(true_positions[i] - est_positions[i])
                        rmse_values.append(error)
            
            rmse_dict[filter_name] = rmse_values
        
        # 绘制RMSE对比图
        perf_visualizer = PerformanceVisualizer()
        
        fig1 = perf_visualizer.plot_rmse_comparison(
            rmse_dict,
            time_steps=np.arange(len(time_steps)),
            title="RMSE Comparison",
            save_path=save_path.replace('.png', '_rmse.png') if save_path else None
        )
        
        # 绘制性能摘要
        fig2 = perf_visualizer.plot_performance_summary(
            metrics_dict,
            title="Performance Summary",
            save_path=save_path.replace('.png', '_summary.png') if save_path else None
        )
        
        plt.show()
    
    print("\n演示完成!")
    
    return {
        'scenario_manager': scenario_manager,
        'results': results,
        'metrics': metrics_dict,
        'true_trajectories': true_trajectories,
        'measurements': measurements_by_time
    }


if __name__ == '__main__':
    # 运行演示
    results = run_performance_comparison_demo(
        duration=50.0,
        n_targets=1,
        show_plot=True
    )

