"""
线性跟踪演示"""
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
from visualize import ScenarioVisualizer, MetricTracker


def run_linear_tracking_demo(duration: float = 100.0,
                              n_targets: int = 1,
                              time_step: float = 1.0,
                              process_noise_std: float = 0.1,
                              measurement_noise_std: float = 1.0,
                              clutter_rate: float = 5.0,
                              detection_probability: float = 0.9,
                              random_seed: int = 42,
                              show_plot: bool = True,
                              save_path: Optional[str] = None) -> Dict:
    """运行线性跟踪演示    
    演示使用卡尔曼滤波器进行线性目标跟踪    
    Args:
        duration: 场景时长（秒，        n_targets: 目标数量
        time_step: 时间步长
        process_noise_std: 过程噪声标准巨        measurement_noise_std: 观测噪声标准巨        clutter_rate: 杂波率        detection_probability: 检测概率        random_seed: 随机种子
        show_plot: 是否显示图形
        save_path: 保存路径
        
    Returns:
        结果字典
    """
    print("=" * 60)
    print("线性跟踪演示(Kalman Filter)")
    print("=" * 60)
    
    # 设置随机种子
    np.random.seed(random_seed)
    
    # 创建场景管理器
    scenario_manager = ScenarioManager(
        time_step=time_step,
        process_noise_std=process_noise_std,
        random_seed=random_seed
    )
    
    # 创建线性场景
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
    
    # 生成观测
    print("生成观测数据...")
    measurements_by_time = measurement_simulator.generate_scenario_measurements(
        scenario_data, time_step
    )
    
    # 创建卡尔曼滤波器
    kf = KalmanFilter(
        state_dim=4,
        measurement_dim=2,
        process_noise_std=process_noise_std,
        measurement_noise_std=measurement_noise_std
    )
    
    # 初始化滤波器
    # 使用第一个目标的初始状态
    initial_target = list(scenario_manager.targets)[0]
    initial_state = initial_target.initial_state
    initial_covariance = np.eye(4) * 100.0
    kf.initialize(initial_state, initial_covariance)
    
    # 运行滤波
    print("运行卡尔曼滤波..")
    estimated_states = []
    estimated_positions = []
    
    time_steps = sorted(measurements_by_time.keys())
    
    for t in time_steps:
        measurement_set = measurements_by_time[t]
        
        # 预测
        kf.predict(time_step)
        
        # 更新（使用最近的观测，
        if measurement_set.size > 0:
            # 简单选择最近的观测（这里简化处理）
            meas_array = measurement_set.get_measurements_array()
            
            # 计算预测观测
            predicted_meas = kf.state[[0, 2]]  # KF: z_pred = [x, y]
            
            # 找到最近的观测
            distances = np.linalg.norm(meas_array - predicted_meas, axis=1)
            nearest_idx = np.argmin(distances)
            nearest_meas = meas_array[nearest_idx]
            
            # 门限检某
            if distances[nearest_idx] < 100:  # 简单门限
                kf.update(nearest_meas)
        
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
            # 对齐长度
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
    
    # 可视区
    if show_plot:
        print("\n生成可视区..")
        visualizer = ScenarioVisualizer()
        
        # 准备观测数据
        measurement_list = [measurements_by_time[t].get_measurements_array() 
                           for t in time_steps]
        
        fig = visualizer.plot_scenario(
            true_trajectories=true_trajectories,
            measurements=measurement_list,
            estimated_trajectories=estimated_trajectories,
            title="Linear Tracking Demo (Kalman Filter)",
            show_measurements=True,
            show_estimates=True,
            show_true=True,
            save_path=save_path
        )
        
        plt.show()
    
    print("\n演示完成!")
    
    return {
        'scenario_manager': scenario_manager,
        'kf': kf,
        'estimated_states': estimated_states,
        'estimated_positions': estimated_positions,
        'metrics': metrics_summary,
        'true_trajectories': true_trajectories,
        'measurements': measurements_by_time
    }


if __name__ == '__main__':
    # 运行演示
    results = run_linear_tracking_demo(
        duration=50.0,
        n_targets=1,
        show_plot=True
    )

