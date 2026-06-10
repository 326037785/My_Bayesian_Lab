"""
正确的场景演示
场景1：线性跟踪（笛卡尔坐标）
- 状态：[x, vx, y, vy]
- 测量：z = [x, y] = H @ x
- 滤波器：KF（最优）

场景2：雷达跟踪（极坐标观测）
- 状态：[x, vx, y, vy]
- 测量：z = [range, bearing] = h(x)
- 滤波器：EKF（雅可比）或 UKF（sigma点）
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
from typing import Dict, List, Callable

from filter.tracking_backends import LinearKFBackend, EKFBackend, UKFBackend
from filter.filter_backend import MultiTargetBackend, PredictedState
from data_association import KNearestNeighborAssociation, JPDAFilter
from ground_truth import ScenarioManager
from measurements import (
    MeasurementSimulator, MeasurementType, 
    LinearMeasurementNoise, PolarMeasurementNoise, UniformClutter
)


# ============================================================
# 雷达测量函数（极坐标，# ============================================================

def radar_measurement(x: np.ndarray) -> np.ndarray:
    """雷达测量函数：笛卡尔 →极坐标    
    x = [px, vx, py, vy]
    z = [range, bearing] = [sqrt(px^2 + py^2), atan2(py, px)]
    """
    px, py = x[0], x[2]
    range_meas = np.sqrt(px**2 + py**2)
    bearing_meas = np.arctan2(py, px)
    return np.array([range_meas, bearing_meas])


def radar_jacobian(x: np.ndarray) -> np.ndarray:
    """雷达测量雅可比矩阵    
    dh/dx = [[px/r, 0, py/r, 0],
             [-py/r^2, 0, px/r^2, 0]]
    """
    px, py = x[0], x[2]
    r = np.sqrt(px**2 + py**2)
    r2 = r**2
    
    H = np.array([
        [px/r, 0, py/r, 0],
        [-py/r2, 0, px/r2, 0]
    ])
    return H


# ============================================================
# 多目标管理器
# ============================================================

class MultiTargetManager(MultiTargetBackend):
    """多目标管理器"""
    
    def __init__(self, backend_factory):
        self._factory = backend_factory
        self._targets = {}
    
    def add_target(self, target_id, state, cov=None):
        backend = self._factory()
        backend.initialize(state, cov)
        self._targets[target_id] = backend
    
    def remove_target(self, target_id):
        if target_id in self._targets:
            del self._targets[target_id]
    
    def predict_all(self, dt):
        for b in self._targets.values():
            b.predict(dt)
    
    def get_predicted_states(self):
        states = []
        for tid, b in self._targets.items():
            if b.is_initialized:
                states.append(PredictedState(
                    state=b.get_state(),
                    covariance=b.get_covariance(),
                    predicted_meas=b.get_predicted_measurement(),
                    innovation_cov=b.get_innovation_covariance(),
                    target_id=tid
                ))
        return states
    
    def update_target(self, target_id, z):
        if target_id in self._targets:
            self._targets[target_id].update(z)
    
    def get_target_ids(self):
        return list(self._targets.keys())


# ============================================================
# 场景1：线性跟踪# ============================================================

def run_linear_tracking(n_targets=3, duration=50.0):
    """线性跟踪场景
    - 测量：z = [x, y]（笛卡尔坐标直接观测）
    - 滤波器：KF（线性最优）
    """
    print("=" * 60)
    print("场景1：线性跟踪（KF）")
    print("测量模型：z = H @ x = [x, y]")
    print("=" * 60)
    
    # 创建场景
    np.random.seed(42)
    scenario = ScenarioManager(time_step=1.0, process_noise_std=0.1, random_seed=42)
    scenario.create_linear_scenario(n_targets=n_targets, duration=duration)
    scenario_data = scenario.generate_scenario(duration)
    
    # 线性测量模拟器
    noise_model = LinearMeasurementNoise(x_std=1.0, y_std=1.0)
    clutter = create_uniform_clutter(clutter_rate=5.0)
    meas_sim = MeasurementSimulator(
        measurement_type=MeasurementType.LINEAR,
        noise_model=noise_model,
        clutter_model=clutter,
        detection_probability=0.9,
        random_seed=42
    )
    meas_sim.set_surveillance_region((-600, 600), (-600, 600))
    measurements = meas_sim.generate_scenario_measurements(scenario_data, 1.0)
    
    # KF后端（线性测量）
    manager = MultiTargetManager(lambda: LinearKFBackend(
        state_dim=4, meas_dim=2,
        process_noise_std=0.1, measurement_noise_std=1.0
    ))
    
    # 关联器
    association = KNearestNeighborAssociation(k=3)
    
    # 运行跟踪
    est_trajectories = {}
    time_steps = sorted(measurements.keys())
    
    for t in time_steps:
        meas_set = measurements[t]
        manager.predict_all(1.0)
        
        # 添加新目标
        for target in scenario.targets:
            if target.is_alive(t) and target.target_id not in manager.get_target_ids():
                state = target.get_state_at_time(t)
                if state is not None:
                    manager.add_target(target.target_id, state, np.eye(4) * 100)
        
        # 数据关联
        if meas_set.size > 0 and len(manager.get_target_ids()) > 0:
            meas_array = meas_set.get_measurements_array()
            pred_states = manager.get_predicted_states()
            
            pred_meas = np.array([ps.predicted_meas for ps in pred_states])
            innov_covs = [ps.innovation_cov for ps in pred_states]
            
            result = association.associate(meas_array, pred_meas, innovation_covariances=innov_covs)
            
            target_ids = manager.get_target_ids()
            for mi, ti in result.associations.items():
                if ti < len(target_ids):
                    manager.update_target(target_ids[ti], meas_array[mi])
        
        # 记录
        for ps in manager.get_predicted_states():
            if ps.target_id not in est_trajectories:
                est_trajectories[ps.target_id] = []
            est_trajectories[ps.target_id].append(np.array([ps.state[0], ps.state[2]]))
        
        # 移除消失目标
        for target in scenario.targets:
            if not target.is_alive(t) and target.target_id in manager.get_target_ids():
                manager.remove_target(target.target_id)
    
    # 计算指标
    from visualize import MetricTracker
    tracker = MetricTracker()
    for target in scenario.targets:
        if target.target_id in est_trajectories:
            true_pos = target.get_trajectory_positions()
            est_pos = np.array(est_trajectories[target.target_id])
            min_len = min(len(true_pos), len(est_pos))
            if min_len > 0:
                tracker.update(true_pos[:min_len], est_pos[:min_len])
    
    metrics = tracker.get_summary()
    print(f"\n性能指标:")
    for k, v in metrics.items():
        print(f"  {k}: {v:.4f}")
    
    return metrics


# ============================================================
# 场景2：雷达跟踪（极坐标观测）
# ============================================================

def run_radar_tracking(n_targets=3, duration=50.0, use_ekf=True):
    """雷达跟踪场景
    
    - 测量：z = [range, bearing]（极坐标）
    - 滤波器：EKF（雅可比）或 UKF（sigma点）
    """
    filter_name = "EKF" if use_ekf else "UKF"
    print("=" * 60)
    print(f"场景2：雷达跟踪（{filter_name}）")
    print("测量模型：z = [range, bearing] = h(x)")
    print("=" * 60)
    
    # 创建场景
    np.random.seed(42)
    scenario = ScenarioManager(time_step=1.0, process_noise_std=0.1, random_seed=42)
    scenario.create_linear_scenario(n_targets=n_targets, duration=duration)
    scenario_data = scenario.generate_scenario(duration)
    
    # 极坐标测量模拟器
    noise_model = PolarMeasurementNoise(range_std=10.0, bearing_std=0.01)
    clutter = create_uniform_clutter(clutter_rate=5.0)
    meas_sim = MeasurementSimulator(
        measurement_type=MeasurementType.POLAR,
        noise_model=noise_model,
        clutter_model=clutter,
        detection_probability=0.9,
        random_seed=42
    )
    meas_sim.set_surveillance_region((-600, 600), (-600, 600))
    measurements = meas_sim.generate_scenario_measurements(scenario_data, 1.0)
    
    # 创建后端（非线性测量）
    if use_ekf:
        polar_R = np.diag([10.0**2, 0.01**2])  # range_std=10m, bearing_std=0.01rad
        backend_factory = lambda: EKFBackend(
            state_dim=4, meas_dim=2,
            process_noise_std=0.1, measurement_noise_std=1.0,
            h=radar_measurement,
            H_jacobian=radar_jacobian,
            R=polar_R
        )
    else:
        polar_R = np.diag([10.0**2, 0.01**2])  # range_std=10m, bearing_std=0.01rad
        backend_factory = lambda: UKFBackend(
            state_dim=4, meas_dim=2,
            process_noise_std=0.1, measurement_noise_std=1.0,
            h=radar_measurement,
            R=polar_R
        )
    
    manager = MultiTargetManager(backend_factory)
    association = KNearestNeighborAssociation(k=3)
    
    # 运行跟踪
    est_trajectories = {}
    time_steps = sorted(measurements.keys())
    
    for t in time_steps:
        meas_set = measurements[t]
        manager.predict_all(1.0)
        
        # 添加新目标
        for target in scenario.targets:
            if target.is_alive(t) and target.target_id not in manager.get_target_ids():
                state = target.get_state_at_time(t)
                if state is not None:
                    manager.add_target(target.target_id, state, np.eye(4) * 100)
        
        # 数据关联（在极坐标空间比较）
        if meas_set.size > 0 and len(manager.get_target_ids()) > 0:
            meas_array = meas_set.get_measurements_array()
            pred_states = manager.get_predicted_states()
            
            # 预测观测是极坐标
            pred_meas = np.array([ps.predicted_meas for ps in pred_states])
            innov_covs = [ps.innovation_cov for ps in pred_states]
            
            result = association.associate(meas_array, pred_meas, innovation_covariances=innov_covs)
            
            target_ids = manager.get_target_ids()
            for mi, ti in result.associations.items():
                if ti < len(target_ids):
                    manager.update_target(target_ids[ti], meas_array[mi])
        
        # 记录轨迹（从状态提取位置）
        for ps in manager.get_predicted_states():
            if ps.target_id not in est_trajectories:
                est_trajectories[ps.target_id] = []
            est_trajectories[ps.target_id].append(np.array([ps.state[0], ps.state[2]]))
        
        # 移除消失目标
        for target in scenario.targets:
            if not target.is_alive(t) and target.target_id in manager.get_target_ids():
                manager.remove_target(target.target_id)
    
    # 计算指标
    from visualize import MetricTracker
    tracker = MetricTracker()
    for target in scenario.targets:
        if target.target_id in est_trajectories:
            true_pos = target.get_trajectory_positions()
            est_pos = np.array(est_trajectories[target.target_id])
            min_len = min(len(true_pos), len(est_pos))
            if min_len > 0:
                tracker.update(true_pos[:min_len], est_pos[:min_len])
    
    metrics = tracker.get_summary()
    print(f"\n性能指标:")
    for k, v in metrics.items():
        print(f"  {k}: {v:.4f}")
    
    return metrics


# ============================================================
# 主程序# ============================================================

if __name__ == "__main__":
    # 场景1：线性跟踪（KF最优）
    m1 = run_linear_tracking(n_targets=3, duration=30.0)
    
    print("\n")
    
    # 场景2：雷达跟踪（EKF，
    m2 = run_radar_tracking(n_targets=3, duration=30.0, use_ekf=True)
    
    print("\n")
    
    # 场景2：雷达跟踪（UKF，
    m3 = run_radar_tracking(n_targets=3, duration=30.0, use_ekf=False)

