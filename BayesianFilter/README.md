# Bayesian Filter 示例项目

一个完整的贝叶斯滤波使用示例，包含从目标生成、观测模型、滤波算法到性能评估的全流程�?
## 功能特�?
### 1. 真实场景生成 (ground_truth)
- **运动模型**：匀速直�?CV)、匀加�?CA)、协调转�?CT)、随机游�?RW)
- **场景管理**：多目标、目标出�?消失时间控制

### 2. 观测模型 (measurements)
- **观测类型**：线性观测、极坐标观测
- **噪声模型**：高斯噪声、极坐标噪声、非高斯噪声
- **杂波模型**：泊松杂波、均匀杂波、非均匀杂波

### 3. 滤波算法 (filter)

#### 高斯滤波器族
- **KF** - 卡尔曼滤波器
- **EKF** - 扩展卡尔曼滤波器
- **UKF** - 无迹卡尔曼滤波器
- **CKF** - 容积卡尔曼滤波器

#### 粒子滤波器族
- **PF** - 标准粒子滤波�?- **APF** - 辅助粒子滤波�?- **RBPF** - Rao-Blackwellized粒子滤波�?- **UPF** - 无迹粒子滤波�?
### 4. 数据关联 (data_association)
- **KNN** - K最近邻
- **JPDA** - 联合概率数据关联
- **MHT** - 多假设跟�?- **PHD** - 概率假设密度滤波�?
### 5. 性能评估 (visualize)
- **RMSE** - 均方根误�?- **OSPA** - 最优子模式分配
- **GOSPA** - 广义最优子模式分配

## 安装

```bash
# 克隆项目
git clone <repository-url>
cd BayesianFilter

# 安装依赖
pip install -r requirements.txt
```

## 使用方法

### 命令行运�?
```bash
# 线性跟踪演�?python run_demo.py --demo linear --duration 50 --n_targets 1

# 非线性跟踪演�?(UKF)
python run_demo.py --demo nonlinear --duration 50 --filter UKF

# 多目标跟踪演�?python run_demo.py --demo multi --duration 50 --n_targets 3

# 多目标跟�?+ JPDA 关联
python run_demo.py --demo multi --duration 50 --n_targets 3 --association jpda

# 多目标跟�?+ MHT 关联
python run_demo.py --demo multi --duration 50 --n_targets 3 --association mht

# 多目标跟�?+ PHD 滤波器（完整多目标滤波器，非关联算法�?python run_demo.py --demo multi --duration 50 --n_targets 3 --association phd

# 多目标跟�?+ 最近邻关联
python run_demo.py --demo multi --duration 50 --n_targets 3 --association nn

# 性能对比演示
python run_demo.py --demo compare --duration 50

# 不显示图�?python run_demo.py --demo linear --no_plot

# 保存图形
python run_demo.py --demo linear --save result.png
```

### Python API

#### 基础用法：线性场�?+ 卡尔曼滤�?
```python
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from ground_truth import ScenarioManager
from measurements import MeasurementSimulator, MeasurementType, LinearMeasurementNoise, UniformClutter
from filter import KalmanFilter
from visualize import MetricTracker

# 创建场景
scenario_manager = ScenarioManager(time_step=1.0, process_noise_std=0.1)
scenario_manager.create_linear_scenario(n_targets=1, duration=50.0)
scenario_data = scenario_manager.generate_scenario(50.0)

# 创建观测模拟�?noise_model = LinearMeasurementNoise(x_std=1.0, y_std=1.0)
clutter_model = UniformClutter(clutter_rate=5.0)
measurement_simulator = MeasurementSimulator(
    measurement_type=MeasurementType.LINEAR,
    noise_model=noise_model,
    clutter_model=clutter_model,
    detection_probability=0.9
)

# 生成观测
measurements_by_time = measurement_simulator.generate_scenario_measurements(scenario_data, 1.0)

# 创建卡尔曼滤波器
kf = KalmanFilter(state_dim=4, measurement_dim=2)
kf.initialize(initial_state, initial_covariance)

# 运行滤波
for t, measurement_set in measurements_by_time.items():
    kf.predict(1.0)
    if measurement_set.size > 0:
        kf.update(measurement_set.get_measurements_array())
    position = kf.get_position()
    print(f"t={t}: 估计位置 = {position}")
```

#### 机动目标场景

```python
from ground_truth import ScenarioManager

# 协调转弯场景（所有目标匀速转弯）
scenario_manager = ScenarioManager(time_step=1.0, process_noise_std=0.1)
scenario_manager.create_coordinated_turn_scenario(
    n_targets=3,
    duration=100.0,
    turn_rate=0.1  # 转弯速率 (rad/s)
)
scenario_data = scenario_manager.generate_scenario(100.0)

# 随机游走场景（布朗运动）
scenario_manager = ScenarioManager(time_step=1.0, process_noise_std=0.1)
scenario_manager.create_random_walk_scenario(
    n_targets=2,
    duration=100.0,
    walk_std=5.0  # 游走标准�?)
scenario_data = scenario_manager.generate_scenario(100.0)

# 混合机动场景（CV/CA/CT 随机分配�?scenario_manager = ScenarioManager(time_step=1.0, process_noise_std=0.1)
scenario_manager.create_multi_model_scenario(
    n_targets=5,
    duration=200.0
)
scenario_data = scenario_manager.generate_scenario(200.0)
```

#### 性能对比：线�?vs 非线性场�?
```python
from demo.demo_performance_comparison import run_performance_comparison_demo

# 线性场景（KF 应该最优，所有滤波器使用线性测量模型）
results_linear = run_performance_comparison_demo(
    duration=100.0,
    scenario_type="linear",      # 线性观�?z = H @ x
    show_plot=True
)

# 非线性场景（极坐标观测，EKF/UKF/CKF 应优�?KF�?results_nonlinear = run_performance_comparison_demo(
    duration=100.0,
    scenario_type="nonlinear",   # 极坐标观�?[range, bearing]
    polar_range_std=10.0,        # 距离噪声标准�?    polar_bearing_std=0.01,      # 角度噪声标准�?(rad)
    show_plot=True
)
```

#### 多目标跟踪：选择数据关联算法

```python
from demo.demo_multi_target import run_multi_target_demo

# 默认 KNN 关联
results = run_multi_target_demo(duration=100.0, n_targets=3)

# 使用最近邻关联
results = run_multi_target_demo(
    duration=100.0,
    n_targets=3,
    association_method="nn"
)

# 使用 JPDA 关联并自定义参数
results = run_multi_target_demo(
    duration=100.0,
    n_targets=3,
    association_method="jpda",
    association_params={
        "detection_probability": 0.9,
        "clutter_density": 1e-3
    }
)

# 使用 MHT 关联
results = run_multi_target_demo(
    duration=100.0,
    n_targets=3,
    association_method="mht"
)

# 使用自定义关联算法对�?from data_association import JPDAFilter
my_jpda = JPDAFilter(
    gating_threshold=9.21,
    detection_probability=0.95,
    clutter_density=1e-4
)
results = run_multi_target_demo(
    duration=100.0,
    n_targets=3,
    association=my_jpda
)
```

## 项目结构

```
BayesianFilter/
├── ground_truth/          # 真实场景生成模块
�?  ├── __init__.py
�?  ├── target.py          # 目标基类
�?  ├── motion_models.py   # 运动模型实现
�?  └── scenario_manager.py # 场景管理�?├── measurements/          # 观测模型模块
�?  ├── __init__.py
�?  ├── measurement.py     # 观测数据结构
�?  ├── noise_models.py    # 噪声模型
�?  ├── clutter_models.py  # 杂波模型
�?  └── measurement_simulator.py # 观测模拟�?├── filter/                # 滤波算法模块
�?  ├── __init__.py
�?  ├── base_filter.py     # 滤波器基�?�?  ├── kalman_filter.py   # 卡尔曼滤波器
�?  ├── extended_kalman_filter.py # 扩展卡尔曼滤波器
�?  ├── unscented_kalman_filter.py # 无迹卡尔曼滤波器
�?  ├── cubature_kalman_filter.py # 容积卡尔曼滤波器
�?  ├── particle_filter.py # 粒子滤波�?�?  ├── auxiliary_particle_filter.py # 辅助粒子滤波�?�?  ├── rao_blackwellized_particle_filter.py # RBPF
�?  └── unscented_particle_filter.py # UPF
├── data_association/      # 数据关联模块
�?  ├── __init__.py
�?  ├── base_association.py # 关联基类
�?  ├── nearest_neighbor.py # 最近邻关联
�?  ├── jpda.py            # JPDA
�?  ├── mht.py             # MHT
�?  └── phd_filter.py      # PHD滤波�?├── visualize/             # 可视化模�?�?  ├── __init__.py
�?  ├── metrics.py         # 性能指标
�?  ├── scenario_visualizer.py # 场景可视�?�?  └── performance_visualizer.py # 性能可视�?├── demo/                  # 演示脚本
�?  ├── __init__.py
�?  ├── demo_linear_tracking.py
�?  ├── demo_nonlinear_tracking.py
�?  ├── demo_multi_target.py
�?  └── demo_performance_comparison.py
├── utils/                 # 工具函数
�?  ├── __init__.py
�?  ├── math_utils.py
�?  └── coordinate.py
├── run_demo.py            # 主运行脚�?├── requirements.txt       # 依赖列表
└── README.md             # 项目说明
```

## 算法说明

### 卡尔曼滤波器 (KF)
适用于线性高斯系统，通过预测-更新框架进行状态估计�?
### 扩展卡尔曼滤波器 (EKF)
通过雅可比矩阵对非线性系统进行线性化，适用于弱非线性系统�?
### 无迹卡尔曼滤波器 (UKF)
使用Sigma点捕获非线性变换后的统计特性，无需计算雅可比矩阵�?
### 容积卡尔曼滤波器 (CKF)
使用三阶球面-径向容积规则，数值稳定性更好�?
### 粒子滤波�?(PF)
使用序贯蒙特卡洛方法，适用于任意非线性和非高斯系统�?
### PHD滤波�?使用随机有限集理论，能够自动处理目标数量变化�?
## 性能指标

### RMSE (均方根误�?
衡量估计位置与真实位置的偏差�?
### OSPA (最优子模式分配)
综合考虑定位误差和基数误差的多目标跟踪指标�?
### GOSPA (广义最优子模式分配)
将误差分解为定位误差、漏检误差和虚警误差�?
## 参考文�?
1. Bar-Shalom, Y., Li, X. R., & Kirubarajan, T. (2004). Estimation with applications to tracking and navigation.
2. Vo, B. N., & Ma, W. K. (2006). The Gaussian mixture probability hypothesis density filter.
3. Arulampalam, M. S., et al. (2002). A tutorial on particle filters for online nonlinear/non-Gaussian Bayesian tracking.

## 许可�?
MIT License
