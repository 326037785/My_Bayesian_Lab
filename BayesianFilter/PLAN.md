# Bayesian Filter 示例项目

## 项目概述
创建一个完整的贝叶斯滤波使用示例，包含从目标生成、观测模型、滤波算法到性能评估的全流程。

**状态: 已完成**

## 目录结构
```
BayesianFilter/
├── ground_truth/          # 真实场景生成
├── measurements/          # 观测模型
├── filter/                # 滤波算法
├── data_association/      # 数据关联
├── demo/                  # 测试样例
├── visualize/             # 可视化
├── utils/                 # 工具函数
└── requirements.txt       # 依赖
```

---

## Phase 1: 基础框架与真实场景生成
**目标**: 建立项目基础结构，实现目标运动模型

### 1.1 项目初始化
- [ ] 创建目录结构
- [ ] 创建 requirements.txt
- [ ] 创建 utils/ 工具模块（坐标转换、数学工具等）

### 1.2 ground_truth 模块
- [ ] 定义目标基类 `Target`
- [ ] 实现运动模型：
  - [ ] 匀速直线运动 (CV - Constant Velocity)
  - [ ] 匀加速运动 (CA - Constant Acceleration)
  - [ ] 协调转弯 (CT - Coordinated Turn)
  - [ ] 随机游走 (RW - Random Walk)
- [ ] 实现场景管理器 `ScenarioManager`
  - [ ] 目标出现/消失时间控制
  - [ ] 多目标轨迹生成
- [ ] 实现线性/非线性状态转移

---

## Phase 2: 观测模型与杂波生成
**目标**: 实现完整的观测生成流程

### 2.1 measurements 模块
- [ ] 定义观测基类 `Measurement`
- [ ] 实现观测模型：
  - [ ] 线性观测模型 (直接观测 + 噪声)
  - [ ] 非线性观测模型 (距离/角度观测)
- [ ] 实现噪声生成器：
  - [ ] 高斯噪声
  - [ ] 非高斯噪声选项
- [ ] 实现杂波生成器：
  - [ ] 泊松杂波模型
  - [ ] 均匀分布杂波
  - [ ] 杂波密度可配置

### 2.2 观测管理器
- [ ] `MeasurementSimulator` 类
- [ ] 整合真实目标观测 + 杂波
- [ ] 支持目标检测概率 `Pd`

---

## Phase 3: 高斯滤波器族
**目标**: 实现基于高斯假设的滤波算法

### 3.1 滤波器基类
- [ ] `BaseFilter` 抽象类
- [ ] 标准接口：`predict()`, `update()`, `initialize()`

### 3.2 线性滤波器
- [ ] Kalman Filter (KF)
  - [ ] 状态预测
  - [ ] 量测更新
  - [ ] 协方差传播

### 3.3 非线性滤波器
- [ ] Extended Kalman Filter (EKF)
  - [ ] 雅可比矩阵计算
  - [ ] 线性化误差处理
- [ ] Unscented Kalman Filter (UKF)
  - [ ] Sigma点生成
  - [ ] 无迹变换
- [ ] Cubature Kalman Filter (CKF)
  - [ ] 容积规则
  - [ ] 三阶球面-径向规则

---

## Phase 4: 粒子滤波器族
**目标**: 实现基于序贯蒙特卡洛的滤波算法

### 4.1 SMC基础
- [ ] 粒子表示 `Particle`
- [ ] 重采样算法：
  - [ ] 多项式重采样
  - [ ] 系统重采样
  - [ ] 残差重采样

### 4.2 粒子滤波器变体
- [ ] Standard Particle Filter (PF)
- [ ] Auxiliary Particle Filter (APF)
- [ ] Rao-Blackwellized Particle Filter (RBPF)
- [ ] Unscented Particle Filter (UPF)

### 4.3 自适应机制
- [ ] 粒子数自适应
- [ ] 建议分布选择

---

## Phase 5: 数据关联算法
**目标**: 实现多目标跟踪中的数据关联

### 5.1 关联基类
- [ ] `DataAssociation` 抽象类
- [ ] 关联矩阵计算

### 5.2 关联算法
- [ ] Nearest Neighbor (NN)
- [ ] K-Nearest Neighbor (KNN)
- [ ] Joint Probabilistic Data Association (JPDA)
- [ ] Multiple Hypothesis Tracking (MHT)

### 5.3 概率假设密度
- [ ] PHD滤波器
  - [ ] 高斯混合PHD (GM-PHD)
  - [ ] 粒子PHD (PHD-PF)

---

## Phase 6: 性能评估与可视化
**目标**: 实现评估指标和可视化工具

### 6.1 性能指标
- [ ] RMSE (Root Mean Square Error)
- [ ] OSPA (Optimal Subpattern Assignment)
- [ ] GOSPA (Generalized OSPA)
- [ ] 航迹连续性指标

### 6.2 可视化模块
- [ ] 场景可视化 `ScenarioVisualizer`
  - [ ] 真实轨迹绘制
  - [ ] 观测点绘制
  - [ ] 滤波结果叠加
- [ ] 性能可视化 `PerformanceVisualizer`
  - [ ] RMSE时序图
  - [ ] OSPA/GOSPA对比图
  - [ ] 误差分布图

---

## Phase 7: 集成测试与演示
**目标**: 创建完整的演示样例

### 7.1 基础演示
- [ ] 单目标线性跟踪演示
- [ ] 单目标非线性跟踪演示
- [ ] 多目标跟踪演示

### 7.2 高级演示
- [ ] 目标出现/消失场景
- [ ] 高杂波环境测试
- [ ] 滤波器性能对比

### 7.3 综合演示脚本
- [ ] `demo_linear_tracking.py`
- [ ] `demo_nonlinear_tracking.py`
- [ ] `demo_multi_target.py`
- [ ] `demo_performance_comparison.py`

---

## 技术规范

### 环境要求
- Python 3.10
- NumPy >= 1.21
- Matplotlib >= 3.5
- SciPy >= 1.7

### 代码规范
- 使用类型注解
- 每个模块包含 `__init__.py`
- 关键类和方法添加文档字符串
- 保持模块间低耦合

### 测试策略
- 单元测试覆盖核心算法
- 集成测试验证端到端流程
- 可视化验证结果合理性

---

## 实施顺序
1. **Phase 1-2**: 基础框架（优先级最高）
2. **Phase 3-4**: 滤波算法核心
3. **Phase 5**: 数据关联（依赖Phase 3-4）
4. **Phase 6**: 评估与可视化
5. **Phase 7**: 演示与集成

每个Phase完成后进行测试验证，确保功能正确后再进入下一阶段。
