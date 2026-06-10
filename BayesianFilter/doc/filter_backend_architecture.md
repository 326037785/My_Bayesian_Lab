# FilterBackend 架构设计

## 当前问题

```
┌─────────────┐     ┌─────────────┐     ┌─────────────┐
│  KF Filter  │     │  EKF Filter │     │  UKF Filter │
│  - state    │     │  - state    │     │  - state    │
│  - cov      │     │  - cov      │     │  - cov      │
│  - H, R     │     │  - f, h, F_j│     │  - f, h     │
└──────┬──────┘     └──────┬──────┘     └──────┬──────┘
       │                   │                   │
       │ 直接访问内部属性   │                   │
       ▼                   ▼                   ▼
┌─────────────────────────────────────────────────────┐
│              Demo 代码 (硬编码)                      │
│  z_pred = kf.H @ kf.state                          │
│  S = kf.H @ kf.covariance @ kf.H.T + kf.R         │
└─────────────────────────┬───────────────────────────┘
                          │
                          ▼
┌─────────────────────────────────────────────────────┐
│              数据关联算法                             │
│  - KNN / JPDA / MHT                                │
│  - 需要 predicted_measurements                      │
│  - 需要 innovation_covariances                      │
└─────────────────────────────────────────────────────┘
```

**问题：**
1. 数据关联算法直接依赖滤波器内部结构
2. 预测观测/新息协方差计算散落在demo中
3. 切换滤波器需要修改多处代码
4. 门限计算在关联算法中重复实现

## 新架构：FilterBackend 抽象层

```
┌─────────────────────────────────────────────────────────────┐
│                    FilterBackend 接口                        │
│  - get_predicted_measurement() → z_pred                     │
│  - get_innovation_covariance() → S                          │
│  - compute_likelihood(z) → p(z|x)                          │
│  - compute_mahalanobis(z) → d                               │
│  - gating_test(z, threshold) → bool                         │
└─────────────────────────┬───────────────────────────────────┘
                          │
         ┌────────────────┼────────────────┐
         │                │                │
         ▼                ▼                ▼
┌─────────────┐   ┌─────────────┐   ┌─────────────┐
│  KFBackend  │   │ EKFBackend  │   │ UKFBackend  │
│  (线性)     │   │ (非线性雅可比)│   │ (sigma点)   │
└─────────────┘   └─────────────┘   └─────────────┘
         │                │                │
         └────────────────┼────────────────┘
                          │
                          ▼
┌─────────────────────────────────────────────────────────────┐
│              MultiTargetFilterManager                        │
│  - add_target(id, state)                                    │
│  - predict_all(dt)                                          │
│  - get_predicted_states() → List[PredictedState]           │
│  - update_target(id, z)                                     │
└─────────────────────────┬───────────────────────────────────┘
                          │
                          │ PredictedState {state, cov, z_pred, S}
                          ▼
┌─────────────────────────────────────────────────────────────┐
│              数据关联算法 (只依赖接口)                        │
│  - KNN / JPDA / MHT                                        │
│  - 接收 PredictedState 列表                                 │
│  - 不关心具体滤波器实现                                     │
└─────────────────────────────────────────────────────────────┘
```

## 使用示例

```python
from filter import KFBackend, UKFBackend, MultiTargetFilterManager
from data_association import KNearestNeighborAssociation

# 1. 创建后端（可轻松切换）
backend_factory = lambda: KFBackend(state_dim=4, meas_dim=2)
# 或: backend_factory = lambda: UKFBackend(state_dim=4, meas_dim=2)

# 2. 创建管理器
manager = MultiTargetFilterManager(backend_factory)

# 3. 添加目标
manager.add_target(0, initial_state=np.array([100, 10, 50, 5]))
manager.add_target(1, initial_state=np.array([200, -5, 150, 10]))

# 4. 预测
manager.predict_all(dt=1.0)

# 5. 获取预测状态（统一接口）
predicted_states = manager.get_predicted_states()
# predicted_states[0].predicted_meas  ← 预测观测
# predicted_states[0].innovation_cov  ← 新息协方差

# 6. 数据关联（不依赖具体滤波器）
association = KNearestNeighborAssociation()
result = association.associate(measurements, predicted_measurements, ...)

# 7. 更新
manager.update_target(target_id, measurement)
```

## 优势

| 方面 | 旧架构 | 新架构 |
|------|--------|--------|
| 耦合度 | 高（直接访问内部属性） | 低（只依赖接口） |
| 可切换性 | 需要修改多处代码 | 只需改 backend_factory |
| 代码复用 | 散落在demo中 | 统一在后端中 |
| 可测试性 | 难以mock | 接口易于测试 |
| 扩展性 | 需要修改关联算法 | 只需实现新后端 |
