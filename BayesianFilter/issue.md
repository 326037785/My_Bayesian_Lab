# BayesianFilter 项目代码审查 Issue

> **审查日期**: 2026-06-09  
> **审查范围**: `filter/`, `data_association/`, `demo/` 全部实现  
> **评价基准**: 硬核数学推理正确性 + 对标 GitHub 高星项目 (filterpy ~3.3k★, Stone-Soup ~570★, particles ~390★)

---

## 目录

1. [严重问题 (Critical)](#1-严重问题-critical)
2. [重要问题 (Major)](#2-重要问题-major)
3. [中等问题 (Moderate)](#3-中等问题-moderate)
4. [演示缺陷 (Demo Issues)](#4-演示缺陷-demo-issues)
5. [对标高星项目差距](#5-对标高星项目差距)
6. [总结与优先级建议](#6-总结与优先级建议)

---

## 1. 严重问题 (Critical)

### 1.1 APF 辅助粒子滤波器 — 数学实现错误

**文件**: `filter/auxiliary_particle_filter.py`

**问题描述**: APF 的核心机制是通过引入辅助变量来减少权重方差——在重采样阶段使用当前观测信息选择"有前途"的粒子。但当前实现在 `predict()` 行 91-101 中，辅助权重仅设为先验权重：

```python
# 行 101: 辅助权重 = 先验权重，完全未使用观测信息
auxiliary_weights[i] = particle.weight
```

然后在 `update()` 行 159-164 中，辅助似然被设为等于当前似然：

```python
# 行 164: 辅助似然 = 当前似然 → 最终权重恒为 1.0
auxiliary_likelihood = current_likelihood
```

这导致 `particle.weight = current_likelihood / auxiliary_likelihood = 1.0`，APF 完全退化为标准 PF。

**正确的数学公式**:  
APF 应先计算每个粒子预测均值的似然 `p(z_k | μ_k^{(i)})`，然后辅助权重应为 `w_{k-1}^{(i)} × p(z_k | μ_k^{(i)})`，最终权重为 `p(z_k | x_k^{(i)}) / p(z_k | μ_k^{(i)})`。

**对标**: filterpy 没有 APF；particles 库通过 `GuidedPF` 正确实现了完全适应的辅助粒子滤波器——其 proposal 分布使用了当前观测信息。

---

### 1.2 JPDA — 关联概率计算为占位代码

**文件**: `data_association/jpda.py`, 行 130-162 (`_compute_association_probabilities`)

**问题描述**: JPDA 的核心是联合关联概率计算，但当前实现将所有有效关联对设为相同概率：

```python
# 行 153: 所有有效关联对的概率恒为 1.0！
prob = 1.0  # 这里应该根据实际距离计算
```

正确的 JPDA 需要：
1. 枚举所有联合关联假设 θ
2. 计算每个假设的后验概率：`P{θ|Z^k} ∝ λ^φ · P_D^δ · (1-P_D)^{n-δ} · Π N(z; Hx, S)`
3. 计算边际关联概率：`β_{jt} = Σ_{θ:θ_{jt}=1} P{θ|Z^k}`

虽然 `compute_joint_probabilities()` 方法（行 164-205）提供了完整实现框架，但 **`associate()` 从不调用它**，且其内部的 `_compute_hypothesis_probability` 也简化过度。

**对标**: Stone-Soup 实现了完整的 JPDA，包含正确的假设概率计算、P_D 和杂波密度参数的真实使用。其 `JPDAMixtureReducer` 正确计算了 β 系数。

---

### 1.3 PHD 滤波器 `update()` 中的变量引用错误

**文件**: `data_association/phd_filter.py`, 行 274

```python
# 行 274: self.covariance 在 PHDFilter 中未定义！
K = self.covariance @ self.H.T @ S_inv
```

**问题描述**: `PHDFilter` **不继承** `BaseFilter`，`self.covariance` 从未在 `__init__` 中定义。正确的写法应该是 `comp.covariance`。这会导致 **`AttributeError` 运行时崩溃**。

**对标**: `danstowell/gmphd` (72★) 和 `djape24394/gmphd_filter` (28★) 都使用了正确的 `comp.cov` 引用方式。

---

### 1.4 UPF 无迹粒子滤波器 — 重要性权重计算错误

**文件**: `filter/unscented_particle_filter.py`, 行 263-297

**问题描述**: UPF 使用 UKF 作为 proposal 分布 `q(x_k | x_{k-1}, z_k)`，但权重更新（行 297）仅设为：

```python
particle.weight = likelihood  # 缺少 proposal 密度项
```

正确的重要性权重应为：  
`w_k ∝ p(z_k | x_k) × p(x_k | x_{k-1}) / q(x_k | x_{k-1}, z_k)`

其中 proposal 密度 `q(·)` 是 UKF 后验的多元高斯 PDF。当前实现缺少分母中的 proposal 密度，导致权重系统性偏差。

此外，`SquareRootUnscentedParticleFilter`（行 384-429）的 `_create_ukf()` 方法是占位实现，注释写明"这里仍然使用标准UKF"。

**对标**: particles 库的 `GuidedPF` 正确计算了 `log_target - log_proposal` 的完整权重比。

---

## 2. 重要问题 (Major)

### 2.1 Q 矩阵代码重复 8 次

**影响文件**: `kalman_filter.py:83-114`, `extended_kalman_filter.py:168-197`, `unscented_kalman_filter.py:122-150`, `cubature_kalman_filter.py:108-136`, `particle_filter.py:149-177`, `auxiliary_particle_filter.py` (继承), `rao_blackwellized_particle_filter.py:157-171`, `unscented_particle_filter.py:159-191`

完全相同的 `_get_process_noise_matrix(dt)` 实现被复制了 8 次。任何对 Q 矩阵公式的修正需要修改 8 个文件。应该提取到 `base_filter.py` 或独立的工具模块。

**对标**: filterpy 使用 `Q_discrete_white_noise()` 单一函数生成 Q 矩阵，所有滤波器共享。

---

### 2.2 缺少平方根形式的 UKF/CKF

**文件**: `filter/unscented_kalman_filter.py`, `filter/cubature_kalman_filter.py`

UKF 和 CKF 使用标准协方差传播 `P = P - K·S·K'`，而非 Joseph 形式 `(I-KH)P(I-KH)' + KRK'`。更重要的是，没有实现平方根 UKF/CKF（使用 Cholesky 因子传播），在病态条件下会丢失正定性。

**对标**: filterpy 提供了 `SquareRootKalmanFilter`；particles 库在内部使用 Cholesky 因子传播。

---

### 2.3 RBPF 默认函数不可用

**文件**: `filter/rao_blackwellized_particle_filter.py`

- `_default_nonlinear_transition` (行 120-124): 恒等函数，不产生任何状态变化
- `_default_measurement` (行 148-151): 仅观测线性部分，忽略了非线性部分对观测的贡献
- `_update_state_estimate` (行 255): 协方差设为硬编码的 `np.eye(self.state_dim) * 10.0`

RBPF **没有自定义模型函数就无法正常工作**，这违反了"开箱即用"的原则。

---

### 2.4 缺少 IMM (Interacting Multiple Model)

项目中完全没有 IMM 滤波器的实现。对于机动目标跟踪，IMM 是工业标准方案。Stone-Soup 和 filterpy 都提供了 IMM 实现。

---

### 2.5 MHT 缺乏 N-scan 剪枝和正确的似然比得分

**文件**: `data_association/mht.py`

- 假设管理仅按得分的 top-K 剪枝 (行 218-219)，缺乏标准的 N-scan 剪枝策略
- 假设得分 (行 186-204) 使用 `exp(-0.5·d)` 而非正确的对数似然比 `LLR = log(p(z|H1) / p(z|H0))`
- `GlobalHypothesisMHT._extend_hypothesis` 每步生成 O(n_meas × n_targets) 个新假设，无约束增长

**对标**: Stone-Soup 的 MHT 使用 Murty's k-best 分配和 N-scan 剪枝。`erikliland/pyMHT` (105★) 使用正确的 track score（对数似然比）。

---

### 2.6 默认似然函数中重复矩阵求逆

**文件**: `filter/particle_filter.py`, 行 115-143

`_default_likelihood` 对**每个粒子每次更新**调用 `np.linalg.inv(R)` 和 `np.linalg.det(R)`，但 R 在滤波器生命周期内不变。在大粒子数场景下浪费严重——应当预计算并缓存。

---

### 2.7 缺少自适应噪声估计

所有滤波器使用固定的过程噪声和测量噪声参数。实践中噪声特性可能未知或时变。缺少：
- Sage-Husa 自适应卡尔曼滤波
- 基于新息的自适应噪声协方差估计
- 多模型自适应估计 (MMAE)

---

## 3. 中等问题 (Moderate)

### 3.1 状态向量约定硬编码

整个代码库假定状态为 `[x, vx, y, vy]`（4维）或 `[x, vx, ax, y, vy, ay]`（6维）。位置提取全部使用硬编码索引 `[0, 2]`。这使得以下场景无法处理：
- 3D 跟踪 `[x, vx, y, vy, z, vz]`
- 不同的变量排序（如 `[x, y, vx, vy]`）

**对标**: filterpy 将状态向量结构留给用户定义，不硬编码索引。Stone-Soup 使用 `State` 对象封装状态向量及其映射关系。

---

### 3.2 缺少分层重采样 (Stratified Resampling)

**文件**: `filter/particle_filter.py`

实现了 3 种重采样方法（multinomial, systematic, residual），但缺少 stratified resampling，它通常提供最低的方差。particles 库全部实现了 5 种方法。

---

### 3.3 UKF/UPF 缺少 Spherical Simplex Sigma 点

标准 UKF 使用 2n+1 个 sigma 点。spherical simplex 方法只需要 n+2 个点，可显著减少计算量。在高维状态空间中这是重要的优化。

---

### 3.4 缺少数值雅可比计算

**文件**: `filter/extended_kalman_filter.py`

当用户提供自定义非线性函数但未提供雅可比时，EKF 没有自动数值微分能力（有限差分或复数步长微分）。

**对标**: filterpy 的 EKF 支持通过 `compute_jacobian` 进行数值雅可比计算。

---

### 3.5 PF 缺少 Regularization/MCMC Move Step

标准粒子滤波器在重采样后容易出现样本贫化（sampling impoverishment）。常见的缓解方法：
- 正则化粒子滤波 (RPF): 从连续核密度中重采样
- MCMC move step: 重采样后用 MCMC 移动粒子

当前实现两者皆无。

**对标**: particles 库的 `SMC` 类支持 `move` 步骤，可以进行 MCMC 扰动。

---

### 3.6 缺少信息滤波器形式 (Information Filter)

对于某些传感器配置（如极大测量维度或低过程噪声），信息滤波器形式在数值上更稳定。项目中只有标准的协方差形式。

---

### 3.7 PHD 滤波器缺少 Cardinized PHD (CPHD)

GM-PHD 的基数估计方差较大。CPHD 通过传播完整的基数分布来改善这一点。当前仅有 GM-PHD 的基础形式。

---

### 3.8 代码无向量化优化

UKF 的 `_unscented_transform` 和 CKF 的 `_cubature_transform` 使用显式 Python `for` 循环而非 NumPy 批量操作，对于大量 sigma/cubature 点效率低下。

---

### 3.9 缺少单元测试

项目中**完全没有测试文件**。对于数学库而言，缺少：

- 协方差正定性回归测试
- 已知场景下的 RMSE 基准测试
- 滤波器一致性测试（NEES/NIS chi-squared test）

**对标**: filterpy 有完整的测试套件；Stone-Soup 拥有 CI/CD 和覆盖率报告。

---

### 3.10 GlobalNearestNeighbor 放在错误位置

**文件**: `data_association/base_association.py`, 行 166-238

`GlobalNearestNeighbor` 类逻辑上属于 `nearest_neighbor.py`，但被放在了 `base_association.py` 中——组织不一致。

---

## 4. 演示缺陷 (Demo Issues)

### 4.1 所有 Demo 使用简化测量选择而非真正的数据关联

**文件**: `demo/demo_linear_tracking.py`, `demo/demo_nonlinear_tracking.py`

线性/非线性 demo 使用欧氏距离最近邻选测量（行 142-148），并使用任意门限值（线性用 100，非线性用 5.0）：

```python
# demo_linear_tracking.py:147
if distances[nearest_idx] < 100:  # 简单门限 — 无统计意义
    kf.update(nearest_meas)
```

没有使用马氏距离门限，没有使用项目中已实现的任何数据关联算法。

---

### 4.2 非线性 Demo 在极坐标下使用欧氏距离

**文件**: `demo/demo_nonlinear_tracking.py`, 行 178

```python
distances = np.linalg.norm(meas_array - z_pred, axis=1)
```

距离和角度具有不同的单位和尺度，直接用欧氏距离在极坐标空间中比较是无意义的——一个角度差 0.1 rad ≈ 距离差 100m（取决于实际距离）。

---

### 4.3 多目标 Demo 脆弱的目标索引映射

**文件**: `demo/demo_multi_target.py`, 行 154

```python
target_id = sorted(filters.keys())[target_idx]
```

关联器返回的 `target_idx` 是基于 `predicted_measurements` 数组的行索引，但 `sorted(filters.keys())` 的排序可能与 `predicted_measurements` 的构建顺序不一致。这会导致关联结果错误地映射到滤波器。

---

### 4.4 多目标 Demo 使用 Ground Truth 初始化新目标

**文件**: `demo/demo_multi_target.py`, 行 168-169

```python
initial_state = target.get_state_at_time(t)
```

新出现目标的滤波器使用**真实状态**初始化——在实际应用中这是不可用的。应从未关联的测量中初始化。

---

### 4.5 性能对比 Demo 在线性场景下测试非线性滤波器

**文件**: `demo/demo_performance_comparison.py`

在纯线性场景下对比 KF vs EKF vs UKF vs CKF——此时 KF 总是最优的（MSE 意义上），对比没有信息量。应该在非线性场景（如极坐标观测）下对比。同时缺少粒子滤波器参与对比。

---

### 4.6 EKF/UKF/CKF 在线性 Demo 中错误导入极坐标函数

**文件**: `demo/demo_performance_comparison.py`, 行 166-177

线性观测场景下却导入了 `polar_measurement_function` 等函数用于计算预测测量——线性观测应直接使用 `kf.H @ kf.state`。非 KF 滤波器在行 170-175 中的分支逻辑实际上用的是线性场景，不应导入极坐标函数。

---

## 5. 对标高星项目差距

### 5.1 vs filterpy (⭐ ~3,300)

| 方面 | 本项目 | filterpy |
|------|--------|----------|
| KF 实现正确性 | ✅ Joseph 形式正确 | ✅ Joseph 形式 |
| UKF | ❌ 无平方根形式 | ✅ 有 SquareRootUKF |
| EKF | ❌ 无数值雅可比 | ✅ 自动数值微分 |
| PF | ⚠️ 无stratified/regularization | ⚠️ 同样缺失部分高级特性 |
| APF | ❌ 数学实现错误 | — 未实现 |
| 文档 | ❌ 仅 README | ✅ 完整 RTD 文档 + 配套书籍 |
| 测试 | ❌ 无 | ✅ 有测试套件 |
| 代码复用 | ❌ Q 矩阵重复 8 次 | ✅ 共享工具函数 |

### 5.2 vs Stone-Soup (⭐ ~570)

| 方面 | 本项目 | Stone-Soup |
|------|--------|------------|
| JPDA | ❌ 概率计算为占位 | ✅ 完整正确实现 |
| MHT | ⚠️ 无 N-scan, 错误得分 | ✅ Murty's k-best + 正确 LLR |
| PHD | ⚠️ 有一个 bug (self.covariance) | ✅ 成熟实现 |
| 架构 | ⚠️ 单一继承树 | ✅ 高度模块化组件系统 |
| 综合算法 | ❌ 缺 PMBM, PMB, GLMB, LMB 等 | ✅ 覆盖现代 RFS 全家族 |
| Track 管理 | ❌ 无 | ✅ Initiator/Deleter/Updater 组件 |
| 传感器模型 | ⚠️ 仅线性/极坐标 | ✅ 雷达/红外/声呐等多种传感器 |

### 5.3 vs particles (⭐ ~390)

| 方面 | 本项目 | particles |
|------|--------|------------|
| 重采样方法 | 3 种 | 5 种 (含 stratified, SSP) |
| SMC 理论 | ⚠️ 基础实现 | ✅ SQMC, SMC², PMCMC |
| MCMC move | ❌ 无 | ✅ 支持 |
| 学术质量 | ❌ 多处数学错误 | ✅ 配套 Springer 教材 |

---

## 6. 总结与优先级建议

### 必须立即修复 (P0)

1. **APF**: 重写 predict/update 以正确使用当前观测计算辅助权重 (1.1)
2. **JPDA**: 实现真正的 `_compute_association_probabilities`——当前为占位代码 (1.2)
3. **PHD**: 修复 `self.covariance` → `comp.covariance` (1.3)
4. **UPF**: 修复重要性权重，加入 proposal 密度 (1.4)

### 应尽快修复 (P1)

5. **Q 矩阵去重**: 提取到共享模块，8 处重复 → 1 处 (2.1)
6. **RBPF**: 提供可用的默认模型或明确要求自定义函数 (2.3)
7. **添加 SQUK/SQCKF**: 数值稳定性的平方根形式 (2.2)
8. **MHT**: 实现正确的 LLR 得分和 N-scan 剪枝 (2.5)
9. **Demo 修复**: 使用正确的数据关联 + 马氏距离门限 (4.1-4.4)

### 应规划实现 (P2)

10. **IMM 滤波器**: 机动目标跟踪的关键缺失 (2.4)
11. **自适应噪声估计**: Sage-Husa 或类似方案 (2.7)
12. **Stratified 重采样**: 完成 PF 重采样方法集 (3.2)
13. **单元测试**: 至少覆盖 KF/EKF 的一致性 (NEES/NIS) 测试 (3.9)
14. **CPHD**: PHD 的基数增强版本 (3.7)

### 长期改进 (P3)

15. 解耦状态向量约定 (3.1)
16. 向量化优化 UKF/CKF/PF 循环 (3.8)
17. 添加信息滤波器形式 (3.6)
18. 正则化/MCMC particle filter (3.5)
19. API 文档 (Sphinx/RTD)
20. 跟踪管理系统 (航迹起始/确认/删除)

---

## 参考项目

| 项目 | Stars | 强项 |
|------|-------|------|
| [rlabbe/filterpy](https://github.com/rlabbe/filterpy) | ~3,300 | 教育清晰度 + 平方根滤波 |
| [dstl/Stone-Soup](https://github.com/dstl/Stone-Soup) | ~570 | 综合框架 + 数据关联 + RFS |
| [nchopin/particles](https://github.com/nchopin/particles) | ~390 | SMC 理论深度 + 5种重采样 |
| [erikliland/pyMHT](https://github.com/erikliland/pyMHT) | ~105 | MHT 正确实现 |
| [danstowell/gmphd](https://github.com/danstowell/gmphd) | ~72 | GM-PHD 正确参考实现 |
| [apennisi/jpdaf_tracking](https://github.com/apennisi/jpdaf_tracking) | ~240 | JPDAF 专项实现 |

---

*由 Claude Code 基于硬核数学推理和 GitHub 高星项目对标审查生成*
