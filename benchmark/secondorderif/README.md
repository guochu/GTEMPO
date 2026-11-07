# secondorderif — Weyl 二阶 GTEMPO（完整版）

基于 GTEMPO 框架（Z2 GMPS）实现 **Weyl 二阶算法**（`weyl_tti` 方案），结构
与 `~/Documents/Missile/U1GTEMPO2` 对齐，但全部构建在 GTEMPO 已有机制之上。
对**任意 `ConstImpurityHamiltonian`（任意带数）与一般有界谱浴**成立。

```
src/                       SecondOrderWeyl 模块（完整版算法）
  SecondOrderWeyl.jl         模块入口
  contours.jl                WeylRealTime / WeylImagTime / weyl_lattice /
                             ChronoKeldyshOrdering（实时时序顺序存储）
  kernels.jl                 单元级浴核：虚时 QuAPI Δτ（一般浴）；实时
                             QuAPI Δt 谱路径统一（离散谱精确极点装配；δ 闭式特例已于 2026-11-07 移除）
  tensors.jl                 单元核 → Z2 site tensors（逐变量分级精确 SVD）
  system.jl                  local_weight / system_factor / WeylInsertion
                             （实时 K = ChronoKeldysh 时序直铺，零 permute!）
  influence.jl               IF（复用 GTEMPO hybriddynamics）
  contraction.jl             积分（复用 GTEMPO integrate）
  solver.jl                  Problem / Solution / solve
  measurement.jl             occupation（总/逐带）/ greater / lesser /
                             raw_correlation（等时 0⁺）/ matsubara
  diagnostics.jl             fit_order / convergence_report / maxbond
validation/                共享校验（ED 参考 + 核双路径 + 实时快路径 +
                           等时两点/逐带占据回归 + 连续浴 U1 对齐 +
                           U=0 因子化：K₂ vs K₁∧K₁ 楔积 / 端到端 mixed+rt +
                           IF 带不变性探针 + mixed-U1 对齐报告）
toulouse/                  无相互作用 Toulouse 模型验证（单带 + 单模浴）
singlemodebath/            两带 Anderson 模型（U≠0，每带单模浴）验证
plot_scaling.py            error-scaling 图（matplotlib）
HISTORY.md                 2026-11-04/06 十轮修复史（部件级对齐、跨框架
                           对比教训、假象补偿清算、每带 XOR 补偿）——重要档案
```

## 运行

```bash
julia --project=. benchmark/secondorderif/toulouse/convergence.jl .            # ①②
julia --project=. benchmark/secondorderif/singlemodebath/convergence.jl .      # ③④
julia --project=. benchmark/secondorderif/toulouse/convergence_mixed.jl .      # 单带混合
julia --project=. benchmark/secondorderif/singlemodebath/convergence_mixed.jl . # 两带混合
python3 benchmark/secondorderif/plot_scaling.py                                # 出图
# 共享校验：
julia --project=. benchmark/secondorderif/validation/kernel_paths_check.jl .   # 核双路径 + 双δ浴端到端
julia --project=. benchmark/secondorderif/validation/rt_fastpath_check.jl .    # 实时快路径回归
julia --project=. benchmark/secondorderif/validation/twopoint_band_check.jl .  # 等时 0⁺ + 逐带占据
# U=0 因子化（本征锚，M≥2 装配判决）：
julia --project=. benchmark/secondorderif/validation/k_factor_mixed.jl .       # K₂ vs K₁∧K₁ 逐支撑楔积（0/400）
julia --project=. benchmark/secondorderif/validation/u0_factor_mixed.jl .      # 混合端到端（Z₂/Z₁² + vs 2×单带）
julia --project=. benchmark/secondorderif/validation/u0_factor_rt.jl .         # 实时端到端（同款）
julia --project=. benchmark/secondorderif/validation/if_band_invariance.jl .   # IF₁ 带嵌入探针（比值判据）
julia --project=. benchmark/secondorderif/validation/mixed_twopoint_check.jl . # 混合两点：M=1 全网格 + M=2 等时 + 防守确认
julia --project=. benchmark/secondorderif/validation/twopoint_m2_frontier.jl . # M≥2 跨格奇插入边界判决（已知限制记录）
```

## 已验证结果（ED 参考，2026-11-06 第十轮后）

| 模型 | 观测量 | order 1 斜率 | order 2 斜率 | 最细误差 |
|------|--------|--------------|--------------|----------|
| Toulouse | 虚时 ⟨n⟩ | 1.380 | **2.002** | 5.9e-6 (N=40) |
| Toulouse | 实时 ⟨n(0.4)⟩ | 0.822 | **1.996** | 1.9e-6 (N=8) |
| Toulouse | 混合 ⟨n(0.4)⟩ | 0.932 | **2.005** | 3.3e-5 (N=16) |
| Anderson 双带 | 虚时 ⟨n₁+n₂⟩ | 1.607 | **2.004** | 7.5e-5 (N=16) |
| Anderson 双带 | 实时 ⟨n₁+n₂(0.4)⟩ | 0.987 | **2.001** | 8.3e-8 (N=32) |
| Anderson 双带 | 混合 ⟨n₁+n₂(0.4)⟩ | 0.944 | **2.004** | 6.4e-5 (N=16) |

Anderson 实时/混合的 U≠0 二阶收敛干净（第十轮修复前 rt 斜率 1.927、
mixed 0.641——均带 M≥2 一阶工件，已闭合）。逐带 ⟨n_b⟩、等时两点 0⁺
（⟨dd†⟩ = 1−⟨n⟩）、occupation(j=0) 折叠均 ED 验证通过
（`validation/twopoint_band_check.jl` 四节全绿，④ ≤1.06e-5）。

**U=0 因子化（本征锚，全部精确）**：双带 K₂ = K₁∧K₁（逐支撑楔积 0/400）、
Z₂/Z₁² = 1.0、⟨n₁+n₂⟩ = 2×单带（机器精度）——mixed 与 rt 两路径
（`k_factor_mixed.jl` / `u0_factor_mixed.jl` / `u0_factor_rt.jl`）。

**U1GTEMPO2 端到端对齐（两条独立实现管线）**：
- 连续谱浴两带 Anderson（τ + rt，order 1/2 × 3 网格 = 12 点）：
  **全部机器精度级**（虚时 ≤7.8e-12、实时 ≤1.4e-11，
  `validation/continuous_bath_alignment.md`）；
- 混合轮廓 Anderson 双带（order 1/2 × N=4/8 = 4 点）：全部逐位一致
  （`validation/mixed_u1_alignment.md`）。

## 一般浴支持

- 虚时核：QuAPI `Δτ`，任意有界谱（含 `DiscreteSpectrum` 多 δ 谱）；
- 实时核：全部浴型统一走 QuAPI `Δt`（离散谱精确极点装配；一般谱 quadgk；
  `[1:N,1:N]` 直接截取）——双路径在 Dirac-δ 上互验 ≤4e-17；双 δ（两模）
  浴端到端 vs ED 通过（`validation/kernel_paths_check.jl`）；
- `α_weight(spec::DiscreteSpectrum) = Σ αₙ`（多 δ 谱的总耦合）。

## 连续谱浴 + U1GTEMPO2 对齐（2026-11-05 建；2026-11-06 第十轮复测）

连续带谱 A(ε) = 0.25·cos²(π(ε−0.2)/4)（∫A = 0.5）上的两带 Anderson：
`validation/continuous_bath_{mine,u1}.jl` 两条独立实现管线端到端逐点互比
（逐带 ⟨n₁⟩ ≡ U1 逐自旋 ⟨n↑⟩，τ 与实时、order 1/2、各 N 网格）：

- **虚时：逐点偏差 ≤ 7.9e-12（机器精度级）**；
- **实时：逐点偏差 ≤ 1.4e-11（机器精度级）**——第十轮 rt 每带 XOR 修复前
  为 4.8e-5→6.6e-6 随 N 收敛（当时归因 U1 指数拟合容差，实为 rt 工件）。

报告：`validation/continuous_bath_alignment.md`；逐点数据
`validation/continuous_bath_{u1,mine}.csv`。

## 实时 K 的 ChronoKeldysh 直铺（性能）

实时 K 的单元张量按时序链本征顺序存储（`ChronoKeldyshOrdering`：块序 =
[(+,1)…(+,N),(−,N)…(−,1)]，块内 (a,ā) 奇偶对齐），MPS = 时序链本身，
**零 permute!**。旧实现重排到顶点交错序需 swap 门链，M≥2、N≥4 中间键
组合爆炸（N=4 >9 min；SVD 截断又会丢真秩给出错误值）。ChronoKeldysh
下每格 2M 个位置、(a,ā) 奇偶对齐 ⟹ AC 积分器位置对恰为共轭对且不跨格，
`_ac_integrate` 直接消费；N=8 双带实时 K+观测 1.3 s（约 400×）。
各网格无关——IF/观测量的构造全部经 `lattice[index]` 驱动。

## 完整版移植要点（2026-11-05 第八轮）

1. 模块从 `singlemodebath/src`（单模特例）提升为 `secondorderif/src`
   （`SecondOrderWeyl`），路径/引用全部重定位；旧特例版本已删除。
2. 一般谱实时核（QuAPI 谱路径）+ Dirac-δ 精确快路径双轨。
3. `ChronoKeldyshOrdering` 时序直铺（见上）。
4. 已实证不变式（完整证据链见 `HISTORY.md`）：IF = e^{+Δ} 精确指数、
   U=0 因子化（Z₂=Z₁²、K₂=K₁∧K₁ 逐支撑）、配对结合律、等时 0⁺、
   逐带占据二阶收敛 2.002。

## 第十轮：M≥2 混合/实时工件闭合（2026-11-06）

**症状**：双带 mixed order 2 斜率 0.641（O(h) 工件 c₁≈0.033）、rt 两带
1.927；有限 N vs U1 残差 ~6e-4；U=0 双带亦有一阶衰减分量。

**判决链**（全程本征锚）：IF₁ ≡ IF_single × 全局标量（比值 spread 3.3e-15，
IF 无罪——第九轮"结构性偏移"系探针绝对容差误读）→ U=0 因子化判决
Z₂/Z₁² ≈ 0.996（K 有罪）→ 逐支撑楔积 96/400 纯反号且全部 deg₁≥2∧deg₂≥2
⟹ 总宇称补偿 (−1)^{c_total·p_total} 相对每带 XOR ⊕_b c_b·(p_b⊕F_b) 恰差
交叉项 c₁p₂⊕c₂p₁。

**修复**（`src/system.jl`，mixed+rt 两路径）：补偿改每带 XOR；M=1 精确
退化（单带基线逐位不变）。**结果**：楔积 0/400；Z₂/Z₁² = 1.0；
mixed 斜率 **2.004**（与 U1 逐位一致）、rt 两带 **2.001**（err 8.3e-8）、
连续谱 U1 实时对齐 ≤1.4e-11。工具：`validation/{k_factor_mixed,
u0_factor_mixed,u0_factor_rt}.jl`。完整证据链见 `HISTORY.md` 第十轮。

## 第十轮补遗：两点函数全面闭合（2026-11-06）

**两个重大根因（均已修复，rt/mixed 同修）**：
1. **ysv↔site 带映射互换**：装配循环把 slot code 低位当 site 低带 ⟹ M≥2
   两条带在装配中被互换。带对称性掩蔽（对称观测量全部照常通过）——直到
   带插入 d₁⊗I 打破带对称才暴露。τ 路径（_code_from_ys）无此病。
2. **补偿缺双向跨带交叉项 ⊕_{b≠b'} c_b·p_{b'}**（楔积 inv 的逐元素实现）。
   之前「每带 XOR 修复 0/400」实为 bug 与补偿的巧合抵消。

**闭合（ED 端到端）**：混合 M=2 (i≥j) 全网格 greater/lesser 6.8e-4（离散
水平）、rt 全网格 Nt=8 1.54e-3（二阶收敛）、α→0 双轮廓机器精度、
matsubara M=1 全网格 + M=2 等时 ✅、occupation 四象限。另修复混合 − 支
插入切位 off-by-one（(−,0) 曾被静默丢弃）。单带基线与全部对称观测量
逐位不变。判决链与四次数训见 HISTORY 第十轮补遗 39–43。

**τ 状态（F 区门控跨带补偿后）**：matsubara M=1 全网格 3.67e-5 ✅、M=2
等时 3.73e-5 ✅、j≥1 两点大部 ≤1e-3 ✅；**唯一残留 = j==0（ddag 折入
fold 的 R，~1e-2 待推导）**——matsubara 对 j==0 显式抛错（force 绕过）。
API 承诺区间 i ≥ j（raw_correlation/matsubara）。
工具：validation/{twopoint_m2_check, mixed_twopoint_check,
twopoint_m2_frontier}.jl。
