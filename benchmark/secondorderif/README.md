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
                             Dirac-δ 闭式精确路径 + QuAPI Δt 谱路径（一般浴）
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
                           等时两点/逐带占据回归 + 连续浴 U1 对齐）
toulouse/                  无相互作用 Toulouse 模型验证（单带 + 单模浴）
singlemodebath/            两带 Anderson 模型（U≠0，每带单模浴）验证
plot_scaling.py            error-scaling 图（matplotlib）
HISTORY.md                 2026-11-04/05 七轮修复史（部件级对齐、跨框架
                           对比教训、假象补偿清算）——重要档案
```

## 运行

```bash
julia --project=. benchmark/secondorderif/toulouse/convergence.jl .            # ①②
julia --project=. benchmark/secondorderif/singlemodebath/convergence.jl .      # ③④
python3 benchmark/secondorderif/plot_scaling.py                                # 出图
# 共享校验：
julia --project=. benchmark/secondorderif/validation/kernel_paths_check.jl .   # 核双路径 + 双δ浴端到端
julia --project=. benchmark/secondorderif/validation/rt_fastpath_check.jl .    # 实时快路径回归
julia --project=. benchmark/secondorderif/validation/twopoint_band_check.jl .  # 等时 0⁺ + 逐带占据
```

## 已验证结果（ED 参考，2026-11-05）

| 模型 | 观测量 | order 1 斜率 | order 2 斜率 | 最细误差 |
|------|--------|--------------|--------------|----------|
| Toulouse | 虚时 ⟨n⟩ | 1.380 | **2.002** | 5.9e-6 (N=40) |
| Toulouse | 实时 ⟨n(0.4)⟩ | 0.822 | **1.996** | 1.9e-6 (N=8) |
| Anderson 双带 | 虚时 ⟨n₁+n₂⟩ | 1.607 | **2.004** | 7.5e-5 (N=16) |
| Anderson 双带 | 实时 ⟨n₁+n₂(0.4)⟩ | 0.584 | **1.927** | 3.0e-6 (N=32) |

Anderson 实时的 U≠0 离散误差常数较大（err/h² → 1.9e-2），纯二阶渐近区
N≥8。逐带 ⟨n_b⟩、等时两点 0⁺（⟨dd†⟩ = 1−⟨n⟩）、occupation(j=0) 折叠
均 ED 验证通过（`validation/twopoint_band_check.jl` 四节全绿）。

## 一般浴支持

- 虚时核：QuAPI `Δτ`，任意有界谱（含 `DiscreteSpectrum` 多 δ 谱）；
- 实时核：Dirac-δ 走闭式精确路径；一般谱走 QuAPI `Δt`（单元级积分值，
  `[1:N,1:N]` 直接截取）——双路径在 Dirac-δ 上互验 ≤4e-17；双 δ（两模）
  浴端到端 vs ED 通过（`validation/kernel_paths_check.jl`）；
- `α_weight(spec::DiscreteSpectrum) = Σ αₙ`（多 δ 谱的总耦合）。

## 连续谱浴 + U1GTEMPO2 对齐（2026-11-05）

连续带谱 A(ε) = 0.25·cos²(π(ε−0.2)/4)（∫A = 0.5）上的两带 Anderson：
`validation/continuous_bath_{mine,u1}.jl` 两条独立实现管线端到端逐点互比
（逐带 ⟨n₁⟩ ≡ U1 逐自旋 ⟨n↑⟩，τ 与实时、order 1/2、各 N 网格）：

- **虚时：逐点偏差 ≤ 7.9e-12（机器精度级）**；
- **实时：偏差 ≤ 4.8e-5 且随 N 收敛**（N=2 4.8e-5 → N=8 6.6e-6，趋同
  连续极限；残差来自 U1 的指数拟合容差与双方数值容差）。

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
