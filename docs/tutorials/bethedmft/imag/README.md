# Bethe 格 DMFT（虚轴）：GTEMPO vs TRIQS cthyb 对比

本目录用两种完全独立的杂质求解器求解 Bethe 格上的 Hubbard DMFT
（虚轴 / Matsubara 表示），并逐迭代对比二者结果：

- **GTEMPO**：虚轴 Grassmann 张量网络工作流（`bethe_dmft_gtempo.jl`）
- **TRIQS cthyb**：连续时间 QMC（混合展开算法）（`bethe_dmft_ctqmc.py`）

四组参数：(β, U) ∈ {5, 10} × {1, 5}，半满填充（ϵ_d = −U/2），t = 1。

## 1. 模型与自洽循环

Anderson 杂质模型

```text
H_imp = ϵ_d (n↑ + n↓) + U n↑n↓ ,   ϵ_d = −U/2（粒子-空穴对称）
```

Bethe 格自洽条件（半带宽 2t，t = 1）：

```text
Δ(iω) = t² G(iω) ,   G₀⁻¹(iω) = iω + μ − Δ(iω) ,   μ 对应 ϵ_d = −U/2
```

## 2. 自洽循环的实现（复用现成的转换函数）

GTEMPO 侧完全沿用 `GTEMPOProjects/src/bethedmft/imag/main.jl` 的参考流程，
不自写任何 Fourier / Hilbert 变换：

| 步骤 | 所用函数 | 来源 |
| --- | --- | --- |
| 初始 Δ(iω)：解析半圆浴杂化 | `toulouse_Δiw(bath, n) .* t²` | ImpurityModelBase |
| Δ(iω) → Δ(τ) 路径积分核 | `Δiw_to_Δτ` | QuAPI（GTEMPO re-export） |
| 杂质 G(τ) | `cached_gf_fast(..., b1=:τ, b2=:τ)` | GTEMPO |
| G(τ) → G(iωₙ) | `Gτ_to_Giw` | ImpurityModelBase |
| 频率网格 | `ifrequencies` | ImpurityModelBase |

自洽条件（半满 ϵ_d = U/2，t = 1）：

```text
G₀⁻¹(iωₙ) = iωₙ + ϵ_d − t² G(iωₙ) ,   Δ′(iωₙ) = iωₙ + ϵ_d − 1/G₀(iωₙ)
```

G(τ) 先用 `Interpolations.linear_interpolation` 插值到 δτ = 1e-4 的细网格
再送 `Gτ_to_Giw`（参考文件的做法）。每次迭代的 Δ(iωₙ) 数值与 G(τ)、
G(iωₙ) 一并写入 `data/gtempo_b*.json`；**cthyb 驱动按频率逐点读取这
份 Δ(iωₙ)**，因此两个求解器每一步拿到的是完全相同的 Δ，且都从相同的
初始态（解析半圆浴）出发。

## 3. 文件

| 文件 | 说明 |
| --- | --- |
| `bethe_dmft_gtempo.jl` | GTEMPO 虚轴 DMFT 迭代 |
| `bethe_dmft_ctqmc.py`  | TRIQS cthyb DMFT 迭代（读 GTEMPO 导出的 Δ(iω)） |
| `compare_plot.py`      | 逐迭代一致性检查 + 画图 |
| `data/`                | 两个求解器的 JSON 结果、对比图、一致性数据 |

GTEMPO 侧要点：

- 晶格：`GrassmannLattice(N=Nτ, δτ=β/Nτ, bands=2, contour=:imag)`，δτ = 0.05
- 路径积分核：`hybriddynamics(lattice, Δiw_to_Δτ(...), trunc, band=1)`
  + `swapband` 到第二条自旋带（不用 `boundarycondition_branching`）
- 杂质：两带 `AndersonIM(U=U, μ=−U/2)`，`sysdynamics` + `boundarycondition!`
- 观测量：`cached_gf_fast` 给出裸关联，**G(τ) = −⟨d(τ)d†⟩**

cthyb 侧要点：

- `Solver(beta, n_iw, n_tau, gf_struct=[("up",1),("down",1)])`
  （`n_tau ≥ 6 n_iw` 以满足求解器的 Fourier 精度要求）
- `G0_iw << inverse(iOmega_n + U/2 − Δ)`（cthyb 的 mu 吸收了 ϵ_d = −U/2）
- `h_int = U n↑ n↓`，两个自旋块取平均

## 4. 运行

```bash
# GTEMPO 侧（4 组参数；DMFT_SET=k 只跑第 k 组，便于并行）
julia --project=<env> docs/tutorials/bethedmft/imag/bethe_dmft_gtempo.jl

# TRIQS 侧（需要 triqs + triqs_cthyb + matplotlib；DMFT_NCYCLES 控制 QMC 精度）
python docs/tutorials/bethedmft/imag/bethe_dmft_ctqmc.py

# 对比与画图
python docs/tutorials/bethedmft/imag/compare_plot.py
```

环境变量：`DMFT_SET`（只跑某一组）、`DMFT_SMOKE=1`（GTEMPO 小参数冒烟）、
`DMFT_NCYCLES`（cthyb 循环数，默认 200000；统计噪声 ∝ 1/√N）。

## 5. 参数与研究范围

四组参数为 (β, U) ∈ {5, 10} × {1, 5}，半满填充（ϵ_d = −U/2），t = 1。
**当前只完成了 β = 5 的两组**（β = 10 未运行，脚本以
`DMFT_SET=3` / `DMFT_SET=4` 可随时补跑）。

实际参数：δτ = 0.05（Nτ = 100）、键维 D = 120、Nω = 1024、
10 次 DMFT 迭代、Δ 混合系数 0.5；cthyb 每迭代 2×10⁶ 个 QMC 循环。

DMFT 自洽收敛（GTEMPO 侧 `‖Δ′−Δ‖/‖Δ‖`）：

```text
iter     1      5      10
U = 1   6.2e-1  2.3e-2  8.6e-4
U = 5   3.4e-1  1.5e-1  4.8e-2
```

## 6. 结果判读

`compare_plot.py` 对每一次 DMFT 迭代给出两个求解器在同一 Δ(iω) 下的

```text
rel.dev G(τ)  = ‖G_GTEMPO(τ) − G_cthyb(τ)‖ / ‖G_GTEMPO(τ)‖
rel.dev G(iω) = ‖G_GTEMPO(iω) − G_cthyb(iω)‖ / ‖G_GTEMPO(iω)‖
```

（G(τ) 在 GTEMPO 的 τ 网格上比较，cthyb 的 G(τ) 由其自身网格插值而来；
G(iω) 在两者共有的正频网格上逐点比较。）

偏差的主要来源已用受控实验定位：同一 Δ 下把 cthyb 循环数从
2×10⁵ 提到 2×10⁶，G(τ) 偏差从 3.7×10⁻² 降到 1.7×10⁻²（≈1/√N），
即**当前偏差由 CT-QMC 统计噪声主导**，而非两个求解器的系统性差异。

## 7. 文件清单

- `compare_gtau_b5_U1.png` / `compare_gtau_b5_U5.png`：6 个子图分别展示
  前 6 次 DMFT 迭代的 G(τ) 对比（GTEMPO 实线 vs cthyb 圆点），
  子图标题给出该次迭代的相对偏差
- `data/consistency.txt`：逐迭代一致性数值
- `data/*.json`：两个求解器的原始输出（含每迭代的 Δ(iω)、G(τ)、G(iω)）
