# Bethe 格超导（Nambu）DMFT（虚轴）：GTEMPO vs TRIQS cthyb 对比

本目录用两种完全独立的杂质求解器求解**超导态** Bethe 格 DMFT
（Nambu / Werner–Millis 矩阵表示，虚轴），并逐迭代对比二者结果：

- **GTEMPO**：虚轴 Grassmann 张量网络工作流，BCS 浴（`bethe_dmft_nambu_gtempo.jl`）
- **TRIQS cthyb**：连续时间 QMC，Nambu 矩阵模式 + delta 接口（`bethe_dmft_nambu_ctqmc.py`）

参数：β = 1，U ∈ {1, 5}，半满填充（ϵ_d = −U/2），t = 1。
GTEMPO 侧流程参考 `GTEMPOProjects/src/bcsbath/siam/main.jl`（main_imag_dmft）。

## 1. 模型与自洽循环

Anderson 杂质模型（与正常态脚本一致）：

```text
H_imp = ϵ_d (n↑ + n↓) + U n↑n↓ ,   ϵ_d = −U/2（粒子-空穴对称）
```

超导浴：Bethe 格自洽条件写成 Nambu 四分量形式（参考代码约定）：

```text
Δuu′ = t² G_uu ,   Δdd′ = −t² conj(G_dd) ,
Δud′ = t² G_ud ,   Δdu′ = t² G_du
```

其中 uu/dd 为正常分量、ud/du 为反常分量。初始 Weiss 场取
Δuu = t²·G₀ᴮᵉᵗʰᵉ（解析半圆，注意平方根支选保证 Δ(−iωₙ) = Δ(iωₙ)*），
Δdd = −conj(Δuu)，反常分量取宽度 ω_c = 1、幅度 f₀ = 0.3 的实偶函数。

## 2. GTEMPO 侧要点

| 步骤 | 所用函数 |
| --- | --- |
| 四分量 Δ(iω) → Δ(τ) 核 | `Δiw_to_Δτ`（QuAPI，GTEMPO re-export） |
| BCS 影响泛函（单个 IF，不用 branching） | `hybriddynamics_naive(lattice, BCSCorrelationFunction(...), orbital=1)` |
| 环境缓存 | `environments(lattice, mpsK, mpsI)` |
| 正常分量 G(τ) | `cached_gf_fast(..., band=1/2)`（返回 Nτ+1 点，G = −C） |
| 反常分量 G(τ) | 逐点 `cached_gf` + 跨 band `ContourIndex`（端点 F(β) = −F(0)） |
| G(τ) → G(iωₙ) | 细网格插值 + `Gτ_to_Giw` |
| 频率网格 | `ifrequencies` |

## 3. cthyb 侧要点（Nambu 矩阵模式 + delta 接口）

抽象 flavor 取 c₀ = d↑、c₁ = d↓†（n₁ = 1 − n↓）。抽象 flavor 作用
S_hyb = ∫∫ c̄ᵢ(τ) Δᵢⱼ(τ−τ′) cⱼ(τ′) 表明 **GTEMPO 的四个块直接放入 2×2
矩阵，无需做任何粒子-空穴共轭变换**：

```python
S = Solver(beta=beta, n_iw=..., n_tau=..., gf_struct=[("nm", 2)],
           delta_interface=True)
S.Delta_tau["nm"] << make_gf_from_fourier(Delta_iw)   # [[Δuu, Δud], [Δdu, Δdd]]
S.solve(h_int=U * n("nm",0) * (1 - n("nm",1)),
        h_loc0=(-U/2) * (n("nm",0) - n("nm",1)), ...)  # h_loc0 通过 solve 传入
```

测量的抽象 G 分量换算回物理分量：

```text
G_uu(τ) = G₀₀(τ) ,   G_dd(τ) = +G₁₁(β−τ) ,
G_ud(τ) = −G₀₁(τ) ,  G_du(τ) = −G₁₀(β−τ)
```

以上构造（尤其 Δ₁₁ = Δdd 原样放入）经过 U = 0 精确解交叉验证：
U = 0 时精确解为 Nambu 矩阵求逆 G = [iω·1 − Δ]⁻¹，两个求解器均与之吻合
（`data/gtempo_nambu_b5_U0.json`、`data/ctqmc_nambu_b5_U0.json`、
`compare_gtau_b5_U0.png`）；用错误的 p-h 共轭构造（Δ₁₁ = −conj(Δdd(−iω))）
则会使 cthyb 平均符号崩塌到 ~10⁻³，正确构造下平均符号 = 1（U = 0）。

## 4. 文件

| 文件 | 说明 |
| --- | --- |
| `bethe_dmft_nambu_gtempo.jl` | GTEMPO Nambu 虚轴 DMFT 迭代 |
| `bethe_dmft_nambu_ctqmc.py`  | TRIQS cthyb Nambu DMFT（读 GTEMPO 导出的 Δ(iω)） |
| `compare_plot.py`            | 逐迭代一致性检查 + 画图 |
| `data/`                      | JSON 结果、对比图、一致性数据 |

## 5. 运行

```bash
# GTEMPO 侧（DMFT_SET=k 只跑第 k 组，便于并行）
julia --project=<env> docs/tutorials/bethedmft/nambu_imag/bethe_dmft_nambu_gtempo.jl

# TRIQS 侧（需要 triqs + triqs_cthyb；DMFT_NCYCLES 控制 QMC 精度）
python docs/tutorials/bethedmft/nambu_imag/bethe_dmft_nambu_ctqmc.py

# 对比与画图
python docs/tutorials/bethedmft/nambu_imag/compare_plot.py
```

环境变量：`DMFT_SET`（只跑某一组）、`DMFT_SMOKE=1`（小参数冒烟）、
`DMFT_U`（覆盖 U，U = 0 用于精确解交叉验证）、`DMFT_NCYCLES`
（cthyb 循环数，默认 2×10⁶）。

生产参数：GTEMPO δτ = 0.05（Nτ = 20）、D = 64、Nω = 512、
10 次迭代、Δ 混合系数 0.5；cthyb 每迭代 2×10⁷ 个 QMC 循环。
（U = 0 精确交叉验证在 β = 5 下完成，构造约定与 β 无关。）

## 6. 结果判读

`compare_plot.py` 对每次迭代给出四个物理分量在相同 Δ(iω) 下的

```text
rel.dev G_c(τ) = ‖G_c,GTEMPO(τ) − G_c,cthyb(τ)‖ / ‖G_c,GTEMPO(τ)‖
```

**β = 1 实测结果**（cthyb 每迭代 2×10⁷ 循环，平均符号 = 1）：

- 正常分量（有效信号）：全部迭代 rel.dev ≈ 1.2–2.1×10⁻²
  （2×10⁶ 循环时为 3–7×10⁻²，按 1/√N 改善），由残余 QMC 噪声与
  GTEMPO 离散化误差共同决定；
- 反常分量：β = 1 温度高，F 幅度小且被自洽流压低（U = 1 从 ~0.03 衰减到
  ~10⁻³，U = 5 仅 ~0.004 → 0，排斥 U 下高温不维持超导）。信号可分辨的
  前几次迭代两者吻合（U = 1 rel.dev ~6–14×10⁻²），信号衰减到 10⁻⁴ 量级后
  cthyb 统计噪声主导；反常通道的**约定正确性**由 U = 0 精确解交叉验证
  保证（该测试中反常分量有良好分辨）；
- DMFT 自洽收敛（GTEMPO 侧 ‖Δ′−Δ‖/‖Δ‖，U = 1）：iter 1 → 10：
  3.3×10⁻¹ → 7.2×10⁻⁴。

对比图 `compare_gtau_b1_U{1,5}.png`：每个 U 一张图、6 个子图展示前 6 次
迭代；正常分量（左轴实线）与反常分量（右轴虚线，放大）分别绘制，
GTEMPO 为线、cthyb 为空心圆点。
