# bcsweyl — BCS 浴的 `coupling_norm` 定义：Weyl 标量与二阶配分函数

- 日期：2026-11-17
- 问题：`coupling_norm`（`src/correlationfunction.jl`）对 normal 浴已有定义
  （`DiracDelta→α`、`DiscreteSpectrum→Σαₙ`、连续谱→∫f dε），对 BCS 浴
  **未定义**（反常通道谱强度语义未定）。Weyl 二阶 K 的标量补偿
  `exp(cn/4·Σ step²)` 直接消费该值——BCS 浴下取什么，才能使配分函数 Z
  保持二阶收敛（`plan/weyl2_bcs_plan.md` F5 开放点）？
- 答案（推导 §2 + 数值裁决 §3）：

> **`coupling_norm(bath::AbstractBCSBath) = coupling_norm(bath.spectrum)`**
> ——即裸谱强度 ∫f dε（每带），**Δ、μ、β 都不进**。K 侧全带总量 =
> Σ_bands α_b（一个 BCS 对 = 两带各 α，总 2α）。反常通道（uv 块）对
> Weyl 标量**零贡献**。

## 1. 背景与机器

Weyl 二阶 = `sysdynamics_weyl2order`：对称半步 `S·M·S` + q=1/2 局部权重
（`{qI−n, −d†, d, qI+n}`）+ **Weyl 标量** `lscale = cn/4·Σ step²`
（τ 格 +δτ²、实时格 −δt²，`setlogscale!` 挂链）。标量的来源：接触对从
cell 端点重排到 cell 中点产生的**解析 c-number**
（`QuAPI/weyl/weyl_anderson_singlemode.tex` §logscale：恒等式级别的重整，
normal 浴已由单模/Anderson/Toulouse 六基准斜率 2.00 验证）。

BCS 侧机器：`bcs_partial_if`（四块有向指数
`ā↑Δ↑↑a↑ + ā↑Δ↑↓ā↓ + a↓Δ↓↑a↑ + ā↓Δ↓↓a↓`；normal 块 = u²/v² 谱、
反常块 = uv 谱，Bogoliubov 系数 `get_all_fs`）；浴模型
（`test/util.jl` `bcs_ed`）：

```
H = H_imp + ω₀(n₁+n₂) − (Δc₁†c₂† + Δ*c₂c₁) + √α(a₁†c₁+c₁†a₁+a₂†c₂+c₂†a₂)
```

浴自身带配对项，杂质只通过 normal 通道耦合；反常关联由浴配对态产生。

## 2. 推导：为什么反常道不进标量

**(a) 通道重分解。** Bogoliubov 旋转 `c₁ = uγ₁ + vγ₂†`、
`c₂ = uγ₂ − vγ₁†`（u²+v²=1，E=√(ω₀²+Δ²)）把浴对角化为 normal 准粒子浴。
杂质耦合展开为**四条线性通道**：(a₁↔γ₁: √α·u)、(a₂↔γ₂: √α·u)
（粒子道，宽 αu²）；(a₁↔γ₂†: √α·v)、(a₂↔γ₁†: √α·v)（空穴道，宽 αv²）。
BCS 四块正是这四条通道的 γ-传播子组合（`bcscorr.jl` 的 u²/uv/vu/v² 加权
同构）；准粒子浴对 γ 是 normal 的（⟨γγ⟩=0）。

**(b) 标量只数对角接触宽度。** Weyl 标量的 λ²/4 因子来自**同一对变量
的自接触**（āσaσ 密度道）指数的中点重排。反常接触 `ā↑κā↓ + a↓κa↑`
配对的是**不同变量**（单源槽，eq:weights 中重排不移动的槽位）；其进入
log Z 只能通过 κ² 级收缩——每格 O(λ⁴)、全环 O(λ³)，**污染不到 O(λ) 的
标量系数**。

**(c) 正则性收尾。** 每个准粒子模的影响泛函是杂质对 (a₁,a₂) 上的
2×2 Nambu 高斯：对角接触宽度 αu²（a₁ 槽）与 αv²（a₂ 槽）。
标量 = λ²/4·Σ 对角宽度（×两模）= λ²/4·2α(u²+v²) = **λ²/4·2α**。
u²+v²=1（{c,c†}=1 的正则反对易在 Bogoliubov 下不变）⟹ **cn 与 Δ 无关**，
等于 normal 浴的总通道范数 Σ_b α_b。μ、β 同 normal 情形一样不进
（谱权重 vs 温度核函数的分工不变）。

**(d) 定义。** `coupling_norm(::AbstractBCSBath) = coupling_norm(bath.spectrum)`
（每带裸谱强度；Δ=0 时与两条独立 normal 带的定义严格一致——回归兼容）。
若连续谱：`quadgkwrapper(f)` 同款；语义与 normal 浴的扩展点完全平行。

## 3. 数值裁决

### 3.1 接触缺口提取（`contact_extract_refined.jl`，平凡系统 U=ϵ_d=0）

U = ϵ_d = 0 ⟹ S = I，order1/2 之差纯为 {q=1/2 权重平移, W₄, 标量}——
接触机器被完全隔离。naive2（cn=0）的对数缺口
`g(δτ) = A₀·δτ + B·δτ²` 二次拟合直接给出 `cn* = −4A₀/β`
（N∈{16,24,32,48,64}，IF 截断 D=256/rtol=1e-14）：

| Δ | v² | cn*(Δ) | cn*/2α | B |
|-----|-------|---------|--------|------|
| 0.0 | 0     | 0.999886 | 0.99989 | 0.035 |
| 0.3 | 0.021 | 0.999880 | 0.99988 | 0.037 |
| 0.6 | 0.071 | 0.999863 | 0.99986 | 0.042 |
| 1.0 | 0.146 | 0.999822 | 0.99982 | 0.055 |
| 2.0 | 0.276 | 0.999632 | 0.99963 | 0.107 |

**裁决：cn*(Δ) ≡ 2α（=Σ_b α_b），Δ 无关，残差 ≤ 4×10⁻⁴ 相对**。
机器验证内嵌于 Δ=0 点（normal 浴理论已知精确值 2α，提取复现至 1.1×10⁻⁴）；
残差随 B 同步增长（Δ=2 时 B×3、残差×3.4）= O(δτ³) 拟合系统差，
非物理 Δ 修正——任何正则保持的结构性修正（v⁴ 级等）都会出现在 ≥1e-2 档。
（初版 `contact_extract.jl` 线性拟合被 O(δτ²) 项偏置 ~2%，曾误报 Δ 依赖；
其数据档案保留在 `results/contact_extract.csv` 作对照。）

### 3.2 generic H_S 收敛阶（`z_convergence.jl`，U=1, ϵ_d≠0, β=0.8）

cn = Σα_b = 2α 下 |Z−Z_ed| 的 log-log 斜率（N∈{8,16,24,32,48}）：

| Δ | order1 | order2 (cn=2α) |
|-----|--------|----------------|
| 0.3 | 0.975  | **1.980** |
| 0.6 | 0.974  | **1.977** |
| 1.0 | 0.972  | **1.973** |

order2 + cn=Σα_b 全 Δ 干净二阶；斜率略低于 2.00 与 normal 浴基准同款
（拟合区间含 O(δτ²) 项曲率的系统差）。**naive2（cn=0）对照实测：斜率
0.982–0.986（全 Δ 退化到一阶）**——Weyl 标量携带全部 O(δτ) 修正，
定义取错/不取即失二阶。数据：`results/z_convergence_sum_alpha.csv`、
`results/z_convergence_zero.csv`。

## 4. 结论与落地

- `src/correlationfunction.jl`：`coupling_norm(bath::AbstractBCSBath) =
  coupling_norm(bath.spectrum)` 已落地（每带裸谱强度；docstring 注明反常道
  零贡献与推导出处）；`test/bcs/runtests.jl` 补语义 testset（Δ/μ/β 不进、
  Δ=0 与 normal 一致、DiscreteSpectrum/BCSVacuum 覆盖、带对 2× 总量）。
- 语义对齐：`bcs_hybriddynamics` 一个浴 = 一个带对；K 侧 cn =
  `2·coupling_norm(bath)`（两带）。
- `plan/weyl2_bcs_plan.md` F5 的裁决：naive 组合 + **正确 cn** 二阶收敛；
  反常道不进 Weyl 标量（Phase 2 形式化的"normal 通道 α 和"语义由此证实，
  且 BCS 情形的"normal 通道"= 裸谱 α——正则性使 Δ 不出现）。

## 5. 范围与后续

- 实时/混合轮廓：Weyl 标量符号随轮廓段翻转（`_cell_step2`），通道宽度
  论证轮廓无关——cn 定义不变；实轴 Z（真空存活率）基准未跑（后续项）。
- 连续谱 BCS 浴：定义走 `quadgkwrapper(f)` 同款（u²+v² 正则性逐模成立），
  未单独跑收敛基准。
- 本目录脚本为一次性裁决档案（`results/` 留档），不入 `test/` 套件；
  语义回归由 test/bcs 的 coupling_norm testset 承担。
