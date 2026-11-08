# 虚轴严格 IF：integrate 的算符迹结构（理论笔记，2026-11-08）

状态：管线配方已机器精度闭环（`imag_closure.jl` / `imag_general.jl`，N=2/o2 全绿）；
本文记录 τ2 泛函的解码与闭式识别的进展/障碍，供后续理论工作接续。

## 1. τ2 泛函（integrate 的收缩语义）— 已完全测定

黑盒测量（`imag_trace_oracle.jl`，1024 次 pair-注入 + Möbius，`imag_tau2.csv`）：
`Z(S, I) = Σ_{a,b} S_a I_b τ2[(a,b)]`，τ2[(A_K, A_I)] = ±1 当且仅当
A_K ∪ A_I = 完整切片块的并，且每块的占据 ∈ {全K, 全I, 空, K/I 各一腿}；
"同切片 K 对 ∧ I 对"（4 腿 0）被禁；跨链"各取一腿"带 JW 型符号
（符号 = −1 当两根 K 腿同为 eta 或同为 bar）。

**算符迹解释**：imag-lattice 的 GMPS 是算符链（每切片一个 c/c† 对 = 时间副本），
`integrate(K, I)` = **两条算符链乘积的迹** Tr[ρ_K · ρ_I]，其中
- 每切片块的占据规则 = 同一生成器 θ² 幂零性（K 与 I 共享切片生成器 ⇒ 双占即灭）；
- 分裂情形的符号 = 把乘积重排到正则序的费米子交换符号；
- 跨切片"共享副本"结构给出 Tr[c] = 0 型相消（K={0,3} 单链迹为 0 等）。

这一解释把全部 17 个探针值（`imag_trace_probe.jl`）与 τ2 表统一起来。

## 2. 管线配方（已锁定，机器精度）

- **K 链**（自建，GTerm 注入）：门 exp(GTerm(bar_{g+1}, eta_{g+1}, e^{−τε_g}−1))
  + 杂质跨位分辨率 exp(GTerm(bar_{p+1}, eta_{p+2 环闭合}, −1))
  + BC exp(GTerm(bar_1, eta_k, −1))。
  ⚠ 库门张量（`_propagator_sparsegmps`+`mult!`）的多项式 ≠ (e^{−τε}−1) kernel 形式
  —— combo 判决：门用 ±e^{−τε} 时联合行不一致（0.43），(e^{−τε}−1) 系一致（1e-11）。
- **I 链** = 影响多项式 I*：联合 Z+n 行（h/g 均经 I 侧 Möbius 到单项式基；
  g_raw 经 `GTEMPO._mult_A(GTerm, mpsK)` 直测）min-norm 解 → 序贯消元
  exp(GTerm) 因子（任意偶腿）实现。ν_eff = I* 的 ∅ 系数。
- **读出**：n = 1 − g_end，g_raw/Z 经 `_mult_A` + integrate。
- 结果（N=2, order 2）：Z err ~1e-12、n err ~1e-13（εd 0..1.4）。

## 3. 结构发现

- **缝因子解耦**：I* = (1 + θ_{4}θ_{5}) ⊗ I' 精确成立（seam 切片自分辨率 +1，
  与 K 链 BC 无关地单独出现）。
- **I' 的 pair 投影** = 2×3 影响矩阵（行 = 切片 2,3 的 bar；列 = 切片 1,2,3 的 eta），
  切片 2/3 对称（Δ[2,3] = Δ[3,2]、Δ[2,2] = Δ[3,3]）——环平移不变性的体现。
- **3×3 exp-form 族（复规范）精确且唯一**（`imag_seamfit2.jl` 命中 9e-14；
  24 行联合约束 rank ≈ 11 复 ≥ 19 实参数 ⇒ 解唯一）。实行积族（F2）无解 ⇒
  多项式 I* 不是任何"实行积"的像。
- **h 行秩 = 3**（Z 的 εd 依赖 = x^0..x²，x = e^{−βεd}）；g 通道秩 7-8
  （n̂ 插入通道阶更高）——两者给联立系统的有效约束维数。

## 4. 闭式识别（未完成——障碍与路线）

u(k) 代数封闭性论证：双线性算符 c†ᵢcⱼ 封闭（[c†ᵢcⱼ, c†_kc_l] = δ_jk c†ᵢc_l − δ_il c†_kc_j），
故"分辨率/BC/门的逆 × 影响算符 Ω"都是 exp(bilinear)，I 链 = 单个 exp-form，
Δ̃ = log(ρ_K⁻¹ e^{Ω}) 应为矩阵公式。**但** `imag_matrixid.jl` 的直接矩阵公式扫描
（Ω ∈ {log(I±D̃), ±D̃} × ρ_K 序 × 对角相似变换）不命中（最优 94% 误差）。

障碍 = **img-弦规范**：多项式 → 算符像的映射带 JW 型弦（跨副本），弦使
"多项式 exp-form 的系数矩阵"与"算符 bilinear 矩阵"之间差一个非对角的
（Fock 基下对角、切片指标下非对角的）规范变换。数值唯一解 Δ̃ 含 O(1) 虚部
（如 [3,2] 虚部 2.918），正是弦相位的量级。

**路线**（下一轮）：①从探针值直接钉死 img 映射的弦（~5 个探针即可标定）；
②在弦规范下重推 ρ_K 的矩阵（弦 → Fock 基对角算符，切片指标下=多体对算符，
u(k) 论证在"弦规范矩阵"层面改用 2^k 维显式矩阵：8×8（N=2）—— 直接构造
ρ_K、Ω 的 8×8 矩阵，Δ̃ = log 回拉到 bilinear 块——数值上完全可行）；
③核对 exp-form(Δ̃) = I* 后写出闭式。

## 5. 数据/脚本索引

- 数据：`imag_tau2.csv`（τ2 全表）、`imag_h.csv`/`imag_g.csv`（库 K 链行，历史）、
  `imag_n2o2_*.csv` / `imag_n3o2_*.csv`（v2 K 链行，`imag_general.jl` 产出）
- 判决/求解：`imag_trace_probe.jl`（17 探针）、`imag_trace_oracle.jl`（τ2 全测）、
  `imag_k2test.jl`（combo 一致性）、`imag_seamfit2.jl`（3×3 族唯一解）、
  `imag_closure.jl`（N=2 闭环）、`imag_general.jl`（参数化 N/order 全流程）、
  `imag_identify.jl` / `imag_matrixid.jl`（闭式识别尝试，未完）
