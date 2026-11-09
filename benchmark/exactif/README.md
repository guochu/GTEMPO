# exactif —— IF 严格核方法（exact_Δt）基准

以**严格核** `exact_Δt` 替代 QuAPI 提供的冻结核 `Δt`，在 GTEMPO 实时 Keldysh
管线上构建并验证各阶 Trotter 演化。

参考：N. Ng, G. Park, A. J. Millis, G. K.-L. Chan, D. R. Reichman,
*Real-time evolution of Anderson impurity models via tensor network influence
functionals*, **PRB 107, 125103 (2023)**（`benchmark/NgReichman2023.pdf`）；
外部 IFGTEMPO 包（两自旋杂质 + 有限 bath 的同型实现）。

## 方法

- **冻结核 `Δt`**（QuAPI `fermionic_Δt`，现状）：cell 积分（双胞三角积分）
  把 system-bath 耦合源冻结在 cell 内，顶点幅只有 O(V²) 精度。
- **严格核 `exact_Δt`**（QuAPI `exactfermionic.jl`）：bath + 耦合的二次型
  cell 演化 `B = e^{-δt·Σ(εc†c + V(d†c+c†d))}` 被**精确**处理（单粒子传播子
  显式矩阵指数），杂质相干腿在节点插入，bath 腿经 Grassmann 高斯 Schur
  消元精确积掉。本目录经 `exact_Δt(bath, cells)`（`src/NativeExact.jl`）适配为
  GTEMPO `AbstractRealCorrelationFunction`，直接供 `hybriddynamics!` 的
  partialmpo 行扫消费（无 Prony、无 Toeplitz 近似）。

### Trotter 调度（`trotter_factors`）

杂质演化与二次型 cell 的 Trotter 分割由 `order ∈ (1, 2, 4)` 选择：

| order | 分割 | lattice（pair 数） |
|---|---|---|
| 1 | `e^{-iHδt} ≈ A(δt)·B(δt)`（同 pair 放置） | N |
| 2 | `A(δt/2)·B(δt)·A(δt/2)`（合并链交错放置） | 2N+1 |
| 4 | Suzuki `S2(aδt)²·S2(bδt)`，`2a+b=1`（b<0） | 6N+1 |

order ≥ 2 的因子**必须交错占位**（A 门与 B cell 各占一个 pair；合并后相邻
A 吸收为全步长），零长度 cell 补齐使核节点与格点切片一一对应；每个 pair
都放置传播子门（τ=0 的恒等门携带相干态分辨率重叠 ∏(1+c̄c)，是测度的一部分，
跳过会使 Z 偏 ~1/2）。

## 文件

**目录结构**：`src/` = 核心代码（NativeExact.jl 单核心模块）；`commons/` = 两模型文件夹共享的 helper（realtime 引擎 `scaling_compact_models.jl` + `formb_dense.jl` 稠密参照机器 + 虚轴标度 `imag_commons.jl` + 共享绘图脚本）；`toulouse/`、`singlemodebath/` = 各自模型的实时/虚轴测试脚本 + 数据 + 图；`theory/` = 数学原理文档；`imag_theory_notes.md` = 虚轴严格 IF 全程理论笔记（τ2/探针族/单 MPS/调度阶数）。`validation/` 已清空——工作期脚本与数据在 git 历史（commit 843b92f..c25ea43），结论沉淀于 `imag_theory_notes.md` 与下表。

| 文件 | 内容 |
|---|---|

| `src/NativeExact.jl` | **IF 严格核核心模块**（IFGTEMPO-free，原 common.jl + NativeExact 合并）：①原生管线 `exactif_native_compact`（紧凑消元，orders 1/2/4 主力）/`exactif_native_quapi`（orders 1/2，QuAPI-Δ 数据源）/`compact_delta`/`sysdynamics_exact!`/`gtlt_compact`（双时 gt/lt）/`open_rho_native`（ρ(T) 开放腿读出）；②Z 证书机器 `exact_Δt` 适配器 + `exactif_pipeline`/`frozen_pipeline` 双管线（zidentity 依赖）；③`plot_png` 绘图助手 |
| `theory/MATHEMATICS.tex/.pdf` | **形态 B 数学原理文档（LaTeX 正式版，8 页）**：两扇区路径积分 → IF 变量高斯消元（Δ_eff 推导）→ 收缩语义/迹闭合 → GTEMPO pair 机制对接（反测度/graded 读出符号）（xelatex 编译，`ctexart`） |
| `theory/keldysh.md/.tex/.pdf` | **严格 IF 实时核全流程闭环理论（Keldysh 轮廓）**：①bath 精确积掉（两扇区 Gaussian + Schur）→ ②IF Grassmann 变量积掉（Pfaffian 恒等式与消元核 D）→ ③1/2/4 阶调度与杂质部分完全组装（含回文 Suzuki 必要性）→ 读出约定与闭环验证阶梯 |
| `theory/MATHEMATICS.md` | 数学原理文档的 Markdown 版（与 .tex 同内容；供仓库内直接预览） |
| `validation/zidentity.jl` | **Z 恒等式证书**（ε_d=0 时 Z·prefactor == 1，全部 order 机器精度） |
| `validation/formb_dense.jl` | **形态 B（消元式单扇区）稠密系数级验证 + 活跃验证链的共享稠密机器**：IF 变量从两扇区 lattice 精确积掉，与 Fock 对齐（`FORMB_QUICK=1` 快测档）。**已 IFGTEMPO-free**（Fock 参照链/调度表本地化，鸭子类型 bath 参数） |
| `validation/verify_compact.jl` | 紧凑路线对齐矩阵（离散链 + Z/V=0 双证书；`ARGS=trunc` 截断档） |
| `validation/scaling_compact_models.jl` | **两模型（Toulouse + 单模 bath Anderson U≠0）error scaling + 真阶判定**；结果分别写入 `toulouse/`、`singlemodebath/` |
| `validation/verify_align.jl` | 与 IFGTEMPO GMPS 的直接对齐（两自旋 U=0 解耦恒等式）+ quapi 回归（`ARGS=quapi`） |
| `validation/verify_align_u15.jl` | **IFGTEMPO 直接同模型对齐（U=1.5 Anderson）**：双可观测量 n_up + ⟨n₁n₂⟩，8 组 ≤2.2e-14 |
| `validation/verify_gtlt_align.jl` | **双时 gt/lt Green 函数系列对齐**：①单带 vs IFGTEMPO 解耦 flavor-1（4 组）、⑤双带 U=1.5 直接同模型（3 组），全部 ≤8.6e-14；附 ED 有效时刻 Trotter 收敛 |
| `toulouse/realtime.jl` + `realtime.{png,csv}` | **Toulouse（U=0, ϵ_d=0.7）基于严格 IF 的 1/2/4 阶 Trotter error scaling**（紧凑消元路线；机器级预检 + 真阶判定 PASS；共享引擎 = `validation/scaling_compact_models.jl`） |
| `singlemodebath/realtime.jl` + `realtime.{png,csv}` | **单模 bath Anderson（μ=0.7, U=1.5, bands=2）同上**（16 维 Fock 参照含 U 关联；PASS） |

注：`twosector/`（两扇区权威标度基准，旧调度下的 p=4.05-4.11）已随 2026-11-06 目录清理
移除——其 order-4 数值属旧调度的交叉区假象（见"order-4 调度修正"一节）；如需两扇区
交叉验证，用 IFGTEMPO 的 `build_impurity_path`（已修回文调度）重建。

**2026-11-07 早期代码清理**：以下已被紧凑路线取代/仅具历史价值的代码已删除
（结论保留于下方沿革记录与 git 历史）：`src/FormB.jl`（形态 B GMPS 移植）、
`validation/verify_formb.jl`/`scaling_formb.jl`/`formb_legacy_ifgtempo.jl`
（FormB 验证链，IFGTEMPO 依赖）、`validation/native_probe.jl`/`restrictd_probe.jl`
（判决实验，结论已锁定进实现）、`validation/verify_native.jl`/`calibrate_native.jl`/`validate.jl`
（contour_G 路线与早期验证——contour_G 路线的 `exact_delta`/`exactif_native_pipeline`
已从 NativeExact.jl 移除）、`validation/scaling_compact.jl`（被 scaling_compact_models.jl
取代）、`formb_scaling_*.{csv,png}`（旧调度标度图记录）。

**依赖状态（2026-11-07）**：`src/` 与活跃验证链（zidentity / verify_compact /
scaling_compact_models / 各模型 realtime）**已 IFGTEMPO-free**——Fock 参照链与
调度表本地化（`formb_dense.jl` 的鸭子类型 bath 参数 + `NativeExact.native_schedule`）。
IFGTEMPO 仅在对齐参照脚本出现（verify_align* / verify_gtlt_align——
"用 IFGTEMPO 验证自己"的证据链，保留）。

运行（GTEMPO 根目录）：

```
julia --project=. benchmark/exactif/validation/zidentity.jl               # Z 恒等式证书
julia --project=. benchmark/exactif/validation/verify_compact.jl          # 紧凑路线对齐矩阵 + 双证书
julia --project=. benchmark/exactif/validation/verify_compact.jl trunc    # 截断档（order 4 N=2,3）
julia --project=. benchmark/exactif/validation/scaling_compact_models.jl  # 两模型 1/2/4 阶标度（出图）
julia --project=. benchmark/exactif/validation/verify_align.jl            # IFGTEMPO 直接对齐（U=0 解耦）
julia --project=. benchmark/exactif/validation/verify_align.jl quapi      # quapi 路线回归
julia --project=. benchmark/exactif/validation/verify_align_u15.jl        # IFGTEMPO 同模型 U=1.5 对齐
julia --project=. benchmark/exactif/validation/formb_dense.jl             # 形态 B 消元验证（稠密，全量）
FORMB_QUICK=1 julia --project=. benchmark/exactif/validation/formb_dense.jl   # 快测档
julia --project=. benchmark/exactif/toulouse/realtime.jl                  # Toulouse 严格 IF 1/2/4 阶标度（出图）
julia --project=. benchmark/exactif/singlemodebath/realtime.jl            # Anderson 严格 IF 1/2/4 阶标度（出图）
```

**注意**：`toulouse/` 与 `singlemodebath/` 的 `realtime.jl`（2026-11-06 起）是
**紧凑消元路线的严格 IF 1/2/4 阶标度**——旧的单扇区 gt/lt 诊断脚本（受闭合约定
差异限制、order 4 不收敛）已删除。两扇区交叉验证看 IFGTEMPO 的
`test_trotter_orders.jl`（已修回文调度）。

## 已验证的结论

1. **Z 恒等式证书**（`validation/zidentity.jl`）：ε_d = 0 时
   `Z·prefactor == Tr[ρ₀] = 1` 对全部 order、全部 N 成立到机器精度
   （最差 4.5e-13）。冻结核同一恒等式偏差 5e-3 ~ 1e-2（O(δt)）。这是区分
   严格核与冻结核的最锐利证书，也逐位锁定了适配器的全部符号/测度约定
   （读序节点取值 + 跨支块内嵌 Keldysh 负号 + 每 pair 恒等门）。

2. **IFGTEMPO 交叉验证**（外部包，本轮已修复其对当前 GTEMPO API 的适配并
   新增 `test/test_trotter_orders.jl`）：两自旋杂质 + 有限 bath 上，同一
   严格核方案对**精确连续 Fock 演化**的密度矩阵误差：
   - order 1（全步长放置）：**全局二阶**（ratio 4.00）——精确顶点幅已重求和
     `[A,B]δt²` 首项（QuAPI exactfermionic.jl 的 O(V²δt³)/cell 结论）；
   - order 2（Strang）：**全局二阶**（ratio 4.00）；
   - order 4（Suzuki）：**四阶**（ratio 16.5/17.3）。**勘误（2026-11-06）**：该比值
     出自非回文序 S2(a)²S2(b) 的交叉区——修正为回文序 S2(a)S2(b)S2(a) 后，同一测试
     给出干净比值 15.96/15.93（9/9 仍过；见"order-4 调度修正"一节）。

3. **GTEMPO 冻结核基线**：`frozen Δt` 的 gt/lt 系列误差标度 p ≈ 1.03
   （Toulouse）/ 0.95（Anderson）——与 `secondorder/convergence.jl` 的一阶
   结论一致。

## 与 GTEMPO 现有管线的差异（对齐 IFGTEMPO 的修复记录，备将来集成）

### 根源：接口测度的耦合结构

严格核（exact_Δt）链在 GTEMPO 现有单扇区管线上，ε_d ≠ 0 时可观测量
（gf / occupation）只表现 O(δt) 收敛，order 4 的 gt 系列不随 δt 收敛，
而同一链在 QuAPI 稠密表示与 Fock 直算下是干净的二阶/四阶。两轮排查把
根源精确定位到**接口测度的耦合结构**（见 IFGTEMPO/docs/
IF_DERIVATION_AND_PACKAGE_LOGIC.md 第 3、9 节）：

- **IFGTEMPO/稠密约定（正确）**：每个 U_imp–U_bath 接口是一个独立恒等
  resolution，其测度是**交叉耦合**——`M = exp[−ā_IF·x_K − x̄_K·a_IF]`，
  即 IF 扇区（bath 响应）与 K 扇区（杂质轨迹）使用**分离的生成器**
  （`IF(t,±,bar/eta)` 与 `IMP(t,±,bar/eta)`），由该测度逐接口互绑。
  bath 响应在解析上是对 K 轨迹的有效双线性（IF 场高斯消元后
  `exp(−x̄·Δ·x)`），且**轨迹变量上没有额外自耦合测度**。
- **GTEMPO 单扇区约定（现状）**：每切片仅一套生成器，闭合 = 显式折返
  双线性 `exp(−ā⁻a⁺)` + 每切片隐含自耦合测度。该结构在 ε_d = 0 时给出
  与稠密表示逐位相同的 Z（`validation/zidentity.jl` 证书），但在 ε_d ≠ 0 时相位
  落位不同：迹偏离 1 达 O(δt³)（N=1: 0.99952），跨节点共轭腿对的两点
  函数存在真实结构差异（穷举 N=1 全部 8 腿两点函数证实）。

### 修复：`twosector/realtime.jl`（与 IFGTEMPO 对齐的权威基准；已随 2026-11-06 目录清理移除）

`twosector/realtime.jl` 直接采用 IFGTEMPO 的两扇区构造（path 引入姊妹
仓库；其全部底层为 GTEMPO 的 GrassmannMPS 机制），配合本轮加入的
Trotter 调度（`trotter_schedule(order, N, dt)`），对精确连续 Fock 演化
的密度矩阵误差实测：

| 模型 | order 1 | order 2 | order 4 |
|---|---|---|---|
| 两自旋 (ε=0.35,−0.22, U=1.0) | p = 2.00 | p = 2.00 | **p = 4.05** |
| Anderson 型 (ε=0.7,0.7, U=1.5) | p = 2.00 | p = 2.00 | **p = 4.10** |

（两自旋情形 order 1 与 order 2 误差逐位相同——对角初态下两链的约化
密度矩阵恰一致，见 IFGTEMPO `test_trotter_orders.jl` 的注记。）

### 将来集成到 GTEMPO 需要的改动清单

1. **接口变量分离**：实时 Grassmann lattice 的每接口需容纳 IF/K 两套
   生成器（IFGTEMPO 的 `IF(t,±,bar/eta)` + `IMP(t,±,bar/eta)`，每接口
   4 生成器），或等价地在现有 4 腿切片外增加 bath 响应腿。
2. **逐对测度语义**：GTEMPO 的 pair `g_trace` 当前按切片实现自耦合
   测度；需支持/切换为上述交叉耦合测度（IFGTEMPO 的 lattice 排列
   `(ā_IF, x_K, x̄_K, a_IF)` 使相邻 site 对恰为交叉对）。
3. **闭合方式**：显式 `boundarycondition!` 折返（`exp(−ā⁻a⁺)`）应替换为
   恒等核闭合（IFGTEMPO 的零时长恒等核：N+1 个核、无额外折返双线性；
   迹在该约定下逐位精确）。
4. **可观测提取**：`insert_n!`（清零位点张量偶宇称块）带 O(δt) 有效时刻
   偏移；折返跨越 2 腿 `GTerm` 投影（`ntT` 约定，已并入读出实现）是与稠密参考
   一致的约定，应作为 ⟨n(T)⟩/dd-型可观测的标准实现。
5. **多带**：`boundarycondition!` 对 bands > 1 必须逐 band 施加（本轮
   在 benchmark 侧修复，库侧 `sysdynamics` 路径已正确）。

### 关键实验发现（第三轮）：轨迹变量上的过程张量不是自由 Gaussian

`IFDiscretization/debug/gaussianity.jl`（`if_path_values` 的过程张量，
N=2、单模 bath）：穷举全部 4 腿/6 腿构型，`I[S] ≠ c·Pf(Θ[S])`——
114 个构型 I 精确为零（踢算符在占据基对角，强制 α_m=β_m 的数守恒约束；
自由 Gaussian 不满足），另有若干构型 Pfaffian 偏差。**结论**：
"Δ_eff = 轨迹变量上的纯双线性"作为 `exact_fermionic_Δt` 的输出**不可行**
——严格离散 IF 的正确形态是把踢算符**分布到子胞边界的顶点扩展构造**
（与 IFGTEMPO 两扇区/交叉测度在物理上等价，见 strictkernel.jl 的
`_split_cell` 先例：虚时同款踢分布已验证 order 1/2/4/6）。

### 形态 B（消元式单扇区）的定位（经用户指正后的修正结论）

**修正**：此前一轮曾基于 Gaussianity 判决实验给出"轨迹变量上的纯双线性
Δ_eff 数学上不可行"的结论——**该结论有误，现予撤回**。判决实验
（`IFDiscretization/debug/gaussianity.jl`）实际证明的是：**占据基的
path values（`if_path_values` 的 I[S]）与相干态多项式的 Pfaffian 系数
不是同一个数**（path values 携带踢算符的对角 δ 约束：114 个构型精确为
零）——这是两类表示的差别，不构成对 Δ_eff 双线性形式的否定。

**正确结论**：消去 IF 变量是精确的高斯积分，不损失信息：

    Z = Pf(G) · ∫Dx exp(xᵀΔ_eff x) · K[x] · ρ₀[x],   Δ_eff = −½PᵀG⁻¹P

合并收缩中，踢门的对角结构使 Δ_eff 的"额外分量"自然无贡献（与用户的
合并逻辑一致）；IF 侧包含比两扇区约束 IF 更多的分量完全没有问题。
**形态 B 在原理上成立**；剩余的未解项是稠密接线的测度/闭合簿记
（详见实施规范一节的当前状态），属于可定位的工程问题。

### 形态 B（消元式单扇区）实施规范（稠密端验证已完成，2026-11-04）

> **数学原理的完整推导见 [`theory/MATHEMATICS.md`](theory/MATHEMATICS.md)**
>（两扇区多项式 → 源嵌入与高斯恒等式 → Δ_eff 的反对称形式 → 投影收缩
> 与反周期迹闭合 → GTEMPO 隐式 resolution/反测度/graded 读出符号的
> 对接，每条公式附验证入口）。本节保留规范、锁定约定与结果。

目标：在只有一套杂质轨迹变量的 lattice 上，用"消元后的有效核"替换
"接口交叉测度 + IF 因子"，达到与 IFGTEMPO 两扇区 / Fock 离散链逐位一致，
为将来集成到 GTEMPO 单扇区 lattice 做准备。稠密验证在
`validation/formb_dense.jl`（保留，兼作活跃验证链的共享稠密机器）；
GMPS 移植（原 `src/FormB.jl`）已于 2026-11-07 随早期代码清理删除。
`FORMB_QUICK=1` 为快测档。

**公式（单自旋记号；消元是精确恒等式）**：

    两扇区：P(a,x) = ν·exp(½aᵀGa) · ∏_t exp(−ā_t x_t − x̄_t a_t) · K[x] · M₀
    消元后：P_elim(x) = ν·Pf(G) · exp(Σ_{λ<λ'} D[λ,λ'] θ_λ θ_λ') · K[x] · M₀
    D[λ,λ'] = −½(G⁻¹[i,k] − G⁻¹[k,i])   （i,k = λ,λ' 的 IF 槽位，源 j = −partner）

即 README 早先的 `Δ_eff = −½PᵀG⁻¹P` 的反对称多项式形式。D 只作用在接口
轨迹腿上（层 0 与开放腿无源）；归一化 = ν·Pf(G)（ν = ZN/Z bath 配分）。
高斯积分恒等式 ∫exp(½aᵀGa + jᵀa) = Pf(G)·exp(−½jᵀG⁻¹j) 已用 n=2 直算
核验，暴力 Berezin 消元与解析式逐系数一致（≤2e-16，全部 order/N）。

**锁定约定**（稠密表示下的经验标定，全部由 Z=1 + ρ=Fock 双重锁定）：

1. **收缩 = 纯系数投影**：逐生成器提取系数（∫dθ θf → f），**无定向
   符号**（每 monomial 中生成元至多出现一次，投影与次序无关）；
2. **接口测度在消元后彻底删除**（交叉测度已被 Δ_eff 吸收）；layer 0 的
   M₀ resolution 保留；
3. **迹闭合 = 每 flavor 独立 resolution `(1 + x̄⁺_{N+1}x⁻_{N+1})`**——
   系数 **+1**（反周期边界 ⟨−x̄| 的翻转并入闭合符号；−1 给出周期迹，
   Z ≠ 1）。单自旋与两自旋同符号；
4. kernel 单项式符号照搬 `impurity_kernel_gmps` 配方
   （coeff = (−1)^{C(ket,2)}·K[bra+1,ket+1]，local 序 (bar…,eta…)）；
5. lattice 布局 = `timelocallattice`（单自旋取单 flavor 投影；层层数 =
   `length(cells)`，order 4 每宏步 3 个 subcell ⇒ 3N 层）；
6. 同因子内不同恒等/pair 必须拆成**独立因子**再相乘（exp 的 pair/行分解
   ——合并成单因子表会丢失因子间的组合，这是首轮"全零"错误的根源）。

**验证结果**（`formb_dense.jl` 四步，对照 = 同 schedule 的 Fock 离散链
`Tr_bath[U(ρ₀⊗ρ_B)U†]`，U = A(link_{Nc+1})·Q_{Nc}·A(link_{Nc})···Q_1·A_1）：

| Step | 内容 | 结果 |
|---|---|---|
| A | 稠密**两扇区**副本 vs Fock | Z−1 ≈ 2.6e-18，ρ-err ≈ 1.1e-16 |
| B | 暴力 ∫Da vs ν·Pf·exp(D)（系数级） | max|Δcoeff| ≤ 2.2e-16，elim_sign=−1 |
| C | **消元式单扇区**收缩 vs Fock（单自旋） | orders 1/2 × N=1..4 + order 4 × N=1 全部 ≤1e-16（见下） |
| D | 消元式**两自旋** vs Fock vs IFGTEMPO GMPS | ρ 三方一致 ≤1e-15（`gtempo_density_matrix_direct`），Z ≈ 1 |

即**消元后的单套变量 lattice 与两扇区/GMPS/Fock 三方逐位对齐**——形态 B
从原理到稠密接线全部完成。⟨n(T)⟩ 直接取 ρ(T) 的对角（一致到机器精度）。

**规模边界（dense 的结构性限制）**：消元后 L = 4+4·nl+2（单自旋）/
12+8·nl（两自旋），nl = `length(cells)`。order 4（nl=3N）下两自旋 N=1
已是 2^36、单自旋 N=2 是 2^30——dense 验证覆盖 orders 1/2 × N≤4 与
order 4 × N=1（单自旋）；order 4 的完整标度验证以 `twosector/`（GMPS）
为准，消元式的 MPS 端实现不受此限制（系数以 GTerm/row 因子稀疏承载）。

**将来集成到 GTEMPO（MPS 端）的改动清单**（在原有 5 项清单之上的
形态 B 专用增量）：

1. `Δ_eff` 适配器：把 `D[λ,λ']`（G⁻¹ 反对称部分 × 源符号）实现为
   `AbstractRealCorrelationFunction`；行因子化 `1 + θ_λ·Σ D θ_λ'` 与
   `get_if_gmps` 的逐行 partialmpo 同构，可直接复用；
2. 归一化：`ν·Pf(G)` 以 scaling 承载（替代 exact 核 prefactor）；
3. pair 测度语义：接口对**无测度**（纯投影收缩）、layer 0 对保留
   resolution、开放对用 `+` 号 resolution——GTEMPO 的 pair `g_trace`
   若强制 `exp(−x̄x)`，需以反测度 GTerm（`exp(+x̄x)` 型）补偿开放对；
4. 可观测：ρ(T) 从开放腿系数提取（`impurity_kernel_gmps` 逆配方）；
   ⟨n⟩ 用 ntT 的开放腿投影约定；
5. order 4：MPS 端按 3N 层铺设（同 twosector），无 dense 的规模问题。

### 形态 B 的 GTEMPO GMPS 移植（已完成 2026-11-04；代码已于 2026-11-07 清理删除）

> 历史记录：本节的实现 `src/FormB.jl` 与验证 `validation/verify_formb.jl`/
> `scaling_formb.jl` 已删除（被紧凑消元路线 `compact_delta`/`exactif_native_compact`
> 取代）；数学原理文档 `theory/MATHEMATICS.*` 与稠密验证 `validation/formb_dense.jl`
> 保留。标度图 `formb_scaling_*` 属旧调度（非回文序）记录，随脚本删除，
> 勘误见"order-4 调度修正"一节。

稠密验证锁定的装配原样移植到 GTEMPO 的 GrassmannMPS 机制
（原 `src/FormB.jl` + 对齐验证 `validation/verify_formb.jl`）：

- **lattice**：IMP-only 紧致格（timelocallattice 变量序滤掉 :IF，相邻
  pair 结构保持；order 4 时层数 = 3N）。
- **因子**：(Δ_eff+反测度, bwd, ρ₀+前向) 三个非平凡 GMPS——Δ_eff 行
  因子用 `partialmpo × mult`（get_if_gmps 同款），K 核两 flavor 用
  `impurity_kernel_gmps`（单 flavor 对角 kernel 以双线性 partialmpo
  直构）。
- **收缩**：Z = `l_LL × GrassmannTransferMatrix × scalar`（GTEMPO
  pair transfer 的隐式 g_trace 即 dense 里的 resolution 语义）；ρ(T)
  复用 `gtempo_open_kernel_direct`/`gtempo_density_matrix_direct`
  （单 flavor 需自写 2 开放腿版 `open_kernel_single`——gtempo 版硬编码
  4 条开放腿）。
- **GTEMPO 隐式 g_trace 与 dense 约定的对接**（标定结果）：
  ① 接口对的隐式 resolution 已被 Δ_eff 吸收 → 逐接口对乘幂零反测度
  `(1 − x_t·x̄_t)`（系数 **−1**，(1−x̄x)(1+x̄x)=1 精确抵消）；layer 0 /
  开放对的隐式 resolution 恰为 dense 的 M₀ / 反周期闭合，不显式放因子；
  ② 因子排列 (Δ, bwd, ρ₀+前向)（bwd 置中）与开放系数提取的 graded
  符号联合标定：单 flavor 直接映射；两 flavor 乘
  `diag(1,1,1,−1)`（仅 |↑↓⟩ 对角扇区差一个费米子符号）。
- **对齐结果**（`verify_formb.jl`，对照 = 同 schedule Fock 离散链 +
  IFGTEMPO GMPS；分块运行 `ARGS = calib | single_exact | single_trunc |
  two_exact | two_trunc` 限单进程峰值内存）：
  单自旋精确档（orders 1/2 × N=1..4 + order 4 × N=1，notrunc）
  Z ≈ 1、ρ 逐位一致 ≤1.5e-14；单自旋截断档（orders 1/2 × N=6，
  `truncdimcutoff(D=512, ϵ=1e-12)` 同 twosector 生产用法）≤2.8e-14，
  bond 压至 8；两自旋精确档（orders 1/2 × N=1,2）三方 ≤7.2e-14；
  两自旋截断档（orders 1/2 × N=3）三方 ≤1.7e-12。
- **Trotter 阶数标度**（`validation/scaling_formb.jl`，图在 `validation/`）：对精确
  连续 Fock 演化的密度矩阵误差——单自旋（T=0.4）orders 1/2 → p=2.00、
  order 4 → p=3.99；两自旋（T=0.34 + twosector 同款占据初态
  diag(0,0.3,0.2,0.5)）orders 1/2 → p=2.00、order 4 → p=4.16
  （比值 16.5/17.3 与 twosector 发表值一致；order 4 平台行为依赖初态，
  |0⟩⟨0| 初态在固定 T 下比值 ~10-12）。图上叠加 IFGTEMPO 两扇区参照点
  ——与 FormB 曲线逐位重合（对齐的直接可视化）。注意：本目录的对齐
  全部在**实时（Keldysh）轴**；虚轴严格核属 IFDiscretization 工作线。
- **移植中的两个教训**（防复发）：
  ① ν·Pf(G) 必须用 **Parlett–Reid 数值 Pfaffian**（O(nl³)）：早期
  pair 因子稠密展开既丢 monomial 合并的排列符号（order 4 偏 8.3e-6，
  解析完美匹配求和判定），又按 2^{4nl} 耗内存（N=8 截断档 = 2³²×16B
  ≈ 54GB → **内存爆机**）；
  ② GTEMPO 的测度（相干态重叠 ∏(1+c̄c)）在 **sysdynamics**（
  `generalimpurity/sysdynamics.jl` 的 overlap-removal GTerms）里施加，
  不在 pair-integration 层——bench 移植走 IFGTEMPO 收缩路径（g_trace
  隐式 resolution，反测度系数 −1）；将来集成到 GTEMPO 单扇区管线时，
  Δ_eff 适配器与反测度补偿应接在 sysdynamics 的测度语义上。
- **已知边界**：①初态密度须为对角（偶宇称）——`impurity_kernel_gmps`
  的电荷零分级构造不携带奇宇称 kernel 单项式（IFGTEMPO 全部测试同为
  对角初态类）；②单 flavor 仅支持对角 kernel（踢门/对角初态，已覆盖
  U=0 类验证）；③精确档 bond 按 4^N 增长（Δ_eff 的 G⁻¹ 非局域性），
  大 N 走截断档。

### 管线现状（2026-11-07 清理后的最终接口）

- **生产路线 = 紧凑消元管线**：`exactif_native_compact`（`compact_delta` +
  `sysdynamics_exact!`，orders 1/2/4）——机器精度对齐 IFGTEMPO 与 Fock 离散链，
  真阶标度 2/2/4（见下方"已验证"与"order-4 调度修正"）；双时 gt/lt 走
  `gtlt_compact`，ρ(T) 走 `open_rho_native`。
- **`exactif_native_quapi`**（orders 1/2，QuAPI-Δ 数据源）保留作回归
  （verify_align.jl quapi 档；缓存 bug 修复后逐位/机器精度）。
- **`exact_Δt` 适配器 + `exactif_pipeline`/`frozen_pipeline`** 仅保留作
  **Z 证书对照**（`validation/zidentity.jl`）：ε_d = 0 时严格核
  Z·prefactor == 1 全 order 机器精度，冻结核偏差 O(δt)——区分严格/冻结的
  最锐利证书。
- 早期路线（contour_G 路线的 `exact_delta`/`exactif_native_pipeline`、形态 B
  GMPS 移植 `FormB.jl`）已于 2026-11-07 删除；两扇区交叉验证走 IFGTEMPO
  `build_impurity_path`。将来的 GTEMPO 库内集成在两扇区构造（接口变量分离 +
  交叉测度，lattice 层改动大）与消元式单扇区（一套变量 + Δ_eff 核，
  与现有单扇区管线结构最接近，即本目录紧凑路线）之间选择。


## GTEMPO 原生管线（NativeExact）——开发记录（2026-11-05~06）

> 本节为开发过程记录。其中提到的 contour_G 路线（`exact_delta`/
> `exactif_native_pipeline`）与判决实验脚本（`native_probe.jl`/
> `restrictd_probe.jl`/`verify_native.jl`/`calibrate_native.jl`）已于
> 2026-11-07 删除，最终接口以文件表与"管线现状"一节为准。

把形态 B 的消元式构造搬到 GTEMPO 原生接口（`src/NativeExact.jl`，IFGTEMPO-free）：

1. **lattice**：`GrassmannLattice(N = nl+1, δt, contour=:real, bands)`——只含杂质
   Grassmann 变量；node 1 = t=0（ρ₀ 层）、node 2..nl+1 = 接口（D 源所在）、
   node k = 开放端（物理 T：折返 + ρ(T) 读出）。
2. **CorrelationFunction**：`ExactDelta <: QuAPI.AbstractRealCorrelationFunction`
   封装消元核 D = −½(G⁻¹[i,k]−G⁻¹[k,i])（两扇区 bath Gaussian 的 Ginv 反对称化，
   per-(branch,node) pair 布局 + Schur；`contour_G`），ExactIFKernel 读序约定 +
   Keldysh 跨支符号，**直接由 GTEMPO `hybriddynamics`/`partialif_hybrid` 行扫
   消费**。判决实验（`validation/restrictd_probe.jl`）：**消元核 D 本身就只含
   bar→eta 有向对（0 个同 kind 单项式）——行扫完整承载，机制成立**。
3. **`sysdynamics_exact!`**：调度门（`sysdynamics_factors!` 同款全窗口，
   τ=0 恒等门保留）+ 接口分辨率反测度（(1+cc·x̄x)，FormB 同款幂零抵消）+
   `boundarycondition!` + `systhermalstate!`。
4. **观测量**：`integrate`/Z 证书原生可用；ρ(T) 读出见下。

### 已验证

- **V=0（自由杂质）Z 证书：机器精度通过**（真实 GMPS 机器，
  `integrate(lat, mpsK, mpsIs…)·norm == 1`，3.3e-15）——骨架、调度、
  折叠、ρ₀ 全部正确。
- **稠密判决**（`validation/native_probe.jl`）：复合多项式（ExactDelta IF ×
  native 骨架：门核 + bulk 链接 + ρ₀ + node-1 分辨率 + 折返）在
  **开放腿提取约定**下与 Fock 离散链逐位一致（V=0 与 ε_d≠0 双重通过；
  `res=boundary + fold=neg + nodek 分辨率 + norm=Pf` 组合，Z 与 ρ 双证书）。

### 接口矛盾（进行中）

- GTEMPO 原生观测量（`occupation`/`gf` 的 GTerm 插入 + 全收缩）读的是
  **另一种 T 端约定**：同一复合多项式，插入式读出在 V=0 自由情形给出
  0.287（应为 0）；开放腿对 `(bar⁺(k), eta⁻(k))` 插入给出 0.713 = 1−0.287。
  两种约定不能腿对腿互迁（原生 node k 4 腿 = 折返+开放对+额外分辨率结构，
  compact 开放端 2 腿——FormB 清单因此规划"ρ(T) 从开放腿系数提取"）。
- `src/NativeExact.jl` 已含单带开放腿 ρ(T) 读出（右扫 `update_pair_right`
  + 边界块吸收 + 折返 `g_trace_phy` 收缩 + sector 直读），全链可运行、
  ρ 呈正确对角结构，但约定常数未标定（V=0 给 diag(0.452, 0.302)，
  比例不随 sc/符号改变 → 腿/扇区映射层面需系统性标定，同 FormB 当年
  Step A/C 的流程）。下一阶段工作：完成该标定 → 两带推广 → Trotter 标度。

### ★ 打通（2026-11-06）：QuAPI-Δ 数据源配方全绿

`exactif_native_quapi`（src/NativeExact.jl）：**QuAPI 严格核**
（`exact_fermionic_Δt_cells`，IFGTEMPO-free）作 Δ 数据源 + 补零 pair 网格表
（`trotter_factors` 移植：每 gap 一个 Trotter 因子，点数 = gap 数+1）+
`sysdynamics_exact!`（**无反测度**：QuAPI-Δ 相对测度+单位链骨架定义，
隐式 g_trace 正是其所要）+ 折返腿 GTerm 插入读出。

一致性（εd=0.7 单模 bath，热态 ρ₀）：

| order | N | n(T)（管线） | n(T)（Fock） | Δn |
|---|---|---|---|---|
| 2 | 1 | 0.28787724 | 0.28788771 | 1.1e-5 |
| 1 | 2 | 0.28789854 | 0.28788834 | 1.0e-5 |
| 2 | 2 | 0.28789213 | 0.28788834 | 3.8e-6 |

Δn ~ 1e-5 = ntT 插入约定的已知地板（与旧 exactif_pipeline 同款）。
εd=0 时 Z 恒等式机器精度（1.0000000000）。

**标定教训（防复发）**：①Δ 的原始块 = partialif 行扫的物理指数系数，
**勿加跨支符号**（加了会引入 O(δt) 系统偏差——ExactIFKernel 的内嵌符号
与其 Δ 布局配对，换布局后净符号约定不同）；②cell/gate 表必须用**补零
pair 网格**（native_schedule 的紧凑表与 QuAPI 点基不配对 → O(δt) 假失败）；
③对照的 Fock 参照与管线必须**同初态**（热态 vs 纯态错配曾造成一系列
假阴性）、ρ₀ 只能插一次（sysdynamics_exact! 已含）。

**Order 4 的 Δ 数据源限制（已定位）**：QuAPI 的点基 Δ 无法承载 order 4——
Suzuki S4 的杂质门插在 bath cell **之间**（gap 3、5），而 QuAPI-Δ 的 cell 链
点对点直连（内部杂质线假设 cell 间无门）。Order 2 能对齐是 Strang 结构特例
（单 cell、门只在两端、外段补零 cell 对齐）。Order 4 的原生管线需改用
**紧凑消元 Δ**（FormB 路线，contour_G 的紧凑布线版）——其正确性已由基线的
twosector/FormB 标度独立验证（p = 4.05-4.16）；order 4 N=1 在 D=256 下
Δn = 5.4e-4（D 无关 = 结构性，V=0 机器精度 5.6e-16 确认骨架无错）。
另：N=4 以上出现 Julia 原生崩溃（深层栈，非内存）——待查。

**紧凑消元 Δ 的进展（2026-11-06）**：`compact_delta`（src/NativeExact.jl）
= IFGTEMPO `_analytic_contour_G` 的逐字移植（compact 布线：共享 trace pair +
trace_bra 翻转 + 热边 + pair 测度；Schur 消 bath）。**G 装配已逐位验证**：
vs IFGTEMPO influence_matrix_G，order 2（4×4）与 order 4 含负 cell（12×12）
max|ΔG| = 0.000e+00。期间修复 `_coupled_h` 的**共轭约定 bug**：
h[imp,bath] = conj(V)、h[bath,imp] = V（V 复数时旧代码非厄米——这也解释了
contour_G 路线 V≠0 的全部失败：V=0/实 V 不受影响）。
剩余工作 = 紧凑 Δ 在原生 lattice 骨架上的标定：FormB 紧凑格的 pair 序为
(bar, eta) 而原生 lattice 为 (eta, bar)，隐式测度费米符号不同（接口需
cc = −1 类反测度 + node-1 修正 (1 − 2x̄x) 量级），且 far-end 结构不同
（FormB 交叉对 vs 原生全节点+折返）——需按 FormB Step A/C 流程重标定。
Order 4 物理结果已由基线 FormB 路线独立验证（p = 4.05-4.16，eae55d4）。

**紧凑路线标定补充（2026-11-06，commit c4c8a6f+）**：裁决网格（修复 cc 静默失效 bug 后）命中
**counter_c = −1、sgn = −1、无 siteflip、无 node-1 修正**：
- order 2 N=1：n vs Fock **3e-7**（essentially exact，优于 QuAPI 路线的 1e-5 地板）
- order 2 N=2：**3.8e-6** PASS
- order 4 N=1：**1.57e-4**（0.05% 相对；D=384/512 逐位同 = 非截断；V=0 机器精度 = 骨架无错）

Order 4 残余的定位实验：全配对 Δ（≡ FormB 已验证 D，含同 kind 块）× 我的骨架 → 8.4e-3（更差）；
(bar,eta) 块版 × 我的骨架 → 1.57e-4。两者都未达地板 ⇒ 紧凑 Δ 与骨架的约定**成对锁定**：
FormB D ↔ FormB 骨架（M₀ + 交叉闭合、无折返）；我的 (bar,eta) 块 Δ ↔ 我的骨架（分辨率+折返）。
Order 4 的剩余差异 = 同 kind 块在原生骨架上的正确配对方式（fold/交叉闭合与 (eta⁺,eta⁻)/(bar⁺,bar⁻)
对的符号耦合）——需要稠密-vs-GMPS 多项式逐项 diff 定位。Order 4 物理结论由基线 FormB 路线独立背书
（p = 4.05-4.16）。另：N≥2 部分案例有 Julia 原生崩溃（深层栈，非内存），待查。

### ★ 打通（2026-11-06 晚）：order 4 残余根因 = 门传播子缓存键不含时长——全部对齐

刀口实验链（`verify_compact.jl` / 判决记录）：
1. 残余 ∝ **V²**（V→V/10 时 ×100）——bath 响应侧；
2. **εd = 0（门恒等）时残余完全消失**（order 2 N=2 与 order 4 N=1 都是机器精度）；
   结合 V=0 机器精度 ⇒ 错误只在 **D×门干涉**——而 V=0 时 n(T)=n₀ 对任意门时长不变、
   εd=0 时门为恒等，两条既有证书对"门时长正确性"**双双失明**；
3. 同 kind GA 项在全部 order **恒为零**（compact G 的 U(1) 结构：G⁻¹ 反对称化后只有
   (bar,eta) 块；数值验证 0.000e+00）——"同 kind 配对"假说排除（此前 8.4e-3 的
   "全配对"实验系 delta_eff_factor 的 bar→bar 直连接线错误所致混淆，该探针已移除）；
4. order 2 N=1（links = [δt/2, δt/2] 全同）机器精确而 N=2（[δt/2, δt, δt/2]）出错
   ⇒ 嫌疑锁定在**门→gap 的映射**：`_propagator_sparsegmps` 的缓存键
   `(branch, t, brel, krel)` 不含门时长（其契约是均匀 δt 管线"每步门矩阵相同"），
   `_scheduled_gates!` 对调度门（各 gap 时长不同）共享缓存 → 同支后续 gap 静默复用
   第一个 gap 的传播子张量（时长错、位置对）。

判决（无缓存重跑）：order 2 N=2 7.8e-7 → **3.4e-15**；order 4 N=1 1.57e-4 → **2.4e-15**。
修复 = `_scheduled_gates!` 传 `cache=nothing`（含防复发注释）。连带效应：
- **order 4 N=2,3 全绿**（truncdimcutoff D=512, ϵ=1e-12：1.9e-15 / 3.8e-14）——
  原"N≥2 Julia 原生崩溃"随缓存修复一并消失（NoTruncation 下 order 4 N≥2 的键
  4^nl 增长本就会爆，截断档是正确用法）；
- **quapi 路线的 1e-5 地板也是同一 bug**（恒等门 gap 被上一 gap 的门张量替换）：
  order 2 N=1 1.1e-5 → **0.0**（逐位）、N=2 3.8e-6 → **2.9e-15**。

**最终验证状态（全部机器精度级）**：

| 判据 | 结果 |
|---|---|
| 离散链对齐（n vs Fock 离散 Trotter 链） | orders 1/2 × N=1..3 + order 4 × N=1（notrunc）≤3.8e-15；order 4 × N=2,3（D=512）≤3.8e-14 |
| Z 证书（εd=0，Z = integrate·ν·Pf） | ≤1.6e-14（norm 已含 ν = bath_norm，`normdelta`） |
| V=0 证书（n = n₀ 热态期望） | ≤5.6e-17 |
| **IFGTEMPO 直接对齐**（两自旋 U=0 解耦恒等式：相同 bath + ρth⊗ρth ⇒ 双自旋 = 单自旋⊗自身） | orders 2/4 × N=1,2 全部 **≤2.4e-14** |
| **Trotter 标度 vs 精确连续演化**（单模 4×4 精确对角化） | order 1: p=2.01、order 2: p=2.01（对角初态下与 order 1 逐位同链）、**order 4: p=3.99**（比值 16.0，同 twosector/FormB/IFGTEMPO 发表值） |

验证/标度/对齐入口（GTEMPO 根目录）：

```
julia --project=. benchmark/exactif/validation/verify_compact.jl          # 精确档矩阵 + 双证书
julia --project=. benchmark/exactif/validation/verify_compact.jl trunc    # 截断档（order 4 N=2,3）
julia --project=. benchmark/exactif/validation/scaling_compact_models.jl  # p = 2/2/4 标度
julia --project=. benchmark/exactif/validation/verify_align.jl            # IFGTEMPO 直接对齐
julia --project=. benchmark/exactif/validation/verify_align.jl quapi      # quapi 路线回归
```

读出约定：`occupation`（insert_n）带 O(δt) 有效时刻偏移，不可用于此管线；标准读出 =
`ntT` 同款**折返点自由腿对 (a⁺_k, ā⁻_k) 的 2 腿 GTerm 投影**（见 verify_compact.jl）。

### ★★ order-4 调度修正（2026-11-06 深夜）：真四阶 = 回文组合 S2(a)·S2(b)·S2(a)

两模型 error scaling（`validation/scaling_compact_models.jl`，Toulouse + 单模 bath Anderson
U=1.5 bands=2，诊断基准同款参数 T=0.4/β=1.0/ω=1.0/α=0.5）暴露：**实现中的 order-4
调度是非回文组合 `S2(a)²·S2(b)`**（cells=[a,a,b]、links=[a/2,a,c,b/2]）——算符序
不对称 ⇒ 奇幂误差项不恒零 ⇒ 携带 O(δt³)。判决链：
1. Anderson 强耦合（α=0.5）下 order-4 比值 12.7→12.0→10.8→9.8→9.0 **递减趋向 8**
   （p≈3.17），N=16/32 纯 16 维链直算同型——排除管线/截断（D=512/1024/2048 逐位同）；
2. Suzuki 权重 a=1/(2−∛2)、b=1−2a 的解条件 2a³+b³=0 是**回文序** S2(a)S2(b)S2(a)
   的三阶消失条件——非回文序不满足；
3. 回文调度（合并链 links=[a/2,c,c,(a,c,c)×(N−1),a/2]、cells=[a,b,a]×N）纯链验证：
   比值 14.4→15.6→15.9→16.0→16.00，强耦合 Anderson 与 U=0 对照**双恢复 16×**。

修正（共享 bug，两处同步）：`NativeExact.native_schedule`/`trotter_factors` 与
**IFGTEMPO `trotter_schedule`**（ImpurityPath.jl）的 order-4 全部改为回文组合。
修正后（紧凑消元路线，`scaling_compact_models.jl`）：

| 模型 | order 1 | order 2 | order 4 | order-4 每减半比值 |
|---|---|---|---|---|
| Toulouse（bands=1, U=0） | p=2.01 | p=2.01 | **p=3.94** | 14.5→15.6→15.9→16 |
| Anderson（bands=2, U=1.5, μ=0.7） | p=1.98 | p=1.98 | **p=3.94** | 14.4→15.6→15.9→16 |

（order 1 与 order 2 逐位同链 = 对角 H_imp + 对角热态初态的已知恒等式；order-4 的
p 拟合 3.94 是 N=1..8 四点拟合仍含预渐近段，比值收敛到 16 为真四阶判据，纯链
N=32 达 16.00。机器级预检：两模型 × orders 2/4 × N=1,2 管线 vs 离散链 ≤1.2e-14。）

**渐近定点验证（N=16）**：order 4 比值 vs N=8 = Anderson **15.97**、Toulouse **16.12**，
渐近段（N=2..16）拟合 p = 3.984 / 3.989——真四阶。粗端说明：Anderson 图右端
（δt=0.4）order 4 只比 order 1 好 1.5×（4.17e-5 vs 6.24e-5）——Suzuki 负 cell
误差常数大，粗端预渐近属预期；优势随 δt 减半 1.5×→5.6×→22×→88×，由斜率兑现。

**IFGTEMPO 直接同模型对齐升级（`validation/verify_align_u15.jl`）**：不再依赖 U=0 解耦
恒等式——IFGTEMPO 两自旋路线直接跑同一交互模型（TwoSpinImpurityModel(0.7,0.7,1.5)
+ 相同单模 bath ×2 + 相同热态 4×4 初态），双可观测量 n_up 与**双占据 ⟨n₁n₂⟩**
（band-1+band-2 双折返腿 GTerm 乘积投影，先用 Fock 链机器级校验 ≤1.1e-14）：
orders 2/4 × N=1,2 共 8 组全部 ≤2.2e-14 PASS。

**波及面（重要）**：twosector/`scaling_formb` 此前发表的 "order 4 → p=4.05-4.16、
比值 16.5/17.3" 是非回文调度在弱耦合参数下的**交叉区假象**（δt³ 项系数小但非零，
比值越过 16 后终将滑向 8；IFGTEMPO `test_trotter_orders.jl` 的 16.5/17.3 同源——
修正后该测试比值 15.96/15.93，9/9 仍过）。相关标度图（twosector、FormB、
singlemodebath/toulouse 诊断）建议用修正后调度重跑更新。验证全绿：verify_compact
两档 + IFGTEMPO 直接对齐（orders 2/4 × N=1,2 ≤2.7e-14）。

图：`toulouse/realtime.{png,csv}` 与 `singlemodebath/realtime.{png,csv}`（含 O(δt)/O(δt²)/O(δt⁴)
参考线与拟合斜率标注）。

### 判决实验（结论已锁定进实现，脚本于 2026-11-07 删除）

- Q1（消元核 D 的 bar→eta 块充分性）：**PASS**——D 只含 bar→eta 有向对
  （0 个同 kind 单项式），行扫完整承载。
- Q2（native 腿集两扇区消元 × 原生骨架的稠密标定，旋钮
  trace/fold/res/nodek/norm）：判定开放腿约定下的精确组合
  （`res=boundary + fold=neg + nodek 分辨率 + norm=Pf`；注意其对 GTEMPO
  g_trace 语义的稠密仿真与真实机器有 δt³ 级偏差，GMPS 判决以真实机器为准）。

## 相关修改（本轮）


- `GTEMPO/benchmark/exactif/validation/formb_dense.jl`（本轮新增）：形态 B
  稠密系数级验证——IF 变量从两扇区 lattice 的精确高斯消元
  （D = −½(G⁻¹[i,k]−G⁻¹[k,i])，暴力 Berezin 逐系数对照），消元式单扇区
  收缩与 Fock 离散链 / IFGTEMPO GMPS（`gtempo_density_matrix_direct`）
  三方逐位对齐（≤1e-15）；锁定全部收缩约定（投影式积分无定向符号、
  迹闭合 `(1+x̄⁺x⁻)` 系数 +1、逐恒等/pair 独立因子化）。
- `IFGTEMPO`（外部包）：适配当前 GTEMPO API（`GrassmannTensorMap` 移除、
  `mult` 返回元组、`stable_tsvd` → `Z2Tensors.tsvd`、`bonddim` 重命名）；
  `get_k_gmps`/`build_impurity_path`/`influence_matrix_G` 支持 Trotter 调度
  （`trotter_schedule(order, N, dt)`）；新增 Trotter 阶数验证测试。
- `QuAPI`（外部包）：`exact_fermionic_Δt_cells`——严格核的逐 cell 调度版
  （均匀网格逐位一致；支持 Suzuki 负 cell）。
- `GTEMPO`：`test/runtests.jl` 的 `branches` 消歧（ImpurityModelBase 上游
  新导出导致的名称冲突）；`benchmark/exactif/` 全部内容（含本轮
  `twosector/`、`ntT`、多带折返修复、恒等门测度修复与 README）。
