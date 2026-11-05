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
  消元精确积掉。本目录经 `exact_Δt(bath, cells)`（`common.jl`）适配为
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

| 文件 | 内容 |
|---|---|
| `common.jl` | exact_Δt 适配器、Trotter 调度、调度版 sysdynamics、完整管线、ED 参考、标度驱动 |
| `zidentity.jl` | **Z 恒等式证书**（ε_d=0 时 Z·prefactor == 1，全部 order 机器精度） |
| `twosector/realtime.jl` | **两扇区对齐基准（权威标度）**：order 1/2 → 二阶、order 4 → 四阶 |
| `twosector/realtime_*.png` | 对齐基准的标度图（两自旋 / Anderson 型各一张，p=2.00/2.00/4.05-4.11） |
| `form_b/MATHEMATICS.tex/.pdf` | **形态 B 数学原理文档（LaTeX 正式版，8 页）**：两扇区路径积分 → IF 变量高斯消元（Δ_eff 推导）→ 收缩语义/迹闭合 → GTEMPO pair 机制对接（反测度/graded 读出符号），附标度图与验证入口；PDF 亦在 `results/MATHEMATICS.pdf`（xelatex 编译，`ctexart`） |
| `form_b/MATHEMATICS.md` | 数学原理文档的 Markdown 版（与 .tex 同内容；供仓库内直接预览） |
| `form_b/formb_dense.jl` | **形态 B（消元式单扇区）稠密系数级验证**：IF 变量从两扇区 lattice 精确积掉，与 Fock/GMPS 三方对齐（`FORMB_QUICK=1` 快测档） |
| `src/FormB.jl` | **形态 B 的 GTEMPO GMPS 移植**：紧凑 IMP-only lattice 上的 Δ_eff/反测度/K 装配（GTEMPO GrassmannMPS 因子 + pair transfer 收缩） |
| `src/verify_formb.jl` | FormB 移植的对齐验证：单自旋 vs Fock、两自旋 vs Fock vs IFGTEMPO GMPS 三方 |
| `src/scaling_formb.jl` | FormB 的 Trotter 阶数标度（vs 精确连续演化，叠加 IFGTEMPO 参照点；`ARGS = single | twospin | plot`） |
| `results/formb_scaling_*.png/.csv` | FormB 标度图：orders 1/2 → p=2.00、order 4 → p=3.99（单自旋 T=0.4）/ 4.16（两自旋 T=0.34，同 twosector 基准设置） |
| `form_b/validate.jl` | 形态 B 早期框架验证（G⁻¹ bath 响应核；已被 formb_dense.jl 取代，保留作记录） |
| `toulouse/realtime.jl` | Toulouse 模型（U=0, ϵ_d=0.7）的 δt 误差标度（单扇区，见差异记录） |
| `singlemodebath/realtime.jl` | 单模 bath Anderson 模型（μ=0.7, U=1.5）的 δt 误差标度（同上） |

运行（GTEMPO 根目录）：

```
julia --project=. benchmark/exactif/zidentity.jl        # Z 恒等式证书
julia --project=. benchmark/exactif/twosector/realtime.jl   # 两扇区对齐基准（权威，出图）
julia --project=. benchmark/exactif/form_b/formb_dense.jl   # 形态 B 消元验证（稠密，全量）
FORMB_QUICK=1 julia --project=. benchmark/exactif/form_b/formb_dense.jl   # 快测档
julia --project=. benchmark/exactif/src/verify_formb.jl     # FormB GMPS 移植对齐（全量）
FORMB_QUICK=1 julia --project=. benchmark/exactif/src/verify_formb.jl   # 快测档
julia --project=. benchmark/exactif/toulouse/realtime.jl    # 单扇区诊断（见差异记录）
julia --project=. benchmark/exactif/singlemodebath/realtime.jl  # 单扇区诊断（同上）
```

**注意**：`toulouse/` 与 `singlemodebath/` 是 GTEMPO 单扇区管线的**诊断基准**，
受上述闭合约定差异限制，其 gt/lt 系列只表现 O(δt) 收敛（order 1/2）或
不收敛（order 4），**不**是对严格核 2/4 阶的验证；真正的 2/2/4 阶验证
看 `twosector/`（图）与 IFGTEMPO 的 `test_trotter_orders.jl`（测试）。

## 已验证的结论

1. **Z 恒等式证书**（`zidentity.jl`）：ε_d = 0 时
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
   - order 4（Suzuki）：**四阶**（ratio 16.5/17.3）。

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
  与稠密表示逐位相同的 Z（`zidentity.jl` 证书），但在 ε_d ≠ 0 时相位
  落位不同：迹偏离 1 达 O(δt³)（N=1: 0.99952），跨节点共轭腿对的两点
  函数存在真实结构差异（穷举 N=1 全部 8 腿两点函数证实）。

### 修复：`twosector/realtime.jl`（与 IFGTEMPO 对齐的权威基准）

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
   偏移；`ntT`（`common.jl`）的折返跨越 2 腿 `GTerm` 投影是与稠密参考
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

> **数学原理的完整推导见 [`form_b/MATHEMATICS.md`](form_b/MATHEMATICS.md)**
>（两扇区多项式 → 源嵌入与高斯恒等式 → Δ_eff 的反对称形式 → 投影收缩
> 与反周期迹闭合 → GTEMPO 隐式 resolution/反测度/graded 读出符号的
> 对接，每条公式附验证入口）。本节保留规范、锁定约定与结果。

目标：在只有一套杂质轨迹变量的 lattice 上，用"消元后的有效核"替换
"接口交叉测度 + IF 因子"，达到与 IFGTEMPO 两扇区 / Fock 离散链逐位一致，
为将来集成到 GTEMPO 单扇区 lattice 做准备。实现与全部验证在
`form_b/formb_dense.jl`（稠密系数级）与 `src/FormB.jl`（GMPS 移植）；
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

### 形态 B 的 GTEMPO GMPS 移植（`src/FormB.jl`，已完成 2026-11-04）

稠密验证锁定的装配原样移植到 GTEMPO 的 GrassmannMPS 机制
（`src/FormB.jl` + 对齐验证 `src/verify_formb.jl`）：

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
- **Trotter 阶数标度**（`scaling_formb.jl`，图在 `results/`）：对精确
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

### 单扇区管线的现状（保留并标注）

单扇区 exactif 管线（`common.jl` 的 `exactif_pipeline`）保留作为对照：
- Z 恒等式证书（ε_d = 0，全 order 机器精度）仍是最锐利的严格性证书；
- ⟨n(T)⟩（`ntT`）与 Fock 一致到 ~1e-5 地板；
- gt/lt 系列（`gtlt_exactif`）受上述闭合差异限制，只表现 O(δt)
  （order 1/2）与不收敛（order 4 的 gt）；
- 与 IFGTEMPO 完全对齐的**两条已验证路径**：①两扇区构造
  （`twosector/`，GMPS，已实测 2/2/4 阶标度）；②**消元式单扇区**
  （`form_b/formb_dense.jl`，稠密系数级，三方逐位对齐——见
  "形态 B 实施规范"一节）。将来的 GTEMPO 库内集成可在两者之间选择：
  两扇区=接口变量分离+交叉测度（lattice 层改动大），消元式=一套变量+
  Δ_eff 核（`AbstractRealCorrelationFunction` 适配器 + 开放对反测度，
  与现有单扇区管线结构最接近）。

## 相关修改（本轮）


- `GTEMPO/benchmark/exactif/form_b/formb_dense.jl`（本轮新增）：形态 B
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
