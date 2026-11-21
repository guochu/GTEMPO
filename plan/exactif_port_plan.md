# exactif 移植计划（v2.3）：extensions/exactif → 当前 GTEMPO（复用 HybridizationKernel + 统一常数契约 + 类型层级）

> **v2（2026-11-20 用户裁决：统一常数封装）**：Weyl 2 阶 QuAPI 的 Weyl 标量与严格 IF 的
> normdelta 都是"乘在整体收缩上的 Z 常数"——统一存储于 HybridizationKernel，下游
> `sysdynamics` 分派消费。E4 的元组配对返回废除，normdelta 单一存储进 kernel。
> **v2.2（2026-11-20 用户裁决：类型层级取代 scheme 字段）**：`AbstractHybridizationKernel`
> + `FrozenHybridizationKernel`（现行 QuAPI 核改名）/`ExactHybridizationKernel`（对接 IFD
> 严格 IF）两具体类型；**合并接口与否 = Exact 行为完全摸清后再裁**（显式推迟点）。
> **v2.3（2026-11-20 用户裁决：hybriddynamics 算法支持集）**：Exact 核天然破坏时间平移
> 不变性 ⟹ `hybriddynamics` 对 Exact 核只支持 `alg ∈ (PartialIF, PolarIF)`（类型级守卫）。
> **v2.4（2026-11-20 用户定标：对比金标准）**：主流路线物理锚两个——G1 AndersonIM × 单模
> bath、G2 Toulouse × 连续 bath（直接 Trotter vs 严格核装配，order ∈ {1,2} 严格对齐）；
> **范围 = normal bath，BCS 不做**。本修订与 `fockmpo_dress_selfread_plan.md`（合并计划：
> fock_sysdynamics 穿戴 / logscale E2 语义 / fock_initstate! 删除）对齐。

- 日期：2026-11-17（v1）；2026-11-20（v2 统一常数契约 → v2.2 类型层级 → v2.3 hybriddynamics 算法支持集 → v2.4 对比金标准 G1/G2）
- 基线：GTEMPO HEAD `a6672a84`；移植源 = `extensions/exactif/src/NativeExact.jl`（1698 行）+ `validation/`（18 个 CSV + imag_identify/matrixid/opfit 脚本 + theory notes）
- 现状：extension 针对旧栈编写（`Z2Tensors` 截断、旧 `sysdynamics!`/`boundarycondition!`/`systhermalstate!`、`GTEMPO.Defaults.alg_trunc_sys`、旧 `cached_greater/lesser`）——框架大改后失效
- 关联：`if_expmult_backend_plan.md` Phase 1（**前置已落地**——2026-11-19 十七批，`HybridizationKernel{T}` + `ContourBlockIndex` 在产）；`fockmpo_dress_selfread_plan.md`（合并计划：穿戴语义 P1/P3、logscale E2 语义、`fock_initstate!` 删除——F3 映射表与 §2.3 依赖之）；`ifd_backend_switch_plan.md` R5（imag + Exact 已接线，real/mixed 拒绝面 = 本计划关闭）

## 0. 已裁决

- **E1**：ExactIF 核心内核机制移植——`native_schedule`/`trotter_factors`（order ∈ 1,2 的 Trotter 调度）/`compact_delta`/`ExactDelta`/`ExactIFKernel`/`exact_Δt`/`quapi_delta`；核封装改用**新 `HybridizationKernel{T}`**（ContourBlockIndex 键）——**PartialIF（现行或新 partialif.jl）直接消费**。
- **E2**：管线移植——`exactif_pipeline`/`frozen_pipeline`/`exactif_native_compact`/`exactif_native_quapi`/`gtlt_compact`/`open_rho_native` 到新 API（GM `truncate!`/`mult`、MAK 策略、新 sysdynamics/integrate）。
- **E3**：imag_* 研究工具机械移植（偶宇称子集基/探针族/Möbius/design/OMP/CSV I/O，~700 行）+ `validation/` CSV 回归（18 个档案 = 移植正确性锚点）。
- **E4′（v2 修订，取代 v1-E4 元组配对）**：**normdelta 单一存储进 `ExactHybridizationKernel.normdelta`**（带号，v2.2 类型化），不折进 blocks、不元组配对返回。移植管线（研究工具链）与主流管线**同一取数点** = kernel 常数按类型取（Frozen `.coupling_norm` = Σα / Exact `.normdelta`）；`Z = normdelta · det(装配)` 契约不变，常数的家搬进 kernel。imag + Exact 路线（`correlationfunction` 元组解包 → normdelta 透传）已按同一契约在产——v2.2 起该路线容器切 `ExactHybridizationKernel`，两构造器一致。
- **E5**（2026-11-17 用户裁决）：**只考虑等距离散化**（非等距 cells 支持不移植）；**exactif 不支持 order 4 Trotter 分解**——S4 门序列涉及负时间演化，移植范围 = order ∈ {1, 2}，与新 `fock_sysdynamics` 的 order 集对齐。
- **E6（v2.2 定稿，类型层级取代 scheme 字段）**：

  ```julia
  abstract type AbstractHybridizationKernel end

  struct FrozenHybridizationKernel{T<:Number} <: AbstractHybridizationKernel
      blocks::Dict{Tuple{ContourBlockIndex,ContourBlockIndex},Union{Matrix{T},CorrelationMatrix{T}}}
      coupling_norm::Float64     # Σα ≥ 0（Weyl 标量输入；现行守卫原样保留）
      bands::Vector{Int}
      anomalous::Bool
  end   # = 现行 HybridizationKernel 改名（QuAPI 离散化）

  struct ExactHybridizationKernel{T<:Number} <: AbstractHybridizationKernel
      blocks::Dict{Tuple{ContourBlockIndex,ContourBlockIndex},Union{Matrix{T},CorrelationMatrix{T}}}
      normdelta::Float64         # 带号 Z 常数（守卫 = isfinite + imag 尾 ~1e-15 取实）
      bands::Vector{Int}
      anomalous::Bool
  end   # 对接 IFD 严格 IF（exactif_kernel 与 imag+Exact correlationfunction 两构造点）
  ```

  - **改名面**：现行 `HybridizationKernel` 全 call site 改 `FrozenHybridizationKernel`（in-repo 全可控，无 alias）；`HybridizationKernel` 名字退役——**合并接口若裁合一，统一类型回收此名**；export 增 Abstract/两具体型；
  - **共享泛型**：`changebands`/`_band_kernel`/anomalous 键派生（返回同型）；`Base.:+` 按类型分派：Frozen+Frozen = 相加（Σα，现行）；Exact+Exact = 因子化裁决门（见下）；混合 = 抛错；`TK.scalartype` 各一行；
  - **下游消费签名**：`correlationfunction`（返回类型按构造路线）/`fullfactor`/`hybriddynamics`/PartialIF 系/measurements kernel kwarg → `AbstractHybridizationKernel`；
  - **hybriddynamics 算法支持集（v2.3 用户裁决）**：`hybriddynamics(lattice, k::ExactHybridizationKernel, alg)` **只支持 `alg ∈ (PartialIF, PolarIF)`**——Exact 核天然破坏时间平移不变性（逐时刻独立的非平移不变块），ExactTTI 的"核块平移不变指数拟合"前提不成立；`alg = ExactTTI` → `ArgumentError`（**类型级守卫**：hybriddynamics 入口显式拒绝，不依赖运行时 `is_tti` 门顺带抛错；错误信息注明"Exact 核非平移不变，ExactTTI 不适用"）。Frozen 核支持集不变：PartialIF/PolarIF/ExactTTI 三算法照旧；
  - **常数消费分派（类型分派，替代 v2 的 scheme kwarg 设计）**：
    ```julia
    _constant_lscale(k::FrozenHybridizationKernel, lattice, order) =
        order == 2 ? k.coupling_norm/4 * Σ step² : 0.0        # 现行 _fock_weyl_lscale
    _constant_lscale(k::ExactHybridizationKernel, lattice, order) =
        log(abs(k.normdelta))                                  # 号折进链幅：首站 × sign(normdelta)
    ```
  - **sysdynamics 接口**：`kernel::Union{Nothing,AbstractHybridizationKernel}=nothing` kwarg——`kernel !== nothing` → 常数按类型分派；`kernel === nothing` → 现行 `coupling_norm` → Weyl（Frozen 缺省语义）；两者同给且耦合不一致 → 抛错。测量层 8 入口转发 kernel（替换 v2 的 scheme 转发设计）；
  - **构造守卫**：Frozen 保持 `coupling_norm ≥ 0`；Exact = `isfinite(normdelta)`（带号——Δτ 无耦合极限 (−1)^{NK}Z_bath、Δm 显式 (−1)^{n(K+1)}·det(BB)/Z_bath；号约定锚 = IFD 金标准 `test/exact/fermionic_{imagtime,realtime,mixedtime}.jl`）；
  - **多浴加法 = 因子化裁决门（v2.1 修订，取代 v2 初稿的"一律抛错"）**：blocks 相加在 Exact 下**恰好正确**——结构恒等式：联合浴 bath 块对角 ⟹ Schur 补线性叠加，Δ_comb = Δ₁ + Δ₂（imag/real/mixed 同理），IF/比值观测量加法无损；不服从加法的只有 **normdelta**（行列式型常数）。联合构造 vs 单独构造的关系**是否恰好相乘（因子化）待数值裁决**：EDToulouseSolver 的 h₀ = 杂质+浴联合本征求解（λ 谱 ≠ 单浴并集），bath 块经共享杂质侧/源链的交叉项大概率存在 ⟹ 因子化不成立的概率高，但以对拍定谳。**Phase B 裁决步骤**：exactif_kernel 双浴联合构造 vs（blocks 相加 + normdelta 相乘）对拍 IFD 金标准精度。**成立** ⟹ `Base.:+(::Exact, ::Exact)` = blocks 相加 + `normdelta` **相乘**（带号自然处理；双 lead 加法免费获得，无需联合构造接口）；**不成立** ⟹ `Base.:+(::Exact, ::Exact)` 抛错 + `exactif_kernel` 接浴列表联合构造（同一机器喂两浴，接口成本小，裁决前不预建）。恒错禁止项："blocks 相加 + normdelta 相加"（Frozen 语义套用到 Exact）；"允许加法但绝对 Z 静默降级"（比值对、Z 错的 footgun）；
  - **消费语义**：常数进 K 链幅（logscale，E2 累加语义）⟹ Z 与分子同携带——比值观测量严格相消、绝对 Z 与真空形式 n̄_vac 需要（IFD normdelta 语义定谳一致）。

## 1. 探查定谳

**F1【ExactDelta 的本质与新形态】** 旧栈它伪装成 `QuAPI.AbstractRealCorrelationFunction`（getindex 带跨支 Keldysh 负号）喂给旧 `hybriddynamics!` 的 partialif 行扫。字段：`pos`（GrassmannContourIndex → 列号）、`M`（(ā,a) 读矩阵，仅接口节点 2..nl+1 有源）、`G`（两扇区 bath Gaussian）、`ν`/`pf`（归一化 = ν·Pf）。**新形态 = `ExactHybridizationKernel`（v2.2）**：`M` 按 (branch, band) 分块 → `blocks[(ContourBlockIndex(b₁,band₁), ContourBlockIndex(b₂,band₂))] = Matrix`（(ā,a) 味对的子矩阵；conj 相同的块 = 零/缺键）；**`ν·Pf` → `normdelta` 字段（E4′）**。精确消元核**非平移不变** → real/mixed 值表示取 **dense `Matrix{T}`**（不入 CorrelationMatrix）；imag 严格核环上 Toeplitz → 走 `correlationfunction` 现行路线（CorrelationMatrix 块），容器同为 `ExactHybridizationKernel`。

**F2【符号约定的收敛】** 旧 ExactDelta.getindex 内嵌跨支负号（"partialif 行扫读数处补偿回物理块"）；IF 后端计划 W2 已裁「符号处理收进构造期、blocks 存物理终值」——移植时负号**折进 blocks 构造**，ExactDelta.getindex 的补号逻辑随之删除（同一约定收敛，两计划共享）。

**F3【旧 → 新 API 映射表】**（v2：`fock_initstate!` 已随合并计划删除，行 4 更新）

| 旧（extension 依赖） | 新 |
|---|---|
| `Z2Tensors.TruncationScheme`/`TruncateDimCutoff(D,ϵ,addback)` | MAK 策略：`truncdimcutoff(D, ϵ)`（addback 语义并入）/`truncrelerr` |
| `GTEMPO.Defaults.alg_trunc_sys()` | `truncdimcutoff(…)` 显式参数 |
| 旧 `sysdynamics!`/`sysdynamics_exact!`（counter_c=+1, FormB 标定） | 移植自身机制（Trotter 调度自管，**order ∈ 1,2、等距离散**——与 `fock_sysdynamics` 的 order 集对齐；E5） |
| 旧 `boundarycondition!`/`systhermalstate!` | 管线内折叠，语义同合并计划 P3（缺省 `fock_vacuumstate`、real `initial_state` 覆盖；`fock_initstate!` 已删） |
| 旧 `hybriddynamics!` partialif 行扫 | `ExactDelta → HybridizationKernel` → partialif（E1/E4′） |
| 旧 `cached_greater/lesser` | measurements/operators 路线（快照扫掠） |
| 旧 `integrate(lattice, K, Is…)` | 新 `integrate` |
| `mult!(…, TruncateDimCutoff)` | `mult!(…; trunc=truncdimcutoff(…))` |
| 旧管线 Z 常数（pipeline 内乘 norm） | kernel 常数字段按类型直取（E4′——Frozen `.coupling_norm` / Exact `.normdelta`，单一取数点） |

**F4【imag_* 工具族的规模与形态】** `imag_lattice`/`imag_inject`/`ImagProbes`/`imag_mobius`/`imag_zeta`/`imag_h_row`/`imag_g_row`/`imag_design`/`imag_fit_svd`/`imag_omp_path`/`imag_read_rows` 系——全部消费 GMPS + 旧 GM API；机械移植 = API 映射表套用，算法零改动。

**F5【回归锚点】** `validation/` 的 18 个 CSV + `imag_identify.jl`/`imag_matrixid.jl`/`imag_opfit.jl` 验证脚本 = 移植后的逐数值回归基准；`imag_theory_notes.md` = 推导文档随迁。

**F6【PartialIF 直接支持的落点】** ExactIF 核 → `ExactHybridizationKernel` 后：现行 `partial_if`（行扫）与新 `partialif.jl`（W3 `fullfactor`×PartialExp）**都能消费**（PolarIF 同样适用——逐块 SVD 秩一分解任意核；ExactTTI 排除，E6 v2.3 算法支持集）——ExactIF 从"独立 IF 算法"定位转为"**严格核构造器**"（kernel 层），IF 算法层零特例。文件落点 = `src/influencefunctional/exactif.jl`（W1 文件族第三文件）。

**F7【常数消费现状定谳（v2 探查）】** ①imag + Exact：`correlationfunction` 元组解包已把 normdelta 透传进 kernel 常数字段（T5 锚）——消费端**未接线**：`_system_factor_fock` 无条件套 Weyl 公式，拿 normdelta 当 cn 传 = 无意义值（本计划 §2.2 关闭）；②real/mixed + Exact：correlationfunction `ArgumentError` 拒绝面 = 本计划 Phase B 关闭；③BCS + Exact：`exact_bcs_*` 无 normdelta（IFD 刻意留白）——BCS 核保持 Frozen 型 Σα 语义，留白待 BCS 后端立项（不在本计划）；④现行 Frozen 构造守卫 `coupling_norm ≥ 0` 与带号 normdelta 冲突——v2.2 类型层级下自动消解（`ExactHybridizationKernel` 自带 `isfinite` 守卫，Frozen 守卫不动）。

## 2. 设计

### 2.1 核封装（Phase A 核心）

```julia
exactif_kernel(bath, model, N, δt; order, nl, kwargs...)
    -> H::ExactHybridizationKernel
    # compact_delta/native_quapi 两条构造线 → ExactDelta/ExactIFKernel →
    # M 分块（(ā,a) 味对 → ContourBlockIndex 对，跨支负号折入）→ dense Matrix blocks；
    # H.normdelta = ν·Pf（E4′：单一存储，无元组，带号）；
    # bands = kernel 覆盖带号表；anomalous 从键派生
```

- `native_schedule`/`trotter_factors`/`imag_schedule` 移植，**裁剪至 order ∈ {1, 2}**（E5）；`exact_Δt` 的非等距 cells 支持不移植（等距离散，E5）。
- `sysdynamics_exact!`（FormB 标定接口分辨率吸收）按新链格式重接——K 链构造保持在 exactif 机制内（Trotter 调度自管，不复用 `fock_sysdynamics` 的 order 1/2 调度）；链内 ρ₀ 折叠语义同 P3。
- 构造守卫（E6）：`ExactHybridizationKernel` 构造器 = `isfinite(normdelta)` + 实值尾取实（带号允许）；Frozen 型守卫不变。

### 2.2 类型层级 + 统一常数消费落地面（Phase A/B，E6 落地面）

落点（主流侧，非 exactif 机制内）：

1. **类型层级重构**（correlationfunction.jl）：`AbstractHybridizationKernel` + `FrozenHybridizationKernel`（现行 HybridizationKernel 改名，全部 call site 同步）+ `ExactHybridizationKernel`；`changebands`/`_band_kernel`/anomalous 派生改共享泛型（返回同型）；`Base.:+` 按类型分派（Frozen+Frozen 相加 / Exact+Exact 裁决门 / 混合抛错）；`TK.scalartype` 补 Exact 一行；export 三符号；
2. **构造路线归位**：`correlationfunction` 的 imag + Exact（元组解包 normdelta 透传）容器切 `ExactHybridizationKernel`；real/mixed + Exact 拒绝面暂保持（Phase B 关闭改路由）；BCS 路线保持 Frozen 型（normdelta 留白，F7③）；
3. **常数消费分派**（sysdynamics/`_system_factor_fock`）：`kernel::Union{Nothing,AbstractHybridizationKernel}=nothing` kwarg + `_constant_lscale` 双方法（§0 E6）——替换现行无条件 `_fock_weyl_lscale`；`kernel === nothing` 时 `coupling_norm` → Weyl 公式（现行路径原样）；**logscale 累加语义（E2）不变**；号折进链幅（Exact）；
4. **测量层**：8 入口转发 `kernel`（替换 v2 的 scheme 转发设计）；`kernel_coupling_norm` 分派（Frozen → Σα / Exact → normdelta）；
5. **hybriddynamics 算法守卫**（E6 v2.3）：Exact 核 × ExactTTI → `ArgumentError`（类型级，入口显式）；Exact 核 × {PartialIF, PolarIF} 全支持——PartialIF 行扫天然无平移不变前提、PolarIF 逐块 SVD 秩一分解任意核适用（现行"PartialIF/PolarIF 任意 kernel 适用"语义的显式化）。

验收锚：**imag + Exact 端到端 Z**（`correlationfunction` → `sysdynamics(kernel=…)` → `integrate`）vs IFD 金标准 `normdelta · Pf(B)` 机器精度；**frozen 不回归**（Weyl 公式路径 = 现行 203 闸门 + observables golden 逐位）；比值观测量带/不带常数同值（相消实证）；**多浴因子化裁决**（双浴联合构造 vs blocks+ / normdelta× 对拍，E6 裁决门——成立则 `Base.:+` Exact 语义落地为乘法，不成立则拒绝面单测）；**Exact 核 × {PartialIF, PolarIF} 端到端观测量**（vs IFD/ED 参照）+ × ExactTTI 拒绝面（v2.3 支持集实证）。

### 2.3 管线（Phase B 后半）

`exactif_pipeline`/`frozen_pipeline`/`gtlt_compact`/`open_rho_native`/`_right_sweep1`/`_absorb_site2` 按映射表（F3）机械移植；`bath_norm`（ν = ZN/Z 迹比）与 `pfaffian`（Parlett–Reid）原样；Z 常数从 kernel 常数字段按类型取（E4′：Frozen `.coupling_norm` / Exact `.normdelta`）。

### 2.4 imag 工具族（Phase C）

按映射表机械移植（F4 清单全量）；`ImagProbes` 记忆化与族种子机制原样；CSV I/O 原样（`validation/` 档案路径不变，回归直接可比）。

### 2.5 对比金标准（v2.4 用户定标；Phase B 主验收、Phase E 固化）

主流路线的物理定义级锚，两个（**范围 = normal bath；BCS 金标准不做**——见不做清单）。与既有两轴互补：imag Z 证书轴（§2.2，`normdelta·Pf(B)`）+ 研究管线 CSV 回归轴（F5，18 档案）之外的第三条腿——**ExactHybridizationKernel → hybriddynamics(PartialIF/PolarIF) → integrate → 观测量**全链路的物理正确性锚。

- **G1【AndersonIM × 单模 bath——定义级锚】**：AndersonIM 杂质（U ≠ 0 照常——U 走系统因子 K，bath 部分经严格核精确处理，εd 经 kicks 进核，这正是 IF 机器的分工；装配契约 = IFD 金标准已定谳的带踢 Schur，Z ≡ 1 机器精确）+ 离散单模 bath（DiscreteSpectrum 单模）。参照 = **联合 Fock ED 直接 Trotter 门序列演化**（杂质 + 1 bath 模的小维空间，order 1 = Lie / order 2 = Strang）；被检 = 严格核装配同一观测量（最低集：杂质占据 n̄；虚时轴补 Z、n̄(β)）。**判据 = 相同 (N, δt, order) 下机器精度对齐，order 1/2 两档逐对**。
- **G2【Toulouse × 连续 bath】**：ConstImpurityHamiltonian（U = 0 单能级 εd）+ 连续谱 bath（半圆/wide-band 有界谱），GL-M 求积离散化（与 `EDToulouseSolver(nbath)` 的 M 节点**同源**）。参照 = 离散化后 M 模联合 ED 门序列直接 Trotter；被检 = 严格核装配同观测量。判据同 G1。

**对齐前提（硬约束）**：参照门序列与严格核的 cells/kicks 调度**逐 gap 同源**（IFD 金标准 `pair_grid`/`trotter_factors` 同款拆分）——不同调度只保证阶数一致性（δt² 收敛斜率），无逐位对齐。**参照机纪律**：联合 ED 门序列演化 = **测试侧 helper**（`src` 零新增符号——IFD spectral 计划同款纪律；IFD `test/exact/_anchor_utils.jl` 的 ED 门序列机可对照移植；extensions 区不入库，R7 定谳的 ED 参照机按需重建）。

## 3. 阶段与验收

- **Phase 0（前置）**：~~`if_expmult_backend_plan` Phase 1~~ **已落地（2026-11-19 十七批）**。
- **Phase A（类型层级重构 + 常数消费落地面）**：§2.2 全部（改名面/分派/kernel kwarg/测量转发/imag+Exact 容器切换/hybriddynamics 算法守卫）。验收：**frozen 全套不回归**（203 闸门 + observables golden + 测量电池逐位）；imag + Exact 容器 = `ExactHybridizationKernel`（带号 normdelta 守卫单测）；`Base.:+` 分派与混合抛错单测；**hybriddynamics 算法守卫单测**（Exact × ExactTTI → ArgumentError；Frozen × 三算法不回归）。
- **Phase B（exactif 核心 + real/mixed 路由 + 因子化裁决门 + 金标准）**：§2.1 + §2.5；`ifd_backend_switch` 拒绝面（real/mixed + Exact ArgumentError）关闭改路由 `exactif_kernel`。验收：ExactDelta/ExactIFKernel 读数与旧栈（backup 可跑）逐位或机器精度对拍；`H.normdelta` vs IFD 金标准逐位；**imag + Exact 端到端 Z**（`correlationfunction`/`exactif_kernel` → `sysdynamics(kernel=…)` → `integrate`）vs IFD 金标准 `normdelta · Pf(B)` 机器精度；比值观测量带/不带常数同值（相消实证）；**多浴因子化裁决**（双浴联合构造 vs blocks+ / normdelta× 对拍——成立则 `Base.:+` Exact 语义落地为乘法，不成立则拒绝面单测）；**金标准 G1/G2 主验收**（§2.5：直接 Trotter vs 严格核装配，order ∈ {1,2} 逐对机器精度——不过 = 停）。
- **Phase C（管线）**：§2.3。验收：`exactif_pipeline`/`frozen_pipeline` 端到端 gt/lt 对 `validation/` CSV 机器精度回归（dsol 系列档案）；Z 常数取数 = kernel 字段按类型。
- **Phase D（imag 工具族）**：§2.4。验收：三验证脚本对 `validation/imag_*.csv` 全量数值回归。
- **Phase E（清理 + 测试 + 合并接口裁决）**：`test/exactif/` 套件（核心对拍 + CSV 回归精简版 + 类型分派/守卫/裁决门单测 + **金标准 G1/G2 固化**）；`changes.md`；**合并接口裁决点**（用户预设推迟：Exact 行为经 Phase B–D 完全摸清后，裁 Frozen/Exact 合一与否——若合一，统一类型回收 `HybridizationKernel` 名并内化分派；本计划不预设结论）。

## 4. 不做清单

- FockMPO 迁移（exactif 的 Trotter 链暂留现行 GMPS 路径；可选迁移另行）；
- imag 理论的再推导（`imag_theory_notes.md` 冻结随迁——工具移植 ≠ 理论重验）；
- ExactIF 的 `fullfactor` 化（ExactIF = 核构造器定位，不进 W3 因子路线）；
- 旧栈 API 兼容层（一次性切净）；
- order 4 Trotter 分解（S4 序列含负时间演化，用户裁决不支持）与非等距离散化（E5）；
- **BCS strict normdelta**（IFD 刻意留白——待 BCS 后端立项，届时 E4′/E6 契约直接套用）；**BCS bath 金标准**（本计划金标准范围 = normal bath，§2.5）；
- **双 Exact 浴联合构造接口**（按 §2.2/Phase B 因子化裁决：成立 ⟹ 加法（normdelta 乘）覆盖、无需本接口；不成立 ⟹ `exactif_kernel` 接浴列表联合构造——同一机器喂两浴，裁决前不预建）。

## 5. 工作量与顺序

Phase 0 已完成；Phase A ≈ 1 天（类型层级改名面 + 分派 + 测量转发 + hybriddynamics 守卫）；Phase B ≈ 2 天（exactif 核心 + 路由 + 裁决门 + 金标准 G1/G2——ED 参照 helper ≈ +0.5）；Phase C ≈ 2 天；Phase D ≈ 2 天；Phase E ≈ 0.5 天。总 ≈ 7.5 天。顺序：A → B → C → D → E；CSV 回归与金标准不过全停；合并接口裁决 = Phase E 显式用户裁决点，不预设。

## 6. 落地偏差记录（Phase A，2026-11-21）

- **E6 "两者同给且耦合不一致 → 抛错"未按文实施（Frozen 分支）**：落地语义 =
  显式 `coupling_norm` 优先（0 = 未给哨兵 → `kernel.coupling_norm`），Frozen
  不做一致性抛错。原因：two-band ED 锚（fewmodesbath/realtime.jl 两带 U≠0）
  的既有惯例是 kernel 全带 cn = 2α、K 侧显式传 α——比值观测量下常数分子分母
  同因子严格相消、"不一致同给"是合法用法；强行抛错/强用 kernel 值都会破坏
  frozen 基线（实测：强用 kernel 值 → 4e-3 系统偏差）。Exact 分支保留拒绝：
  显式非零 coupling_norm 响抛（normdelta 当耦合范数用 = 恒错）。
- **计划外前置修复（IFDiscretization）**：exact 路线 `_quad_nodes` 此前只有
  GL-M 连续谱方法——DiracDelta 派发后 `f(ϵ)` 抛错、DiscreteSpectrum ±Inf 撑界
  静默产 Inf/NaN。补离散两方法（物理模式直取、nbath 不消费——docstring
  「离散谱 = 谱线数目」的实装），IFD test T5 回归（离散直取 == T2 尖峰路线 ==
  DiracDelta≡单线 DiscreteSpectrum；六金标准全绿）。这是 Phase B 的硬前置：
  exactif extension 验证锚（zidentity 等）全用 DiracDelta 单模 bath。
- 测量层 kernel 转发的耦合面：8 入口 `coupling_norm` 缺省从
  `kernel_coupling_norm(kernel)` 改 0.0（常数改经 kernel 转发按型分派；
  Frozen 数值逐位不变——Weyl 公式输入同源）；`kernel_coupling_norm` 保留为
  E4′ 单一取数点（Frozen → Σα / Exact → normdelta / nothing → 0），供管线
  Z 常数取数。

## 7. 实施进展（2026-11-21 会话）

**已完成**：
- **Phase A 全部落地并验收**：类型层级（`AbstractHybridizationKernel` +
  `FrozenHybridizationKernel` + `ExactHybridizationKernel`，构造期 canonical
  化抽共用 `_canon_kernel_blocks`；`HybridizationKernel` 名退役无 alias）+
  常数消费分派（`_constant_lscale` 双方法 + sysdynamics/`_system_factor_fock`
  `kernel` kwarg + Exact 号折首站 `lmul!`）+ 测量层 8 入口 kernel 转发 +
  hybriddynamics 算法守卫（Exact × ExactTTI → ArgumentError）+ imag+Exact
  构造路线归位（元组解包 → Exact 容器）+ exports。验收 = 58 testsets + 3
  fixtures（weyl 4.5e-16 / BCS 12 / ExactTTI 6）全绿。
- **Phase B 核心（部分）**：`src/influencefunctional/exactif.jl`——
  `native_schedule`/`trotter_factors`（order ∈ {1,2}，E5 裁剪）+ bath Gaussian
  机器（`pfaffian`/`bath_norm` 等 extension 逐字）+ `StrictDelta`（去 QuAPI
  容器化中间载体）+ `quapi_strictdelta`（IFD `exact_fermionic_Δt_cells` 四分支
  块 → legpos 无损重构）+ `exactif_kernel`（跨支负号折入 blockification →
  `ExactHybridizationKernel`）+ `exactif_pipeline`/`exactif_Z`/
  `exactif_expectation`（调度门管线骨架）+ 折叠 `fock_sysdynamics` `gates`
  调度 kwarg（real/mixed、order=1 组合强制、Td 累计时间）。

**Phase B 停点（金标准未过全停——§3 纪律）**：zidentity 证书
（Z = ∫·normdelta ≡ 1，Toulouse εd=0）在 GTEMPO 新栈上 **order1/N=2 机器
精度（4.4e-16）、其余配置 O(δt) 残差（~5e-3）**。取证定谳：

1. **GTEMPO 收缩侧已严格**：∫ 逐位等于 IFD 契约装配 det(A)（紧凑 (ā,a) 子块
   + 测度 −1 + 单位链 + 闭合边 + 真空 ρ₀ 的行列式）——含 N=2 逐位一致。
2. **根因 = K 链机器**：新栈 `fock_sysdynamics` 胞 = S·W（W = anti-Wick/Weyl
   局部权重）——这是冻结核路线的 O(δt)/O(δt²) 测度象（冻结路线由 Weyl 标量
   补偿、ED 锚 1e-4 容差看不见）；严格契约要的是旧栈
   `_propagator_sparsegmps` 的**裸 Berezin 机器**：逐门**跨片二站点传播子
   因子**（bra 变量在节点 idx+1、ket 在节点 idx，Fock 矩阵元的精确相干态
   单项式展开、无 W 权重）+ `boundarycondition!`（exp(GTerm(ā⁻(k),a⁺(k);−1))
   Keldysh 折返）+ `systhermalstate!`（ρ₀）。新栈落点 = GM 因子机器
   （`MonomialSum` 4^B 单项式 / Pairing 二站点）+ expmult!/mult! 链。
3. **节点-片映射定谳**：新格 Weyl 胞布局（a = 片左端 t_{j−1}、ā = 片右端
   t_j）⟹ 片 j 隐含节点 j（a）与 j+1（ā）；旧格 Nt 片 + k = Nt+1 节点一一
   对应新格 Nt 片（**当前 exactif.jl 的 Nt = npairs+1 + 尾窗与 Δ 整片映射
   需按此重做**：Δ 点 p ↔ {a@片p, ā@片p−1}、门 = 同片因子 exp(ā_j S a_j)）。
4. **会话内临时状态**：exactif.jl 为未接线 WIP（不进 runtests），全量回归
   不受影响；`gates` kwarg 缺省 nothing = 现行路径逐位不变。

**下一步（Phase B 续）**：①裸 K 链机器移植（MonomialSum 门因子 + 折返因子
+ ρ₀ 因子，符号约定由 zidentity 证书机器精度钉定）；②按节点映射重做
`_strict_pos`/管线 lattice；③G1/G2 金标准；④因子化裁决门。

## 8. 关联修复（IFDiscretization，2026-11-21，两处实证 bug）

1. **exact 路线谱支持**：`_quad_nodes` 原只有 GL 连续谱方法——DiracDelta
   派发后 `f(ϵ)` 抛错、DiscreteSpectrum ±Inf 撑界静默产 Inf/NaN。补离散两
   方法（物理模式直取、nbath 不消费）+ T5 回归（离散直取 == T2 尖峰路线 ==
   DiracDelta ≡ 单线 DiscreteSpectrum；六金标准绿）。
2. **realtime prefactor 真空权重归一**：IFD 移植把旧 QuAPI 的
   `1/(1+exp(−βε))`（真空权重 1−n̄）改写成物理占据 `1/(exp(βε)+1)`——差每模
   e^{βε}，**对称谱（Σε=0：金标准半圆/T2 双线）下两种公式简并、掩盖至今**；
   非对称离散谱 Z 证书直接破（Z = e^{−βΣε}）。修复 = Π(1−n̄)·detBB（实时
   结构化 + 稠密参照两处；混合路径经 Zb=Π(1+e^{−βε}) 已是真空归一、未动）。
   新增「Z 证书扩展：非对称离散谱」回归（DiracDelta 单模/DiscreteSpectrum
   双线 × orders 1/2 × μ 平移，Z ≡ 1 机器；实时金标准 74/74 绿）。

## 7.1 裸 K 链机器校准状态（2026-11-21 第二轮，sysdynamics_exact 落地）

**已落地**（`sysdynamics_exact`，exactif.jl）：真空 + 逐门跨片传播子因子
（`_fockop_factor`：Fock 矩阵精确相干态展开 Σ M[I,J]·c̄(bra 带升序)·c(ket
带降序)，旧 `_fockpropagator_sites` 序约定逐字；`lattice_monomial` + `togmps`
物化 + `mult!` 压缩上链，truncrelerr(1e-14)）+ Keldysh 折返（带子集展开
(−1)^{|S|}）+ ρ₀（bra ā⁺@(+,1)、ket a⁻@(−,1)）。荷守恒门守卫、恒等门不
短路（单位链项在门因子内）。lattice = Nt = npairs+1 片、Δ 点 p ↔ 片 p。

**证书现状**：|Z−1| ~ 1e-2（εd=0 即全部恒等门也差）。∫/det(A_ifd) 指纹：
- order 1（npairs = 2,3）：0.98935, 0.98121——按 gap 衰减 f ≈ 0.9919/gap
  + 基底偏移 c ≈ 1.0077（fold/ρ₀ 通道）；
- order 2（npairs = 5,7）：**0.958663 完全饱和**（npairs 无关！）——
  额外 gap 零贡献，指向 order 2 调度的零门 gap 上有结构性抵消/缺失。

**约定台账（下一步校准用）**：
1. 旧栈 K 链实为**四件套**：vacuumstate + **bulkconnection!**（逐片逐支
   exp(GTerm(ā@(j+1), a@(j); +1))——+1 系数的跨点分辨率，2026-11-21 深挖
   发现；新栈未移植）+ 门因子 + 折返 + ρ₀。我的机器缺 bulkconnection 等价物。
2. 契约装配目标分布：逐点同支测度 −ā_p a_p（×2 支）+ 单位链 +1 + 闭合 −1。
   我尝试的显式测度块（同支配对 coeff −1）→ 爆炸（×4/点）；跨支配对 → 整体
   偏移——**two-chain 收缩引擎的站点配对语义**（是否内建 Berezin pairing）
   未定谳，是测度块该不该显式上链的先决问题。
3. 恒等门因子符号：旧 `_fockpropagator_sites` 的反序奇偶（seq = [bra 位置,
   ket 位置]）⟹ 恒等门 = 1 − ā@g+1 a@g（跨点 −1）；契约单位链 = +1——
   两者差 bulkconnection(+1) 与测度(−1) 的分解归属，未定谳。
4. 校准正路（下一会话）：旧栈不可运行 ⟹ 以**契约装配 det(A) 为 oracle**
   做单点扫描——先扩 engine 语义探针（两链收缩 vs 显式测度块的四种组合），
   再扫 (bulk ±1、门对角 ±1、折返 ±1) × 测度(无/−1/+1) 的离散空间，
   order 1 npairs=2 单点机器精度命中后全配置扫掠。

## 7.2 裸机器校准第二轮结果 + ρ₀/读出阻塞点（2026-11-21 深夜）

**已定谳并落地**（zidentity 证书双路线机器精度，2026-11-21）：
1. **折返系数 = +1**（lattice_monomial 给定序语义；旧 GTerm −1 的换号）。
2. **反测度 = 1 + η̄η**（lattice_monomial 语义下 counter_c=−1 的换号；腿序
   (η̄,η)；等价变体 = 翻转腿序 × 翻转符号——同一算符）。
3. **compact_strictdelta 移植**（compact_delta 逐字 + blockification：接口
   c=1..nl ↔ 片 c+1、消元负号 sgn=−1、normdelta = ν·Pf(G)）；sysdynamics_exact
   增 counter_c kwarg（接口片 2..Nt−1 逐支逐带幂零因子）。
4. **证书**：QuAPI-Δ 路线与 compact 路线 εd=0 全 orders/N 机器精度
   （quapi ≤1.6e-14 全扫；compact ≤2.4e-15）。
5. **kicks 管线撤销**（用户纠偏：IF 组装不依赖杂质哈密顿量——核 = 纯浴
   对象、杂质全项走 K 链门；IFD kicks 接口只是其金标准测试的效率技巧）。

**阻塞点（重大不明确——待用户裁决/指点）**：观测量读出。判别实验：
- εd=0（恒等门）下插入**平凡热态 ρ₀ = I/2**（因子 = 0.5·vac + 0.5·横支对
  ā⁺@1 a⁻@1）⟹ **Z = 0.25 ≠ 1**——横支 ρ₀ 因子结构性破坏收缩（8 种支归
  属/取向全不满足幺正迹 Z(ρ₀)=1 oracle；非符号级偏差）。
- n(T) 的 2 腿读出（全片×全荷中性腿对×符号穷举，含 verify_compact 的折返点
  (a⁺,ā⁻) 投影）在两路线上均只命中 normdelta 常数通道，杂质内容不可达。

**证据链**：证书（εd=0、真空 ρ₀）机器精度 ⟹ 因子式构造/测度分布/引擎收缩
自洽；但任何非真空 ρ₀ 或观测量插入因子都塌缩到常数通道 ⟹ K 链的杂质线程
在我的构造里没有被放在可被横支因子/读出因子正确楔入的位置。怀疑方向：
旧栈 GTerm→PartialMPO 的链上语义（apply! 的站点楔积与 Z2 扇区键的荷穿线）
与 GM mult! 的站点楔积不等价——横支因子（ρ₀/读出）在旧栈承载跨站荷结构，
在 GM 栈需要不同的表达（可能 = open_rho_native 的 fusion-tree ρ(T) 直读
机制，而非插入因子）。

**下一步选项（待裁决）**：
(a) 移植 open_rho_native 的 ρ(T) 直读（旧 g_fuse/g_trace_phy 的 GM 对应物
    重建——中等工作量）；
(b) ρ₀/读出改为**IF 链侧承载**（把杂质热态折进 IF 或核常数——需推导）；
(c) 用户指点旧栈 GTerm/apply! 站点语义或读出挂链位置的正确理解。

当前可交付：双路线严格核构造 + εd=0 证书机器精度（zidentity 完整复现）。
