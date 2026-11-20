# exactif 移植计划（v2.3）：extensions/exactif → 当前 GTEMPO（复用 HybridizationKernel + 统一常数契约 + 类型层级）

> **v2（2026-11-20 用户裁决：统一常数封装）**：Weyl 2 阶 QuAPI 的 Weyl 标量与严格 IF 的
> normdelta 都是"乘在整体收缩上的 Z 常数"——统一存储于 HybridizationKernel，下游
> `sysdynamics` 分派消费。E4 的元组配对返回废除，normdelta 单一存储进 kernel。
> **v2.2（2026-11-20 用户裁决：类型层级取代 scheme 字段）**：`AbstractHybridizationKernel`
> + `FrozenHybridizationKernel`（现行 QuAPI 核改名）/`ExactHybridizationKernel`（对接 IFD
> 严格 IF）两具体类型；**合并接口与否 = Exact 行为完全摸清后再裁**（显式推迟点）。
> **v2.3（2026-11-20 用户裁决：hybriddynamics 算法支持集）**：Exact 核天然破坏时间平移
> 不变性 ⟹ `hybriddynamics` 对 Exact 核只支持 `alg ∈ (PartialIF, PolarIF)`（类型级守卫）。
> 本修订与 `fockmpo_dress_selfread_plan.md`（合并计划：fock_sysdynamics 穿戴 /
> logscale E2 语义 / fock_initstate! 删除）对齐。

- 日期：2026-11-17（v1）；2026-11-20（v2 统一常数契约修订）
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

## 3. 阶段与验收

- **Phase 0（前置）**：~~`if_expmult_backend_plan` Phase 1~~ **已落地（2026-11-19 十七批）**。
- **Phase A（类型层级重构 + 常数消费落地面）**：§2.2 全部（改名面/分派/kernel kwarg/测量转发/imag+Exact 容器切换/hybriddynamics 算法守卫）。验收：**frozen 全套不回归**（203 闸门 + observables golden + 测量电池逐位）；imag + Exact 容器 = `ExactHybridizationKernel`（带号 normdelta 守卫单测）；`Base.:+` 分派与混合抛错单测；**hybriddynamics 算法守卫单测**（Exact × ExactTTI → ArgumentError；Frozen × 三算法不回归）。
- **Phase B（exactif 核心 + real/mixed 路由 + 因子化裁决门）**：§2.1；`ifd_backend_switch` 拒绝面（real/mixed + Exact ArgumentError）关闭改路由 `exactif_kernel`。验收：ExactDelta/ExactIFKernel 读数与旧栈（backup 可跑）逐位或机器精度对拍；`H.normdelta` vs IFD 金标准逐位；**imag + Exact 端到端 Z**（`correlationfunction`/`exactif_kernel` → `sysdynamics(kernel=…)` → `integrate`）vs IFD 金标准 `normdelta · Pf(B)` 机器精度；比值观测量带/不带常数同值（相消实证）；**多浴因子化裁决**（双浴联合构造 vs blocks+ / normdelta× 对拍——成立则 `Base.:+` Exact 语义落地为乘法，不成立则拒绝面单测）。
- **Phase C（管线）**：§2.3。验收：`exactif_pipeline`/`frozen_pipeline` 端到端 gt/lt 对 `validation/` CSV 机器精度回归（dsol 系列档案）；Z 常数取数 = kernel 字段按类型。
- **Phase D（imag 工具族）**：§2.4。验收：三验证脚本对 `validation/imag_*.csv` 全量数值回归。
- **Phase E（清理 + 测试 + 合并接口裁决）**：`test/exactif/` 套件（核心对拍 + CSV 回归精简版 + 类型分派/守卫/裁决门单测）；`changes.md`；**合并接口裁决点**（用户预设推迟：Exact 行为经 Phase B–D 完全摸清后，裁 Frozen/Exact 合一与否——若合一，统一类型回收 `HybridizationKernel` 名并内化分派；本计划不预设结论）。

## 4. 不做清单

- FockMPO 迁移（exactif 的 Trotter 链暂留现行 GMPS 路径；可选迁移另行）；
- imag 理论的再推导（`imag_theory_notes.md` 冻结随迁——工具移植 ≠ 理论重验）；
- ExactIF 的 `fullfactor` 化（ExactIF = 核构造器定位，不进 W3 因子路线）；
- 旧栈 API 兼容层（一次性切净）；
- order 4 Trotter 分解（S4 序列含负时间演化，用户裁决不支持）与非等距离散化（E5）；
- **BCS strict normdelta**（IFD 刻意留白——待 BCS 后端立项，届时 E4′/E6 契约直接套用）；
- **双 Exact 浴联合构造接口**（按 §2.2/Phase B 因子化裁决：成立 ⟹ 加法（normdelta 乘）覆盖、无需本接口；不成立 ⟹ `exactif_kernel` 接浴列表联合构造——同一机器喂两浴，裁决前不预建）。

## 5. 工作量与顺序

Phase 0 已完成；Phase A ≈ 1 天（类型层级改名面 + 分派 + 测量转发）；Phase B ≈ 1.5 天（exactif 核心 + 路由 + 裁决门）；Phase C ≈ 2 天；Phase D ≈ 2 天；Phase E ≈ 0.5 天。总 ≈ 7 天。顺序：A → B → C → D → E；CSV 回归不过全停；合并接口裁决 = Phase E 显式用户裁决点，不预设。
