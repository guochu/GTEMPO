# IF 装配后端统一计划（v2）：influencefunctional/ 拆分 + HybridizationKernel{T} 重定义 + expmult! 后端

- 日期：2026-11-17（v2：用户五项裁决——文件夹拆分、HybridizationKernel{T} 重定义、fullfactor/PartialExp、tti_* 因子/PolarExp）
- 基线：GTEMPO HEAD `f3ee05c4` 系工作区
- 关联：**GM `polarexp-svd-rankone` 计划 = 本计划的 GM 侧前置**（PolarExp 支持 Pairing/AnomalousPairing SVD 秩一精确指数，旧实现更名 TTIExp——已定稿可实施）；`sysdynamics_fock_mpo_plan.md`（K 侧，独立并行）；kernel/fit Part1（kernel 存 corr 视图已实施，本计划的 blocks 键在其上推广）

## 0. 已裁决（2026-11-17 用户五项）

- **W1**：新建 `src/influencefunctional/` 文件夹，两个文件分别实现 ExactTTI 与 PartialIF。
- **W2**：删除现行 `HybridizationKernel` 定义；`correlationfunction.jl` 定义 `HybridizationKernel{T<:Number}`，仅三 field：`blocks::Dict{Tuple{ContourBlockIndex,ContourBlockIndex}, Union{Matrix{T},CorrelationMatrix{T}}}`、`coupling_norm::Float64`、`bands::Vector{Int}`；配套定义 `ContourBlockIndex`（`branch::Symbol, band::Int`）与 `scalartype`（仅 Type 上定义）。**统一 normal bath 与 BCS bath**（双构造方法 + `changebands` + 符号处理收进构造期，见 §2.1）。
- **W3**：`influencefunctional/partialif.jl` 定义 `fullfactor(HybridizationKernel) -> AnomalousPairing`，调用 GM 后端 `expmult!(gmps, fullfactor(H), PartialExp)` 实现 PartialIF。
- **W4**：`influencefunctional/exacttti.jl` 定义 `tti_imag_factors -> ExpDecaySum`、`tti_real_factors -> ExpDecaySumSum`、`tti_cross_factors -> Pairing`（与现行 imag_/real_/cross_factors 内容一致），分别调用 `expmult!(gmps, factors, PolarExp)` 实现 ExactTTI。
- **W5**（继承 v1）：golden = 机器精度新锚，前提断言 = **每因子 dense 对拍**（各 tti_*/fullfactor 产物 vs 现行因子链）。

## 1. 探查定谳

**F1【GM 因子语义核对（factors.jl 实读）】**
- `AnomalousPairing{T}`：C 为 **2L×2L 一般二次型系数矩阵**（L = 链长，2L = 全变量空间）——正规块与反常块统一容器，`fullfactor` 的天然目标（W3）；
- `Pairing{T}`：`O = Σ C[i,j]·a_{asites[i]}·JW·b_{bsites[j]}`，矩形站点表（asites/bsites 独立、可不等长），`C` 为 |asites|×|bsites| 一般耦合矩阵——cross 通道 SVD 秩一分量的天然容器（W4）；
- `ExpDecaySum{T}`：固定轨道对（cell/asite/bsite）上的多极点级联池（αs/λs 等长）——imag 热通道（锚点列拟合）天然容器；`ExpDecaySumSum{T}` = ExpDecaySum 之和（多轨道/多带池）——real 双通道（greater/lesser 两侧列）天然容器。
- `expmult!` 三算法对因子的支持：PolarExp 现仅 ExpDecay 族（Pairing/AnomalousPairing 分派抛错）——**W4 的 `tti_cross_factors -> Pairing × PolarExp` 依赖 GM `polarexp-svd-rankone` 计划实施**（本计划 GM 侧前置）；PartialExp 行表机制对 AnomalousPairing 的支持同样以前置计划为准（Phase 0 核对）。

**F2【现行 HybridizationKernel（influence.jl:24-45）与新定义的差异】**现行 field：`lattice`（→ 删，`bands::Vector{Int}` 覆盖带号表替代）、`coupling_norm::Union{Nothing,Float64}`（→ 具体化 `Float64`，"order=2 必填"守卫上移 `prepare_system`/调用方——已在位）、`diagnostics::Dict`（→ 删，质量键迁 `InfluenceFunctional.diagnostics`——tti_* 构造器返回 audit 值）、blocks 键 `(branch,branch)` 二元组（→ `ContourBlockIndex` 对，键值可 `CorrelationMatrix{T}` 直接挂 QuAPI Toeplitz——kernel/fit Part1 的 corr 语义推广到全键）。

**F3【blocks 值的两种表示】** `Union{Matrix{T}, CorrelationMatrix{T}}`：`Matrix` = dense 已物化；`CorrelationMatrix` = QuAPI Toeplitz（O(N) 向量、|i−j| 查表）——fullfactor/tti_* 的列读取走 `getindex`（两表示同值），物化仅在构造端需要时发生。

**F4【消费方清单（blocks/field 触点）】** 构造端 `hybridization`（influence.jl:61-83）退役/更名 `correlationfunction`（签名 per §2.1）、`kernel_coupling_norm`（measurements/inserted_k.jl:231——新类型无 nothing 分支，简化为一行）、blocks 索引（partial_if:147、exacttti 通道构造、currents:354/388）、kernel/fit Part1 的 branch 视图 accessor（键 `(branch,branch)` → `ContourBlockIndex` 对）＋ currents 符号补号退役（:354-355 `ηpm = -blocks` → 直读）＋ `hybriddynamics` 签名收窄（删 bath/band/orbital 输入，收 `kernel`）、scalartype（GTEMPO.jl:73 已 re-export GM scalartype——补 Type 方法一处，instance 走泛型回退）。

**F5【BCS 统一的天然性】** 现行 BCS 走 corr 直通（`bcs_partial_if(lattice, corr)`）+ `BCSCorrelationFunction` 四味块；`ContourBlockIndex` 化后 BCS 块 = 同一 Dict 的 (branch, band) 键条目，且 `AnomalousPairing` 的 2L×2L C **天然容纳反常象限**——BCS 的 fullfactor 与 normal 的 fullfactor 同一函数（BCS = C 的反常象限非零特例）。`bcs_symmetry_guard`（仅 :z2）保留。

**F6【算法语义换代声明】** PartialIF 内核从「逐行秩一因子积 + 乘间压缩」→「整形式 AnomalousPairing 幂零指数（PartialExp 截断）+ 压缩」——名称保留、内核换代；两者是同一 IF 算符的两种精确表示（PartialExp 幂零级数无配置限制），互换对拍闸门保留。ExactTTI 的 tti_* 因子内容与现行 imag_/real_/cross_factors 一致（拟合照旧），装配从 memory_cell/emit 换成 expmult!+PolarExp。

**F7【文件/退役清单】** 现行三通道构造器 `imag_factors`/`real_factors`/`cross_factors` 与 `partial_if`/`partial_row`/`memory_cell`/`embed_cells`/`rankone_factor`（influence.jl:99-509）→ 迁入/退役；`influence_functional` 入口与 `InfluenceFunctional` 结构保留（算法分派重接）。

## 2. 设计

### 2.1 `ContourBlockIndex` + `HybridizationKernel{T}`（correlationfunction.jl）

```julia
struct ContourBlockIndex
    branch::Symbol      # :τ / :+ / :-（BCS 同三分支）
    band::Int
end

struct HybridizationKernel{T<:Number}
    blocks::Dict{Tuple{ContourBlockIndex,ContourBlockIndex}, Union{Matrix{T},CorrelationMatrix{T}}}
    coupling_norm::Float64
    bands::Vector{Int}              # 本 kernel 覆盖的带号有序表（单带 → [b]；全带 → [1..B]）
end
scalartype(::Type{HybridizationKernel{T}}) where {T} = T
    # 仅 Type 上定义——TensorKit 泛型回退 scalartype(x) = scalartype(typeof(x))
    # 自动覆盖实例，无需 instance 方法
```

**双构造方法（constructor = correlationfunction，与文件同名）**：

```julia
correlationfunction(normbath::AbstractNormalBath, lattice; band) -> HybridizationKernel
    # normal bath：band kwarg = 单带构造（blocks 仅该带的 ContourBlockIndex 键、bands = [band]）；
    # 缺省全带。内部走 QuAPI Δτ/Δt/Δm → branch 视图挂键

correlationfunction(bcsbath::AbstractBCSBath, lattice; orbital) -> HybridizationKernel
    # BCS bath：orbital kwarg（Nambu 对指标，旧 lattice orbitals = div(bands,2) 语义），
    # bcs_Δτ/Δt/Δm 四味块按 (branch, band) 键入同一 Dict——normal/BCS 从此同一类型
```

**符号处理收进构造期（用户裁决，存储即物理终值）**：现行约定是「存储携带 QuAPI 符号、读点就地补号」（correlationfunction.jl:24-34 注记：currents 的 η⁺⁻ = −blocks、混合 τ↔rt 交叉负号等）——**退役**。`correlationfunction` 构造期一次性把全部跨支相对符号/Keldysh 约定折进存储矩阵：**blocks 存的 = 物理终值**，下游（tti_* 列读取、fullfactor 组装、currents η 系数）直读不再补号（currents :354-355 的 `ηpm = -blocks` 负号删除，F4 同步）。约定文档由「各读点就地处理」改为「构造期统一处理」一段注释。

**`changebands`：带重排**：

```julia
changebands(H::HybridizationKernel, x::Vector{Pair{Int,Int}}) -> HybridizationKernel
    # 按 x = [old => new, …] 重排 band：blocks 键的 ContourBlockIndex.band 逐对映射、
    # bands::Vector{Int} 同步重排；标量/矩阵数据零改动
    # 用途：单带核建一份、重排出多带装配所需的全部带号（免重建 corr）
```

**`Base.:+`：多 bath 耦合同一杂质（用户裁决）**：

```julia
Base.:+(A::HybridizationKernel, B::HybridizationKernel) -> HybridizationKernel
    # blocks 键 = 并集：同键块相加（CorrelationMatrix+CorrelationMatrix = Toeplitz 向量相加；
    # 混合表示 = CorrelationMatrix 侧物化到 dense 再加——任一方 dense 即 dense 结果）、
    # 单侧缺键 = 零块透传；coupling_norm 相加；bands = 并集
    # 形状校验：同键块 size 不匹配 throw
    # 用途：双 lead（异 μ/β 无法合成单 bath）逐 lead 建 kernel 后相加——总耦合标量自动正确
```

**`hybriddynamics` 签名收窄（用户裁决：band/orbital 输入删除——信息在 kernel 中；alg 位置/关键字双形）**：

```julia
hybriddynamics(lattice, kernel::HybridizationKernel, alg) -> 逐带 IF        # 位置形（主形）
hybriddynamics(lattice, kernel::HybridizationKernel; alg) = hybriddynamics(lattice, kernel, alg)
    # 关键字形 = 位置形委托；缺省 alg = PartialIF()（主形缺省实现）
    # band 覆盖 = kernel.bands；多带 kernel 返回逐带 IF 序（Tuple/Vector，包装实现期定）
    # 单带 kernel（bands = [b]）返回单 IF——与现行 Is... 变参消费兼容
```

### 2.2 `src/influencefunctional/` 文件布局

```
src/influencefunctional/
├── influencefunctional.jl   # 头注释 + include + export（IF 算法入口/InfluenceFunctional 分派重接）
├── partialif.jl             # W3：fullfactor + PartialIF（expmult!×PartialExp）
└── exacttti.jl              # W4：tti_* 三构造器 + ExactTTI（expmult!×PolarExp）
```

GTEMPO.jl：`include("influence.jl")` 拆分——IF 装配段迁入，核块/通道词典保留在 correlationfunction/influence 的 kernel 构造段；`partial_if`/`imag_factors`/`real_factors`/`cross_factors`/`partial_row`/`memory_cell`/`embed_cells`/`rankone_factor` 退役（golden 交叉验证期 deprecated 保留）。

### 2.3 `partialif.jl`：`fullfactor` + PartialExp（W3）

```julia
fullfactor(H::HybridizationKernel, lattice; band) -> GM.AnomalousPairing
    # 跨 (ContourBlockIndex, ContourBlockIndex) 组装 2L×2L 一般二次型 C：
    # 变量序 = Chrono 全链变量序（lattice 决定，L = 总变量数）；块路由 = 键荷守恒
    # （CorrelationMatrix 值逐元素 getindex 物化；dense 值直抄）
PartialIF = 组装链：vacuumstate(lattice) →
    expmult!(gmps, fullfactor(H, lattice; band), PartialExp(SVDCompression(trunc)))
```

- AnomalousPairing C 的块布局：对角象限 = 各 (branch, band) 的同指数块（measure/传播子已折叠于 blocks 值内——与现行 cell mats 同源）、非对角象限 = 跨支/跨带耦合块（现行核块的跨支项）。**块布局与变量序的规范表显式成文进代码注释**。
- 语义换代（F6）：逐行精确性由 PartialExp 幂零级数的完备性承接（同一算符）；截断 = 压缩（SVDCompression），与现行乘间压缩同档。

### 2.4 `exacttti.jl`：tti_* 三构造器 + PolarExp（W4）

```julia
tti_imag_factors(H, lattice; band)  -> ExpDecaySum       # 热通道：锚点列拟合池（现行 imag_factors 因子内容）
tti_real_factors(H, lattice; band)  -> ExpDecaySumSum    # 实轴双通道：greater/lesser 两侧列（现行 real_factors 内容）
tti_cross_factors(H, lattice; band) -> Pairing            # 混合交叉：SVD 秩一（现行 cross_factors/rankone 内容）
ExactTTI = 组装链：vacuumstate(lattice) →
    expmult!(gmps, tti_imag_factors(...),  PolarExp(SVDCompression(trunc)))
    expmult!(gmps, tti_real_factors(...),  PolarExp(...))
    expmult!(gmps, tti_cross_factors(...), PolarExp(...))
```

- 拟合层不动（fit_exponentials 供给，anchored pencil）；tti_* 从 blocks（CorrelationMatrix 视图或 dense）读列 → 拟合 → GM 因子。
- **Keldysh/热配对审计保留**（现行 :365-378/:411-412）：tti_* 返回 `(factor, audit::NamedTuple)`，audit 并入 `InfluenceFunctional.diagnostics`（kernel diagnostics field 已删，F2）。
- `tti_cross_factors` 的 Pairing×PolarExp 依赖 GM 前置（F1）；前置未实施期间 cross 通道暂以现行 rankone 发射混装——**裁决点**。
- contact/对角项：现行 contact_factor 语义并入对应 tti_* 的池（对角折叠，ExpDecay 族原生支持）。

## 3. 阶段与验收

- **Phase 0（前置核对）**：GM polarexp-svd-rankone 实施状态确认（Pairing/AnomalousPairing×PolarExp/PartialExp 全通）；未就绪则本计划排队于其后。
- **Phase 1（HybridizationKernel 重定义）**：§2.1 + 消费方迁移（F4 全清单）。验收：全测试绿；kernel/fit corr 视图逐位对拍；BCS/normal 同类型实证。
- **Phase 2（influencefunctional/ 拆分 + 新装配）**：§2.2/2.3/2.4。验收闸门：① **每因子 dense 对拍**（fullfactor 的 AnomalousPairing todense vs 现行 partial_if 链 todense；tti_* 三因子 vs 现行 imag/real/cross 因子链——机器精度，全配置矩阵：三轮廓 × bands 1,2 × spec × 插入 × order）；② 互换对拍 ExactTTI ≡ PartialIF（≤ 现行 1e-10 档）；③ goldens 按值重锚（V3 语义，>1e-10 记录理由）；④ ED 对拍维持。
- **Phase 3（入口切换 + 旧径封存）**：默认切新；旧装配函数 deprecated 保留至交叉验证期满。
- **Phase 4（清理）**：influence.jl 装配段退役；bcs 直通路径收编（guard 保留）；性能基准（expmult! 单调用 vs 逐因子发射）。

## 4. 不做清单

- 缓存/效率（同 V2 口径——FockMPO 计划统一立项）；
- 拟合层改动（anchored pencil 归 ExpExp 计划，已交付）；
- TDVP2/XTRG 后端（GM 有、不引入）；
- 逐位 golden（机器精度契约）。

## 5. 工作量与顺序

Phase 0 ≈ 半天（GM 状态核对）；Phase 1 ≈ 1.5 天（类型重定义 + 消费方矩阵）；Phase 2 ≈ 2.5 天（fullfactor/tti_* 组装 + 对拍矩阵——本计划科学核心）；Phase 3 ≈ 半天；Phase 4 ≈ 半天。总 ≈ 5 天。顺序：GM 前置就绪 → 1 → 2 → 3 → 4；Phase 2 闸门不过全停。
