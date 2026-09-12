# 接口更改说明

本轮重构涉及的接口更改汇总。所有更改均已通过全量测试验证（行为不变或数值等价）。

## 第三批更新

### 杂质哈密顿量类型层级与接口

- 新增抽象类型层级（`src/sysdynamics/sysdynamics.jl`）：
  - `AbstractImpurityHamiltonian`：所有杂质哈密顿量的根类型
  - `ConstImpurityHamiltonian <: AbstractImpurityHamiltonian`：常数哈密顿量（quench 模型属于此类）；**只有 `ConstImpurityHamiltonian` 支持 `sysdynamics_fast`**
  - `AbstractTdImpurityHamiltonian <: AbstractImpurityHamiltonian`：含时哈密顿量（原 `GeneralTdImpurityHamiltonian` 更名）
- 接口函数（`AbstractImpurityHamiltonian` 的两个必须实现）：
  - `fock_propagator(model, branch, dt) -> FockMatrix`：`ConstImpurityHamiltonian` 的接口，`branch ∈ {:+, :-, :τ}`
  - `fock_thermalstate(model, β) -> FockMatrix`：设置 real 轮廓初态时需要；时变模型返回 `hτ` 的热态
  - `fock_propagator(model, branch, dt, t) -> FockMatrix`：`AbstractTdImpurityHamiltonian` 的接口，`t` 为该步所在物理时间区间的左端点（forward: `(idx-1)δt`；backward: `(Nt-idx)δt`）；`Const` 模型提供同签名委托方法（忽略 `t`），内部统一按 4 参调用
  - `num_bands(model) -> Int`：默认读取 `model.bands` 字段；无该字段的模型单独定义（`AndersonIM` → 2，`ToulouseIM` → 1）
- 移除显式传 bands 的旧方法：`fock_propagator(model, branch, dt, bands)`（4 参泛型）与 `fock_thermalstate(model, β, bands)`（3 参泛型）；传播子缓存键加入分支时间 `t`（时变模型每步不复用）。

### 预定义模型

- `AndersonIM`：**固定为两带**（每自旋一带，`H = μ(n₁+n₂) + U n₁n₂`），移除 `bands` 字段与 `bands` 关键字；用于两带 lattice 时不再需要（也不能再传）`bands=2`。
- 新增 `ToulouseIM(μ)`：单带、`U = 0` 的 Anderson 情形（原 `AndersonIM(U=0)` 在单带 lattice 上的角色），`num_bands = 1`。
- 在单带 lattice 上做 `U = 0` 计算的代码需改用 `ToulouseIM`；`AndersonIM(U=0)` 只能用于两带 lattice。

### Quench 与时变模型

- `QuenchImpurityHamiltonian` 更名为 `QuenchedImpurityHamiltonian`，字段 `h0/h1` 更名为 `hτ/ht`（语义不变：τ 分支用 `hτ`，实分支用 `ht`，热态取 `hτ`）；仅能通过构造器构建，无 `push!`。
- 新增 `TdImpurityOp(data, f; bands)`：系数随时间变化的 `NormalTerm` 集合，`t` 时刻贡献 `f(t)·data`。
- 新增 `TdImpurityHamiltonian(hτ, ht, htt; bands)`：时变杂质哈密顿量。
  - τ 分支与热态：常量 `hτ`
  - 实分支 `t` 时刻：`ht + Σ op(t)`，即调用 `model(t)` 得到该时刻的 `ImpurityHamiltonian`
  - `sysdynamics`/`baresysdynamics` 在全部三种轮廓上可用；`sysdynamics_fast` 不支持（仅 Const）
- `systhermalstate!`/`systhermalstate` 的 model 签名放宽为 `AbstractImpurityHamiltonian`；`sysdynamics`/`baresysdynamics` 同样放宽（支持 Const 与 Td），`sysdynamics_fast`/`baresysdynamics_fast` 保持 `ConstImpurityHamiltonian`。

### 测试与文档

- 单模 bath 的 quench（`QuenchedImpurityHamiltonian`）与时变（`TdImpurityHamiltonian`）测试覆盖 imag/real/mixed 三种 lattice，ED 参考采用与 GTEMPO 相同的逐步常数哈密顿量离散化。
- 原单模 bath 测试中 `AndersonIM(U=0)`（单带 lattice）全部改为 `ToulouseIM`；两带 lattice 的 `AndersonIM` 去掉 `bands=2` 关键字。

## 第四批更新

### ExactTTIIF 原地构建接口（in-place）

- **删除** `influencefunctional` 与 `influencefunctional!`（`src/influencefunctional/hybridization/exact_ttiif/`）：与 `hybriddynamics`/`hybriddynamics!` 完全等价，导出去重，导出列表移除 `influencefunctional, influencefunctional!`。
- **唯一入口**：`hybriddynamics(lattice, corr, alg::ExactTTIIF; band)`（从真空态构建完整 IF）与 `hybriddynamics!(gmps, lattice, corr, alg::ExactTTIIF; band)`（把 IF 逐 term 原地乘入已有 `GrassmannMPS`，原 `influencefunctional!` 的 `mult!` 循环内联至 `hybriddynamics!`）。
- **推荐工作流**：`gmps = sysdynamics(lat, model, trunc=trunc)` 得到杂质动力学后，`hybriddynamics!(gmps, lat, corr, ExactTTIIF())` 直接在同一个状态上合并 bath 影响；后续观测只需单状态（`gτ_series`/`gtlt_series` 等不再需要单独的 `mpsI` 参数）。
- **多带处理移入 `_influencefunctional`**：`lattice.bands > 1` 时先 `similar(lattice, bands=1)` 建单带格点（`_influencefunctional_util` 仅支持 `bands == 1`），再对每个 term 用 `fillband(lattice, mps, band=band)` 扩充回全格点。
- 测试固化：
  - `test/api/hybriddynamics.jl`：in-place 与分开构造一致（`hybriddynamics!(vacuumstate(...), ...) == hybriddynamics(...)` 位级一致；`sysdynamics → hybriddynamics!` 与 `mult(sysdynamics, hybriddynamics(...))` 一致，实时容差 1e-3 为截断顺序噪声）；多带 lattice 的内部 `fillband` 路径校验（1e-12）。
  - `test/normalbath/fewmodesbath/imagtime.jl` 与 `realtime.jl`：新增 "in-place ExactTTIIF hybriddynamics!" testset，完整工作流 `sysdynamics → boundarycondition! →（实时加 systhermalstate!）→ hybriddynamics!` 的 GF 与 ED 参考对比（imag 1e-2 / real 3e-2），并与标准分开构造流程一致（1e-3）。

## 第二批更新

### 指数展开算法外包至 ExpExp.jl

- 删除 `src/mpo/mpohamiltonian/schurmpo/exponentialexpansion.jl` 与 `exponentialexpansion2.jl`（Prony/最小二乘拟合、自动步长选择、`first_period`、`cut` 等全部实现），改为依赖独立包 **ExpExp.jl**（path 依赖）。
- 算法类型重命名：

| 旧名（已删除） | 新名（ExpExp 导出，GTEMPO re-export） | 说明 |
|---|---|---|
| `PronyExpansion` | `OverDeterminedProny` | 超定最小二乘 Prony；`stepsize` 为 `Int` 或 `nothing`（自动步长选择） |
| `DeterminedPronyExpansion` | `DeterminedProny` | 方程数恰好确定的经典 Prony |
| `PronyExpansion2` | 由 `OverDeterminedProny(stepsize=nothing)` 覆盖 | 自动步长级联 + 剪枝（`cut`） |
| `LsqExpansion2` | `LeastSquareProny` | Prony 初值 + (norm, phase) 实参数化阻尼 Gauss–Newton 精修 |
| — | `MatrixPencil`（新增） | 矩阵束算法 |

- `ExponentialExpansionAlgorithm`、`AbstractPronyExpansion`、`exponential_expansion`、`expansion_error` 由 ExpExp 提供；GTEMPO 不再扩展 ExpExp 的 `exponential_expansion`，对 `GenericDecayTerm` 的转换函数重命名为 `expand_decayterm`。
- 依赖变化：新增 `ExpExp`；移除 `LsqFit`、`Polynomials`（仅被删除的实现使用）。
- `XTRGIF`/`ExactTTIIF`/`TDVPIF` 的默认 `algexpan` 与各 `influenceoperators`/`hybriddynamics` 入口的默认展开算法同步改为 `OverDeterminedProny`。

### 观测函数接口精简

- 删除 `Gτ`、`Gt`、`Gm` 及其全格点批量版本、`parallel_Gτ`；一律改用 `gf(lattice, (ContourIndex(...), ContourIndex(...)), A, B...; Z, alg)`。
- 删除 `cached_Gτ`、`cached_Gt`、`cached_Gm`、`cached_contour_ordered_Gm`；改用 `cached_gf`。
- 删除 `cached_Gτ_fast`、`cached_Gt_fast`、`cached_Gm_fast`；改用 `cached_gf_fast`（等时端点修正 `r[end] = 1 - r[1]` 在 imag 轮廓由 `cached_gf_fast` 内部处理，mixed 轮廓需调用方显式执行）。
- `cached_greater`/`cached_lesser` 重写为直接构造 `ContourIndex` 调用 `cached_gf`：不再接受 `c1/c2/b1/b2` 关键字（语义固定为 greater/lesser），仅保留 `band` 与缓存相关关键字。
- `parallel_run`、`parallel_integrate`（与 `parallel_Gτ` 一起）整文件注释停用（`src/integration/parallelintegrate.jl`）。
- `occupation2`（等时 Green 函数路径的占据数实现）已注释停用：与 Toulouse ED 严格解对比，Keldysh 轮廓略优于 `occupation`（0.35% vs 0.53%，均为离散化/截断噪声量级），但虚时轮廓在 i ≥ 2 存在约 6% 的错误跳变（仅 i = 1 恰好正确，边界 Grassmann 迹贡献未正确计入）。保留 `occupation`（insert_n 路径）与 `cached_occupation`。

### 其它

- `_normalize!(psi::GrassmannMPS)` 改为重载 `LinearAlgebra.normalize!`，不再导出私有名。
- 项目根目录生成 `Manifest.toml`（此前缺失导致缓存失效后 path 依赖解析失败）；`.gitignore` 加入 `Manifest.toml` 与 `CVTemporaryFiles`，并删除遗留的 `CVTemporaryFiles/` 运行时缓存目录。
- `tutorials/` 移至 `docs/tutorials/`；新增 Documenter 文档结构（`docs/make.jl`、`docs/Project.toml`、`docs/src/`）与 `README.md`。
- 为 `ContourIndex`、`environments`、`hybriddynamics`、`occupation`、`sysdynamicsstepper!` 补充英文 docstring。

## 删除

### CachedVector 磁盘缓存机制
- 删除 `src/auxiliary/CachedVectors.jl`：`CachedVector` 类型及其磁盘换入换出机制（`all2disk!`、`clear!`、`change_sweeporder!`、`destory_copy!`、`@withlock` 等）与 `DefaultUseCache` 常量全部移除，不允许再使用该机制。
- 以下函数签名中的 `useHCache` 关键字参数已移除，内部一律使用普通 `Vector`：
  - `mult!`、`mult_cache`、`parint_cache`、`TwosideExpectationCache`、`TwosideExpectationCache2`
  - CUDA 扩展：`cumult!`、`cu_mult_cache`、`cu_parint_cache`、`cuTwosideExpectationCache2`
- 原 `destory_copy!(dst::Vector, src::Vector)` 调用点改为 `copy!`。

### 依赖
- 移除 `Statistics`：`iterative_error_2` 改为独立实现（样本标准差/均值的显式公式，数值等价）。
- 移除 `Permutations`：`CoxeterDecomposition(Permutation(perm))` 的相邻对换分解由 `TK.permutation2swaps` 替代（语义相同：按序施加对换后 `new[k] = old[perm[k]]`）；置换符号 `sign(p)` 由 `(-1)^length(swaps)` 替代。
- 移除 `Serialization`：仅 CachedVectors 使用；tutorials 各自独立 `using Serialization`，不受影响。

### mult2! / mult2
- 原 `mult!`（res 中转 + copy back 的写法，为 CachedVector 引入）在去掉缓存机制后与 `mult2!` 行为完全一致（scaling、svectors、物理态均相同，仅个别 site 张量相差 bond 规范自由度），已删除；`mult2!` 更名为 `mult!`，`mult2` 包装函数删除（由既有的 `mult(x, y) = mult!(copy(x), y)` 覆盖）。

## 重命名

| 旧名 | 新名 | 说明 |
|---|---|---|
| `compensate_twist!(t, i, j)` | 并入 `compensate_twists!(t, (i, j))` | 单对情况由多变体版本覆盖，单次遍历融合树 |
| `occupation` / `occupation2` | 互换 | 原 `occupation2`（通用 `insert_n` + `integrate` 实现，`AbstractGrassmannLattice`）现名为 `occupation`；原 `occupation`（基于 `Gt`/`Gτ` 的 Green 函数实现）现名为 `occupation2` |
| `DMRGMult1` | `DMRG1` | 单点 DMRG 压缩算法 |
| `DMRGMult2` | `DMRG2` | 两点 DMRG 压缩算法 |
| `CuDMRGMult1` | `CuDMRG1` | CUDA 扩展中的对应包装类型 |
| `TranslationInvariantIF` | `XTRGIF` | 平移不变（XTRG 式加倍）影响泛函算法 |
| `ExactTranslationInvariantIF` | `ExactTTIIF` | 精确时间平移不变影响泛函算法 |

`DMRGMultAlgorithm` / `CuDMRGMultAlgorithm` 抽象类型名保持不变。

### 影响算符/影响泛函构建接口重命名

围绕"影响算符 → 单步指数 → 单步影响泛函 → 完整影响泛函"统一命名：

| 旧名 | 新名 | 说明 |
|---|---|---|
| `influenceoperator` | `influenceoperators` | 影响算符，输出为 MPO 组成的 tuple（虚时为 1 元组，实时为 4 分支元组） |
| `influenceoperatorexponential` | `influenceoperatorsteppers` | 影响算符的单步指数 e^{dt·h}（MPO tuple；虚时 FirstOrderStepper 版由单 MPO 改为 1 元组） |
| `differentialinfluencefunctional` | `influenceoperatorstepper` | 微分影响泛函（单步 IF，MPS）；私有 `_differentialinfluencefunctional` 同步改为 `_influenceoperatorstepper` |
| ExactTTIIF 的 `differentialinfluencefunctional` | `influencefunctional` | ExactTTIIF 的完整 IF 构建入口；私有 `_influenceoperatorexponential`（旧 `_differentialinfluencefunctional`）同步改为 `_influencefunctional` |

导出更新为 `influenceoperators, influenceoperatorsteppers, influenceoperatorstepper, influencefunctional`。

### 删除

- `test/models/bmps_gf.jl` 及其在 `test/runtests.jl` 中的 include。
- `band_boundary` 的 4 个死 method（realtime.jl 中 `A1Ā1B1B̄1b̄1B̄1ā1Ā1`、`A1B1ā1b̄1Ā1B̄1a1b1` 与两个 `A2...` ordering 版本）：这些 ordering 不在 `_AllowedRealGrassmannOrdering` 内或非 TimeLocal，`_fit_to_lattice` 永远不会走到；其余 4 个活跃 method（虚时 2 个、实时 2 个）保留。

### 其它

- `GrassmannMPS` 去掉第 3 个类型参数 `VA<:AbstractVector{A}`，`data` 字段实体化为 `Vector{A}`；三参构造函数中对非 Vector 输入做 `convert`。

## 新增

### TDVPIF 算法
- 新增 `TDVPIF <: InfluenceFunctionalAlgorithm`（`src/influencefunctional/hybridization/tdvpif/`）：
  把影响泛函视为影响算符 H 的"平衡态" IF = exp(H)，用二阶单点 TDVP 虚时流 dz/dτ = H·z（τ: 0 → 1）从恒等影响泛函出发演化得到。
  - `hybriddynamics(lattice, corr, alg::TDVPIF; band)` / `hybriddynamics!(gmps, lattice, corr, alg::TDVPIF; band)`
  - `hybriddynamics!` 支持把影响泛函在单次流中合并进任意初始 `GrassmannMPS`（z(1) = e^H·z(0)），如 `sysdynamics` 的输出
  - 实时格子上驱动流的是 4 个分支影响算符之和（逐个相加并以紧的 `DefaultITruncation` 压缩）
  - 费米符号全部继承自 `mult` 机制（`get_left_xy`、`updatemultleft/right`），键矩阵的有效映射实现为 AC 映射的投影 AL†·(H·z)
  - 测试见 `test/tempo/tdvpif.jl`（含 δ 收敛性验证与 in-place 合并验证）

## 导出更改
- `GrassmannBackend` 不再导出（仅模块内部使用；`@grassmann` 宏仍导出）。
- 其余导出名随上述重命名同步更新。
