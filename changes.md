# 接口更改说明

本轮重构涉及的接口更改汇总。所有更改均已通过全量测试验证（行为不变或数值等价）。

## 新增实轴 Grassmann ordering：`a1ā1A1Ā1b1b̄1B1B̄1a2ā2A2Ā2b2b̄2B2B̄2`（2026-09-17）

- 新增**时间升序**（`RealTimeOrderingStyle = TimeAscending()`）的实轴 time-local ordering：每个时间块内按 `a⁻ā⁻a⁺Ā⁺b⁻b̄⁻b⁺B̄⁺` 排列，块顺序为 `j=1,…,k`（既有实轴 ordering 均为时间降序）。它是 mixed ordering `A1Ā1B1B̄1_a1ā1A1Ā1b1b̄1B1B̄1a2ā2A2Ā2b2b̄2B2B̄2` 的纯实轴对应版本，已导出。按命名约定，只写一个时间点的名字表示该模式按时间降序重复（默认 TimeDscending）；本 ordering 为时间升序，名字显式写出时间点 1 和 2。
- 影响泛函（TTIIF / ExactTTIIF / TDVPIF）对该 ordering **直接构造**（加入 `_AllowedRealGrassmannOrdering`），`_fit_to_lattice` 按 `RealTimeOrderingStyle` 决定时间块布放方向（升序时 `u_left` 端紧邻边界块的 `j=1`），不再走 `changeordering`/permute fallback。
- `makestep` 支持时间升序格点：新时间块追加在右端，旧块位点不变（降序时旧块右移的原行为不变）。

## Grassmann lattice 的 `show`/`print`（2026-09-17）

为所有 `AbstractGrassmannLattice`（imag/real/mixed 三种轮廓、全部 ordering）新增统一打印，基于 `indexmappings` 逐位点反查符号，无需按 ordering 特化：

- 紧凑形式 `show(io, lat)`（`print`/`string`/数组内）只输出符号链；REPL 多行显示（`text/plain`）额外给出类型名、参数（`bands`/`N`/`δt`/`Nτ`/`δτ`）和总位点数，类型名用 `nameof` 不带模块前缀。
- 符号约定全小写：带用 `a,b,c,…`；共轭为字母 + 组合长音符（U+0304）；时间点下标 ₀₁₂…；实轴 branch 上标 `⁺`/`⁻`；虚轴与 `i=0` 边界位点无上标；mixed 格点虚/实段之间以 `_` 连接。
- 行为由 `test/api/lattice.jl` 中 "API: lattice show / printing" 固化（102 项断言，含 7 条完整符号链 fixture、`text/plain` 精确文本及全部导出 ordering 的通用不变量）。

## 第七批更新（2026-09-17）：截断方案与 MPS 算法接口统一（同步 TEMPO 2026-09-17）

同步 TEMPO 同日的接口调整：压缩算法的截断参数统一为 `TruncationScheme` 对象，默认截断常量收敛。Z2Tensors 依赖包同步把截断类型 **`TruncationDimCutoff` 更名为 `TruncateDimCutoff`**（`truncdimcutoff(D, ϵ[, add_back])` 构造不变），GTEMPO 全库跟进。

### `SVDCompression`：`D`/`tol` 字段 → 参数化 `trunc` 字段

| 旧接口（已删除） | 新接口 |
|---|---|
| `SVDCompression(; D=Defaults.D, tol=Defaults.tol, verbosity=0)` | `SVDCompression(; trunc=truncdimcutoff(D=Defaults.D, ϵ=Defaults.tol, add_back=0), verbosity=0)` |
| `SVDCompression(trunc::TruncationDimCutoff; verbosity=0)` | `SVDCompression(trunc::TruncationScheme; verbosity=0)`（接受**任意** `TruncationScheme`，含 `truncdim`/`truncerr`） |
| `alg.trunc`（getproperty 合成）/ `get_trunc(alg)` / `alg.D` / `alg.ϵ` | `alg.trunc`（真实字段） |

- 结构体参数化为 `SVDCompression{T<:TruncationScheme}`；`Base.similar` 同步改为 `trunc` 参数。

### `DMRG1`/`DMRG2`：`trunc` 必须携带键维 `D`

- `trunc` 字段与构造器签名收紧为 **`TruncationWithD = Union{TruncationDimension, TruncationDimCutoff}`**（即 `truncdim(D)` 或 `truncdimcutoff(D, ϵ)`）；原因：迭代乘法的初始猜测（`:svd`/`:rand`/`:pre`）需要 `D` 信息，实现改用 `alg.trunc.D`。
- **删除** `Base.getproperty(::DMRGAlgorithm, :D/:ϵ)` 访问器；`Base.similar` 同步。

### 默认截断常量收敛（`src/defaults.jl`）

| 常量 | 变更 |
|---|---|
| `DefaultKTruncation` | `truncdimcutoff(D=1000, ϵ=1e-10)` → **`truncrelerr(Defaults.tolgauge)`**（只按相对 ϵ 截断，等价 TEMPO 的 `trunccutoff`） |
| `DefaultIntegrationTruncation` | **已删除**，原用点（`bulkconnection!`、`boundarycondition!`、`_permute!`、`BMPSIntegrate`、AC-BMPS 积分）改用 `DefaultKTruncation` |
| `DefaultTruncation` | **已删除**，原用点（`mult!`、`swap!`、`canonicalize!`、`parint_mult`/`_parint_svd_guess` 及 CUDA 扩展对应版本）改用 `DefaultITruncation` |
| `DefaultITruncation` | `truncdimcutoff(D=200, ϵ=1e-10)` → **`truncdimcutoff(D=Defaults.D, ϵ=Defaults.tolgauge)`**；`DefaultMultAlg = DMRG1(DefaultITruncation)` 随之 |
| TDVPIF 的 H 压缩 | `DefaultIntegrationTruncation` → `DefaultKTruncation` |

### 其它接口清理

- 删除 2 参的 `_swap_gate(m1, m2; trunc)`（死代码；4 参的 svector 版本保留）。
- `easy_swap!` 更名为 **`swap!`**（含 FockMPS 路径的调用点），默认截断由 `DefaultTruncation` 改为 `DefaultITruncation`。
- **`GrassmannMPS` 的字段 `svectors` 更名为 `s`**（删除合成 `getproperty` 别名；`svectors_uninitialized` / `unset_svectors!` 接口不变）。与 `FockMPS` 的字段命名对齐。
- **`ToulouseIM` / `AndersonIM` 的字段与关键字 `μ` 更名为 `ϵ_d`**。语义确认：`H = ϵ_d·n̂`，与 `ImpurityModelBase.Toulouse` 的 `ϵ_d` 完全一致（同为 on-site 能量直接乘 n̂，**无符号差**；BCS 测试中的 `ϵ_d=-ϵ_d` 只是参数取值）。`IRLM` 的 `μ` 关键字不受影响。
- 修复 `contour_ordered_gf` 的 `Z::Real` 类型约束 → **`Z::Number`**（实时间格点的配分函数是复数，原约束使该接口不可用，属 bug）。
- 删除 `test/api/prony.jl`：Prony 展开由 ExpExp 包负责测试，GTEMPO 直接使用且无重载（`AbstractPronyExpansion` 的 re-export 修复已在上批提交）。
- 新增 `test/api/grassmanntensor.jl`：系统对比 `@grassmann` 与 `@tensor`——偶宇称张量上两个宏完全一致（收缩/迹/原地加法/多链/order kwargs）；奇宇称下普通收缩仍一致，闭合 U-turn 圈（trace）按费米圈规则差 −1（用 parity 分解逐一验证）。

### 迁移指南

- `SVDCompression(D=χ)` → `SVDCompression(truncdimcutoff(D=χ, ϵ=Defaults.tol))`（或按需 `truncdimcutoff(D=χ, ϵ=...)`）。
- `SVDCompression(D=χ, tol=ε)` → `SVDCompression(truncdimcutoff(D=χ, ϵ=ε))`。
- `DMRG1(trunc)` / `DMRG1(trunc=...)`：`trunc` 需为 `truncdim(D)` 或 `truncdimcutoff(D, ϵ)`；`DMRG1(truncerr(ϵ))` 现在会 `MethodError`。
- `alg.D` / `alg.ϵ` 访问改为 `alg.trunc.D` / `alg.trunc.ϵ`（仅限含 `D` 的方案）。
- `DefaultIntegrationTruncation` / `DefaultTruncation` → `DefaultKTruncation` / `DefaultITruncation`。
- `ToulouseIM(μ=x)` / `AndersonIM(μ=x)` → `...ϵ_d=x`；`psi.svectors` → `psi.s`；`easy_swap!` → `swap!`。

## 第六批更新

### Grassmann ordering 导出与缩写清理

- 删除 `src/lattices/grassmannordering.jl` 中全部 7 个缩写别名（`const AĀBB̄ = ...` 等），所有使用处改用完整 ordering 名（含 `MixedGrassmannLattice1Order` 默认 ordering 参数）。
- 只导出满足 **AdjacentConjugation** 的 ordering（6 个）：`A1Ā1B1B̄1`、`A1Ā1B1B̄1a1ā1b1b̄1`、`A1Ā1a1ā1B1B̄1b1b̄1`、`A2Ā2B2B̄2A1Ā1B1B̄1a1ā1b1b̄1a2ā2b2b̄2` 及两个 mixed-time ordering；8 个 GeneralConjugation ordering 不再导出（内部仍可用，测试入口通过 `using GTEMPO: ...` 显式导入）。

### GrassmannOrdering Unicode 规范统一

- 全仓库统一为组合宏形式（`A` + U+0304），消除同一 ordering 名在不同文件中"预组合（U+0100）/组合宏"混用导致的符号不一致：
  - 修正 `src/influencefunctional/hybridization/ttiif/imaginarytime.jl`、`src/influencefunctional/hybridization/exact_ttiif/imaginarytime.jl` 中因编码不一致而失联的 `index` 等方法特化
  - 统一 `src/grassmanntensor/` 中混用的 `ā` 变量名编码
  - 测试与文档中 8 个文件的宏字符一并规范化

## 第五批更新

### partialintegrate 重构与 DMRG1 缩放修复

- `src/partialintegrate` 扁平化：`integrateband`、`multintegrateband`、`utils`、`svdmult`、`itermult`、`multpartialintegrate` 全部由 `partialintegrate/partialintegrate.jl` 统一 include（原 `generalmultpartialintegrate/` 子目录删除）。
- `partialintegrate` 的 `alg` 改为关键字参数（默认 `DefaultMultAlg`），docstring 标注 **experimental interface**（将来可能大改）。
- 修复 DMRG1 路径的范数丢失：DMRG sweep 与 SVD 路径不同、没有逐步的缩放核算，finalize 后以精确积分校准整体缩放（`target/current` 乘入 `z[1]`，`target` 为 partialintegrate 定义上必须保持的 `∫(xs...)`）。
- API 测试同时覆盖 `SVDCompression` 与 `DMRG1`（`test/api/integration.jl`）。

### 删除三个 Grassmann ordering

连同定义、`ConjugationStyle`/`LayoutStyle` trait、专用 `index` 方法、导出与全部测试项一并移除：

- `A2Ā2A1Ā1a2ā2a1ā1B2B̄2B1B̄1b2b̄2b1b̄1`（real，band-local）
- `A2Ā2A1Ā1B2B̄2B1B̄1`（imaginary，band-local）
- `A1Ā1B1B̄1b̄1B̄1ā1Ā1`（real，含别名 `AaBbb̄B̄āĀ`）

### 算法类型层级整理

- 移除抽象类 `DMRGMultAlgorithm`：`DMRG1`/`DMRG2` 直接继承 `DMRGAlgorithm`，相关方法签名（`mult`/`mult!`/`iterativemult`/`compute!`/`sweep!`/`finalize!`/`parint_iterativemult`/`multintegrateband`/`_ac_partialintegrate`）全部改用 `DMRGAlgorithm`（依赖 Julia 方法按特异性分发，`SVDCompression` 等特化分支不受影响）。
- `DMRG1`/`DMRG2` 的定义移到 `src/algorithms.jl`（与 `SVDCompression` 同处）。
- `src/auxiliary/orth.jl` 的内容（`MatrixProductOrthogonalAlgorithm`、`Orthogonalize`）并入 `src/algorithms.jl`，原文件删除。
- `const DefaultMultAlg = DMRG1(DefaultITruncation)` 移到 `src/defaults.jl`；include 顺序调整为 `algorithms.jl` 先于 `defaults.jl`（`DefaultMultAlg` 的构造依赖 `DMRG1`）。

### 其它

- 移除 `ExactTTIIF` 未使用的 `algmult2` 字段。
- 新增 `DefaultExpansionAlg` 常量，统一 `XTRGIF`/`ExactTTIIF`/`TDVPIF` 及 ttiif、decayterm 相关接口的默认 Prony 展开参数；删除 `DefaultMPOTruncation`（原使用处改用 `DefaultITruncation` / `DefaultIntegrationTruncation`）。
- 注释掉 `partialintegrate` 中未使用的 `my_mult`/`my_mult2`。

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
- 旧 `Gτ`/`Gt`/`Gm`（含 cached、fast 版本）与 `gf`/`cached_gf`/`cached_gf_fast` 的逐条对应公式：

  记 `C(i; c, b, d) = ContourIndex(i, conj=c, branch=b, band=d)`（全部为关键字参数）；`band` 关键字为 `Int` 时 `band₁ = band₂ = band`，为 `(b₁, b₂)` 元组时分别取两带。

  **逐点版本 → `gf`**（旧默认 `alg = ExactIntegrate()`、`Z = integrate(lattice, A, B...; alg)`）：

  | 旧接口 | 等价新写法 |
  |---|---|
  | `Gτ(lat::ImagGrassmannLattice, i, j, A, B...; band=1, c1=false, c2=true, alg, Z)` | `gf(lat, (C(i; c1, :τ, band₁), C(j; c2, :τ, band₂)), A, B...; alg, Z)` |
  | `Gt(lat::RealGrassmannLattice, i, j, A, B...; b1, b2, band=1, c1=true, c2=false, alg, Z)` | `gf(lat, (C(i; c1, b1, band₁), C(j; c2, b2, band₂)), A, B...; alg, Z)` |
  | `Gm(lat::MixedGrassmannLattice, i, j, A, B...; b1, b2, band=1, c1=true, c2=false, alg, Z)` | 同上（`b1, b2 ∈ {:+, :-, :τ}`） |

  与旧版一致的约定：`Gτ(lat, i, A, B...; kwargs...) = Gτ(lat, i, 1, A, B...; kwargs...)`（`j` 默认 1）；Mixed 轮廓上 `Gt(lat, i, j, ...; b1, b2 ∈ {:+, :-}, ...)` 委托 `Gm`，`Gτ(lat, i, j, ...; c1=false, c2=true, ...)` 委托 `Gm(...; b1=:τ, b2=:τ)`。

  **全格点批量版本**（`k = lattice.k`；`parallel_Gτ` 公式与此相同，仅把 `integrate` 换成 `parallel_integrate`，已随 `parallel_run` 一并停用）：

  ```julia
  # Gτ(lat::ImagGrassmannLattice1Order, A, B...; band)：
  g = [Gτ(lat, i, A, B...; band) for i in 1:k-1]   # Matsubara G(τᵢ) = ⟨aᵢ b₁⟩
  push!(g, 1 - g[1])                               # 末点由等时关系补全
  ```

  **cached 逐点版本 → `cached_gf`**（`cache = environments(lat, A, B...)`，不再有 `alg/Z`）：

  | 旧接口 | 等价新写法 |
  |---|---|
  | `cached_Gτ(lat, i, j, A, B...; cache, band=1, c1=false, c2=true)` | `cached_gf(lat, (C(i; c1, :τ, band₁), C(j; c2, :τ, band₂)), A, B...; cache)` |
  | `cached_Gt(lat::RealGrassmannLattice, i, j, A, B...; cache, b1, b2, band=1, c1=true, c2=false)` | `cached_gf(lat, (C(i; c1, b1, band₁), C(j; c2, b2, band₂)), A, B...; cache)` |
  | `cached_Gm(lat::MixedGrassmannLattice, i, j, A, B...; cache, b1, b2, band=1, c1=true, c2=false)` | 同上（`b1, b2 ∈ {:+, :-, :τ}`；Mixed 上 `cached_Gt`/`cached_Gτ` 的委托关系同逐点表） |
  | `cached_contour_ordered_Gm(lat, i, j, A, B...; cache, b1, b2, band=1)` | `cached_contour_ordered_gf(lat, C(i; false, b1, band₁), C(j; true, b2, band₂), A, B...; cache)`，其中 `cached_contour_ordered_gf(lat, a, b, …) = a < b ? -cached_gf(lat, (b, a), …) : cached_gf(lat, (a, b), …)` |

  批量版本 `cached_Gτ(lat::Union{ImagGrassmannLattice, MixedGrassmannLattice}, A, B...; cache, band=1)`：`g[i] = cached_Gτ(lat, i, A, B...; cache, band)`（i = 1 … kτ−1），`g[kτ] = 1 - g[1]`。

  **fast 版本 → `cached_gf_fast`**（一次给出 `⟨aᵢ b₁⟩, i = 1 … N`，`N` 由轮廓与分支自动确定）：

  | 旧接口 | 等价新写法 |
  |---|---|
  | `cached_Gτ_fast(lat::ImagGrassmannLattice, A, B...; c1=false, c2=true)` | `cached_gf_fast(lat, A, B...; b1=:τ, b2=:τ, c1=c1, c2=c2)`（末点修正 `r[end] = 1 - r[1]` 内部处理） |
  | `cached_Gτ_fast(lat::MixedGrassmannLattice, A, B...; c1=false, c2=true)` | 上式之后再执行 `r[end] = 1 - r[1]` |
  | `cached_Gt_fast(lat::RealGrassmannLattice, A, B...; b1, b2, c1=true, c2=false)` | `cached_gf_fast(lat, A, B...; b1=b1, b2=b2, c1=c1, c2=c2)` |
  | `cached_Gm_fast(lat::MixedGrassmannLattice, A, B...; b1, b2, c1=true, c2=false)` | 同上（`b1 ∈ {:+, :-}` 时长度 `kt`、`:τ` 时长度 `kτ`；末点修正需调用方执行） |
  | `cached_greater_fast(lat, A, B...)` | `cached_gf_fast(lat, A, B...; b1=:+, b2=:+, c1=false, c2=true)` |
  | `cached_lesser_fast(lat, A, B...)` | `-cached_gf_fast(lat, A, B...; b1=:+, b2=:-, c1=false, c2=true)`（注意前导负号） |
- `parallel_run`、`parallel_integrate`（与 `parallel_Gτ` 一起）整文件注释停用（`src/integration/parallelintegrate.jl`）。
- `occupation2`（等时 Green 函数路径的占据数实现）已注释停用：与 Toulouse ED 严格解对比，Keldysh 轮廓略优于 `occupation`（0.35% vs 0.53%，均为离散化/截断噪声量级），但虚时轮廓在 i ≥ 2 存在约 6% 的错误跳变（仅 i = 1 恰好正确，边界 Grassmann 迹贡献未正确计入）。保留 `occupation`（insert_n 路径）与 `cached_occupation`。

### mult 迭代收敛判据（iterative_compute!）

- `iterative_compute!` 的收敛判据从"每轮 sweep 残差的相对标准差（`iterative_error_2`，即 σ/|μ|）"改为"相邻两轮 sweep 损失的相对变化"（与 ITensor/TeNPy/quimb/block2 的 DMRG 判据一致；第一轮总是执行，对应 `delta = 2*tol`）。
- 每轮 sweep 的损失定义为该轮 sweep 输出的**最后一个**残差 ‖mpsj‖（整轮 left+right 扫描结束、site 2 处的残差）；`iterative_compute!` 的返回值 `kvals` 相应变为每轮 sweep 的损失序列。理论上 sweep 内 ‖mpsj‖ 沿扫描方向单调递增，收敛后达到平台，故相邻轮损失差趋于零。
- 删除不再使用的 `iterative_error_2`。

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
- `band_boundary` 的 4 个死 method（realtime.jl 中 `A1Ā1B1B̄1b̄1B̄1ā1Ā1`、`A1B1ā1b̄1Ā1B̄1a1b1` 与两个 `A2...` ordering 版本）：这些 ordering 不在 `_AllowedRealGrassmannOrdering` 内或非 TimeLocal，`_fit_to_lattice` 永远不会走到；其余 4 个活跃 method（虚时 2 个、实时 2 个）保留。

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
