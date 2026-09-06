# 接口更改说明

本轮重构涉及的接口更改汇总。所有更改均已通过全量测试验证（行为不变或数值等价）。

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
  - 实时格子上驱动流的是 4 个分支影响算符之和（逐个相加并以紧的 `DefaultMPOTruncation` 压缩）
  - 费米符号全部继承自 `mult` 机制（`get_left_xy`、`updatemultleft/right`），键矩阵的有效映射实现为 AC 映射的投影 AL†·(H·z)
  - 测试见 `test/tempo/tdvpif.jl`（含 δ 收敛性验证与 in-place 合并验证）

## 导出更改
- `GrassmannBackend` 不再导出（仅模块内部使用；`@grassmann` 宏仍导出）。
- 其余导出名随上述重命名同步更新。
