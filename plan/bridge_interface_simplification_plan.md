# GTEMPO 桥层接口简化计划：对称性/标签/截断对齐 GrassmannMPS 原生接口

- 日期：2026-11-16
- 基线：GTEMPO HEAD `c01f446c` + 工作区未提交修正（bridge.jl / influence.jl / lattices.jl 的 2026-11-16 bands_guard 放宽与 var_charge 轨道交错修正——本计划在其上叠加）
- 参考实现：GrassmannMPS（`grassmannspaces.jl` 注册表 + `truncation.jl` 策略 + `canonical.jl` truncate!）与 TensorKit/TensorKitSectors 扇区代数
- 范围：**仅 GTEMPO**。U1GTEMPO2 / backup/（参考实现与历史备份，含各自 BoundedCache/LegacyTruncation 拷贝）零改动；extensions/exactif 用的是 Z2Tensors.TruncationScheme，不受影响。

## 0. 已裁决（2026-11-16 会话）

- **A1** 荷标签排序口径 = **选项 A（零扰动）**：标签全改 TK 扇区类型，但 qspace/split_cell 等所有 `sort` 走 `sortkey` 复现现序（Int 升序 / 元组字典序 / Z2Label 按 n）——键扇区插入序逐位不变，golden todense 与标量锚原样通过。
- **A2** 导出 API 返回类型变更**接受**：`var_charge`/`cell_charges`/`fock_charges`/`charge_dimensions` 返回值从 Int/NTuple/Z2Label 变为 TK 扇区；`symmetry_spec` 返回 `Type{<:Sector}` 而非单例。
- **A3** 删除 `Truncation`、`CompressionReport`、`LegacyTruncation`（bridge.jl 自造截断机器整体退役）；`compress!(x; trunc)` 调用点 → GM `truncate!`；`multiply` 调用点 → GM `mult`。
- **A4** 删除 `BoundedCache`（GTEMPO 内零调用点，非导出，删除零影响）。

## 1. 目标

桥层（`src/bridge.jl`）从「自造对称性单例 + 自造荷标签代数 + 自造截断机器」收缩为「TK 扇区类型即 spec 即标签 + GM 原生截断/乘法」，净删约 130–160 行自造机器；构造层、层几何、pairize/from_core/split_cell/monomial、integrate/EnvironmentCache 等领域语义全部保留。

## 2. 探查定谳（计划赖以成立的关键事实）

**F1【TK 扇区代数完备——自造标签层的前提已消失】**
四种扇区（`Z2Irrep`/`U1Irrep`/`Z2Z2`/`ChargeCharge`）自带全部所需代数：
- `dual`：Z2 自逆（`dual(c) = typeof(c)(N − c.n)`，znirrep.jl:82）、U1 取负荷（u1irrep.jl:30）、ProductSector 逐分量（product.jl:74）⇒ `qneg = dual`、`qsub(a,b) = qadd(a, dual(b))`；
- 融合：`first(⊗(a, b))` 即 qadd——GM 自身惯用法（grassmannbackend.jl:16、factors.jl:424、schurgmps.jl:188）；
- `==`/`hash`/`getindex`（ProductSector 逐分量，product.jl:20,224-225）标准定义——`Z2Label == Int 恒 false` 一类的包装类型陷阱整类消失。

**F2【TK 扇区 isless 与现标签序不同——必须 sortkey（A1 的依据）】**
`U1Irrep.isless` 按 |charge| 再符号（序为 0,1,−1,2,…，u1irrep.jl:46）；`ProductSector.isless` 按 findindex 距离和再字典序（product.jl:228）——都不是现在的 Int 升序/元组字典序。直接 `sort(sector)` 会改变 qspace/split_cell 的扇区插入序 → 张量块布局重排 → todense 数组级 golden 破。故所有排序走 `sortkey`（§4.1）复现现序。`ZNIrrep.isless` 按 `.n`（znirrep.jl:97）与 Z2Label 同序，但统一走 sortkey 不例外。

**F3【MAK 策略组合不支持 SectorVector——交集求值须由 GM 补一个方法】**
MAK 的 `TruncationIntersection`（`&`，truncation.jl:205）泛型实现返回**扁平索引**，而 TensorKit 的 svd_trunc 消费端要求 sector-aware inds（`SectorDict{Sector,UnitRange}`，TensorKit factorizations/truncation.jl:100-125 注释明言）——这正是 GM 与 bridge 各自定义 `findtruncated_svd(::SectorVector, …)` 的原因。`SectorVector <: AbstractVector`（sectorvector.jl:9），在 GM 补一个 SectorVector 交集方法（~8 行）后，`truncdimcutoff(D, 0) & truncrelerr(rtol)` 即可在 GTEMPO 直接用，GTEMPO 侧零自造策略。

**F4【GM truncate!/mult 语义与返回值】**
`GM.truncate!(x; trunc)` = `canonicalize!` = `leftorth!`（QR 前扫 + 入口 `_rescale_sites!`）+ `rightorth!`（SVD 回扫按 trunc 截断、Σs² 归一吸收 logscale），返回 `(x, truncerr)`，truncerr = 逐键相对截断误差之和（canonical.jl:128-139）。`GM.mult(x, y; trunc, verbosity)` = 逐站楔积 + truncate! + logscale(y) 吸收，返回 `(x, truncerr)`（mult.jl:34-45）。

**F5【截断参数映射——cutoff 无调用方，随退役删除】**
`Truncation(maxdim, rtol, cutoff)` 三条件合取：top-maxdim 全局预算 ∧ 相对尾容差（tails ≤ rtol²·total 断）∧ 绝对阈值（带 logscale 漂移补偿）。全库扫描 `cutoff=` 非零调用方**零个**（默认 0，PartialIF.cutoff=0.0）→ cutoff 旋钮随 LegacyTruncation 删除。映射：
- `Truncation(maxdim=D, rtol=r)` → `truncdimcutoff(D, 0.0) & truncrelerr(r)`；
  - `TruncateDimCutoff(D, 0)`：严格 top-D（v=0 条目与旧 cutoff_eff=0 同断）；
  - `TruncateRelError(r)`：从最小端丢弃至相对尾能 > r²·total——与 LegacyTruncation 的 tails 前缀判据逐条等价（单条目谱保 1 条、全零谱两侧都趋空，见 R2）；
- `Truncation(maxdim=typemax(Int), rtol=r)`（cleanup 路径）→ `truncrelerr(r)`。

**F6【reports 无 src 消费方——降为标量记账零风险】**
`CompressionReport`（discarded_norms/relative_errors/bond_dimensions/saturated/discarded_logabs）在 src 内**没有任何读取方**（grep 空）；唯一读者是 test/api/integration.jl 的记账断言。`InfluenceFunctional.reports::Vector{CompressionReport}`（influence.jl:533）字段改 `Vector{Float64}`（逐步 truncerr）。丢失逐键分解/saturated 旗/_logabs 记账，换得 GM 原生标量误差。

**F7【BoundedCache 零影响】**
scaling.jl:10-43 定义 + Base.empty! 扩展，GTEMPO 全树（src/test/validation/docs/extensions）**零调用点**、非导出。整段删除即可；头注释同步。

**F8【docs/tutorials 已 stale——换轨前即破】**
教程 12 处 `trunc::TruncationScheme`（Z2Tensors 类型）传入 GTEMPO `hybriddynamics` 等（trunc kwarg 声明 `::Truncation`）——换轨前传 Z2Tensors 实例即 TypeError，本就跑不通。本次一并机械迁移到新策略写法，不再单独维护 stale 面。

**F9【精确调用点清单】**
- `multiply`（3 处 src）：influence.jl:148、559、654；
- `compress!`（2 处 src）：influence.jl:557（cleanup，rtol-only）、system.jl:276（K 尾部压缩）；
- `Truncation` 签名面：influence.jl 129/136/503（_trunc 映射）/537/543/583/621/631/667；system.jl 150/160（`trunc=nothing` 无类型，流经 276）；
- `CompressionReport`/`reports`：influence.jl 136/533/545/553/632/672；system.jl 160/277；
- 导出面：GTEMPO.jl:89（Truncation/CompressionReport/LegacyTruncation）、:90（multiply/compress!）、:33（GM compress! import——“桥层扩展目标”，随 bridge compress! 删除一并撤）；
- spec/标签代数：bridge.jl 23-154（spec 单例 + Z2Label + qadd 塔 + sector/charge_label/sector_type/sector_spec + var_charge + station_space + qspace/pspace + cell_charges）；system.jl 103-106（parity_label，唯一调用 227）、165（bands_guard(symmetry_symbol…))、150/161（spec=Z2Symmetry() 默认）；测试 test/api/sysdynamics.jl:28（var_charge）、bridge.jl:11（charge_dimensions）；
- `symmetry_spec`/`symmetry_symbol` 调用：influence.jl 137/548/633、measurement.jl 51/90/135/167/198/333、lattices.jl 47/100/149/271、system.jl 348/361。

## 3. 改动设计

### 3.1 spec：四个单例 → 扇区类型本身

```julia
# spec = Type{<:Sector}（Julia 类型即标签惯用法）
const SYMMETRY_SPECS = IdDict{Symbol,Type}(
    :z2 => TK.Z2Irrep, :z2z2 => TK.Z2Z2, :u1 => TK.U1Irrep, :u1u1 => GM.ChargeCharge)
symmetry_spec(s::Symbol)::Type = get(SYMMETRY_SPECS, s, nothing) === nothing ?
    throw(ArgumentError("unknown symmetry :$s (allowed: :z2, :z2z2, :u1, :u1u1)")) :
    SYMMETRY_SPECS[s]
```

- 删除：`AbstractGrassmannSymmetry`、`Z2Symmetry`/`Z2Z2Symmetry`/`U1Symmetry`/`U1U1Symmetry`、`symmetry_symbol`（唯一消费方 bands_guard 改按类型分派）、`sector_type`、`sector_spec`、`chain_spec`（=`TK.sectortype(x[1])` 恒等，内联）。
- D1 裁决的 Symbol 手选 API（`symmetry=:` 关键字 + 未知 Symbol 显式拒绝）不变。
- `station_space` 收敛为注册表一查（GM 已导出 `canonical_grassmann_spaces`，GrassmannMPS.jl:40）：

```julia
"站形：band b 的物理站空间（单荷形全站同形；双荷形按带奇偶交替朝向）"
station_space(::Type{G}, band::Int) where {G<:Sector} =
    GM.canonical_grassmann_spaces(G)[isodd(band) ? 1 : 2]
```

- `bands_guard(::Type{G}, bands)`：双荷形（`G === TK.Z2Z2 || G === GM.ChargeCharge`，或按 `G <: TK.ProductSector`）要求偶数 bands，文案不变。

### 3.2 标签：Z2Label/qadd 塔 → TK 扇区 + sortkey（A1）

```julia
@inline qadd(a::Sector, b::Sector) = first(⊗(a, b))
@inline qneg(a::Sector) = dual(a)
@inline qsub(a::Sector, b::Sector) = qadd(a, dual(b))
@inline zero_charge(::Type{G}) where {G<:Sector} = one(G)

"荷排序键（A1 裁决：复现旧标签序——Int 升序 / 元组字典序 / Z2Label 按 n）"
sortkey(s::TK.Z2Irrep) = Int(s.n)
sortkey(s::TK.U1Irrep) = Int(s.charge)
sortkey(s::TK.Z2Z2) = (Int(s[1].n), Int(s[2].n))
sortkey(s::GM.ChargeCharge) = (Int(s[1].charge), Int(s[2].charge))
```

- 删除：`Z2Label`（含 4 个 Base 方法）、qadd/qsub/qneg 的 14 个分型方法、`zero_charge` 4 分型、`sector`/`charge_label` 8 个双向转换（标签即扇区，转换层消失）。
- 排序点全部 `sort(unique(labels); by = sortkey)`：qspace（bridge.jl:139）、split_cell（bridge.jl:417）。
- `sector(spec, q)`/`charge_label(spec, s)` 调用点（qspace/pspace/from_core/split_cell/monomial/_convert_site/charge_dimensions）一律去包装直用扇区。
- `charge_dimensions(x)`（导出保名）：Dict 键直接为扇区（`Dict(q => TK.dim(s, q) for q in TK.sectors(s))`）。
- system.jl `parity_label`（103-106）按扇区重签（4 个一行方法）；备选：给 GM grassmannspaces.jl 补 `fermionparity(::Z2Irrep) = isodd(s.n)`（注册表收尾一行）后系统层直用 `GM.fermionparity`——**主路线选 parity_label 重签**（不动 GM，零耦合），GM 补法记为后续可选。
- `var_charge`/`cell_charges`/`fock_charges`/`station_flavor` 重签为扇区返回值（逻辑逐行不变，`:z2z2` 的 `q[2].n` 拆分照旧）。

### 3.3 截断：LegacyTruncation 退役 → GM 策略 + GM 交集（A3）

**GM 侧（前置小改，~8 行）**：truncation.jl 补 sector-aware 交集求值：

```julia
function MAK.findtruncated_svd(values::SectorVector, s::TruncationIntersection)
    keep = MAK.findtruncated_svd(values, first(s.components))
    for t in Base.tail(s.components)
        other = MAK.findtruncated_svd(values, t)
        keep = SectorDict(c => (1:min(length(get(keep, c, 1:0)), length(get(other, c, 1:0))))
                          for c in keys(values))
    end
    keep
end
```

（组件均返回 SectorDict 前缀区间，交集 = 逐扇区 min 长度；更具体的 SectorVector 方法覆盖 MAK 泛型无歧义。）

**GTEMPO 侧**：
- 删除 bridge.jl:179-319（Truncation、CompressionReport、LegacyTruncation、compress! 全段）与 multiply（325-332）。
- bridge 层保留一个默认策略常量：`const DEFAULT_TRUNC = truncdimcutoff(128) & truncrelerr(1e-12)`（= 旧 `Truncation()` 默认）。
- `_trunc(alg::PartialIF)` → `truncdimcutoff(alg.maxdim) & truncrelerr(alg.rtol)`（alg.cutoff 字段删除——F5 无调用方；PartialIF 结构体删 cutoff 字段）。
- 调用点替换（签名 `trunc::Truncation=Truncation()` → `trunc=DEFAULT_TRUNC`，类型放宽为 `MAK.TruncationStrategy`）：
  - influence.jl:148/559/654 `out = multiply(...)` → `out, err = mult(out, …; trunc)` + `push!(truncerrs, err)`；
  - influence.jl:557 cleanup `compress!(f; trunc=Truncation(maxdim=typemax(Int), rtol=…))` → `truncate!(f; trunc=truncrelerr(…))`（取 err 忽略态或记账）；
  - system.jl:276 `r = compress!(out; trunc)` → `_, err = truncate!(out; trunc)`，`push!(truncerrs, err)`；
  - reports 传参全线改 `truncerrs::Vector{Float64} = Float64[]`；`InfluenceFunctional.reports::Vector{Float64}`（F6）。
- GTEMPO.jl：import 撤 `compress!`（:33）、增 `truncate!`、`TruncateDimCutoff`、`TruncateRelError`；export 删 `Truncation/CompressionReport/LegacyTruncation/multiply/compress!`，增 `truncate!`、`TruncateDimCutoff`、`TruncateRelError`（truncdimcutoff/truncrelerr/NoTruncation 已在 :70）。`import MatrixAlgebraKit as MAK` 保留（策略类型标注用）。

### 3.4 BoundedCache 删除（A4）

scaling.jl:10-43 整段删除；文件头注释删「与有界缓存」。无其它改动（F7）。

### 3.5 机械跟随面

- system.jl：`spec=Z2Symmetry()` 默认 → `spec=TK.Z2Irrep`（150/161 两处）；`bands_guard(symmetry_symbol(spec), bands)` → `bands_guard(spec, bands)`（165）；`parity_label`（103-106）；`fock_charges`/`local_weights` 重签随标签类型。
- influence.jl / measurement.jl / lattices.jl 的 `symmetry_spec(lattice.symmetry)` 调用：调用形不变，返回类型变 Type——下游签名（`spec::AbstractGrassmannSymmetry` → `spec::Type{<:Sector}`）机械重签。
- 测试：test/api/bridge.jl:28-32（multiply+Truncation → mult+策略）、test/api/integration.jl:16-35（“multiply / compress! / 截断记账” testset 重写为 truncate!/mult + truncerr 断言，报告类型断言删除）、sysdynamics.jl:28（var_charge 返回值断言随扇区）、lattice.jl / port_smoke.jl / entropy.jl 涉 spec/标签处机械跟随。
- docs/tutorials 12 文件（F8）：`TruncationScheme` → GM 策略写法，机械迁移。

## 4. 语义变化与风险清单

| # | 变化 | 性质 | 处置 |
|---|---|---|---|
| S1 | 截断 tie-break：LegacyTruncation 按 (−v, charge_label, j)，GM 策略按 (−\|v\|, 扇区规范序, j) | 仅跨扇区 \|v\| **精确相等**时选择不同（浮点测度零）；逐扇区前缀语义使非 ties 完全一致 | 记录即可；golden 回归确认 |
| S2 | 绝对 cutoff（logscale 漂移补偿版）能力消失 | 全库无调用方（F5） | 随退役删除；未来需要时在 GM 加 TruncationByValue 的 SectorVector 方法 |
| S3 | reports（逐键分解/saturated/_logabs）→ 标量 truncerr | src 无消费方（F6） | 字段改 `Vector{Float64}`；测试断言改写 |
| S4 | compress! 逐步 rescale（QR 每步重归一）→ GM leftorth! 入口一次性 `_rescale_sites!` | 每步 QR 的 carry r 为 O(1)（上一步张量已被入口归一），幅值不级联——入口一次 + 每步 O(1) 增量在效果上等价；且新路径下每次 multiply 都过 truncate!，站点幅值常态有界 | longchain_benchmark 对照确认无溢出/精度回归（R1） |
| S5 | 全零链：旧 compress! 有零链早退（保键维填零）；GM truncate! 对零谱可能缩到 0 维键 | partial_row/cell_maps 的 coeff=0 零链是真实可达路径 | R2：Phase 2 加零链回归测试；若 GM 行为破坏后续楔积，在 3 个调用点加 `iszero(norm)` 跳过守卫 |
| S6 | 导出 API 返回类型（A2 已接受） | symmetry_spec→Type、var_charge/cell_charges/fock_charges/charge_dimensions→扇区、multiply/compress!/Truncation/CompressionReport/LegacyTruncation export 撤销 | changes.md 记录 breaking 变更清单 |

## 5. 实施阶段与验收闸门

- **Phase 1（spec/标签）**：§3.1 + §3.2 + §3.5 的 spec/标签跟随（不含截断）。验收：全测试套件绿；golden（:z2 通道 todense 逐位 + 标量机器精度）原样通过——布局逐位不变是本阶段断言（A1）。
- **Phase 2（截断/乘法退役）**：GM 交集方法（§3.3 前置，GM 独立提交 + GM 自身测试绿）→ GTEMPO 截断机器删除 + 调用点替换 + reports 降标量 + 测试重写。验收：全测试绿；golden 通过（截断参数不变时 truncation 决策应与 S1 之外逐位一致）；longchain_benchmark 对照无回归（R1）；零链回归（R2）。
- **Phase 3（清尾）**：BoundedCache 删除 + docs/tutorials 迁移 + changes.md breaking 清单。
- 每阶段独立提交，U1GTEMPO2 / backup / extensions/exactif 全程零改动。

## 6. 净效果估计

- bridge.jl：523 行 → 约 360–380 行（spec/label 段 ~120 行压到 ~40，截断段 ~140 行整体撤）；
- GTEMPO 自造截断/缓存机器归零：截断 = GM 策略组合，乘法 = GM.mult，正则化 = GM.truncate!；
- 保留不动：var_charge 交错约定、bands_guard、bcs_symmetry_guard、cell_charges/station_of/site_spaces/pspace、pairize/_convert_site/_fuse_adjacent、from_core/split_cell/monomial、_occupation_index、integrate_scaled/EnvironmentCache、缩放标量（GM.ScaledScalar 皮肤）。
