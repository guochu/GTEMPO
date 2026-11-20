# 合并计划（v2）：插入算符符号化 + fock_sysdynamics 穿戴 + fockmpo2gmps 荷账自读

> **本文件合并取代两份 v1 计划**（`fockmpo_dress_selfread_plan.md` v1 + 
> `insertion_localop_plan.md`，均已废弃）。合并裁决（2026-11-20 用户）：
> **无独立 dressing 层**——v1 提案的 `_fockmpo_dress` 没有意义，insertions/
> initial_state 直接加进 `fock_sysdynamics` 签名，装配与穿戴一个函数完成。
> **计划成文待实施；本会话只定计划未改源码。**

- 日期：2026-11-20（两会话五轮调整收敛后合并）
- 基线：GTEMPO HEAD `a6672a84`（if-expmult-backend；R7 定谳后）
- 波及面：`src/sysdynamics/sysdynamics.jl`（OperatorInsertion/fock_sysdynamics）、
  `src/sysdynamics/fockmpo.jl`（fockmpo2gmps 重写 + `_fockmpo_with_insertions` 删除）、
  **新文件 `src/sysdynamics/localops.jl`**、`src/measurements/inserted_k.jl`（17 构造点）、
  `src/GTEMPO.jl`（export）、`test/api/{fockmpo,sysdynamics,lattice}.jl`（38 构造点）
- 不碰：IF 侧 `correlationfunction.jl`/`exactif_port`、`gmps2fockmpo` stub、
  测量缓存/摊销、ρ₀/`FockMatrix` 路线
- 命名冲突已查：GTEMPO/ImpurityModelBase src 均无 `LocalOp`/`FockOp`

## 0. 已裁决与裁决点

**已裁决（2026-11-20 用户多轮提案）：**

*符号算符侧：*

- **V1（j 改名）**：field `boundary` → `j`，直接代表虚时/实时**时间步数**；实时不论
  branch，第 j 时间步输入 j。**解读钉定**：j = 物理时间步口径（(+,j) 与 (−,j) 同一
  物理时间点），`insertion_cut_index` 的 − 支反向映射**不变**。
- **V2（符号算符）**：`matrix` field → `op::FockOp`。
- **V3（charge 删除）**：荷从符号按位派生，无字段、无手工标签——"prefix 真凶"类
  标签错误从源头消灭。
- **V4（FockOp 无乘法）**：FockOp 之间（及与 LocalOp 混合）不定义乘法——被动容器。
- **V5（LocalOp 乘法封闭）**：`LocalOp * LocalOp → LocalOp` 符号级乘法表封闭；
  新增 `:empty` = 零算符（a·a = 0 等封闭性要求）。
- **V6（位置即带）**：**LocalOp 不携带 band 信息**（纯 `s::Symbol`）；**FockOp.data
  长度 = bands，data[b] = 第 b 带因子**（1-based）。带间积算符序 = 带序（band 1 最左，
  `local_weights` fold 同款）。
- **V7（装配点物化）**：不保留 raw-matrix 构造器；FockOp 物化发生在**生成 FockMPO 时**
  （`sparsefockmpo`/三明治装配点），产物直接写 TensorMap 站张量块，插入对象全程符号态。

*管线侧（合并裁决）：*

- **P1（无 dressing 层）**：`fock_sysdynamics(lattice, model; order=1,
  insertions=OperatorInsertion[], initial_state=nothing)` 一次产出**穿戴完毕**的链
  （开链装配 → 切位 ≥1 稀疏乘法 → ρ₀/切位 0 三明治折首站 → oddlater 逐站号）。
  v1 的 `_fockmpo_dress`/现行 `_fockmpo_with_insertions` 均不存在；
  `fock_initstate!` 删除（无独立折叠路径）。
- **P2（converter 纯翻译 + 荷账自读）**：`fockmpo2gmps(f, lattice; trunc, truncerrs)`
  ——槽荷表 + 站读数 δ 自读 prefix（§2.2），不消费任何插入/初态元数据。
- **P3（边界态恒折 Fock 层）**：`fock_vacuumstate` 缺省（imag/mixed = I、real = 真空
  投影 P）、`initial_state` 仅实时轮廓覆盖——**恒折**，converter 首格直读不折任何态。
- **P4（oddlater 折进链）**：逐站 slot 对角 Koszul 号在 `fock_sysdynamics` 内折入，
  奇偶源 = 符号荷；converter 无 sgn。
- **P5（符号荷 vs 链 δ 交叉自检）**：`fock_sysdynamics` 末尾把符号荷累计与穿戴链的
  δ 读数逐项核对（不一致 = 物化 bug——站内一致性检查抓不住的"一致但错"缺口，
  如 JW 符号错折）。converter 侧独立再读，双路径互证。

**裁决点（唯一待裁；推荐即缺省）：**

- **A1′（电荷非对角初态语义）**：现行对 initial_state 的荷静默忽略；自读路线的
  converter 站内一致性校验会**抛错**。**推荐：抛错**（荷流必须闭合；现行支持面
  I/真空投影子/对角态全部电荷对角、行为不变）。

*v1 遗留裁决的归宿*：dress v1-A1（初态抛错）= A1′；dress v1-A2（删 `fock_initstate!`）
= P1 结构性强制；dress v1-A3（`_fockmpo_dress` 更名）= P1 溶解；localop-A1（位置即带）
= V6；localop-A2（纯符号）= V7；localop-A3（编排顺序）= 合并溶解。

## 1. 目标形状

```
OperatorInsertion(branch, j, FockOp(...); side)         # 全程符号态
    ↓
fock_sysdynamics(lattice, model; order, insertions, initial_state)
    = 开链装配 → 切位≥1 稀疏乘法 → ρ₀/切位0 三明治折首站（在稀疏乘法之后！）
      → oddlater 逐站 slot 号 → 符号荷 vs 链 δ 交叉自检（P5）
    → 穿戴完毕 FockMPO（ρ₀ 已折、仍开链——迹闭合归 converter）
    ↓
fockmpo2gmps(f, lattice; trunc, truncerrs)
    = 槽荷表 + 站读数 δ + 站内一致性 + Σδ=0 环闭合守卫 + labels
      + 首格直读（无 ρ 收缩、无 sgn）→ split_cell/pairize/reorder → GMPS
    ↓
_system_factor_fock：setlogscale!(out, logscale + _fock_weyl_lscale(...))   # 不动
```

**接口 diff（导出面）：**

| 符号 | 现行 | 目标 |
|---|---|---|
| `fock_sysdynamics`（export） | `(lattice, model; order)` 开链 | `(lattice, model; order, insertions, initial_state)` 穿戴链 |
| `fockmpo2gmps`（export） | `(f, lattice; rho0, insertions, trunc, truncerrs)` | `(f, lattice; trunc, truncerrs)`，**breaking** |
| `fock_initstate!`（export） | 原地 ρ 折叠 + 切位 0 警告 | **删除** |
| `fock_vacuumstate`（export） | — | 保留（缺省边界态） |
| `OperatorInsertion`（export） | `(branch, boundary, matrix, charge, side)` | `(branch, j, op::FockOp, side)` |
| `LocalOp`/`FockOp` | — | **新增 export** |
| `sparsefockmpo`/`SparseFockMPO`/`Base.:*` | — | 保留（fock_sysdynamics 内部机制 + 闸门手工路线） |
| `sysdynamics` 公开皮肤 | kwargs 透传 | **签名不变**（测量层零改动调用面） |
| `_fockmpo_with_insertions` | 切位 ≥1 乘法 | **删除**（吸收进 fock_sysdynamics） |

v1 三红线归宿：红线 1（三明治在切位 ≥1 乘法**之后**折——切位 1 插入作用在站 1，
`(Ol₁·M₁·Or₁)·(Ol₀·ρ·Or₀)` 语义）**保留**；红线 2（real 缺省投影子整体随迁、
converter 彻底不折）= P3 **保留**；红线 3（deepcopy 防别名穿透）**消失**——
装配与穿戴同函数，链为函数内新建，无外部输入可污染。

## 2. 语义钉定（必须按此实现）

### 2.1 符号算符

**j/branch 定位**：`insertion_cut_index`（`sysdynamics.jl:181-191`）逐字保留，字段改
`ins.j`：τ j → cut j；+ j → ni + j；− j → ni + 2Nt − j（反向）。(τ,0)/(+,0) = ρ₀
切位（三明治）、(−,0) = 链末切位（合法物理插入）——过滤按切位，不按 j 字段。
`side`：`:left` = O·M（轮廓更晚读法）、`:right` = M·O。构造守卫增 `j ≥ 0`；上限
（j ≤ 该分支格数）在 `insertion_cut_index` 补校验（越界提前抛）。

**LocalOp 六符号（纯 Symbol，无 band——V6）**：

| s | 单带矩阵 | 荷 |
|---|---|---|
| `:a` | d_b | −q_b |
| `:adag` | d_b† | +q_b |
| `:n` | d_b† d_b | 0 |
| `:aadg` | d_b d_b†（= I − n_b，物化为积保留语义） | 0 |
| `:I` | I | 0 |
| `:empty` | 0 | 0 |

单带矩阵引擎 = `_fock_annihilators(bands)[b]`（band 1 = MSB、JW 穿越更低索引带）。
`:aadg`/`:n`/`:I` 两两矩阵互异，六符号无一冗余。

**LocalOp 乘法表（V5；符号级，行 = 左因子、列 = 右因子，右先作用）**——单带
span{I, d, d†, n} = M₂ 代数 + 零吸收：

| · | I | a | adag | n | aadg | empty |
|---|---|---|---|---|---|---|
| **I** | I | a | adag | n | aadg | empty |
| **a** | a | empty | aadg | **a** | empty | empty |
| **adag** | adag | n | empty | adag | empty | empty |
| **n** | n | empty | adag | n | empty | empty |
| **aadg** | aadg | a | empty | aadg | empty | empty |
| **empty** | empty | empty | empty | empty | empty | empty |

非对称注记：`a·n = a` 而 `n·a = empty`（d·d†·d = d vs d†·d·d = 0）——按行列序查。
结合律成立（单测）。自动化简**不做**；跨带组合走 FockOp 位置排布（不是乘法）。

**FockOp（V4/V6）**：`FockOp(data::Vector{LocalOp})`，`length(data) == bands`
（物化点校验，构造不自校验），`data[b]` = 第 b 带因子。物化（V7）= 逐带矩阵 fold
**band 1 最左**（`O = O * mat_b, b = 1..bands`，`local_weights` 同款）；每带全空间
矩阵自带 JW 符号，普通矩阵积符号天然正确——**不加任何 graded⊗ 符号**（
`local_weights` 的 graded⊗ 因子是槽编码装置、现行 `graded_sign_factor = 1` 关闭，
与算符代数无关）。`data` 含 `:empty` ⟹ 整积零矩阵（物化短路、荷 0）。
跨带同积算符序 = 带序（设计属性；其它序 = 未来扩展点，现行无需求）。

**荷派生（V3）**：`_fockop_charge(op, G, bands) = Σ_b charge(data[b], b)`；
奇偶 = `GM.fermionparity`。

**total-n̂ 归宿（V6 推论，采纳）**：Σ_b n_b 不可表达为单 FockOp——
`occupation(band=nothing)` 改逐带求和 Σ_b occupation(band=b)（期望线性、数值恒等）。
单带因子 helper `_single_band_op(s, band, bands)`（目标带放符号、其余 `:I` 填充）
= 测量 API 主构造路径。

### 2.2 管线：穿戴 + 荷账自读

**穿戴（fock_sysdynamics 内，P1）**——顺序即红线：

1. 开链装配（现行主体不动：`chronological_cells` 逐格 S/W 块对角张量）；
2. 切位 ≥1 稀疏乘法（`sparsefockmpo` 过滤 cut ≥ 1、左右各一、`Wr*(Wl*K)`；
   过滤按 `insertion_cut_index` 非 j 字段）；
3. ρ₀/切位 0 三明治折首站（**在步骤 2 之后**）：ρ = `fock_vacuumstate(lattice)`、
   `initial_state` 仅实时覆盖（现行 `_system_factor_fock:460` 门控随迁）；三明治序
   foldr/foldl（现行 `fockmpo.jl:311-316` 逐字）；右乘进首站全部 slot 块；
   现行首格核 `(M₁·Ol₀·ρ·Or₀)[s,t]` 与链上预折逐位恒等（所有 slot 右乘同一矩阵）；
4. oddlater 逐站 slot 对角号（现行 `fockmpo.jl:336-342` 谓词逐字搬运：`cut ≥ k`、
   奇偶源 = 符号荷）；标量因子对线性核构造可提出，折进链与折进核逐位恒等
   （L==1/首/末/中四分支同）；
5. P5 交叉自检：符号荷累计（prefix(键 k) = Σ_{j≤k} δⱼ^symbol）vs 链 δ 读数
   （§2.2 下）逐项一致 + Σδ = 0，不一致 `ArgumentError`。

**converter 荷账自读（P2）**——机制三件：

- **槽名义荷表（先验）**：槽编码 m 每带 2bit 组 `slots[b] = (m >> 2(b−1)) & 3`
  （band 1 = 低位组）：`0 = q·I−n`（0）、`1 = −d†`（+q_b）、`2 = d`（−q_b）、
  `3 = W₄`（0）——`q_m = Σ_b ([1]−[2])·var_charge(G, b, bands)`，纯函数。
  带间 graded⊗ 符号与 JW ±1 均纯标量，不影响稀疏性。
- **站读数**：站张量腿序 = (键L ⊗ 单项式出) ← (单项式入 ⊗ 键R)，槽对角。未 dressing
  链 `blk[s, m+1, m+1, u] ≠ 0 ⟺ f_q(s) − f_q(u) = q_m`；穿戴后移位 `q_m + δ_k`：
  **δ_k = f_q(out) − f_q(in) − q_m**（任一非零条目；容差 = 1e-12 × 站内 max|blk|，
  全条目低于容差 → 零链早退，保持现行 `isempty(column)` 语义与链长结构）。
  ρ₀ 不进 δ（现行 prefix 同样只让插入荷进；电荷对角 ρ 对首站读数零污染——
  对角 ρ 下 `(B·ρ)[s,t] = B[s,t]·ρ[t,t]`，s−t 差不变）。
- **prefix 恒等式**：现行键 k 处 prefix = 切位 0..k 插入荷之和
  （`fockmpo.jl:319,331-335`）；穿戴后切位 0 → 站 1、切位 k → 站 k，故
  **prefix(键 k) = Σ_{j≤k} δ_j**，逐位等价。
  验证例（闸门"多插(偶+奇)"：n̂@cut0 + d@cut1(−q) + ddag@cut2(+q)）：
  键 1 = δ₁ = −q ✓、键 2 = δ₁+δ₂ = 0 ✓（与现行 cut 累计一致）。

**converter 站内自验证三件套**：

1. 站内一致性：同站不同条目/槽的 δ 读数全等，不一致 = 装配电荷不相干
   （电荷非对角初态，A1′）→ `ArgumentError`；
2. Σ_k δ_k = 0 环闭合守卫 ≡ 现行总荷守卫（:303），源 = 矩阵；
3. slot 对角断言：离对角 slot 条目 ≈ 0（读数前置检查）。

**开发期锚（Phase 3 先立，改造前后都绿）**：未 dressing 链（现行 `fock_sysdynamics`
直出）全站 δ ≡ 0——槽荷表符号约定与 `local_weights` 逐位核对的机器可检不变量。

**labels 公式不变**：首/末格 [z]、内部 `qsub(qsub(fqs1[s], fqs1[t]), prefix)`
（prefix 改 δ 累计）；首格分支直读 `blk'[s, m+1, m+1, t]`；`split_cell`/`pairize`/
`reorder!`/trunc 透传全不动（`fockmpo.jl:373-391`）。

## 3. 设计

```julia
# src/sysdynamics/localops.jl（新文件；include 于 OperatorInsertion 定义之前）
struct LocalOp
    s::Symbol                    # :a :adag :n :aadg :I :empty（无 band——V6）
end
struct FockOp
    data::Vector{LocalOp}        # data[b] = 第 b 带因子；length == bands（物化点校验）
end
Base.:(*)(x::LocalOp, y::LocalOp)::LocalOp             # 2.1 表；唯一乘法
_single_band_op(s::Symbol, band::Int, bands::Int)::FockOp
_fockop_matrix(op::FockOp, bands::Int)                 # V7：装配点物化 → TensorMap 块
_fockop_charge(op::FockOp, ::Type{G}, bands::Int) where G<:Sector

struct OperatorInsertion
    branch::Symbol
    j::Int                       # 原 boundary：时间步数
    op::FockOp
    side::Symbol
    # 守卫：j ≥ 0；side ∈ (:left,:right)；charge 字段删除（V3）
end
```

```julia
# fock_sysdynamics（P1）：签名扩展 + 穿戴内联
function fock_sysdynamics(lattice, model; order=1, insertions=OperatorInsertion[],
                          initial_state=nothing)
    # 1. 开链装配（现行主体）
    # 2. 切位≥1 稀疏乘法（现行 _fockmpo_with_insertions 主体迁入）
    # 3. ρ₀/切位0 三明治折首站（initial_state 仅实时覆盖；步骤 2 之后！）
    # 4. oddlater 逐站 slot 号（符号荷）
    # 5. P5 交叉自检
end

# converter（P2）：签名收敛 + 自读
function fockmpo2gmps(f::FockMPO, lattice; trunc=nothing, truncerrs=nothing)
    # 删：rho0/insertions kwargs、cuts dict、总荷守卫、三明治折叠、
    #     prefix-自-insertions、oddlater、四分支 sgn（现行 :296-342）
    # 增：_slot_charge 表 + _site_charge_shift 读数（~20-25 行）+ 三件套
end

# _system_factor_fock 收敛为三行：
Kf = fock_sysdynamics(lattice, model; order, insertions, initial_state)
out = fockmpo2gmps(Kf, lattice; trunc, truncerrs)
GM.setlogscale!(out, GM.logscale(out) + _fock_weyl_lscale(lattice, order, coupling_norm))
```

`==`/`hash`：OperatorInsertion 结构化（branch、j、side、op 符号逐位）；物化同矩阵
不同符号形视为不同对象——缓存多算一份、数值不受影响。`sysdynamics` 公开皮肤签名
不变（kwargs 透传）。`FockMPO` struct docstring 更新：链 = 穿戴态（ρ₀ 已折）、
仍开链（迹闭合归 converter）。

## 4. 调用面改造清单

**src/measurements/inserted_k.jl（6 函数 17 构造点；`qc`/`var_charge` 局部全删）：**

| 函数 | 现行 | 目标 |
|---|---|---|
| `occupation` (:72) | `_number_op(bands, band)`；nothing = Σn_b | band=Int：`_single_band_op(:n, band, bands)`；**nothing：逐带求和**（线性恒等） |
| `raw_correlation` (:103-108) | d/ddag + charge= | `_single_band_op(:a/:adag, 1, bands)` |
| `_normalized_two_point` (:146-151) | 同上 | 同上 |
| `matsubara` (:173-174) | τ 支 d/ddag + charge | 同上 |
| `nn` (:195-198) | `_number_op(bands, band₁/₂)` | `_single_band_op(:n, band1/2, bands)` |
| `raw_lead_correlator` (:235-236) | ddag/d + charge | `_single_band_op(:adag/:a, band, bands)` |

`_number_op` 删；`_fock_ops` 在 inserted_k 的用途删（`_fock_annihilators` 留包内作物化
引擎）；`electriccurrent`/`heatcurrent`/`_inserted_expectation` 不碰。

**src/sysdynamics/sysdynamics.jl**：OperatorInsertion 定义（:113-134）、`==`/`hash`
（:141-153）、`insertion_cut_index`（:181-191）+ docstring；fock_sysdynamics 签名 +
穿戴内联。

**src/sysdynamics/fockmpo.jl**：converter 重写（删 :296-342 插入/初态 plumbing、增
δ 读数）、`_fockmpo_with_insertions`（:419-425）与 `fock_initstate!`（:164-175）删除、
`fock_vacuumstate` 保留、FockMPO/文件头 docstring 更新。

**src/GTEMPO.jl**：export 增 `LocalOp`/`FockOp`，删 `fock_initstate!`。

**tests（38 构造点）**：`fockmpo.jl` 21、`sysdynamics.jl` 15、`lattice.jl` 2——
`GTEMPO._fock_ops`/`GTEMPO._number_op` 换 `_single_band_op`/符号；charge= 全删；
`fock_initstate!` 双路测试删除；闸门手工闭环重写（§5 Phase 3）；total-n̂ 配置换逐带
（对拍双侧同换）。

## 5. 阶段与验收

**Phase 1（纯增量：localops）：** 类型 + 乘法表 + 物化 + 荷派生 + 单测：
- 乘法表 36 项 vs 手工 2×2 矩阵乘逐位；结合律；`:empty` 短路；
- 物化：单带因子 FockOp（`:I` 填充）vs 现行 `_fock_ops`/`_number_op` 矩阵**逐位 ==**
  （bands{1,2} × 六符号 × 每带位置）；fold 序锚 `[:a, :adag]` ≡ `d₁*d₂†` 逐位；
  JW 实证（d₁·d₂ = −d₂·d₁）；荷表 vs `var_charge`；push 序 ≡ 拼接 FockOp 等价。

**Phase 2（符号切换）：** OperatorInsertion（j/op/charge 删）+ inserted_k/tests 重写 +
converter kwargs 路线荷源换 `_fockop_charge`（装配管线其余不动）。
**硬锚 = 物化矩阵 vs 现行手工矩阵逐位 ==** ⟹ 下游数值逐位不变；闸门矩阵 +
inserted_k 电池全绿（total-n̂ 逐带恒等另测）。

**Phase 3（管线）：** fock_sysdynamics 穿戴内联 + converter 自读重写 +
`_fockmpo_with_insertions`/`fock_initstate!` 删除 + export/docstring。
先立锚（不改语义即可测）：未 dressing 链 δ ≡ 0（全轮廓 × bands × spec × order）；
单带荷插入 δ = 符号荷、他站零；prefix 恒等式（手工 cut 累计 vs δ 累计）。
然后闸门重写：手工闭环 = fock_sysdynamics 裸调用（真空穿戴）→ 手工稀疏乘法
（公开 `sparsefockmpo *`）→ 手工三明治折首站 → 手工 oddlater 号 → `fockmpo2gmps`——
对拍公开入口（穿戴 + 自读两条独立簿记路径互验）。抛错面：电荷不相干 dressing →
ArgumentError；(−,0) 不误入切位 0 支。

**Phase 4：** 全套件（终点 = Part4 Kanamori testset）+ observables golden
（U1GTEMPO2 机器精度）+ ED 电池。**任何位级偏差 = 停。**

## 6. 不做清单

- FockOp 乘法、LocalOp↔FockOp 混合乘法（V4）；自动规范化/化简（2.1）。
- raw-matrix 插入构造器（V7 定谳：纯符号）。
- ρ₀/initial_state/`FockMatrix` 路线本体不动（仅物化点随穿戴迁入 fock_sysdynamics）。
- `gmps2fockmpo` stub、测量缓存/摊销、IF 侧 `correlationfunction.jl`/`exactif_port`
  不碰。
- 跨带非带序算符序的 FockOp 扩展（2.1 设计属性，现行无需求）。
- 装配的 split_cell/pairize/reorder 机制不动（converter 只换荷账来源与首格分支）。

## 7. 工作量与顺序

Phase 1 ~1.5h → Phase 2 ~2h → Phase 3 ~3h → Phase 4 ~1.5h（单会话可完或分两会话）。
硬约束：Phase 1 单测全绿才动 Phase 2；Phase 2 硬锚（物化逐位 ==）不过 = 停；
Phase 3 锚先行（δ ≡ 0 未 dressing 链）再重写 converter。

## 8. 实施状态

**✅ 已实施（2026-11-20，本会话；A1′ 按推荐 = 抛错落地）。** 四阶段全部完成，
Part 1 API 闸门全绿（含 fockmpo 对拍矩阵 300 项 + localops 单测 313 项 +
golden U1GTEMPO2 机器精度锚 + ED 电池）；Part 2–4 全套件见套件运行记录。

实施偏差（相对计划文本，均有据）：

- **乘法表 4 项修正**：§2.1 表中 (adag,n)/(adag,aadg)/(aadg,n)/(aadg,aadg)
  四项与 2×2 矩阵代数矛盾（两行的 n↔aadg 列互换——代数：d†·n = 0、
  d†·(I−n) = d†、(I−n)·n = 0、(I−n)² = I−n）。实现按矩阵代数落地，
  Phase 1 的 36 项对拍闸门（计划自定验收判据）全绿。
- **闸门配方修正**：Phase 3 的"裸调用（真空穿戴）→ 手工稀疏乘法 → 手工
  三明治"与 P3 冲突——裸调用已恒折边界态（真空投影子），其上再折 ρ₀ =
  双重折叠且乘法序不可复原。闸门改为复刻步骤 ①原始装配（未折 ρ）再走
  手工稀疏乘法/三明治/oddlater，独立簿记路径语义不变。
- **P5 补充钉定**：切位 0 插入（三明治）计入站 1 的 expected 荷
  （`cut == k || (cut == 0 && k == 1)`）——prefix 恒等式 δ₁ 含切位 0/1。
- **A1′**：电荷非对角初态 → `_site_charge_shift` 站内一致性抛 ArgumentError
  （fock_sysdynamics 内先于 converter 触发）；对角态行为不变（含
  sysdynamics 旧热化对角初态测试）。

**实施后追加裁决（2026-11-20 用户）**：`fock_vacuumstate` 删除（计划原定保留，
用户推翻）——边界态构造内联进 `fock_sysdynamics` 步骤 ③（imag/mixed = 热迹 I、
real = 真空投影 P、`initial_state` 仅实时覆盖）；export 随删，闸门手工路径同式
内联复刻。附带发现：旧闸门 `_wf_build` 里 `initial_state isa FockMatrix` 分支
是从未执行的死代码（旧闸门从不传非 nothing 初态，三元惰性求值绕过名字解析）
——新热化初态用例首次踩到，测试侧改用 `GTEMPO.FockMatrix` 限定名。全套件复跑
绿。
