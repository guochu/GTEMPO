# fock_sysdynamics 计划（v2）：Fock 表象直积 MPO + fockmpo2gmps/gmps2fockmpo 互逆转换

> **✅ 已完结（2026-11-17，四阶段全落地，见 §7 实施状态）。**

- 日期：2026-11-17（v2：按用户锚点修订——TEMPO NonDiagonalHyb 先例、格点数确认、转换数学原理重推）
- 基线：GTEMPO HEAD `f3ee05c4` 系工作区（measurement cache / 接口简化已实施）
- 关联：`if_expmult_backend_plan.md`（IF 装配，独立）；`measurement_cache_refactor_plan.md`（缓存**本计划推迟**，V2）
- 先例：TEMPO `NonDiagonalHyb` 直接输出 Fock space MPO（玻色、无 bands 概念——所有 bands 合并进一个 site）；本计划 = 同一架构的费米 + 多带推广

## 0. 已裁决

- **V1**（2026-11-17 用户提案）：K 链装配重构三段——`fock_sysdynamics` 产出 **Fock 表象直积 MPO** → 观测量插入 = **插入算符 MPO 作用到 K-Fock**（SparseFockMPO 相乘）→ **一对互逆转换** `fockmpo2gmps` / `gmps2fockmpo`。prefix/pair-labels/切位 foldr·foldl/split_cell-for-K 退役。
- **V2**（用户裁决）：**缓存与效率整体推迟**——新路径无缓存裸跑，金标准对齐后另立缓存/效率计划；本计划期间 `preparation` kwarg 不接新路径。
- **V3**：golden = **机器精度重锚**（数学值一致、张量 gauge 不同）；现行路径保留至 Phase 4 末作交叉验证。

## 1. 目标形状（2026-11-17 实测定谳，非推算）

**实测记录**（ julia --project 各自包，ToulouseIM(ϵ_d=0.5)、δτ=0.25、N=4 / δt=0.1、N=3）：

| 对象 | 站数（bands=1） | 迹 tr(ρ₀·Π O) |
|---|---|---|
| 旧栈变量链（backup，开链+boundarycondition! 闭合） | 2(Nτ+1)=10 / 4(Nt+1)=16（变量粒度：每切片 2/4 个变量站；k=N+1 / Nt+1，lattices/imaginarytime.jl `k=N+1`） | 1.6065306597126**423** |
| 旧链 pair 化站数（用户口径） | **Nτ+1 = 5 / 2(Nt+1) = 8** | 同上（同一链） |
| 新包现行 GMPS | Nτ=4 / 2Nt=6 | 1.6065306597126**323** |

**Δ ≈ 1e-14 = 机器精度等效**：新链把旧链的 β 闭合边界片（变量腿内容恒等的平凡站）折叠进末格矩阵——站数差 1/分支，算符恒等。用户口径 Nτ+1 / 2(Nt+1) = 旧生产基准（golden 参照）的 pair 化站数。

| 层 | 格点/site 数（bands=1） | 说明 |
|---|---|---|
| FockMPO site | **Nτ / 2·Nt**（每 cell 一 site，bands 合并进 site 内部：Fock 维 2^bands、monomial 腿 4^bands） | 与现行新包 cell 数一致；Nτ 个 cell 携带 Nτ+1 个插入切位（0..kτ，含 ρ₀ 边界与 β 闭合边界），插入独立 site 插在切位上 |
| 转换产物 GMPS | **Nτ·bands / 2Nt·bands**（pair 布局，与现行新包同形） | 值 = 旧基准（机器精度，V3 重锚）；形状 = 新形（非旧 Nτ+1 形——边界片已折叠） |

先例：TEMPO `NonDiagonalHyb` 直接输出 Fock space MPO（玻色、无 bands——所有 bands 合并进一个 site）；本计划 = 同一架构的费米 + 多带推广。

## 2. 数学原理：`fockmpo2gmps` 的变换本体（本计划核心，必须按此实现）

### 2.1 对象

Fock MPO site 张量：`T[a_out, m, a_in] = ⟨a_out| W_m |a_in⟩`，其中
`W_m = (Grassmann 单项式 m 的系数算符) = graded 符号 × Π_{b=1..B} O^{(b)}_{slot_b(m)}`，
`O^{(b)}_{slot} ∈ {qI−n_b, −d_b†, d_b, W₄}`（local_weights 词典，物理不动）。
K = Σ_{m₁..m_L} G_{m₁}···G_{m_L} · tr(ρ₀ W_{m₁}···W_{m_L})。

### 2.2 转换 = 相干态单位分解的两次应用（无 SVD）

**第一层（格内分解）——Fock 单位分解插入**：monomial 腿按带因子化 m = (slot₁..slot_B)，乘积算符的矩阵元经 **Fock 恒等分解**精确拆成链：

```
⟨a_out|W_m|a_in⟩ = Σ_{a'_1..a'_{B−1}} ⟨a_out|Õ^{(B)}|a'_{B−1}⟩···⟨a'_1|Õ^{(1)}|a_in⟩
```

中间键 = **Fock 指标**（维 2^bands）——这就是格内自动机空间本身，**不需要 SVD 去发现**（现行 split_cell 的 SVD 是在 pair-荷差编码下寻找同一空间，编码换掉后分解变为显式）。** graded 符号分配（bands≥2）**：标量 sign_m 依赖全部 slot，不能挂单因子；按增量 Koszul 分配——放置第 b 带时获得因子 Π_{b'<b}(−1)^{popcount(slot_{b'})·popcount(slot_b)}，等价于中间 Fock 键携带的费米宇并相位。M=1 时无跨带对、符号恒 1（现行 `graded_sign_factor(1)=1` 一致）。**bands≥2 的分配正确性由 Phase 1 对拍裁决（显式覆盖 bands=2）。**

**第二层（键转换）——Fock 指标 → 荷扇区**：`Fock 指标 a ↔ (扇区 q = fqs[a], 简并内局部序 = 枚举位置)`；键空间 `Vect[G](q => 计数)`。信息无损（fqs 表 + 枚举序是互逆映射）；数据按守恒律 `q_out = q_in + Σ 变量荷` 重排进块。**纯重排，零算术。**

**物理腿（关键同构）**：Fock 单位分解直接产出**每带一个站张量**（不再有变量粒度中间层）：子 site 物理腿 = slot 腿，4 态 {none, ā, a, āā} 荷 {0, −q, +q, 0} → `Vect[G](0=>2, ±q=>1)`——**与 GM `grassmannspace(G)`（vac/ηη̄/η/η̄，扇区 0=>2、±q=>1）空间同构、逐态一一对应**（none↔vac、ā↔η̄、a↔η、āā↔ηη̄）。故现行管线的「split_cell 变量列 + pairize」两步在转换中**溶解为腿再标记**：每 cell 直接产出 bands 个站（bands=1 → 每 cell 1 站，与实测新包 4 站/N=4 一致）。

### 2.3 ρ₀/迹闭合与插入的时序

- **ρ₀ 折叠 = `fock_initstate!`（Fock 层、原地）**：ρ₀ 收缩进首 site 入腿（现行首格 `(M·ρ)[s,t]` 同语义）；**迹闭合 = `fockmpo2gmps` 的职责**（末格对角求和）——转换的闭形职责，Fock 链保持开链（可复用于其它收缩口径）。产出与现行结构同形的迹闭合 GMPS（真空端标签）。
- **插入 = SparseFockMPO 算符链乘法（左/右作用）**：插入算符是 (dim,dim) 两腿算符（荷跳跃 q_W 无法被 2 维 pspace 变量腿承载，必须作为独立算符乘在链上、不与 measure 权重混叠）——左作用 = `:left`（cell 前）、右作用 = `:right`（cell 后，§3.0 映射表）；同切位 0⁺ 序 = 同支撑位的列乘序。
- **oddlater Koszul 符号不盲目搬运**：现行 `(-1)^(oddlater·popcount(m))` 是 pair-label 路由约定下的补偿；新架构从「全局变量链序」第一性重推——预期被变量槽位规范吸收或以等价公式挂在插位后的 site。**Phase 1 全配置对拍为唯一裁决，不预设结论。**

## 3. 设计

### 3.0 `FockMPO` 类型与 convenience 层（新文件 `src/fockmpo.jl`）

**设计约束（用户锚点）：FockMPO 的对称性空间与 GrassmannMPS 完全一样**——直接复用 TensorKit/GrassmannMPS 的 `Vect[G]` 分级空间机制（扇区分块、fusion tree、ComplexF64 块存储），不发明 dense 中间层。Fock 侧与 GM 侧的差异只在**腿语义**：

| 腿 | GM（GrassmannMPS） | FockMPO | 对应 |
|---|---|---|---|
| 键 | Vect[G] 键差荷扇区 | Vect[G] **Fock 占据荷**扇区（`fqs[a]` = 占据带 var_charge 和；简并 = 同荷计数，局部序 = Fock 枚举位置） | 空间机制同一（Vect[G]），扇区内容不同；转换 = 键重封块 |
| 物理 | grassmannspace 4 维（vac/ηη̄/η/η̄，扇区 0=>2、±q=>1） | **slot 腿** Vect[G](0=>2, ±q=>1)（none/āā/ā/a 四态） | **空间同构、逐态一一对应**（none↔vac、ā↔η̄、a↔η、āā↔ηη̄）——转换 = 同形腿再标记 |

```julia
# GM 侧前置（grassmannspaces.jl，与 MPSTensor:52 并排，一行别名 + export）：
const MPOTensor{T,S} = AbstractGrassmannTensor{T,S,2,2}
    # 位型：(bond_L ⊗ slot_out) ← (slot_in ⊗ bond_R)；FockMPO site 的 slot 结构块对角
    # （slot 指标过站不变，每 slot 一个 (dim,dim) 算符块），2,2 参数化精确匹配

# GTEMPO 侧（fockmpo.jl）——纯容器（用户裁决：无 sites/logscale field）：
struct FockMPO{T<:MPOTensor}
    data::Vector{T}                 # 子 site 张量——与 GMPS{A<:MPSTensor} 同构的类型参数化
end
# 站位解释由 lattice 携带：FockMPO 的 site k ↔ chronological_cells(lattice)[k]
# （所有 FockMPO 函数都配套 lattice 参数）；logscale 不挂 FockMPO——Weyl 标量在
# fockmpo2gmps 产出 GMPS 后 setlogscale!（现行收尾语义原样，sysdynamics.jl:370-372）

# 稀疏算符链（用户锚点：算符生成的 MPO 用 Sparse 定义，限制作用范围、为缓存优化留口
# ——镜像 GM SparseGMPS gmps.jl:218）：
struct SparseFockMPO{T<:MPOTensor}
    positions::Vector{Int}          # 支撑 site（链位 = 插入切位）
    data::Vector{T}                 # 支撑位算符（(dim,dim) 两腿 MPOTensor，slot 腿 1 维平凡）
end
```

**Convenience 函数（闭环必需）**：

```julia
fock_vacuumstate(lattice; initial_state=nothing) -> FockMPO
    # Fock 侧初态/ρ₀ 对象（现行 ρ₀ 逻辑 sysdynamics.jl:174-179 的镜像：imag = I/给定态、
    # real = 真空投影，扇区分块；GM vacuumstate 的 Fock 镜像）
    # fock_initstate! 的缺省 ρ₀ 来源

prodfockmpo(spaces, sites, ops) -> FockMPO
    # GM prodmps 同款构造器：spaces = 键空间序列、sites = (branch,j) 元数据表、ops = 每 site 算符
    # fock_sysdynamics 的装配出口

sparsefockmpo(lattice, positions, ops) -> SparseFockMPO
    # 观测量算符 → SparseFockMPO：positions = 插入切位（lattice 映射到链位）、ops = 每位算符
    # （同位多算符 = 列表序，列乘序 = 0⁺ 序）；greater/lesser 双位置 = 一个 SparseFockMPO 带两位

Base.:*(A::SparseFockMPO, K::FockMPO) -> FockMPO    # 左作用（现行 :left：A_k·O_k = 插在 cell k 前）
Base.:*(K::FockMPO, A::SparseFockMPO) -> FockMPO    # 右作用（现行 :right：O_k·A_k = 插在 cell k 后）
    # 稀疏应用：支撑位 = 键收缩的算符乘积、其余 site 恒等透传——作用范围受限、缓存优化锚点（V2）
    # 纯玻色式张量收缩（无 graded 符号：费米符号全部在 d/d† 矩阵元内；
    # Koszul 只在 fockmpo2gmps 转换边界出现一次）

fock_initstate!(lattice, fockmpo, ρ₀=nothing) -> FockMPO
    # ρ₀ 折叠（原地）：ρ₀ 收缩进首 site 入腿（=(M·ρ)[s,t] 同语义）；ρ₀ 缺省 =
    # fock_vacuumstate(lattice) 的轮廓缺省（imag/mixed 热迹 I、real 真空投影）
    # 迹闭合不在本函数——留在 fockmpo2gmps（末格对角求和，转换的闭形职责）

fockmpo2gmps(f::FockMPO, lattice) -> GMPS
gmps2fockmpo(g::GMPS, lattice) -> FockMPO
    # 互逆转换（同键扇区约定下精确互逆）；gmps2fockmpo 仅适用于 fockmpo2gmps 产出的
    # 新形 GMPS（键 = Fock 占据荷扇区），对现行旧形（pair-荷差键）不适用


闭环工作流（用户语义：插入算符自己也走 FockMPO、相乘即可；两侧插入 = 左右各乘一个）：

```julia
Kf = fock_sysdynamics(lattice, model; order, coupling_norm)    # prodfockmpo 装配（开链，无 ρ₀/闭合）
Kf = fock_initstate!(lattice, Kf, ρ₀)                          # ρ₀ 折叠（缺省 = fock_vacuumstate 缺省）
Wl = sparsefockmpo(lattice, cuts_left, ops_left)               # 左侧插入算符链
Wr = sparsefockmpo(lattice, cuts_right, ops_right)             # 右侧插入算符链（可省）
K  = fockmpo2gmps(Wr * (Wl * Kf), lattice)                     # 迹闭合 + 腿再标记 + 键重封块 + 增量 Koszul + setlogscale!
Z  = integrate(K, Is...)
```

**插入位置语义映射（现行 (cut, side) → MPO 位置，实现在注释成文）**：`:left@cut k`（现行 `Ol·M_k`）= 算符插在 cell k **之前** = 左作用 sparse 支撑位 k；`:right@cut k`（`M_k·Or`）= 插在 cell k **之后** = 右作用支撑位 k。cut 0 = 链最左（ρ₀ 边界前）、cut kτ = 链最右（β 闭合前）。全部现行 (cut, side, 0⁺ 序) 组合在此映射下无损表达；Phase 1 对拍矩阵裁决。

### 3.1 `fock_sysdynamics`（装配，纯 Fock 无 Grassmann 内容；开链）

- 逐格算符 `O_{k,m} = S·M_m·S`（order 2）/ `S·M_m`、`M_m·S`（order 1）→ **Fock 单位分解插入**拆成 bands 个子 site（每带一个，slot 腿 = grassmannspace 形）→ `prodfockmpo` 组链（纯系统因子链，不含插入——插入由 `sparsefockmpo` 算符链 `*` 作用）。
- 开链输出（无 ρ₀/迹闭合——ρ₀ 归 `fock_initstate!`、闭合归 `fockmpo2gmps`，§2.3）。
- 输出 `FockMPO`（纯容器）。

### 3.2 `fockmpo2gmps` / `gmps2fockmpo`（互逆转换对）

**fockmpo2gmps**：末格对角求和闭迹 → per-子-site **slot 腿再标记**（none/ā/a/āā → vac/ηη̄/η/η̄ 四态一一对应，同形空间）+ **键重封块**（Fock 占据荷 → 现 GMPS 键差荷扇区，从 fqs 表重算扇区计数与局部序，纯重排）+ **增量 Koszul**（bands≥2，§2.2）→ 直出每 cell bands 站（**无变量列、无 pairize 步**）→ reorder!（Chrono 存储序）+ setlogscale!（Weyl 标量）→ GMPS。分块填充工具（abelian 四形通用，键腿/物理腿复用，~40 行）。插入 site（1 维 slot 腿）转换时吸收进邻站或保留 1 维物理腿站——实现期定。

**gmps2fockmpo**：fockmpo2gmps 的逆映射（同键扇区约定：扇区 → Fock 指标 = 简并内局部序反查枚举、slot 再标记反向）；用途 = 对新形 GMPS 做 Fock 层再操作（追加插入、跨口径调试对拍）。对现行旧形（pair-荷差键）GMPS 不适用。

### 3.3 与现行路径的关系

- `sysdynamics`/`sysdynamics_weyl2order` 签名不变；`preparation` kwarg 不接新路径（V2，给 prep 的调用暂走旧路径，两径并存至 Phase 4）。
- `split_cell`/`from_core` 本体保留（partial_if:272 / imag_factors:428 仍消费）；`local_weights` 词典原样；**现行管线的 pairize 在 K 路径不再消费**（转换直出站形）。
- **评估项（Phase 3 后可选）**：Z 类用途下单站/格形式（4^bands 物理腿按 monomial 总荷分扇区，免带分解）可能足够——若对拍通过可进一步减站；先不做，保持 per-band 站 drop-in。

## 4. 阶段与验收

- **Phase 1（类型层 + 转换函数 + 全配置对拍闸门）**：GM `MPOTensor` 别名 + GTEMPO `FockMPO{T<:MPOTensor}`（prodfockmpo/fock_vacuumstate/fock_initstate! 最小实现）+ `fockmpo2gmps`/`gmps2fockmpo` + 分块填充工具；对拍 = 新链 `todense` vs 现行链 `todense`，机器精度。配置矩阵：{Imag, Real, Mixed} × bands ∈ {1, 2}（bands=2 显式覆盖 graded 符号分配与扇区简并块）× spec ∈ {:z2, :u1, :u1u1}（:z2z2 抽查）× 插入 ∈ {无、偶（n̂，cut 0/内部/末端）、奇（d/ddag、同切位 0⁺ 双插）、多插入（现行 greater/lesser 测试同款）} × order ∈ {1,2}。**不过全停。**
- **Phase 2（装配 + 乘法 + 端到端）**：`fock_sysdynamics` + `fock_initstate!` + `SparseFockMPO` 乘法（`sparsefockmpo` 构造 + 左/右作用）+ 新旧并存开关。验收：矩阵端到端机器精度；测试套件切新开关下全绿（goldens 按值重锚——V3，>1e-10 记录理由）；ED 对拍维持。
- **Phase 3（入口切换 + 旧径封存）**：默认切新；旧路径 deprecated 保留（交叉验证基准）。
- **Phase 4（清理）**：prefix/pair-labels/foldr·foldl 切位/split_cell-for-K 与 pairize-for-K 调用点退役；**缓存/效率另立计划**（V2：site 键不含切位算符、转换按 (lattice, spec) 缓存、性能基准）。

## 5. 不做清单

缓存与效率（V2）；partial_if/ExactTTI 装配（归 IF 后端计划）；local_weights/split_cell/from_core 本体删除；逐位 golden（V3）；单站/格形式切换（§3.3 评估项，不在本计划）。

## 6. 工作量与顺序

Phase 1 ≈ 2 天（类型层 + 转换数学 + 全配置矩阵，本计划科学核心；含 GM `MPOTensor` 别名一行前置）；Phase 2 ≈ 1 天（`FockMPO * FockMPO` 收缩 ~50 行 + 端到端）；Phase 3 ≈ 半天；Phase 4 ≈ 半天。总 ≈ 4 天。顺序严格 1→2→3→4；Phase 1 不过全停。

## 7. 实施状态（2026-11-17 完结）

- **Phase 1 ✅（前会话，GTEMPO 68e90cb0 + GM 789a1c6）**：类型层 + `fockmpo2gmps` + 对拍闸门（163 项）。关键 pivot（原型实证，`debug/fockmpo_proto/`）：无 SVD 直构被 GM 楔积规范相位 σ(occ) 定谳否决 → 分站装置复用基线 split_cell+pairize（符号规范逐位继承）；FockMPO 键腿 = 非分格 Vect[G]（带荷对角 MPO 腿不可分格表达）。§2 的「无 SVD 两层变换」数学方案由此被取代。
- **Phase 2 残余 ✅（本会话）**：`SYSFACTOR_BACKEND` 开关（`:fock`/`:reference`，`set_sysfactor_backend!`）+ `_system_factor_fock`（lattice/model 与 prep 双版）+ prep `fockmpo` 单槽缓存（FockMPO 与插入无关——扫描 automaton 成本 O(L) 一次 + 每测量仅切位邻域稀疏乘法）。**期间定谳并修复 V1 切位 bug**：`sparsefockmpo`/`fockmpo2gmps` 三明治原按 `boundary==0` 过滤，而 `:-` 支 boundary 0 = **链末切位**（cut = ni+2Nt）非 ρ₀——lesser(j=0) 的 ddag 被错折进链头 → 链归零（real 轮廓被真空 ρ₀ 投影湮灭而侥幸不可见，mixed ρ₀=I 露馅）；改按 `insertion_cut_index` 过滤（src 与闸门同口径），闸门补 U1×bands=2 与 (:-,0) 链末用例（163→203 项）。**全套测试切 `:fock` 绿；双后端 todense 逐位 Δ=0**（抽样 bands∈{1,2} × order∈{1,2} × 三轮廓 × 链末切位）——**golden 零漂移，V3 重锚不需要**。
- **Phase 3 ✅（本会话）**：缺省 `backend=:fock`；`_system_factor_reference` deprecated 保留 = 交叉验证基准（对拍闸门锚定侧；`backend=:reference` 显式可达）。
- **Phase 4 ✅（本会话，按 V1 pivot 重释）**：调用点审计——生产入口（sysdynamics/sysdynamics_weyl2order ×4 + 测量层 `_inserted_expectation`/`electriccurrent_fast`）全走新径。plan 原文的退役清单被 pivot 取代：split_cell/pairize-for-K 即 `fockmpo2gmps` 分站装置本体（保留）；foldr·foldl 切位折叠与 prefix 在新径同样需要（保留）；pair-labels 机制不存在于 src。旧径整体删除推迟至 V2 裁决（闸门锚定侧）。
- **V2（另立计划，未立项）**：缓存/效率（FockMPO 按 (lattice,spec) 缓存、§2 无 SVD 直构路线以 σ 谱为起点）、`gmps2fockmpo`（楔积基底映射 + 合并秩）、单站/格形式评估（§3.3）。
