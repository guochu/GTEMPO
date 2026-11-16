# GTEMPO 测量层重构计划：measurements/ 拆分 + inserted-K 路线缓存

- 日期：2026-11-16
- 基线：GTEMPO HEAD `c01f446c` + 工作区（行号以当前布局为准，Phase A 拆分后引用以新文件为准）
- 设计先例：U1GTEMPO2 `SystemPreparation.local_maps`（cell_maps 键结构，prepared_system.jl:39-41）——注意其 `:reference` 路径（所有实时/混合默认路径）零缓存、传播子无独立缓存（2026-11-16 已核实），本计划的两层缓存均为 **GTEMPO 新建**，非移植
- 关联计划：`bridge_interface_simplification_plan.md` 独立并行；本计划 Phase A 先做（纯移动），其后两计划的行号引用互相按新布局换算

## 0. 已裁决（2026-11-16 会话）

- **M1** 新建 `src/measurements/` 文件夹，把两种并行测量方式分开写（D5 两轨各自成文件）。
- **M2** 缓存 = **传播子缓存 + 格张量缓存** 两层；**结果记忆化不立项**（不假设用户反复请求同一观测量——U1GTEMPO2 `Solution.measurements` 通道在 GTEMPO 不重建）。
- 路线 1 的边界语义矩阵表示维持不变（前次会话已论证：`OperatorInsertion` 不改扇区/变量表示，ρ₀ 折叠与任意矩阵表达力是变量路线不可替代的）。

## 1. 目标

1. `measurement.jl`（407 行，两轨现混排且 `nn` 交错在路线 2 段落中）拆分为 measurements/ 三个文件，两轨各自独立、公共面最小；
2. inserted-K 路线引入 `SystemFactorPreparation` 制备对象，携带传播子缓存 + 格张量缓存——扫描型测量（如 greater(i,0) i=1..N）的 automaton 成本从 O(N·L) 摊销到 O(L + N·切位邻域)，单点逼近"一次分子收缩 O(L)"的路线下限；
3. **零数值语义变化**：缓存命中返回与重算逐位一致的张量，golden/双栈对拍无需重锚。

## 2. 探查定谳（计划赖以成立的关键事实）

**F1【两轨现状与拆分清单】** measurement.jl 当前布局：
- 路线 1（inserted-K 边界语义）：头注释 1-26（含 2026-11-15 五处口径差定案①-⑤）、`_fock_ops`:33、`_number_op`:36、`_inserted_expectation`:48、`occupation`:63、`default_branch`:73、`raw_correlation`:84、`greater`:112、`lesser`:121、`_normalized_two_point`:132、`matsubara`:159、kernel 元数据助手 176-181、`nn`:277（**交错在路线 2 段落中**）、`raw_lead_correlator`:328、`electriccurrent`:350、`heatcurrent`:384、export 块 402-407；
- 路线 2（算符/胞坐标语义）：段注释 183-187、`lattice_monomial`:196、`environments_cache`:232、`cached_Z`:240、`cached_gf`:243、`cached_correlation`:299；
- `GrassmannContourIndex` 在 lattices.jl（路线 2 消费，不动）；`OperatorInsertion` 在 sysdynamics.jl:121（K 构造 API，测量层只消费，**留原地**）。

**F2【路线 1 无缓存且重算是大头】** `_inserted_expectation` → `sysdynamics` → `_system_factor_reference` 每次调用全量重建：逐格矩阵指数（Td 模型带时间参数）、`local_weights`、split_cell 逐带 SVD、pairize、reorder!；K 重建占单点成本 ~88%（measurement.jl:249 实测注释）。分母 Z 已由调用方算一次传入（共享 ✓）。

**F3【格张量的函数依赖——缓存键的完备性依据】** 格 k 的张量列只依赖：
- `mats` = S（传播子，由 branch/j/order/model 派生）· weights（order 绑定）· Ol·Or（该切位插入组）+ oddlater 符号（sysdynamics.jl:217-224，`(-1)^{oddlater·popcount(m)}`）；
- `rho`（仅 k==1/L==1 进 core：tr(M·ρ)/(M·ρ)；k==L 的 core=M 元素不含 ρ）；
- first/last 分支（core 装配三态，sysdynamics.jl:241-255）；
- `leftlabels/rightlabels` = 由 `prefix`（cut < k 的累计插入荷，sysdynamics.jl:203-236）与 k 派生；
- `charges_per_cell`（spec+bands 绑定）。
扫描 G(tᵢ,0) 时格 k 在 K′(i) 对所有 i ≥ k 逐位相同（oddlater=1、无切位算符）→ distinct 格数 ≈ 2L + N 而非 N·L。

**F4【OperatorInsertion 需结构化相等】** 用户每次测量调用新构造插入对象（occupation/measurement.jl:67 等），对象恒等（`===`）不成立——缓存键必须配结构化 `==`/`hash`（字段逐项 + 矩阵元素相等；dim ≤ 16，成本可忽略）。

**F5【sysdynamics 公共签名】** `sysdynamics(lattice, model; order, kwargs…)`（sysdynamics.jl:344）与 `sysdynamics_weyl2order(lattice, model; coupling_norm, kwargs…)`（:358）都转发 `_system_factor_reference(…; spec=TK.sectortype(lattice), kwargs…)`——spec 由 lattice 派生，制备对象绑定 (lattice, model, order) 即可，spec 不进键。

## 3. Phase A：measurements/ 拆分（纯移动，零语义）

```
src/measurements/
├── measurements.jl    # 头注释（D5 两轨总纲 + 两轨语义差异指针）+ include 两个子文件 + export 块（原 402-407 移入）
├── inserted_k.jl      # 路线 1：口径差①-⑤定案注释、_fock_ops、_number_op、_inserted_expectation、
│                      # occupation、default_branch、raw_correlation、greater、lesser、
│                      # _normalized_two_point、matsubara、nn（从路线 2 段落挪回）、kernel 元数据助手、
│                      # raw_lead_correlator、electriccurrent、heatcurrent
└── operators.jl       # 路线 2：lattice_monomial、environments_cache、cached_Z、cached_gf、cached_correlation
```

- GTEMPO.jl:63 `include("measurement.jl")` → `include("measurements/measurements.jl")`；导出面、函数签名、行为零变化；删除旧 measurement.jl。
- 文件头各自携带本轨语义声明：inserted_k.jl = "边界时刻语义，golden 锚定路径"；operators.jl = "胞坐标 0⁺ 语义，差 O(胞内接触项)（原 cached_gf D5 注原样搬运）"。
- 单独提交；验收 = 全测试绿（纯移动无新语义）。

## 4. Phase B：SystemFactorPreparation + 两层缓存

### 4.1 制备对象（sysdynamics.jl 侧，K 构造层自有）

```julia
struct SystemFactorPreparation
    lattice::AbstractGrassmannLattice
    model::AbstractImpurityHamiltonian
    order::Int
    coupling_norm::Float64
    initial_state            # nothing / FockMatrix / 矩阵
    # 一次性预计算
    cells; weights; fqs; charges_per_cell; dim
    # 两层缓存（私有 Dict；无淘汰——distinct 格数被扫描形态天然有界）
    propagators::Dict{Tuple{Symbol,Int},Matrix{ComplexF64}}   # (branch, j) → S
    cell_tensors::Dict{Tuple,Vector}                           # §4.2 键 → split_cell 张量列
end
prepare_system(lattice, model; order=1, initial_state=nothing, coupling_norm=0.0) -> SystemFactorPreparation
clear_cache!(prep)   # 可选：清空两层 Dict（预计算字段保留）
```

- 命名对齐 U1GTEMPO2 `prepare_system`/`SystemPreparation`（port 纪律）；**不抄** 其 `thermal_model` 字段（GTEMPO 的 τ 支热迹走 ρ₀=I，模型单一）。
- 所有权约定（scaling.jl 老原则）：缓存严格归属本制备对象，**永不做进程级全局**；对象不可变字段一经构造不变，无失效问题——输入变了就新建对象。

### 4.2 缓存键（完备性 = F3 逐项对映）

- **传播子**：`(branch, j)`——dt/t 由 (lattice, branch, j) 派生、模型与 order 绑定制备对象，不进键。
- **格张量**：`(branch, j, Ol组, Or组, oddlater::Int, prefix荷, first::Bool, last::Bool, cut0组)`：
  - Ol组/Or组/cut0组 = 该切位 `:left`/`:right` 过滤后的 `OperatorInsertion` 元组（F4 的结构化 `==`/`hash` 先决；cut0 组经 rho = cutleft(0)·ρ·cutright(0) 只影响 first/last 格，但统一进键最稳）；
  - prefix荷 = 切位 < k 的累计插入荷（sector 值）——leftlabels/rightlabels 的唯一自由度；
  - oddlater = 该格之后奇插入计数（Koszul 符号指数）。

### 4.3 等价性论证与接线

- 键完备 ⟺ 命中返回张量与重算逐位一致（F3 的函数依赖全覆盖）——数值零漂移，golden 无需重锚。
- `OperatorInsertion` 补 `Base.:(==)`/`Base.hash`（字段逐项 + `matrix` 元素相等）。
- 新增方法：`sysdynamics(prep::SystemFactorPreparation; insertions=OperatorInsertion[], trunc=nothing, truncerrs=nothing)` 与 `sysdynamics_weyl2order(prep; coupling_norm, …)`（order=2 同一制备机器）——内部走 `_system_factor_reference` 的制备对象版本（逐格改为 `get!(prep.cell_tensors, key) do … end`，传播子经 `prep.propagators`）。
- 原签名 `sysdynamics(lattice, model; …)` 保持为便利包装（内部一次性 prepare_system 即抛）——**现有调用方与 golden 路径行为逐位不变**。
- 测量层（inserted_k.jl）：`_inserted_expectation` 及全部公开函数（occupation/greater/lesser/matsubara/nn/raw_lead_correlator/electriccurrent/heatcurrent）增 `preparation=nothing` kwarg 透传；`nothing` = 现状路径。扫描用户流：

```julia
prep = prepare_system(lattice, model; order=2, coupling_norm=cn)
G = [greater(lattice, model, kernel, Is...; i=i, Z=Z, preparation=prep) for i in 1:Nt]
```

### 4.4 复杂度核算

- 现状：每点 = L×(exp + 权重积 + split_cell SVD) + pairize + reorder! + O(L) 收缩。
- 制备对象下 N 点扫描：传播子 L 个各一次；distinct 格 ≈ 2L + N；每点增量 = 收缩 O(L) + （新切位格变体 × split_cell）。
- 与路线 2 的差距缩至"全链收缩 vs 稀疏楔积"；路线 2 不动（其缓存已具备：environments_cache + 快照扫掠双模式）。

## 5. 验收闸门

- **Phase A**：全测试绿（纯移动）。
- **Phase B**：
  1. 全测试绿；
  2. 双栈对拍维持：greater/lesser 3 点、matsubara 4 点、occupation 4 点、lead 关联 6 点（含同切位）Δ ≤ 5e-15（measurement.jl 头注释同口径）；
  3. **逐位性测试**（新增）：同一 prep 下，带缓存 K′ 与无缓存重建的 todense 逐位相等（扫描首点/中点/末点）；
  4. **键欠定义回归**（新增）：同格换不同插入（如不同 band 的 d、n̂、空切位）产出不同且各自正确的 K′——防缓存键漏字段；
  5. **扫描基准**：performance/ 下新增（M=2 实时 Nt=32，greater(i,0) 全点），报告 before/after 计时与内存。

## 6. 不做清单

- 结果记忆化（M2 裁决；U1GTEMPO2 `Solution.measurements` 通道不重建）；
- `OperatorInsertion` 改扇区/变量表示（维持边界语义矩阵表示）；
- 路线 2 任何改动；
- U1GTEMPO2 侧回移（其 `:reference` 路径仍无缓存；未来若要，参照本计划 §4）。

## 7. 净效果

- measurement.jl 407 行单文件 → measurements/ 三文件（两轨各自独立、公共面缩到 include+export）；
- inserted-K 路线获得 U1GTEMPO2 直接路径同款的格张量缓存 + 它也没有的传播子缓存；扫描型测量单点成本逼近 O(L) 收缩下限；
- 全程零语义变化，golden 不重锚。
