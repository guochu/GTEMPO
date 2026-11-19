# GTEMPO IF 离散化后端切换：QuAPI → IFDiscretization（scheme 入口）计划

> **状态：立项待实施**（2026-11-19 成文）。用户提案：把 GTEMPO 的 IF 离散化
> 后端从 QuAPI 切换到 IFDiscretization，经
> `correlationfunction(bath, lattice, scheme)` 参数化实现，并为严格 IF 算法
> 留入口。**可行性定谳：可行**（§0）——IFD 侧 scheme 层/统一入口/逐位对齐
> 三要素均已就绪，GTEMPO 侧触碰面收敛于 correlationfunction。

- 日期：2026-11-19
- 基线：GTEMPO 分支 `if-expmult-backend` HEAD `0f64cf78` + 工作区并行 WIP；
  IFDiscretization main HEAD `76e2d79`（已推送，2026-11-19 三金标准
  54+58+34 绿 + IFD↔QuAPI 位级对拍契约）
- GM / 其余包：零改动
- 关联计划：`exactif_port_plan.md`（严格 IF 消费侧——本计划 Part 2 只开口子，
  块级消费契约归该计划，两计划接口对齐点见 R5）；`if_expmult_backend_plan.md`
  W2（kernel 存储约定——本计划不动）；`polarif_plan.md` / `bcs_exacttti_plan.md`
  （文件相邻无冲突）

## 0. 可行性定谳（用户两问）

**Q1：GTEMPO 的 QuAPI 入口是否都收敛在 correlationfunction？——是。**

- src 调用点 **6 处全在 `correlationfunction.jl`**：normal
  `QuAPI.Δτ/Δt/Δm`（:175-179）+ BCS `QuAPI.bcs_Δτ/Δt/Δm`（:264-268），
  全部限定名、无非限定调用（全 src 反查零命中）；
- 结构性触碰 **3 处**：`GTEMPO.jl:22` `@reexport using ImpurityModelBase, QuAPI`
  （容器类型进 GTEMPO 命名空间）、`GTEMPO.jl:23` `import QuAPI: branch`
  （correlationfunction.jl:186/:282 两处 branch 窗口化视图）、
  `correlationfunction.jl:105-108` 第八批本地定义的容器 `scalartype` 四方法；
- 测试性触碰 **2 文件**：`test/api/quapi_crosscheck.jl`（头注释自证
  "T7 退役后唯一保留的 QuAPI 面"）+ `test/api/correlation.jl`（容器引用，
  经 reexport 透明）。

**Q2：经 `correlationfunction(bath, lattice, scheme)` 实现切换是否可行？——可行，IFD 已备好全部前提。**

- **IFD scheme 层现成**（`schemes.jl`）：`DiscretizationScheme` 抽象 ←
  `FrozenKernel` ← `QuAPI()`（一阶 Trotter 冻结核）+ `ExactKernel` ←
  `ExactEqualDistantKernel(solver)`（EDToulouseSolver/AnaToulouseSolver）；
- **IFD 统一入口已按 scheme 分派实装**（fermions/bosons/bcs 三族）：
  `Δτ(bath; N, δτ, scheme=QuAPI())` → `_fermions_Δτ(scheme, bath, N, δτ)`，
  `::QuAPI` 与 `::ExactEqualDistantKernel` 双方法齐全，BCS 同构
  （`_bcs_Δ*_dispatch(scheme, ...)`）——**正是本提案的入口形状**；
- **切换按构造值中性**：IFD↔QuAPI 位级对齐已建立（2026-11-18，对拍 21 逐位
  + 5725 逐元素全绿；IFD 同步过 QuAPI 的 BCS getindex 语义/branch/位移访问器），
  且 IFD 容器同构（`CorrelationMatrix{T}` ηⱼₖ/ηₖⱼ 同字段同约定、
  `GrassmannContourIndex` 访问器含 `_vs_rt_acc` 顶点移位）⟹ GTEMPO 消费端
  （branch 视图、`_window_block`、`_bcs_dense_block` 访问器读取）零适配；
- **scheme 类型 = IFD 的**，GTEMPO 零新类型、零新分派层——GTEMPO 只是把
  "写死的 QuAPI 离散化" 换成 "透传 scheme"。

## 1. 事实清单（设计与风险推导的依据）

- **F1【IFD 无 QuAPI 依赖】** IFD deps = ImpurityModelBase + LinearAlgebra
  （uuid `29079917-c15e-41a7-ab41-1df48b32450b`）——单向依赖，切换后 GTEMPO
  的 QuAPI 依赖可降级为测试专用（对拍锚保留）。
- **F2【exact API 返回形状】** `exact_fermionic_Δτ` →
  `(ImagCorrelationFunction, det(BB))`——**标准容器**元组，GTEMPO 现行消费端
  （branch/窗口/fullfactor 直读）原样可用；`exact_fermionic_Δt/Δm` →
  `(StrictReal/MixedCorrelationFunction, normdelta)`——**严格容器**
  （装配级语义：有向指数系数、骨架+踢契约），branch/窗口化不适用 ⟹
  real/mixed 严格路线 ≠ 本计划范围（exactif_port_plan 承接）。
- **F3【GTEMPO scalartype 面】** `correlationfunction.jl:105-108` 四方法
  定义在 QuAPI 容器类型上（Imag{T} 真参数 / Real、Mixed 结构性 ComplexF64 /
  BCS 委托味块）——QuAPI 退主依赖后这组方法不可加载，须在 IFD 类型上同码
  重定义（IFD 容器参数化同构，方法体逐字平移）。
- **F4【BCS + strict 的就绪度存疑】** IFD 导出 `exact_bcs_Δτ/Δt/Δm` 且
  `bcs/exact/bcs.jl` 有实现，但 IFD README TODO 标 "bosonic/BCS strict
  kernels pending"——实施时核验 `_bcs_Δ*_dispatch(::ExactEqualDistantKernel)`
  接线；未接线则 GTEMPO 侧 BCS+exact 明确 ArgumentError。
- **F5【二阶非阻塞】** GTEMPO 现行 6 调用点全为一阶冻结核（二阶语义在
  lattice/Weyl 层）；IFD 的二阶装配（`QuAPI2 <: FrozenKernel`）是注释扩展位
  未实装——未来二阶 IF 核经同一 scheme 面承接，非本计划范围。
- **F6【crosscheck 自比对陷阱】** `quapi_crosscheck.jl` 现用非限定
  `Δt/Δτ/Δm`（经 reexport）——切换后这些名字绑定 IFD，测试退化为自比对
  恒真；须显式 `QuAPI.` 限定才保住"IFD 核 vs QuAPI 原块逐位 =="的安全网语义。
- **F7【exact 的 μ/εd 面】** IFD exact 签名带 `μ`（纯浴核零杂质信息，
  `Δ(x)=x·Δ(1)`）；带踢/算符忠实装配（ε_d kicks）走
  `exact_fermionic_Δt_cells` 的 cells/kicks 接口——超出 correlationfunction
  的 bath-only 签名面，属 exactif_port 的装配契约（R5）。

## 2. Part 1 设计：后端切换 + scheme 入口（机械迁移，闸门严）

### 2.1 签名（`correlationfunction.jl`）

```julia
correlationfunction(bath::AbstractNormalBath, lattice, scheme=IFD.QuAPI(); band=nothing)
correlationfunction(bath::AbstractBCSBath, lattice, scheme=IFD.QuAPI(); orbital=nothing)
```

位置第三参（house style 同 `hybriddynamics(lattice, kernel, alg)` 第三位
算法对象）；现行两参调用透明（缺省 = IFD `QuAPI()` = 现行值）。scheme 类型
直接用 IFD 的 `DiscretizationScheme` 族——GTEMPO 不定义 scheme 类型。

### 2.2 六调用点替换（同文件内）

- normal 三处 → `IFD.Δτ(bath; N=lattice.Nτ, δτ=lattice.δτ, scheme)` /
  `IFD.Δt(bath; N=lattice.Nt, t=lattice.Nt*lattice.δt, scheme)` /
  `IFD.Δm(bath; Nτ, t, Nt, δτ, scheme)`——bath 对象统一分派
  （`bath.f`/`β`/`μ` 拆包代码删除，IFD 从 bath 直读）；
- BCS 三处 → `IFD.Δτ(bath; N, δτ, scheme)` 等（同上；现
  `QuAPI.bcs_Δτ(bath.f; β=…, μ=…, Δ=…, …)` 的手工拆包退役）。

下游容器消费（branch 视图、`_window_block`、`_bcs_dense_block` 的
GrassmannContourIndex 访问器读取、blocks 装配循环）**零改动**——IFD 容器
同构（Q2）。

### 2.3 命名空间迁移（`GTEMPO.jl` + `correlationfunction.jl`）

- `GTEMPO.jl:22`：`@reexport using ImpurityModelBase, QuAPI` →
  `@reexport using ImpurityModelBase, IFDiscretization`（IFD 导出面与 QuAPI
  同名同形：CorrelationMatrix、branch/branches、Δτ/Δt/Δm、fermionic_*/bcs_*、
  四容器、GrassmannContourIndex——用户侧命名空间透明；实施时跑一次
  names 冲突清单核对，见 R2）；
- `GTEMPO.jl:23`：`import QuAPI: branch` 删（IFD 直接导出 branch）；
- `correlationfunction.jl:105-108`：scalartype 四方法改定义在 IFD 容器
  类型上（方法体逐字平移，F3）；
- `GTEMPO.jl:58` 注释「QuAPI 薄壳」→「IFD scheme 薄壳」类表述。

### 2.4 Project.toml

- `[deps]` 增 `IFDiscretization = "29079917-c15e-41a7-ab41-1df48b32450b"`；
- QuAPI uuid（`49288508-…`）从 `[deps]` 移到 `[extras]` + `[targets].test`
  （quapi_crosscheck 对拍锚保留，F1）——注意测试须以
  `julia --project=. test/runtests.jl` 直跑（运维锚）验证 test target
  依赖解析（IFD↔QuAPI 教训：默认 env 掩盖 test target 缺依赖）。

### 2.5 测试面

- `quapi_crosscheck.jl`：QuAPI 侧显式限定（`QuAPI.Δt(...)` 等，F6），
  断言保持逐位 `==`——升级后语义 = **IFD 核 vs QuAPI 原块位级对拍的常驻
  闸门**（迁移安全网 + 未来 IFD 演化的回归锚）；
- `correlation.jl`：容器引用经 reexport 透明，例行核对。

### 2.6 验收闸门（Part 1）

1. `quapi_crosscheck` 逐位 `==` 全绿（位级中性证明——IFD↔QuAPI 对拍契约
   在 GTEMPO 侧的复验）；
2. 双 golden 按值保持（weyl 8/8、exactttiif 6/6）；
3. 全套绿、零容差改动（现行 49 testsets/739 项口径——以 WIP 后实际为准，
   终点 = Part4 testset）;
4. **任何位级偏差 = 停**（对拍契约破坏，回退诊断；预期零偏差，偏差即
   IFD/QuAPI 对拍本身的 bug）。

## 3. Part 2 设计：严格 IF 入口（imag 端到端 + 路由面）

### 3.1 scheme 面接受 `ExactEqualDistantKernel`

`correlationfunction` 内按 (lattice 轮廓 × scheme) 分派：

- **imag lattice + Exact** → `IFD.Δτ(bath; N, δτ, scheme)` 返回元组
  `(ImagCorrelationFunction, det(BB))`（F2 标准容器）——解包：blocks 装配
  照常（下游零改动），`normdelta` 进 kernel（§3.2）；
- **real/mixed lattice + Exact** → `ArgumentError`（指向
  `exactif_port_plan`：StrictReal/MixedCorrelationFunction 的装配级消费
  ——骨架+踢契约、`Z = normdelta·det(装配矩阵)`——归该计划的
  `ExactDelta → HybridizationKernel` 块级路线，两计划接口对齐见 R5）；
- **BCS + Exact** → 按 F4 核验结果：IFD dispatch 已接线则同 imag 路线；
  未接线则明确 ArgumentError（消息注明 IFD 侧待实装）。

### 3.2 normdelta 归宿（待裁 A1）

`HybridizationKernel` 增第四字段：

```julia
normdelta::Union{Nothing,T} = nothing    # QuAPI scheme = nothing；exact = IFD 元组第二元
```

- 先例 = `coupling_norm`（同为核级标量属性）；T 泛型随核（det(BB) 复核型
  走 ComplexF64 核）；
- W2「纯数据核三字段」扩为四字段——**需用户重批**（A1）；备选案 = 严格核
  走独立构造器返回 `(kernel, normdelta)` 元组、kernel 结构不动（但
  `correlationfunction` 返回形状随 scheme 分裂，更丑）；
- normdelta 的 Z 消费契约（`Z = normdelta·Pf(B)` / RT 版）**不在本计划**
  ——exactif_port_plan 领域；本计划只保证字段忠实透传 IFD 返回值。

### 3.3 守卫与文档

scheme × 轮廓 × bath 族支持矩阵进 `correlationfunction` docstring；
`changes.md` 记录（后端切换 + scheme 面 + 依赖降级三事合一）。

## 4. 测试计划（验收闸门）

- **T1 crosscheck 升级**（§2.5）：QuAPI 显式限定 + 逐位 `==`——三轮廓
  （real 四块/imag τ块/mixed 九块含折叠号）现状断言全保。
- **T2 缺省回归**（Part 1 主闸门）：全套跑绿 + golden 按值（§2.6）——
  缺省 `QuAPI()` scheme 下一切逐位不变。
- **T3 scheme 透传**（Part 2）：`correlationfunction(bath, lat, IFD.QuAPI())`
  ≡ `correlationfunction(bath, lat)` 逐位 ==；imag + Exact →
  kernel.blocks ≡ IFD `exact_fermionic_Δτ` 直读元组第一位逐位 ==、
  `kernel.normdelta` ≡ 元组第二位（det(BB)）逐值 ==。
- **T4 拒绝面**：real/mixed + Exact → ArgumentError（消息含 exactif_port
  指引）；BCS + Exact 按 F4 核验结果锚定（接线则进 T3、未接线则断言抛）；
  非法 scheme 类型 → MethodError 自然冒泡（不加守卫）。
- **T5 物理锚（轻量）**：U=0 Anderson 单模 imag + Exact scheme 经 PartialIF
  装配 → `integrate` 值 vs `normdelta · (现行 QuAPI scheme 同格点的
  bath 归一 Z)` 的一致性（验证 normdelta 因子语义方向正确——完整
  `Z = normdelta·Pf` 消费闭环仍归 exactif_port，此为入口级 smoke）。
- **T6 全套回归**：`julia --project=. --heap-size-hint=8G test/runtests.jl`
  直跑（运维锚不变）；GM 套件不重跑；**test target 依赖解析验证**
  （QuAPI 降 extras 后 quapi_crosscheck 仍可加载——Pkg.test 教训的
  正向验证点）。

## 5. 实施顺序与提交切分

1. **feat(deps)!: IF 离散化后端切换 QuAPI → IFDiscretization（scheme 入口）**
   ——§2 全部 + T1/T2（单提交；位级闸门同笔验证）；
2. **feat(if): 严格 IF scheme 入口（imag + normdelta 透传）**——§3 + T3–T5
   （依赖 A1 裁决）；
3. `changes.md` 追加；§7 回填。

预估：Part 1 ~半天–1 天（机械但闸门必须全绿）；Part 2 ~半天–1 天。

## 6. 风险与边界

- **R1【IFD↔QuAPI 位级漂移】**：IFD 未来演化破坏对拍契约 ⟹ GTEMPO 全套
  值变——crosscheck 常驻闸门就是保险（任何漂移在 CI 层拦截，不进 golden）。
- **R2【reexport 命名冲突】**：IFD 导出面与 GTEMPO 自身导出/ImpurityModelBase
  导出的现有与未来冲突——实施时 `names(GTEMPO)` 冲突清单核对（当前无碰撞，
  静态核对一次即可）；今后 IFD 新导出名需在 GTEMPO 侧复检（写进 changes.md
  提示）。
- **R3【容器结构漂移】**：IFD 容器字段/branch 语义与 QuAPI 分叉 ⟹
  GTEMPO 消费端静默错位——T1 的 branch 直读逐位对拍覆盖（结构异 ⇒ 值异）。
- **R4【WIP 基线】**：同 polarif/bcs_exacttti R3——工作区未提交修改
  （correlationfunction.jl 在列）实施前确认基线。
- **R5【strict 消费缺口是暂时的】**：real/mixed + Exact 的 ArgumentError
  面由 exactif_port_plan 关闭——**两计划的接口对齐点** = 本计划的 scheme
  路由目标 == exactif_port 的 `ExactDelta → HybridizationKernel` 构造器
  签名与 normdelta 字段语义（该计划实施时以本计划 §3.2 的字段契约为准，
  避免二次裁决）。
- **R6【BCS strict 就绪度】**：F4 的 README TODO 与导出面矛盾——实施
  Part 2 前跑一次 IFD `exact_bcs_Δτ` 冒烟定谳。

## 7. 实施状态

**✅ Part 1 + Part 2 已落地（2026-11-19，单会话两批；基线 = 十一批后 HEAD
cbac36e9，R8 无冲突——计划行号基于 0f64cf78 已过时，触点面实测一致）**。

- **位级对拍探针（先行）**：IFD↔QuAPI **全轮廓逐位 0.0——含 BCS-RT 三味
  （uu/ud/du × 四支对）**。顶点移位只存在于 QuAPI 的 GrassmannContourIndex
  访问器层（branch 支视图从不经过它）——GTEMPO Phase R 停用访问器后，两
  后端在 GTEMPO 全部消费面上逐位一致；用户前提（IFD 已移除移位 ⟹ 与稠密
  装配对齐）实证成立。计划 Q2 的「BCS-RT 消费端访问器读取」已随 bcs_exacttti
  Phase R 删除，迁移面进一步收敛。
- **Part 1**：签名（第三位 scheme，缺省 IFD.QuAPI()）+ 六调用点 + 命名空间
  迁移（@reexport 换 IFD；R2 names 清单核对无碰撞）+ scalartype 四方法
  （名称重绑定即迁移，方法体零改动）+ deps 降级（QuAPI → extras+targets）。
  **偏差**：直跑 runtests 的 QuAPI 加载依赖共享环境栈（本机 @v1.10 已 dev）；
  干净机器需 dev QuAPI 或 Pkg.test（离线沙箱损坏的运维锚仍适用）。
- **Part 2**：imag + Exact 元组解包 + `normdelta::Union{Nothing,Float64}`
  字段（A1 按推荐；**F3 修正**——IFD det(BB) 随容器标量型返回 ComplexF64
  带 ~1e-15 虚部尾，构造器守卫后取实）；real/mixed + Exact ArgumentError ✓；
  **F4 定谳：BCS + Exact 三轮廓 dispatch 已接线**（README TODO 过时），
  correlationfunction 透传——δ 谱自然报错自 IFD 内部（"can not call to a
  delta/discrete spectrum function"，求解器要求可调用谱）。
- **T5 偏差**：计划预设的「normdelta·Z_exact ≈ Z_quapi」物理一致性**不成立**
  （normdelta ~1e11、语义 = det(BB) 全环闭式——`Z = normdelta·Pf(B)` 的
  Pf 非 port trace）——按计划边界该消费契约归 exactif_port_plan，本批降为
  透传 + 装配 smoke（exact IF bath 归一 Z 与 QuAPI scheme 同量级）。
- **闸门**：T1 crosscheck 逐位 14/14 ✓（QuAPI 显式限定，F6 规避）；T2 全套
  804/804（53 testsets）+ 三 golden 按值（weyl 8/8、exactttiif 6/6、BCS
  12/12）；T3 scheme 透传逐位 ✓；T4 拒绝面 ✓；T6 直跑验证 ✓。
