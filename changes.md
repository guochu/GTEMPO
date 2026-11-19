# 接口更改说明

## IF 离散化后端切换 QuAPI → IFDiscretization 批次（2026-11-19，plan/ifd_backend_switch_plan.md）

- **后端切换 + scheme 入口（Part 1，Breaking）**：
  `correlationfunction(bath, lattice, scheme=IFD.QuAPI(); band/orbital)` 第三位
  scheme 参数（缺省 = IFD `QuAPI()`，两参调用透明）——六个离散化调用点
  （normal/BCS × imag/real/mixed）全部改走 `IFD.Δτ/Δt/Δm(bath; …, scheme)`
  bath 对象直读分派（QuAPI `bath.f/β/μ/Δ` 手工拆包退役）。命名空间：
  `@reexport using QuAPI` → `IFDiscretization`（容器/branch/Δ* 面同名透明，
  实施时 names 冲突清单核对通过）；`import QuAPI: branch` 删。
  **位级中性实证**：IFD↔QuAPI 全轮廓逐位 0.0（含 BCS-RT 三味——顶点移位只
  存在于 QuAPI GrassmannContourIndex 访问器层，branch 支视图本就 raw，两
  后端逐位一致）。
- **依赖降级**：QuAPI 从 `[deps]` 移至 `[extras]`+`[targets].test`
  （crosscheck 对拍锚保留）；⚠️ 直跑 `julia --project=. test/runtests.jl` 的
  QuAPI 加载依赖共享环境栈（本机 `@v1.10` 已 dev QuAPI）——干净机器上须
  Pkg.test 或 dev QuAPI 后直跑。
- **严格 IF 入口（Part 2）**：`HybridizationKernel` 增 `normdelta::
  Union{Nothing,Float64}` 字段（A1 按推荐落地；frozen scheme = nothing，
  changebands/Base.:+ 透传）——imag + `ExactEqualDistantKernel` → IFD 元组
  `(ImagCorrelationFunction, det(BB))` 解包，normdelta 忠实透传（IFD 容器
  标量型带 ~1e-15 虚部尾，构造器守卫后取实）；**real/mixed + Exact →
  ArgumentError**（strict 容器 = 装配级骨架+踢契约，非朴素指数系数——消费
  路线归 exactif_port_plan）；BCS + Exact = IFD 三轮廓 dispatch 已接线，
  correlationfunction 透传（δ 谱自然报错自 IFD 内部）。normdelta 的
  `Z = normdelta·Pf(B)` 消费契约不在本计划（exactif_port_plan 领域）。
- **测试**：crosscheck 升级常驻位级闸门（QuAPI 侧显式限定——F6 自比对陷阱
  规避）；scheme 面 testset（缺省≡显式 QuAPI 逐位、imag+Exact 元组解包、
  real/mixed 拒绝、MethodError 冒泡、exact IF 装配 smoke）。

## BCS-RT raw 口径 + Weyl2×BCS 批次（2026-11-19，plan/bcs_exacttti_plan.md Phase R/W）

- **BCS-RT 去顶点移位（Phase R，Breaking）**：`_bcs_dense_block` 访问器物化
  退役（顶点移位内嵌口径），BCS 三轮廓统一 raw branch 视图直存——与 normal
  存储约定同构，GTEMPO 不再消费任何移位读取（QuAPI 包内 `_vs_rt_acc` 不动）。
  **T6 探针否决计划预留的「跨支构造期取负」**：bcs_Δt 味块支视图自带与
  normal Δt 相同的指数系数 Keldysh 号（Δ=0 四支对逐位 0.0），取反即错。
  结构锚：Δ=0 BCS 核逐位等于 normal 核；T3 BCS 实时 ED 锚翻绿（max abs
  gt 1.4e-4 / lt 6.5e-5 = 一阶标定档，δt 减半减半 ✓）——去移位前的 1.2e-2
  δt 无关偏差（polarif 批 @test_broken 记录的缺陷）即移位口径错位，就此闭环。
- **Weyl 二阶 × BCS（Phase W）**：`sysdynamics(...; order=2)` × BCS IF 端到端
  实证——K 链 bath 无关、`coupling_norm` = normal 通道 α 和（kernel 携带）。
  套件级闸门（test/bcs）：order=2 对 bcs_ed 偏差较 order=1 收缩 ≥10×
  （实测 ~40×，Nt=8：1.2e-4 → 1.2e-5）。
- **`validation/bcs_realtime_scaling.jl` 斜率标定**（N ∈ {4,8,16,32}）：
  order=1 斜率 1.07–1.09 ✓；order=2 改善 ~40× 但 Nt≥8 撞 ~1e-5 δt 无关地板
  （紧截断不变、随 Nt 线性缓增 = 每格常数残差特征）——渐近斜率 ≈ 2.00 被
  地板掩盖，疑二阶接触门 × BCS 交互（plan R7 排查面；正常 bath 九组基准
  1.996–2.005 无此地板）。**陷阱记录：ED 参照须按时刻对齐切片**
  （`gt_ed[(1:Nt).*stride]`）——首版 `gt_ed[1:Nt]` 半时刻错位伪造了
  order 无关大误差；relerr 分母随窗口变化也伪造斜率，斜率口径 = 绝对偏差。
- **BCS goldens 新档**：`test/fixtures/bcs_realtime_golden.json` +
  `check_golden_bcs.jl`（raw 口径四组合 PartialIF/PolarIF × order1/2 的
  Z/gt/lt 全时序，1e-6 档；12/12 绿）。
- **`fermionic_setup` 加 `alg` kwarg**（测试辅助面，缺省 PartialIF 不变）。

## PolarIF 装配算法批次（2026-11-19，plan/polarif_plan.md 落地）

- **`PolarIF(; rtol=0.0, trunc=Defaults.alg_trunc_if())`** 新算法：
  `expmult!×PolarExp`——fullfactor 的 SVD 秩一精确指数（Pairing = 对角剥离 +
  余阵 SVD；AnomalousPairing = 味扇区折叠逐扇区 SVD）。适用范围与 PartialIF
  一致（任意轮廓/kernel/对称性），恒就地单链；`rtol` = PolarExp SVD 秩截断
  （0 = 全数值秩缺省；AnomalousPairing 跨扇区全局归一）；`trunc` 随算法对象。
  `bcs_hybriddynamics` 守卫放宽 `Union{PartialIF,PolarIF}`（D7）；export 补。
- **闸门测试（if_gates）**：PolarIF ≡ PartialIF 对拍并入既有算例（复用装配态）
  ——虚时/U1/实时单模机器精度档（~1e-15，缩格后压缩无损）、混合连续谱
  5.3e-15（Nt=2/Nτ=3）、BCS 反常 1.2e-14（N=3）；跨带随机核 dense oracle
  ≤ 3.2e-12（vs 精确 exp）+ 恒单链契约 + rtol smoke（0.1 不掉秩逐位同、
  0.5 可观测差 4e-3）。虚时格 N≥4（N=3 触发热配对审计轮盘，已知）。
- **BCS 实时 ED（test/bcs 新 testset）**：🔴 **发现已知缺陷**——BCS 实时核块
  与 normal 约定错位（Δ=0 时 (+,+) 块首行全零/条目移位一格、(−,+) 块差一个
  Keldysh 号；δt 无关、随 t 线性、Δ 无关的系统偏差 max abs ~1.2e-2）。根因 =
  `_bcs_dense_block` 的 QuAPI 访问器物化保留旧栈顶点移位内嵌口径，新栈
  fullfactor 按标准网格消费。ED 锚以 `@test_broken` 记录（互拍锚有效——两
  算法消费同一核）；imag 侧不受影响（branch 视图）。修复方向：访问器读取
  补偿内嵌移位，待裁决。
- **`fermionic_setup` 加 `alg` kwarg**（测试辅助面，缺省 PartialIF 不变）；
  fewmodes realtime 阶语义复验：order=2 双算法（一阶 PolarIF 锚由 BCS 实时
  测试承接）。
- **套件减时**（不缩覆盖）：闸门对拍复用装配态 + 混合/实时格缩格
  （Nt=2/Nτ=3、real N=3）+ PolarIF oracle 4+1→2+1 算例；两个闸门 testset
  合计 ~5.5 min → ~0.5 min（闸门算例的装配数 12+12 → 15 共享）。

## anomalous 派生旗标 + 恒单链装配批次（2026-11-19，Breaking）

- **`HybridizationKernel.anomalous::Bool`**：是否含异味键（a·a / ā·ā——不满足
  荷守恒的块，BCS 反常通道）。**构造器从核键派生**
  （`any(key -> key[1].conj == key[2].conj)`）——单一事实来源，
  `changebands`/`Base.:+` 重装路径自动一致（normal + 反常相加 = 反常）。
- **`fullfactor` 双路由**：`anomalous = true` → `GM.AnomalousPairing`（紧凑
  变量表 2M×2M，BCS 反常键落同味象限，行为同前）；`anomalous = false` →
  `GM.Pairing`（矩形站点表，C = M×M 全矩阵、行 ā 列 a；带解析站号天然互异
  满足站点唯一性契约，接触项折 hloc 角）。non-anomalous 核若带非 canonical
  键响抛。dense oracle 对拍裁决两路由 ≡ 同一算符。
- **PartialIF 恒单链装配（多带逐带拆分退役）**：`_fullfactor` 单因子覆盖
  `kernel.bands` 全部变量（多带带对角核不再 `_band_kernel` 逐带构建），
  `_partialif_state` 恒 `expmult!(gmps, fullfactor, PartialExp)` 就地单链。
  **`hybriddynamics!/hybriddynamics` 的 PartialIF 返回恒为单 GMPS**（多带
  元组契约退役；ExactTTI 多带仍逐带元组）——消费方对单 GMPS 须自包装
  `(Is,)` 再走 `Is...` 变参（`fermionic_setup` 的 `Is isa GMPS && (Is = (Is,))`
  守卫语义从此恒真、对直传点不再需要）。装配值与旧逐带路径截断档一致
  （实测 5e-15）。

## 动态标量型批次（2026-11-19，Breaking）

- **包标量链全程动态（不再硬编码 ComplexF64）**：`correlationfunction.jl` 定义
  QuAPI 三容器 `scalartype`（`ImagCorrelationFunction{T}` 取 T；Real/Mixed =
  ComplexF64——QuAPI 结构约束；`BCSCorrelationFunction{C}` 委托 C）；
  kernel 构造端 `T = scalartype(corr)`——实谱虚时核全程 Float64，复谱/实时/
  混合核 ComplexF64。`fullfactor` 的 C = `zeros(scalartype(kernel), …)`；
  **虚时 lattice `scalartype` = Float64**（实链 × 实因子全程 Float64，GM 原生
  支持；复谱核 promote ComplexF64）。标量边界：`hybriddynamics` 真空起始 =
  promote(lattice, kernel)（不变）；ExactTTI 拟合层因子恒复 ⟹ 实链在
  `_exacttti_state` 装配入口复化（`_complex_mps`，标量乘提升 TensorMap——
  TensorMap 无广播 length、`complex.` 不可用；返回态 ≠ 输入别名）。
- **`_neg_block` 删除**：混合 τ↔rt 折叠直接 `-W`（CorrelationMatrix 走 QuAPI
  `Base.:-` 向量取负；Matrix = view 取负物化）。
- **核相加去稠密化**：`Base.:+(A::HybridizationKernel, B::HybridizationKernel)`
  同键块 = `va + vb` 直加（CM+CM = QuAPI Toeplitz 向量相加、混合/稠密 =
  LinearAlgebra 泛型——不再双物化 Matrix）；缺键透传经 `_convert_block`
  标量对齐（跨标量相加 Float64+ComplexF64 可用，同型零拷贝）；`_is_cm` 删。
- **BCS `orbital = nothing`（缺省）**：全部 Nambu 对（1:bands÷2）共享同一套
  杂化函数——四味块每支对只装配一次（味代表带 (1, 2) 读取，QuAPI 奇偶带味
  约定），多对键零拷贝别名（normal 多带同款）；`orbital::Int` 显式选对不变。
  `coupling_norm` 统一 = bands 数 × `coupling_norm(bath)`（单对 2×、全轨道
  2·轨道数×）；`bcs_hybriddynamics` 的 orbital 缺省同步 nothing。
- **ExactTTI 现代化**：`@kwdef` + `alg_expan` 升类型参数 `AE <:
  ExpExp.ExponentialExpansionAlgorithm`（`ExactTTI{Tr,AE}`）；`factor_rtol`
  field 删除（`_mult` = `SVDCompression(trunc=alg.trunc)`）；`cross_rtol` 校验
  保留在内构造器。

## ContourBlockIndex 腿味键批次（2026-11-19，Breaking）

- **`ContourBlockIndex(branch, band)` → `ContourBlockIndex(branch, band, conj)`**：
  核键自携腿味（`conj = true` ⟺ ā/η̄、`false` ⟺ a/η，GM 变量号奇偶约定）——
  HybridizationKernel 从此对任意 bath 自完备（branch + band + 腿味）。
- **`bcs` 味旗标退役**：`fullfactor(kernel, lattice; bcs)` → `fullfactor(kernel, lattice)`、
  `hybriddynamics!/hybriddynamics(lattice, kernel, alg; bcs)` → 去 kwarg——腿味
  从核键直读（`key.conj`），BCS 反常块（ā↑ā↓ / a↓a↑）与一般跨带法向耦合
  （Δ_b1b2）在键结构上可分辨，不再依赖调用端旗标；`bcs_hybriddynamics` 同步。
- **构造端写定**：`correlationfunction`（normal/BCS 双方法）构造期把腿味写进键
  （同带 = (ā, a)、BCS band1→band2 = (ā, ā)、band2→band1 = (a, a)，`_bcs_leg_conjs`
  味表 = 键的写入口）；`changebands` 键重建保留 conj；ExactTTI 守卫升级为结构
  判定（`_is_normal_conj`：所有键 = (ā, a)——手工构造的反常同带核不再漏进
  ExactTTI 路径，响抛 ArgumentError 而非 KeyError）。
- **消费端**：`_kernel_block`（exacttti）、`_kernel_eta`/`current_operator`
  （measurements）显式三参键；外部消费者（自建 blocks Dict 的测试/oracle）须
  改三参 `ContourBlockIndex`。

## IF 装配后端统一批次（2026-11-18，plan/if_expmult_backend_plan.md W1–W4）

- **`HybridizationKernel` 重定义（W2）**：`HybridizationKernel{L,C}`（lattice +
  QuAPI corr 本体 + `(Symbol,Symbol)` 视图 + diagnostics）→ `HybridizationKernel{T<:Number}`
  纯数据核——`blocks::Dict{Tuple{ContourBlockIndex,ContourBlockIndex}, Union{Matrix{T},CorrelationMatrix{T}}}`
  + `coupling_norm::Float64` + `bands::Vector{Int}`（`ContourBlockIndex(branch, band)` 新公开类型；
  `scalartype` 仅 Type 定义）。定义与构造迁 `correlationfunction.jl`。
- **构造端 API 更名**：`hybridization(lattice, bath; coupling_norm)` **退役** →
  `correlationfunction(bath, lattice; band)`（normal 单带/全带）与
  `correlationfunction(bcsbath, lattice; orbital)`（BCS 轨道）——normal/BCS 同一类型（F5）。
  `coupling_norm` 自动（单带 = `coupling_norm(bath)`，全带 = bands×，BCS 轨道 = 2×）；
  多浴耦合 = `Base.:+`（Toeplitz 向量相加、缺键透传、总耦合标量自动）；
  `changebands(kernel, [old=>new…])` 带重排/单带核扩多带（零数据改动）。
- **`hybriddynamics` 签名收窄**：`(lattice, kernel, alg; trunc)`（位置形主式）+
  `(lattice, kernel; alg, trunc)`（关键字形）——band/orbital 输入删除（信息在 kernel）。
  多带 kernel 返回逐带 IF 元组、单带返回单 IF。
- **存储约定（构造期符号收编）**：blocks 值 = IF 指数系数（纯实时跨支块保留 QuAPI
  内嵌 Keldysh 号；混合 τ↔rt 块构造期整体取负——旧 `_NegView` 惰性预翻物化退役；
  BCS 实时 = 访问器物化（顶点移位内嵌）+ 跨支号折叠）。**电流 η⁺⁻ = −blocks[(:+,:-)]
  负号保留**（电流算符 η = 无号物理块，与指数系数差 Keldysh 基一个相对号——公式层
  常数，ED 锚定；计划 §2.1 的「ηpm 负号删」按此裁决不实施，见 correlationfunction.jl 注）。
- **测量层**：`kernel_order` 退役（diagnostics 字段删，order 显式 kwarg 转发、缺省 1）；
  `kernel_coupling_norm` 一行化；电流 η 块读取 band 键化（`_kernel_eta`）。
- **`src/influencefunctional/`（W1/W3/W4）**：`partialif.jl` = `fullfactor`
  （AnomalousPairing 整形式、有向单一取向存储）+ PartialIF 经 `expmult!×PartialExp`
  （语义换代：同一 IF 算符的幂零指数表示，旧逐行因子积冻结 influence.jl）；
  `exacttti.jl` = `tti_imag/real/cross_factors` + ExactTTI 装配——虚时 ESS
  （上下三角级联 + 等时修正，极点压缩保留；**幂算法分派**：ExpDecay 族走
  PartialExp——WII 有序指数裁剪交叉对对环核不精确 3.2e-4，dense 对拍定谳）、
  实时/混合 Pairing×PolarExp（SVD 秩一；实时 = 拟合重构 C、混合 = 原始块 C
  ——混合格 ExpDecay cell 平铺无法限制分支子格）。混合 ExactTTI 从此可用
  （旧路径 cross 行腿落在 a 型腿 = 潜在味错位，无测试覆盖，已废弃）。
- **BCS**：`bcs_hybriddynamics` 换轨 fullfactor 路线；`bcs_partial_if` 签名改为
  `(lattice, kernel)`（块值直读，Keldysh 号构造期折叠）——deprecated 冻结档案。
- **band 输入退役（2026-11-18 第二批，band 信息在 kernel）**：`fullfactor(kernel, lattice)`、
  `tti_imag/real/cross_factors(kernel, lattice; alg)`、`influence_functional(lattice, kernel; alg, trunc)`
  全部不再接收 band——装配范围由 `kernel.bands` 决定。单带 normal / BCS 轨道核 = 单 IF；
  多带 normal 核 = 逐带 IF 元组（`InfluenceFunctional.bands` 字段、`state` 为 Tuple）。
  `fullfactor`/`tti_*` 仅接受单带 normal 或 BCS 轨道核（多带核逐带切片在装配内部，
  `_band_kernel` 视图零拷贝）。`hybriddynamics` = `influence_functional(...).state` 薄壳。
- **`src/influence.jl` 删除（2026-11-18 第二批）**：旧装配档案（partial_row/partial_if/
  memory_cell/embed_cells/contact_factor/embed_single/rankone_factor/旧 imag_/real_/
  cross_factors/bcs_partial_if）随入口切换期满退役。迁移存活符号：`_cell_var_index` →
  influencefunctional/partialif.jl、`_site_charges` → bridge.jl、`bcs_hybriddynamics` →
  partialif.jl（走 fullfactor 路径）、`integrate/integrate_scaled` lattice 便捷形 →
  influencefunctional.jl。`bcs_partial_if` 导出随之移除。if_gates.jl 同步改为互换/守卫
  闸门（vs 旧链 dense 对拍随参照删除，物理锚 = ED 套件 + golden）。
- **算法对象现代化（2026-11-18 第三批，用户三裁决）**：①`Base.@kwdef` → `@kwdef`；
  ②`PartialIF{Tr<:TruncationStrategy}`——乘法压缩截断随算法对象携带
  （`trunc` 字段 = `Defaults.alg_trunc_if()` 缺省，替代旧 maxdim/rtol 字段）；
  ③`ExactTTI{Tr<:TruncationStrategy}` 增加 `trunc`（同上）与 `alg_expan`
  （`ExpExp.ExponentialExpansionAlgorithm` 面，缺省 `Defaults.alg_expan()`
  = OverDeterminedProny(15, 1e-4)）字段——tti_imag/real 的核列拟合经
  `fit_exponentials(alg_expan, y)` 分派（ExpExp.exponential_expansion 薄壳，
  y[d] = Σ w·q^(d-1) 约定对齐）。**随动语义**：旧 :left/:right 锚点概念随
  alg_expan 退役（Pairing 装配原生支撑任意模长极点）⟹ 实时 "growing real
  pole" 守卫删除，DiracDelta |q|=1 边界极点核实时 ExactTTI 从此可装配
  （机器精度互换）；旧 rtol/maxterms 字段退役（配置随 alg_expan）。
  **精度注记**：缺省 Prony(15,1e-4) 对精确单极点核机器精度、连续谱短核列
  （L=3）残差 ~1e-1 档（小样本本征极限，ExpExp MatrixPencil 同档）⟹
  连续谱实时互换闸门放宽至 1e-2；高精度配置 =
  `ExactTTI(alg_expan = Defaults.alg_expan(tol=1e-10, n=40))`；锚定铅笔
  保留为独立 `GTEMPO.fit_exponentials(y; rtol, maxterms)`。
  `influence_functional`/`hybriddynamics`/`bcs_hybriddynamics` 的 `trunc`
  kwarg 随之退役（截断随算法对象）。
- **hybriddynamics 入口原地化（2026-11-18 第四批，用户四裁决）**：①关键字便捷形
  删除；②`InfluenceFunctional` 结构与 `influence_functional` 入口删除（导出移除）——
  `hybriddynamics` 收 `gmps`+`diagnostics`（truncerrs 包装进 `diagnostics[:truncerrs]`，
  通道审计键合并；截断随算法对象，`trunc` kwarg 无）；③新增原地版
  `hybriddynamics!(gmps, lattice, kernel, alg; diagnostics)`——单带 normal/BCS 轨道核
  就地装配 `gmps` 恒等返回，多带 normal 核逐带从 `copy(gmps)` 起步返回元组
  （**顺带修正**：ExactTTI 多带此前误把各带因子叠进同一 state，现逐带独立装配）；
  `hybriddynamics(lattice, ...) = hybriddynamics!(vacuumstate(promote_type(
  scalartype(lattice), scalartype(kernel)), lattice), lattice, ...)`——提升的是
  两者的 scalartype，配 `vacuumstate(::Type, lattice)`；④`Defaults.alg_expan` 改返回
  `ExpExp.AnchoredMatrixPencil`（ ExpExp 已有算法，n=40/tol=1e-10/pencil_rows=128
  与既有锚定铅笔路径同参）——ExactTTI 缺省分解行为与铅笔路线逐位一致
  （含 :right 锚点展平契约），OverDeterminedProny 仍可经 `alg_expan` 字段插拔。
- **hybriddynamics 二元组返回 + fullfactor 泛化（2026-11-18 第五批，用户五裁决）**：
  ①`fit_exponentials(y; rtol, maxterms, pencil_rows)` 关键字形删除（有界插值回退
  移植进 alg_expan 面：分解残差 > alg.tol 且 L ≤ alg.n 时触发，值判据驱动）；
  ②`fit_exponentials(alg, y)` 内 `as` → `amps`/`bases`；
  ③`fullfactor` 泛化为**紧凑站点表**构造——`GM.AnomalousPairing(C, sites)`，
  C 覆盖 kernel 作用域的全部变量对（2M×2M，M = 作用域站数），lattice 可含比
  kernel 更多的 bands（核外带不进因子）；单带/多带 normal 与 BCS 轨道核全部
  适用（旧单带/单 orbital 限定删除）；
  ④`hybriddynamics`/`hybriddynamics!` 的 `diagnostics` 从输入改为**第二输出**
  （不再作可变输入 kwarg）：`(state, diagnostics)` 二元组返回，约 40 处测试
  调用点解构迁移；extensions/secondorderif 五脚本同步 `first(...)` 迁移。
- **测试**：25 文件迁移新 API；新增 `test/api/if_gates.jl`（每因子 dense 对拍 +
  互换对拍闸门，12 项）。全套绿；两个 golden fixture 检查器按值保持
  （weyl 8/8、exactttiif 6/6，机器精度无需重锚）。
- 旧 API 恢复路径：git 历史 48fd249b..本次提交前（hybridization/旧 kernel/旧装配）。

## kernel.blocks 物化退役批次（2026-11-16，plan/kernel_and_fit_simplification_plan.md Part 1）

- **`HybridizationKernel` 存储改 = QuAPI 关联函数本体**（`Imag/Real/MixedCorrelationFunction`，
  Toeplitz 向量 O(N)）；`blocks` 稠密物化退役，改为**构造端一次建齐的窗口化 branch
  视图**（普通不可变字段）——`kernel.blocks[(b1,b2)][i,j]` 语法与数值**逐位不变**
  （CorrelationMatrix 查表 + 窗口视图；Δm 的 τ↔rt 交叉块净读数负号走 O(1) 惰性
  `_NegView` 包装），构造时做逐块形状/有限性校验（corr 与 lattice 不匹配即抛）。
- **内存/构造**：实轴 4×O(N²) 稠密 → O(N) 向量（N=2048 驻留 256 MiB → 129 KiB），
  构造端零物化复制（视图按需首建）；性能基准脚本 performance/kernel_memory_benchmark.jl
  （gitignore 区，磁盘交付）。
- **`validate_tti` 删除**：Toeplitz 平移不变是 CorrelationMatrix 存储的结构性保证
  （统一约定/Keldysh 恒等式由 QuAPI 构造器校验），旧校验退化为恒真；ExactTTI 入口
  留注释。BCS 路线不经过 kernel，不受影响。
- currents 两处（electriccurrent/heatcurrent）循环底板一次性物化
  `Matrix(kernel.blocks[...])`（与旧稠密循环同性能）。
- fit_exponentials（Part 2）按裁定维持原样不动（案 B 已在 ExpExp 侧另行立项）。
- 注意：对拍测试的 Δm/Δt 装配参数须与 hybridization 同源（如 `t = Nt*δt`——
  3*0.1 ≠ 0.3 的浮点位型差会让两次独立 quadgk 装配差 ~1e-17）；quapi_crosscheck
  已升级为逐位 `==` 硬闸门。

## 对称性单例化批次（2026-11-16，Breaking）

- **Symbol 字符串手选退役，改为对称性单例类型**：`abstract type GrassmannSymmetry`
  + 四单例 **`Z2Symmetry`/`Z2Z2Symmetry`/`U1Symmetry`/`U1U1Symmetry`**（export，
  Symmetry 后缀避免与 GM 扇区类型撞名；GM 扇区 alias `Z2Z2` 的 re-import 因此恢复），
  `sectortype` 按类型/实例分派（`sectortype(U1U1Symmetry()) === ChargeCharge`，
  扩展自 TK.sectortype）。
- **`symmetry_spec(::Symbol)` 与 SYMMETRY_SPECS 注册表删除**（export 撤销）；
  `bcs_symmetry_guard` 改收单例（`bcs_symmetry_guard(Z2Symmetry())`）。
- **lattice 按单例类型参数化**：`AbstractGrassmannLattice{Sym <: GrassmannSymmetry}`，
  `symmetry::Sym` 字段携带单例；构造器 kwarg
  `symmetry::GrassmannSymmetry = Z2Symmetry()`（Symbol 输入在 kwarg 类型层即
  TypeError）；`similar` 缺省 `symmetry::GrassmannSymmetry = x.symmetry`（保持
  自身单例）；**`TK.sectortype(lattice)` 编译期即得扇区类型**，下游
  （hybriddynamics/sysdynamics/measurement 等）不再重复传 spec——
  `sectortype(lattice)` 一查即得。
- 桥层荷代数/装配函数（var_charge/qspace/pairize/monomial/…）签名不变，
  仍按 TK 扇区类型分派（`sectortype(对称性单例)` 的产物）。

## convenience 清理 + src 布局重排批次（2026-11-16）

- **`unscale` 删除**（导出 API 移除）：取值直用 **`value`**（GM 原名，本包新增
  export；`GM.ScaledScalar` 的标准取值口）。
- **`scaled_ratio` 删除**（导出 API 移除；src 零调用点）：比值物化直用
  `value(a / b)`（GM.ScaledScalar 的 scaled 除法原生支持）。`scaling.jl`
  整文件删除（ScaledScalar/value 透传移入主文件 using+export）。
- **`zero_charge` 删除**（非导出内部助手）：荷单位直用 `one(G)`（扇区类型上
  语义恒等）。
- **`parity_label` 删除**（非导出内部助手）：与 GM `fermionparity` 逐方法
  相同（四对称形 χ 同态注册一致），oddlater Koszul 符号调用点直用
  `GM.fermionparity`。
- **src 布局重排**（公开 API 零变化，git mv 保历史）：`lattices/lattices.jl` →
  `lattices.jl`（单文件不设目录）；杂质模型层 `sysdynamics/` →
  **`impuritymodels/`**（入口文件同步更名 `impuritymodels.jl`——该文件是
  AbstractImpurityHamiltonian 族，非 K 构造）；K 构造宿主 `system.jl` →
  **`sysdynamics.jl`**（sysdynamics/sysdynamics_weyl2order 的宿主文件名实对齐）。

## 桥层接口简化批次：spec = TK 扇区类型 + 截断走 GM 标准通道（2026-11-16，plan/bridge_interface_simplification_plan.md A1–A4 三阶段）

Breaking 变更清单（A2 接受的返回类型变更 + 退役面）：

- **spec：单例 → 扇区类型本身**。`symmetry_spec(:sym)` 现返回 `Type{<:Sector}`
  （`:z2→Z2Irrep`、`:z2z2→Z2Z2`、`:u1→U1Irrep`、`:u1u1→ChargeCharge`）；
  `AbstractGrassmannSymmetry`/`Z2Symmetry`/`Z2Z2Symmetry`/`U1Symmetry`/
  `U1U1Symmetry`/`symmetry_symbol`/`sector_type`/`sector_spec`/`chain_spec` 删除。
- **荷标签 = TK 扇区**：`var_charge`/`cell_charges`/`fock_charges` 返回值从
  Int/NTuple/Z2Label 改为 TK 扇区；`charge_dimensions` 的 Dict 键从标签改扇区；
  `OperatorInsertion(; charge=…)` 接受扇区。Z2Label 与 sector/charge_label
  双向转换层删除；荷代数 = `qadd`/`qsub`/`qneg`/`zero_charge`（扇区泛型，
  语义不变）。
- **截断：`Truncation` 升级为 `MAK.TruncationStrategy`**。构造 kwarg
  `maxdim` → **`D`**（`Truncation(D=128, rtol=1e-12)`）；**`cutoff` 字段删除**
  （全库无调用方；无状态实现无法复现 logscale 漂移补偿，保留即无声错误——
  `PartialIF(; cutoff=…)` 同步删字段）。
- **`multiply`/`compress!`/`CompressionReport`/`LegacyTruncation` 删除**：
  压缩乘法 = **`GM.mult`**（已 export）、正则化压缩 = **`GM.truncate!`**
  （新增 export）；误差唯一出口 = 两者返回的 truncerr 标量（逐键相对截断
  误差之和）。`partial_if`/`bcs_partial_if` 的 `reports` kwarg →
  `truncerrs::Vector{Float64}`；`InfluenceFunctional.reports::Vector{Float64}`
  （逐步 truncerr，逐键分解/saturated/_logabs 记账不再输出）；
  `_system_factor_reference(reports=…)` → `truncerrs=…`；ExactTTI diagnostics
  键 `:factor_cleanup_reports` → `:factor_cleanup_truncerrs`。
- **`BoundedCache` 删除**（GTEMPO 全树零调用点）。
- ** tutorials（firstorderif 家族 14 文件）截断面机械迁移**：
  `Z2Tensors.TruncationScheme` 注解 → `Truncation`，`truncdimcutoff(D, ϵ)` →
  `Truncation(D, rtol)`（文件其余部分仍为换轨前 stale API，不在本批范围）。

语义保持：荷标签排序走 `sortkey` 复现旧序（键扇区插入布局逐位不变，golden
四对称形 43 项 vs 退役前逐位一致）；截断原理 = LegacyTruncation 去 cutoff 去
记账（top-D 全局预算 ∧ 相对尾能量 rtol² 判据，rtol 是相对 L2 能量判据、与
GM TruncateDimCutoff.ϵ 的逐值噪声地板无换算关系）。验收：全套 38 testset 绿；
longchain n=8/16 新旧对拍最差 7.2e-15、n=32 全工况有限无溢出；零链回归
（coeff=0 单项式 truncate!/mult 迹归零、键维不塌 0）。

## D1 修订批次：多轨道对称性支持（2026-11-16，进行中）

- **bands_guard 语义修订（D1 修订，用户裁决）**：旧「:z2z2/:u1u1 要求 bands≤2」
  不合理——多轨道（每轨道 up/down 两带）完全支持荷分辨形；新守卫 = **偶数带**
  （轨道交错：奇带 = 各轨道 up、偶带 = down）。**var_charge 轨道交错修正**：
  旧 `band == 1 ? up : down` 写法在 bands ≥ 3 时把轨道 2 的 up 错标 down，
  改 `isodd(band)`。
- **integrate 链上限 3 → 5**（GM 引擎均匀循环真实上限 K≤5；norb=2 = K+4IF
  满员；norb≥3 走合并形 :z2/:u1）。
- **多轨道端到端验证状态（norb=2 KanamoriIM, bands=4, 对角 1:1 混杂，ED =
  kanamori_ed 原生对角约定）**：①裸 K 链迹 = 原子极限 ED 机器精度
  （:z2/:u1/:u1u1 三形 4e-14）；②occupation(band1) vs kanamori_ed =
  **1.7e-4** ✓（δτ=β/6；Nτ=12 时 5.5e-5 收敛 ✓）；③绝对 Z（正确口径
  raw Z vs ED_abs/Z_B）= **0.935（Nτ=6）→ 0.966（Nτ=12）一阶收敛 ✓，
  与 bands=2 套件同行为——无缺陷**。**勘误（2026-11-16 上游复核，
  GrassmannMPS/bugreports/2026-11-16-upstream-verdict-multiband-z.md）**：
  此前报告的「绝对 Z 乘性抑制 ~4.7×」为对拍口径不一致的假阳性——port 的
  `integrate` 约定返回浴归一化配分（Z_full/Z_B，由本套件 bands=2 测试
  raw Z vs ED/Z_B 通过钉死），而探针脚本对 port 侧又除了一次 Z_B（双除）
  ⟹ 0.2116 假象；正确口径数字在原复现脚本第二组打印中本已给出
  （Z_port×Z_B vs ED_abs = 0.934）。「:u1 合并形 5.58×」同为配置错误
  （partial_if(band=1) 在 bands=4 格只覆盖 band-1 变量槽 = 另一模型；
  多带正确用法 = 逐带建 IF）。④**:z2z2 bands=4 间歇 segfault 未确认、
  暂不可修**（单跑不复现 + 崩溃时本机并行 Julia 负载强混杂；需
  ASAN/valgrind + 无并行负载复现后另行立案）。本批次的守卫放宽/
  var_charge 修正/链上限由 ①② 与既有电池背书，不受勘误影响。
  **测试侧勘误收尾（2026-11-16）**：test/multiorbital/runtests.jl 同源
  假阳性 FIXME（「J≠0 缺 54%、疑 J 四算符槽分解标定」= port 侧双除 Z_B：
  Z_port/Z_B = 1.43 误配 ED/Z_B = 3.08）已改为正确口径叙事，断言原样
  （实测比值 0.974（Nτ=8）→ 0.987（Nτ=16）一阶收敛 ✓）。

## 命名收敛批次（2026-11-16）

- **对称性类型重命名**：`AbstractSymmetrySpec` → **`AbstractGrassmannSymmetry`**、
  `Z2Spec` → **`Z2Symmetry`**、`Z2Z2Spec` → **`Z2Z2Symmetry`**、
  `U1Spec` → **`U1Symmetry`**、`U1U1Spec` → **`U1U1Symmetry`**（src + tests；
  `symmetry_spec(::Symbol)`/`sector_spec` 等函数名不变）。
- **`_spec_symbol` → `symmetry_symbol`**：移到 bridge.jl（与对称性类型定义同址）、
  **export**（bands_guard 与 chain_spec 的 :z2/:z2z2/:u1/:u1u1 反查）。
- **convenience 函数删除**：`bond_dimensions`（= GM.bonddims，内部 CompressionReport
  调用点直换）、`scaling`/`setscaling!`（= exp/logscale 换算皮肤）——消费方
  （含测试）改用已导出的 GM 名 `bonddims`/`logscale`/`setlogscale!`。CompressionReport
  的 `bond_dimensions` 数据字段名保留（描述性字段，非函数）。
- 旧 Phase-0 fixtures（check_golden.jl/dump_golden.jl，引用已删 scaling 且全
  旧 API）退役删除；live 的 check_golden_exactttiif/weyl + golden_json.jl 保留。

## ordering 类型层删除批次（2026-11-16）

- **AbstractChronoOrdering 及三个单例（ChronoImagOrdering/ChronoKeldyshOrdering/
  ChronoMixedOrdering）删除**：三个 lattice 结构体去 ordering 类型参数
  （`ImagGrassmannLattice{O}` → `ImagGrassmannLattice` 等），构造/similar/
  _index 派发签名同步去参；GTEMPO.jl 撤销三导出。**lattice 构造无 ordering
  输入**——排布固定 Chrono，是格点几何本身而非可插拔策略（GrassmannOrdering
  家族删除的余步）。块序描述保留为注释文档。代码触点仅 export 一行 +
  lattices.jl（外部零依赖），全套测试通过。

## 接口收敛批次：imports 集中主文件 + ScaledComplex → GM.ScaledScalar（2026-11-16）

- **using/import 全部集中主文件**（GM 包同款纪律）：bridge.jl / system.jl 头部
  的重复 import 段删除；主文件补 `LinearAlgebra: tr/norm/diag`、
  `TensorKit: SectorVector/SectorDict`、`import MatrixAlgebraKit as MAK`。
- **ScaledComplex 退役，公共 API 改为直接透传后端 `GM.ScaledScalar`**：
  GM.ScaledScalar 是旧包装的严格超集（参数化实/复相位——实链结果保持实型、
  同构造归一化、`value` ≡ 旧 unscale（同 OverflowError 守卫）、完整 scaled
  算术 log 域精确）。公共 API 变化：`ScaledComplex` → **`ScaledScalar`**（构造
  签名不变：`ScaledScalar(logabs, phase)` / `ScaledScalar(z)`）；`unscale`
  保留（= `GM.value` 历史名别名）、`scaled_ratio(a, b) = value(a/b)` 保留
  （零分母同抛 ArgumentError）。measurement.jl 的 `_scaled` 换皮删除、
  bridge.jl 的两处 ScaledScalar→ScaledComplex 转换直取。行为等价，全套测试
  通过。

## BCS occ 平反 + 性能基线 + alignment sweep 批次（2026-11-16，grassmannmps-backend 分支）

- **BCS「occ 0.55 vs 0.32」假 bug 平反**：上一批次 @test_broken 的分诊定位到
  根因 = **测试的 ED 参照算符约定错配**——bcs_ed 自带独立算符约定
  （occupied=首基、自带 JW 串），其 H 矩阵配 fermion_annihilators 的位基算符
  = 不同的物理模型（bcs_ed α→0 收缩给出 n=0.310，同对角的
  singlemode_ed/正确物理给 0.548，port 给 0.5501=正确值+δτ 修正）。用 bcs_ed
  自带算符（Ab/Bb，旧套件同款）重锚后：**Z/occupation/奇 matsubara（τ=β=−⟨n⟩
  恒等式）对 ED 全部 δτ 修正档（~3e-5）吻合**——反常 IF 通道（含奇插入穿过
  a/ā 混合列）被 ED 完整锚定，GM 引擎对混合型楔号的嫌疑亦被实证否定。
  @test_broken 全部转正。
- **matsubara/occupation/nn/greater/lesser/kernel 签名放宽**：
  `kernel::Union{Nothing,HybridizationKernel}`（BCS 路线无 normal kernel，
  kernel 仅用于 order/coupling_norm 元数据缺省）。
- **GM 快照扫掠缺口闭合复验**：`GM.expectation(op, env)` 右扫对数幅值补偿
  （cprime-logscale 守卫）与本精确路线对拍 **12 位吻合**（GTEMPO measurement.jl
  的旧「缺口」注释过时，已更新）；精确路径保持缺省，多算符高频查询可换
  快照路线（值等价）。
- **性能基线**（plan R6）：`validation/longchain_benchmark.jl` +
  `validation/records/longchain-20261116/`（CSV+README）。real/imag ×
  order 1/2 × n∈{8,16,32} × reps=3：IF 构建主导（每胞加倍 ~4.3×，n=32 real
  o2 ≈ 17.5s）；K/integrate/观测量近线性（n=32 全部亚秒）。与 U1GTEMPO2
  基线同量级（imag n=16 IF 双栈均 ~1.0s）；stage 语义两栈不同构，本表为
  本栈纵向回归基线。
- **alignment sweep**（Phase 6 验收项）：`validation/alignment/`（双栈脚本 +
  records + REPORT.md）。S1 虚时单带 / S3 双带 U=1：Z 与 occupation 对旧栈
  **机器精度一致**（1e-15 档）；S2 实时 greater：节点约定差一位（旧 k=N+1
  含 0 端点 vs 新 k=N）+ 热态构造截断差 2e-3，移位后逐点 1e-5 吻合。

## 奇算符 GF 修复 + 测试套件重建 + BCS/multiorbital 重接批次（2026-11-15，grassmannmps-backend 分支）

上一批次的「已知未决」奇算符 GF 缺陷**全案告破并修复**（对拍 U1GTEMPO2 机器精度），同时完成 Phase 6 测试套件重建与 BCS/多带重接：

- **奇插入 GF 五处口径差**（greater/lesser/matsubara/raw_lead_correlator；上一批次
  的 ~e-49「port split/pairize 有罪」裁决**撤回**——该结论基于手动复刻 builder，
  本次实验证明手动版在两栈皆死：①参照手动链（prefix 不含 cut-0 荷）迹恰零；
  ②补上 cut-0 荷后与参照真身机器精度一致 ⟹ **port 的 split_cell/from_core/
  pairize/桥层全数无罪**）。真凶与全部口径差（measurement.jl 文件头有完整记录）：
  ①**builder prefix 未含 cut-0 插入荷**（system.jl：prefix 初值须 = cut-0 插入荷、
  逐胞累计而非重置——偶算符 charge=nothing 不受影响，故 occupation/nn 全绿而
  奇 GF 全垃圾）；②荷方向反转（d/ddag 应取算符自荷 −q/+q，此前反向且
  OperatorInsertion 无荷一致性校验——参照栈有）；③`_normalized_two_point` 忽略
  调用方显式 order/coupling_norm kwargs（hybridization diagnostics 不登记
  order → **greater 以 order=1 链对 order=2 的 Z 归一**，恒定复因子 1.07 假象）；
  ④push 序（cutleft foldr 首 push 最左 = 最后作用）：同切位 0⁺ greater = d·ddag、
  lesser/电流 = ddag·d；⑤matsubara 缺时序因子 (i≥j ? −1 : +1)。
  验证：实时 greater/lesser 3 点、虚时 matsubara 4 点、occupation 4 点、
  lead 关联 6 点（含同切位 a†a 极限）Δ ≤ 5e-15。
- **BCS 反常 IF 重接**（`influence.jl` bcs_partial_if + bcs_hybriddynamics，
  D7 仅 :z2；`correlationfunction.jl` 加 BCS 派发）：行/列腿覆盖指数四类单项式
  (ā↑,a↑)/(ā↓,a↓)/(ā↑,ā↓)/(a↓,a↑)，系数直读 QuAPI BCSCorrelationFunction
  （有向指数系数 + 实时跨支 Keldysh −1 补偿，旧栈 partialif_hybrid 同款）。
  验证：Δ=0 退化为两条独立正常带（Z 比值机器精度）；Δ≠0 Z 对 bcs_ed ED
  机器精度档（**归一化分母 = 配对浴配分**，Δ 项破坏粒子数、自由浴配分不可用）；
  occ 带 @test_broken（0.55 vs 0.32，反常通道对单粒子观测的路径 bug，
  分诊路径在 test/bcs/runtests.jl FIXME）。
- **matsubara 混合轮廓放行**：τ↔τ 插入在 Kadanoff 轮廓合法（虚时支 = 热态腿）。
- **kernel=nothing 修复**：kernel_order/kernel_coupling_norm 加 ::Nothing 方法
  （occupation/nn 的显式 order/coupling_norm 路径此前被 kwarg 默认值打断）。
- **cached_gf 补 cache kwarg**（透传兼容 cached_correlation；D5 精确路径语义不变）。
- **symmetry lattice accessor + symmetry_spec/var_charge/cell_charges/fock_charges
  导出**（lattices.jl：export 已声明但 accessor 缺失）。
- **CUDA ext 删除**（D1 裁决：ext/ 整目录移除、Project.toml weakdeps/ext 条目清理）。
- **测试套件重建**（Phase 6）：api/ 12 文件新 API 全绿（lattice/bridge/integration/
  correlation/hybriddynamics/sysdynamics/entropy/observables/quapi_crosscheck/
  physics_battery/port_smoke）；normalbath fewmodes（imag/real/mixed × 单模/
  两模/两带 U≠0/quench/Td，ED 锚）+ continuous（Toulouse 解析锚；imag GF 5e-3、
  实时二次型精确锚）+ BCS（Δ=0 精确退化、反常 ED、D7 守卫）+ multiorbital
  （Kanamori norb=1，occ/nn 机器精度）。**ED 锚定铁律**：port 虚时 Z = 浴归一化
  值（ED 须除纯浴配分 Z_B；Δ≠0 时用配对浴配分）+ 行因子化 IF 的一阶 O(δτ)
  常数（fewmodes ~0.1/δτ、两带 U≠0 ~0.45/δτ，容差按此标定）；实时离散化精确
  （常数 H 精确指数切片合成），容差 1e-4~1e-5。旧栈 Z2Tensors 测试（grassmannmps/
  grassmanntensor/mpo/mult/electronphonon）随栈退役删除；toulouse_Gτ 的随机宇称
  估计器（同参数跨进程整体 ± 翻转）在解析锚处按符号不定分支处理。

本轮重构涉及的接口更改汇总。所有更改均已通过全量测试验证（行为不变或数值等价）。

## GrassmannMPS 后端换轨：M=2/跨对称性/ExactTTI 修复批次（2026-11-14，grassmannmps-backend 分支）

换轨后新栈（Chrono 格 + GM 引擎）的缺陷修复与验证批次，全部以 ED/golden 机器精度锚定：

- **M=2 带互换**（`system.jl` local_weights）：m 槽解码 band 1 = 高位组 vs
  split_cell charges 序 LSB 剥离（band 1 在前）冲突，bands≥2 时带内容互换
  （M=1 平凡不敏感故 golden 全过）。症状：单带浴 occupation 给 α=0 常数、浴
  落在与测量相反的带。修复 = 解码改 `(m >> 2(b−1))`（U1GTEMPO2 同款）。验证：
  U=0 带因子化 ⟨n₁⟩≡M=1 机器精度；ED（绝对 ω 约定）real/imag × U=0/0.8 全收敛。
- **:z2z2 flavor**（`bridge.jl` station_flavor）：Z2Label 与 Int 0 比较
  恒 false → 荷分辨形全部路由 :down 站 → SectorMismatch。修后四对称形
  （:z2/:u1/:z2z2/:u1u1）M=2 逐位相等。
- **ExactTTI Chrono 布局适配**（`influence.jl`）：①real/imag 因子族轮廓守卫；
  ②imag_factors `z=sqrt(wZ)` 遮蔽荷标签；③embed_single 站号表存储序；
  ④real 通道/memory contact 在虚拟 ± 相邻布局构建后 GM.reorder! 到 Chrono 站序；
  ⑤虚时 compact 通道链与 Chrono [格 Nτ…1] 天然同序，删除多余 reorder!。
  修复后 ExactTTI ≡ PartialIF（IF 表逐位一致 ≤3.3e-15）且 ED 锚定
  （imag 8.7e-4 / real 2.1e-5）。旧栈 exactttiif fixtures 改为仅参考
  （其编码旧栈自身在案的 Prony ~5% 缺陷，新实现已超越）。
- **Phase 5b 算符路线**（`measurement.jl`）：lattice_monomial（SparseGMPS
  单项式、荷沿键穿线、中性闭合）+ environments_cache/cached_Z（O(1) 读左端）
  + cached_gf（op∧K 楔积 + 单发 Zvalue 精确路径）。D5 裁决：边界观测保持
  inserted-K 唯轨；算符路线 = 原始关联层（胞坐标 0⁺ 点，与边界插入差
  O(接触项)）。GrassmannMPS 上游配套修 expectation() 的 cprime logscale 守卫
  （右扫 rightlogscales 补偿缺口仍开放，见 GM 包）。
- **nn 密度关联**：双插入精确版，ED 两时点验证 ≤9e-4（δt 二阶）、
  同切位 nn(i,i)=occupation 精确。
- **已知未决**：greater/lesser/matsubara/currents 的奇算符（d/ddag）插入
  K' 链 label-flow 缺陷——插入必须显式带荷（已改，链可构）但迹仍错
  （~e^{-155}，reorder 吸收 logscale −114）。Phase 5a 起 greater/lesser
  从未被断言（golden 只锚 occupation）。electriccurrent/heatcurrent 公式
  已就位但被此问题门控（measurement.jl 文件头有完整记录）。
- **测试资产**：test/api/physics_battery.jl（23 项：ED/因子化/四形一致/
  带对称）；test/fixtures/check_golden_exactttiif.jl（8 项：ExactTTI≡PartialIF
  + ED 锚）；test/fixtures/check_golden_weyl.jl（8 项，机器精度）。
  ED 参照约定铁律：浴模在绝对 ω 无 μ 移；U≠0 时 band-2 浴经 U 反馈，
  单浴 ED 只匹配 (IF₁, vacuum) 设定。

## BCS 影响泛函 ↓↓ 正常项装配修正（2026-11-03）

- **问题**：`bcsinfluencefunctional/partialif.jl` 的 BCS 单 IF 装配把 ↓↓ 正常通道实现为
  `(a↓, ā↓)·(−Δ)`。对 Grassmann 单项式 −(a↓ᵢ)(ā↓ⱼ) = +(ā↓ⱼ)(a↓ᵢ)，这等价于 ↓↓ 核取
  **转置**（canonical 形式应为 `(ā↓, a↓)·(+Δ)`，与 QuAPI Dense 的 `_delta_terms` 分解
  一致）。U=0 时 band-2 只贡献一个标量因子、在 Z 归约中相消（因此 U=0 的测试与基准
  不暴露）；U≠0 时产生 δt 无关、随演化时间累积的系统偏差——firstorderif bcs_realtime
  基准实测 gt 误差卡 ~1.5e-3（p≈0.06）；截断（D=200/400/800 逐位相同）与 Δ=0/0.6
  对照均排除后定位到装配。
- **修复**：重定义 QuAPI `contourindex.jl` 的 `getindex(::BCSCorrelationFunction,
  i::GrassmannContourIndex, j::GrassmannContourIndex)` 为**指数系数语义**——与 Dense 侧
  `_delta_terms` 同一分解：四类单项式 (ā↑,a↑)/(ā↓,a↓)/(ā↑,ā↓)/(a↓,a↑) 全部 + 系数，
  其余腿对（同 conj 的 uu 误配、跨 flavor 反向对）置 0，块值的分支/顶点位移约定不变。
  与普通 bath 访问器的收缩语义（x[a,ā]=+Δ、x[ā,a]=−Δ）不同，BCS 是有向的指数系数，
  装配处直接 `corr[i, j]` 读取、无需 conj 翻转。GTEMPO 的行/列腿结构相应改为
  rows {ā↑, ā↓, a↓} × cols {a↑, a↓, ā↓}（行腿自身节点从列集合剔除——(band2,·) 既是行
  又是列，ξ² = 0 且位置不可重复）。转置（列驱动）的
  `partialif_hybrid(rows, col, ::BCSCorrelationFunction)` 尚无调用方，先显式 `error`。
- **验证**：Anderson（U=1, Δ=0.6）bcs_realtime 基准从误差底 p≈0.06 恢复一阶收敛
  （gt p=1.14 / lt p=1.04）；U=0 基准不变（p=1.02）；test/bcs 全部 33 项通过，
  normalbath / electronphonon / multiorbital 全量通过。测试配置（T=0.25, Nt=5）的
  ED 偏差从 4.2e-4/2.1e-4（大部分为旧系统底）降至 7.2e-5/1.1e-4（纯离散化），
  test/bcs 实时 ED 容差相应收紧 1e-2 → 1e-3。顺带发现：api/lattice 的 4 个
  UndefVarError 是 ImpurityModelBase 与 GTEMPO 都导出 `branch/branches` 的命名冲突
  （预先存在，与本次修复无关；已通过测试端 `using GTEMPO: branches` 消歧解决）。
- **对照论文（PRB 2604.23301, Eq.31/34/36/37）**：论文正文把 ↓↓ 正常项写成
  `a↓ Δ↓↓ ā↓`（a-first，Eq.36 第 4 项）且块取折叠版（Eq.34：Δ↑↑ = Δuu − Δᵀvv、
  Δ↓↓ = Δvv − Δᵀuu、Δ↑↓ = Δuv + Δᵀuv、Δ↓↑ = Δvu + Δᵀvu）。数值判决（τ 环，
  U=1, Δ=0.6, Nτ=10）：实现版 relerr 2.18e-04，论文 Eq.(34)+(36) 字面版 3.48e-02
  ——论文正文的折叠/次序表述与精确 Wick 展开不符（Grassmann 反换位使 a-first
  写法等价于转置核，即本次修复的 bug）。正确表达式为 canonical 四类
  (ā↑,a↑)/(ā↓,a↓)/(ā↑,ā↓)/(a↓,a↑) → raw uu/ud/du 块、全 + 系数（`getindex` 与
  Dense `_delta_terms` 同款）。论文数值结果来自 ED 验证的代码路径、不受影响；
  正文公式建议勘误（`bcscorr.jl` 的 struct docstring 已改为 canonical 表达式）。

## Fock lattice 类型重命名与 LeftExpectationCache 删除（2026-10-25）

- **类型重命名**（与 Grassmann 侧的 `AbstractImagGrassmannLattice`/`ImagGrassmannLattice` 命名对齐）：
  - `ImagFockLattice`（抽象）→ **`AbstractImagFockLattice`**；`ImagFockLattice1Order`（具体）→ **`ImagFockLattice`**；
  - `RealFockLattice` → **`AbstractRealFockLattice`**；`RealFockLattice1Order` → **`RealFockLattice`**；
  - `MixedFockLattice` → **`AbstractMixedFockLattice`**；`MixedFockLattice1Order` → **`MixedFockLattice`**。
- 具体类型名即构造函数（`ImagFockLattice(; δτ, N, bands, ordering)` 等）；`FockLattice(; contour=..., kwargs...)` 快捷构造不变；原 `Imag/Real/MixedFockLattice(; order=1, kwargs...)` 的 order 关键字删除（二阶从未实现），调用处请直接用具体类型构造函数或 `FockLattice`。
- **`LeftExpectationCache` 删除**（`observables/envs.jl`）：4+ GMPS 环境缓存的占位实现及其 `expectation`/`expectationvalue`/注释代码一并移除；`environments` 只返回 `TwosideExpectationCache`。

## electron-phonon 路径重构与 lazy-sum 接口移除（2026-10-24）

### electron-phonon：`FockContourIndex` 与统一的 partial IF

- **新增 `FockContourIndex(j, branch, band)`**（`src/electronphonon/fockcontourindex.jl`，已导出）：electron-phonon partial IF 的腿指标，比 `BosonicContourIndex` 多一个 `band` 字段；`<: QuAPI.AbstractContourIndex`，`corr[i, j]` 读取玻色 Δ 矩阵元（band 不参与取值）。
- **`partialif_hybrid`（Fock 版）** 改为按腿构造：`partialif_hybrid(lattice, rowind::FockContourIndex, colinds::AbstractVector{<:FockContourIndex}, corr)` 返回 `exp(Σ Δ n̂_row n̂_c)`，每个 `(时间, 分支, band)` 只做一次键维 2 的 density-MPS 构造（内部 `_partialif_densemps`，原 `partialif_densemps` 更名）。旧的按节点构造（合并全部 band 后相乘）已移除。
- **`hybriddynamics` / `hybriddynamics_naive`（Fock 版）统一**：不再按 `Imag/Real/MixedFockLattice` 分别实现（原 `influencefunctional/partialif/` 下 3 个 per-lattice 文件删除），单一签名 `(::FockMPS, ::AbstractFockLattice, ::AbstractCorrelationFunction; trunc)`，contour 分支与 band 全部从 lattice 读出。
- **`retardedinteractdynamics_naive` 统一**：`retardedinteract/` 目录（3 个 per-lattice 文件 + 已删除的 `partialmpo_retardedinteract` 死函数）合并为单文件 `src/electronphonon/retardedinteract.jl`；一个方法覆盖三种 Grassmann lattice，循环展开 `(branch, time, band)`，**不再限制 `bands ∈ (1, 2)`**、不再分 1band/2band 两套实现。
- **FockMPO 删除**：仓库内无任何使用点，`fockmpo.jl`、`export FockMPO` 与 docs 条目一并移除（需要 MPO 时直接用 `FiniteMPSAlgorithms.MPO`）。
- **`randomfockmps` / `iscanonical(::FockMPS)` 委托 FMA 后端**：分别转发到 `FMA.randommps(T, L; d=2, D=D)`（返回 canonical 且已归一化的态，旧行为是未 canonical 化的裸 `randn` 张量）与 `FMA.iscanonical(psi.parent; kwargs...)`。

### lazy-sum（`Vector{<:GrassmannMPS}`）接口移除

- **`integrate`** 不再接受第一个参数为 `Vector{<:GrassmannMPS}`（多个 GMPS 的 lazy 加法）；请先显式相加或分别积分。
- **`gf` / `cached_gf` 及整个 observables 家族**（`contour_ordered_gf`、`greater`/`lesser`、`occupation`、`nn`/`nn2`、`electriccurrent`(_fast)、`cached_*`、`environments`）的 `A` 参数一律收紧为 `A::GrassmannMPS`（原 `Union{GrassmannMPS, Vector}`）。
- **`VectorExpectationCache` 删除**（连同 `_environments` 的 Vector 方法、`_mult_A`/`_insert_n`/`_insert_nn` 的 Vector 方法）。
- 多带求和场景不受影响：继续用 splat（`environments(lattice, adt, Is...)`、`cached_gf(lattice, ..., adt, Is...)`）。

### 验证

- electron-phonon / retarded 定向回归（api hybriddynamics + api electronphonon + electronphonon few-mode/continuous，覆盖三种 contour × bands 1/2/3）全部通过。
- 全量测试（82 个 testset）通过。

## 移除 Grassmann lattice 的辅助 index-0 位点（2026-10-22）

### 动机

`index(lattice, 0)` 不是 contour 离散化产生的真实时间片，而是为了标记配分函数最后的 trace/闭合而额外引入的一对辅助 GV（`a₀, ā₀`）。它在整张 lattice 上只出现在 `boundarycondition!` 中；对这对 GV 做 Berezin 积分后可精确化为两个真实端点 GV 之间的双线性项。因此该辅助位点被整体移除，contour 闭合改为由 `boundarycondition!` 直接连接物理端点。

### 接口变更

- **时间步编号**：所有 lattice 的 `index(lattice, j; ...)` 中 `j` 范围由 `0:k` 改为 **`1:k`**（mixed 的实轴分支为 `1:kt`、虚轴分支为 `1:kτ`）。`index(lattice, 0; ...)` 现在抛 `BoundsError`；`ContourIndex(0, ...)` 的 contour 排序仍按原样抛 `ArgumentError`。
- **`boundarycondition!` / `boundarycondition`**：不再作用于辅助位点，直接连接物理端点（数学上与原实现严格等价）：
  - 虚时：`exp(-ā₁ a_k)`（原为 `exp(-ā₀ a_k) · exp(ā₁ a₀)`）
  - 实轴：`exp(-ā_k⁻ a_k⁺)`（Keldysh 折返点；原为 `exp(-ā₀ a_k⁺) · exp(ā_k⁻ a₀)`）
  - 混合：折返点 `exp(ā_{kt}⁻ a_{kt}⁺)`、`τ=0 ↔ t=0` 的 `exp(ā_{1,τ} a_1⁻)` 不变；闭合项改为 `exp(-ā_1⁺ a_{kτ,τ})`
  - 每个 band 的 `apply! + canonicalize!` 次数减少（虚时/实时 2→1，混合 2→1，但混合仍有 3 个连接项）
- **删除 `boundarycondition_branching`**：虚时版本的分支唯一来源就是辅助对，移除后只剩单支且仓库内无调用点；标量/向量两个方法与文档条目一并删除。
- **`band_boundary(lattice, j)`**：不再接受 `j = 0`；实轴版本中 `j == 0` 的专用分支删除。
- **`length(lattice)` 变化**（每 band）：

  | lattice | 旧 | 新 |
  |---|---|---|
  | `ImagGrassmannLattice` | `2b(k+1)` | `2b·k` |
  | `RealGrassmannLattice` | `4b·k + 2b` | `4b·k` |
  | `MixedGrassmannLattice` | `4b·kt + 2b + 2b(Nτ+1)` | `4b·kt + 2b·kτ` |

  所有 `index` 返回的位点相应整体前移 `2b`（实轴/mixed 的少数 `length(x)` 基准表达式随新长度自动适配）。
- **打印**：site 符号链不再出现 `a₀/ā₀…`；mixed 的虚/实分界由 `2b(kτ+1)` 改为 `2b·kτ`。
- **`timesteps(gmps, ::RealGrassmannLattice)`**：不再扣除辅助块（`L = length(gmps)`）。

### 内部同步调整

- `ttiif` 的 `_fit_to_lattice`（虚时、实轴）删除给辅助块填恒等张量的循环，时间块从位置 1 开始。
- `partialintegrate` 的收缩索引由 `0:lattice.k` 改为 `1:lattice.k`。
- `indexmappings`、`swapbandperm`、`fillband`、`changeordering` 等基于 `indexmappings` 的逻辑自动适配。

### 迁移指南

- 删除用户代码中任何 `index(lattice, 0, ...)` / `ContourIndex(0, ...)` 的使用；trace/闭合已由 `boundarycondition!` 内部处理。
- 若代码硬编码了 `length(lattice)` 或位点编号，按上表减去 `2b`（每 band 的辅助对）。
- 不再使用 `boundarycondition_branching`；需要“分支求和”的场景直接用 `boundarycondition!`。

### 验证

- 基线（`cd125bb`）全量测试通过；仅替换 `boundarycondition!` 后全量测试通过；移除 index 0 后全量测试再次通过（无 Fail/Error）。
- Toulouse 杂质 + 单浴模的端到端对比：虚时 `Z = 1.664543186559…`、实轴 `Z = 1.008242880733…`、混合 `Z = 1.640682444049… − 1.629e-4i`，与基线差 ~1e-14；`gτ / gt / lt` 逐点一致（差异来自有限键维下 canonicalize 次序变化）。

## 新增：GrassmannMPS 的纠缠熵（2026-10-20）

- 新增 `entanglement_entropy(psi::GrassmannMPS; bond, α=1)`、`entanglement_spectrum`、`schmidt_values`，签名与语义对齐 `FiniteMPSAlgorithms` 的同名函数（Rényi 参数 `α`、自然对数、未初始化 Schmidt 值时就地 `canonicalize!`）。
- 实现对 Z2 分块 `DiagonalTensorMap` 按 `fusiontrees` 逐扇区提取对角谱；熵内核复用 `FMA.renyi_entropy`。
- 对由 `hybriddynamics` 构造的影响泛函，该值即相应时间割上的 OSEE；与 QuAPI 侧基于 hybridization function 的高斯态 OSEE 在虚时/实时/混合/BCS 四条路径上一致到 1e-13~1e-10（验证脚本见 `osee/`）。

## 第九批更新（2026-10-07）：算法类型 @kwdef 化与 Defaults 接口函数化

### Defaults 接口函数化

- `Defaults` 模块统一提供默认算法的构造函数：`alg_trunc_if(; D, ϵ, add_back)`（= `truncdimcutoff`）、`alg_trunc_sys(; ϵ)`（= `truncrelerr`）、**`alg_evo()`**（= `WII()`，新增）、`alg_mult(; trunc)`（= `DMRG1(trunc=trunc)`）、`alg_expan(; n, tol, verbosity)`（= `OverDeterminedProny`）。
- **删除**模块外的 4 个默认常数 `DefaultITruncation` / `DefaultKTruncation` / `DefaultExpansionAlg` / `DefaultMultAlg`，全部调用点（src/ext/test/performance/docs）改为对应的 `Defaults.alg_*()`。kwdef 字段默认值与关键字默认值在**构造/调用时**求值，每次返回新实例（原先为共享 const 实例；截断方案等均为 immutable，行为等价）。
- Julia 嵌套模块对外层模块的绑定（含 `using` 导入与本模块定义）**完全不可见**，`Defaults` 内显式 `using Z2Tensors: ...`、`using ExpExp: OverDeterminedProny`、`import ..DMRG1, ..WII`。
- `include("defaults.jl")` 移到 `include("mpo/mpo.jl")` 之后（`import ..WII` 要求 WII 已定义）；mpo 内对 `Defaults.*` 的引用均为关键字默认值（惰性求值），不受影响。

### 压缩/正交化算法类型 @kwdef 化

- `DMRG1`、`DMRG2`、`SVDCompression`、`Orthogonalize` 改为 `@kwdef` 定义（字段默认值内联），**删除**全部位置参数便捷构造器与 `Base.similar`；构造一律使用关键字形式。删除随之无用的 `AllowedInitGuesses`（`initguess` 的合法性改由 `iterativemult` 运行时检查）。
- `Orthogonalize` 的 `trunc` 默认值由 `NoTruncation()` 改为 **`Defaults.alg_trunc_sys()`**。
- 调用形式迁移（全仓库同步）：

  | 旧写法（已失效） | 新写法 |
  |---|---|
  | `SVDCompression(trunc)` | `SVDCompression(trunc=trunc)` |
  | `DMRG1(trunc)` / `DMRG1(trunc, maxiter=30)` | `DMRG1(trunc=trunc)` / `DMRG1(trunc=trunc, maxiter=30)` |
  | `DMRG2(trunc)` | `DMRG2(trunc=trunc)` |
  | `Orthogonalize(SVD(), trunc)` | `Orthogonalize(orth=SVD(), trunc=trunc)` |

  `DMRG1(trunc=truncerr(ϵ))` 等不带键维的方案仍为 `MethodError`（`trunc::TruncationWithD` 类型约束不变）。

### 影响泛函算法类型 @kwdef 化与字段改名

- `XTRGIF`、`ExactTTIIF`、`TDVPIF`、`PartialIF` 全部改为 `@kwdef` 定义，删除手写的关键字构造器。
- 字段与关键字统一改名（含 `influenceoperators` / `influenceoperatorstepper(s)` 的关键字参数及全部调用点、测试、性能脚本与文档）：

  | 旧名 | 新名 |
  |---|---|
  | `algexpan` | `alg_expan` |
  | `algevo` | `alg_evo` |
  | `algmult` | `alg_mult` |

- `TDVPIF` 的 `alg_expan` 由抽象类型字段改为类型参数：**`TDVPIF{T<:ExponentialExpansionAlgorithm}`**；构造器中 δ 的归一化校验移除（`_tdvpif_flow!` 内部本就以 `nsteps = round(Int, 1/δ)` 取整，数值行为不变）。

### 文件清理

- **删除两个 `conversion.jl`**：`src/integration/conversion.jl` 的 `changeordering` / `toadjacentordering` 移至 `src/lattices/changeordering.jl`（include 顺延到 lattice 类型定义之后）；`src/electronphonon/conversion.jl` 的 `reweighting` 系列移至 `src/electronphonon/reweighting.jl`。附带清理 `focktograssmann`、`_abba2aabb` 等注释死代码。
- **删除 `changetrunc`**（仅剩注释引用）：改用 `DMRG1(trunc=..., ...)` 等关键字构造。

### 追加：DMRG2 扫掠复用单站点融合因子

- DMRG2 的两站点扫掠不再使用独立的 `g_ac_prime2` 两站点收缩（其内部构造左右融合因子后，环境更新又经 `updatemultleft`/`updatemultright` 重复构造同样的因子），改为复用单站点融合原语 `get_left_xy`/`get_xy_right`：两点张量按 `twositemps[1,2;3,4] := left_xy[1,2,5,6] * xy_right[6,5,3,4]` 组装，每步的左右因子构造一次、两站点张量与环境更新双用——与 DMRG1 扫掠及 InfiniteGTEMPO 无限 DMRG2 引擎（`_mk_mult_local_AC2`/`_mk_mult_idmrg2_sweep!`）的复用模式完全一致。两站点的补偿差异由收缩网络自然吸收。
- 删除死代码 `g_ac_prime2`；`updatemultleft`/`updatemultright` 保留（TDVPIF 的 hstorage 更新与 InfiniteGTEMPO 的转移映射不动点仍在使用）。
- 验证：与原 `g_ac_prime2` 实现在固定随机种子下逐位一致（z 位点张量、Schmidt 值、scaling、sweep 损失历史全部 `isequal`）。

### 迁移指南

- `DefaultITruncation` → `Defaults.alg_trunc_if()`；`DefaultKTruncation` → `Defaults.alg_trunc_sys()`；`DefaultExpansionAlg` → `Defaults.alg_expan()`；`DefaultMultAlg` → `Defaults.alg_mult()`。
- 压缩/正交化算法一律关键字构造（见上表）；`similar(alg; trunc=...)` → 直接重新构造。
- `XTRGIF(algmult=..., algexpan=...)` → `XTRGIF(alg_mult=..., alg_expan=...)`；`alg.algmult` → `alg.alg_mult`，其余字段同理。

## 第八批更新（2026-09-24 ~ 2026-10-07）：FockMPS/FockMPO 封装 FiniteMPSAlgorithms 与 changebond!/mult!/truncate! 接口统一

### FockMPS / FockMPO 封装 FiniteMPSAlgorithms 后端

- `FockMPS` 的存储改为内嵌 FiniteMPSAlgorithms 的 `CanonicalMPS`（字段 `.parent`；`.data`/`.s`/`.scaling` 经 `getproperty` 委托到 payload）；新增 `FockMPO` 内嵌 FMA 的 `MPO`。位点张量、Schmidt 值与逐位点 scaling 全部由 payload 携带，FMA 算法（`canonicalize!`、`hadamard`、`dot`/`norm`、`permute!`、`sum`、`apply!` 等）原地驱动；稠密张量原语（`tie`/`isometry`/…）同样来自 FMA。
- `FiniteMPSAlgorithms` 改为限定引用 `FMA.xxx`（避免与 Grassmann 路径上 Z2Tensors 的导出名冲突），并加入 path 依赖。
- 删除未使用的 `DenseMPO`/`PartialDenseMPO` 类型及其全部线代实现；`ExpNTerm` 直接经 `apply!` 作用。

### `increase_bond!` → `changebond!`

- 按 FMA 惯例更名；`FockMPS` 与 `GrassmannMPS` 均提供，`noise` 关键字默认为 0。
- `GrassmannMPS` 的 `changebond!` 把键维调整到可行性上限 `min(D, ∏_{j≤i} d_j, ∏_{j>i} d_j)`（两端受边界约束，不再被 D 单独钳制）；需要缩键时先做无截断 `canonicalize!`（避免在非规范链上盲目切 leading 分量）；最后以无截断 `rightorth!` 定稿为右正则形式（正交中心在 site 1，整体范数由 `scaling` 携带），输出可直接作为迭代算法的初始猜测，无需再 canonicalize/rightorth。

### TDVPIF：去尾部规范化和 Schmidt 值重置

- 删除 TDVPIF 虚时流结束处冗余的最终 `canonicalize!(z)`；改为 `unset_svectors!(z)` 重置流后不再匹配的 Schmidt 值。
- 删除 `DMRG1` 与 `TDVPIF` 的 `callback` 字段。

### `mult`/`mult!`/`truncate!` 返回值与默认值约定（对齐 FMA）

- GrassmannMPS 的 `canonicalize!` 与 `truncate!`（新增，`truncate!(psi; trunc) = _canonicalize!(psi; alg=Orthogonalize(SVD(), trunc; normalize=false))`）的 `trunc` 默认为 `DefaultKTruncation`。
- 返回值约定（全部调用点已同步更新，取第一个输出）：
  - `mult(x, y;)`、`mult!(x, y;)`、`mult(x, y, SVDCompression)` → `(x, maxerr)`（`maxerr` 为定稿 SVD 扫描的最大截断相对误差）
  - `iterativemult(x, y, alg)`、`mult(x, y, alg::DMRGAlgorithm)`、`mult!(x, y, alg::DMRGAlgorithm)` → `(x, info)`，`info` 为 FMA 提供的 `ALSConvergenceInfo`（字段 `niter`/`converged`/`losses`/`itererr`/`residual`）；收敛判据为相邻两轮 sweep 损失的相对变化 `|lₙ-lₙ₋₁|/|lₙ₋₁| < alg.tol`
- `mult!(x, y, alg::DMRGAlgorithm)` 现在把结果真正写回 x（`copy!` 位点张量与 Schmidt 值并 `setscaling!`），不再只返回新对象的别名。
- 导出 `leftorth!`、`rightorth!`；删除主模块的 `using Z2Tensors: NoTruncation` / `export NoTruncation`（Z2Tensors 已导出）。

### 键维/物理空间访问器更名

| 旧名 | 新名 |
|---|---|
| `bond_dimension` / `bond_dimensions` | `bonddim` / `bonddims` |
| `physical_space` / `physical_spaces` | `physpace` / `physpaces` |
| `ophysical_space` / `ophysical_spaces` | `ophyspace` / `ophyspaces` |
| `iphysical_space` / `iphysical_spaces` | `iphyspace` / `iphyspaces` |

## 测试容差调整：实轴 `bare + bulkconnection == sysdynamics`（2026-09-17）

`test/api/sysdynamics.jl` 中 "impurity hamiltonians & steppers" 的实轴断言由绝对阈值 `< 1.0e-8` 改为相对阈值 `< distance/norm(K) < 1.0e-6`，与虚轴对应断言一致。原因：Z2Tensors 近期更新（norm/fusiontensor/⊗ 修复）改变了 contraction 内部的浮点求和顺序，使 bare+bulkconnection 与 sysdynamics 两条独立构造路径从"逐位一致"变为存在 ~1e-6 的浮点水平差异（相对误差约 2e-8），原绝对阈值对实现顺序过敏。

## 删除 `src/electronphonon/tensorops.jl`（2026-09-17）

Z2Tensors 已提供等价的稠密张量工具，GTEMPO 不再自带重复实现：

- `permute(::AbstractArray, ...)`、`tie`、`_group_extent`、`_truncate!`、`tsvd!`（矩阵与 (left, right) 张量版）直接使用 Z2Tensors 导出的同名函数；`tsvd!` 的输入数组自身作 workspace（不再传 workspace 参数）。
- `_eye` → `TK.isometry`；`tqr!(a, left, right)` / `tlq!(a, left, right)` → `TK.leftorth!` / `TK.rightorth!`（默认算法 QRpos/LQpos，保证 R/L 因子对角非负）。
- `stable_svd!` 删除（Z2Tensors 的 `tsvd!` 内部已用 MatrixAlgebraKit 的 SafeDivideAndConquer 驱动）。
- `Base.kron` 对同秩稠密数组的扩展同样由 Z2Tensors 提供（`auxiliary/tensoroperations.jl`），本地实现删除。

## 删除 `BranchLocalLayout` 与相关 ordering（2026-09-17）

- 删除 `BranchLocalLayout` 布局类型及全部三个 band-local ordering：`A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2`（real）、`A2Ā2B2B̄2A1Ā1B1B̄1a1ā1b1b̄1a2ā2b2b̄2`（real）、`A1B1B̄1Ā1_A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2`（mixed），连同各自的 `index` 方法与导出。
- `sysdynamics_fast` 的实轴 canonical ordering 统一为 `A1B1ā1b̄1Ā1B̄1a1b1`（原对非 TimeLocal lattice 返回 band-local ordering）。
- 实轴 `toadjacentordering` 统一转换到 `A1Ā1a1ā1B1B̄1b1b̄1`；删除 mixed band-local ordering 的专用方法。
- 删除 `boundarycondition_branching` 的实轴方法（其唯一实现依赖 `BranchLocalLayout`；虚轴方法保留）。

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
