# PolarIF：fullfactor × PolarExp 装配算法计划

> **状态：立项待实施**（2026-11-19 成文）。用户提案：在 `src/influencefunctional/`
> 增加 `PolarIF` 算法，对接 GrassmannMPS 的 `PolarExp` 后端；实现基于现有
> `fullfactor`，适用范围与 `PartialIF` 一致。

- 日期：2026-11-19
- 基线：GTEMPO 分支 `if-expmult-backend` HEAD `0f64cf78` **+ 工作区并行 WIP**
  （fullfactor 泛化 / AnomalousPairing 紧凑站点表等未提交修改——本计划的行号与
  代码形状以工作区状态为准，实施前与用户确认 rebase 基线，见 R3）
- GM 前置：GrassmannMPS 分支 `expdecayexp-same-orbital-fix`（`2d16f8a` 已推送；
  PolarExp SVD 秩一语义 = GM main `6355850` 已入库，与 ExactTTI 同一前置集）
- 关联计划：`if_expmult_backend_plan.md` W3（fullfactor 因子面——PolarIF 是该面
  的第二个消费算法，无文本冲突）；GM polarexp-svd-rankone（后端，已实施完结）

## 0. 已裁决 / 待裁

**已定**（用户提案 2026-11-19）：

- 新增 `PolarIF <: InfluenceFunctionalAlgorithm`，装配 =
  `expmult!(gmps, fullfactor(kernel, lattice), GM.PolarExp(...))`；
- 适用范围与 PartialIF 一致：任意轮廓（imag/real/mixed）、normal / 跨带 /
  BCS 反常核、四对称性面；单链恒就地装配契约不变；
- 复用现有 `fullfactor`，**不新增任何核/因子构造器**。

**待裁**（实施前用户确认）：

- **A1【rtol 缺省值】**：推荐 `0.0`（全数值秩——与 PartialIF 同为
  「精确指数 + 仅乘法截断误差」语义类；ExactTTI real/imag 通道同款先例）。
  备选 = GM 缺省 `1e-12`（剔除 LAPACK 噪声尾，但理论上不精确）。下文按 `0.0` 成文。
- **A2【BCS 守卫文案】**：`bcs_hybriddynamics` 守卫放宽为
  `Union{PartialIF,PolarIF}`，报错文案随之更新（"supports PartialIF/PolarIF
  on the GM stack (D7)"）。

## 1. 语义定位（三算法对照）

| | PartialIF | **PolarIF（本计划）** | ExactTTI |
|---|---|---|---|
| GM 指数算法 | `PartialExp` 朴素幂零级数 | `PolarExp` SVD 秩一因子化 | `PolarExp`（tti 因子逐个） |
| 因子来源 | `fullfactor` 直读 | `fullfactor` 直读（同一函数） | `tti_imag/real/cross_factors` 三构造器 |
| 精确类 | 精确 exp(O) | `rtol=0` ⟹ 精确 exp(O)（机器精度） | 拟合重构档（单极点核机器精度） |
| 误差源 | 乘法压缩 `trunc` | SVD 秩截 `rtol` + 乘法压缩 `trunc` | 拟合残差 + 压缩 |
| kernel 范围 | 任意（normal/跨带/BCS） | **任意（同 PartialIF）** | 单带 normal (ā,a) 键（D7） |
| 标量型 | 链×核 promote（实链实谱核全程 Float64） | 同左 | 恒 ComplexF64（拟合层），实链入口复化 |
| 多带返回 | 恒单链 GMPS | 恒单链 GMPS | 多带核逐带 Tuple |

关键事实：

- **F1【同一算符的两种精确表示】** `e^{Σ ā Δ a}` 的 PartialExp 表示（幂零级数）
  与 PolarExp 表示（对角剥离/味扇区折叠 → 秩一分量 Xₖ²=0 ⟹ e^{Xₖ}=1+Xₖ 精确、
  偶元超交换 ⟹ Π(1+Xₖ) 精确，GM `exponentials.jl:631-657`）数学上恒等；
  两者消费**同一个 raw fullfactor** ⟹ 任意轮廓上互换对拍合法且达机器精度档
  （实时轴不像 ExactTTI 需要拟合，无残差档）。
- **F2【GM 后端双因子路线现成】** `_polar_factors(::Pairing)`（对角剥离 +
  余阵 SVD，`exponentials.jl:684`）与 `_polar_factors(::AnomalousPairing)`
  （味扇区折叠逐扇区 SVD，`exponentials.jl:740`）均有 GM 全套测试背书
  （GM 7561 绿）。AnomalousPairing 的 ηη̄ 折叠式
  `D[i,j] = C[2i−1,2j] − C[2j,2i−1]` 与 GTEMPO fullfactor BCS 的有向单取向
  存储（partialif.jl 头注释）**逐项对齐**——BCS 反常核零适配。
- **F3【消费面零新因子】** `fullfactor` 现有两返回路线（工作区 WIP 后）：
  normal 同带键 → `GM.Pairing`（矩形站点表 `C[st(ā), st(a)]`）；跨带/BCS
  反常键 → `GM.AnomalousPairing`（紧凑站点表 2M×2M）。PolarIF 原样消费。
- **F4【同味扇区守卫已被上游钉死】** PolarExp 的 ηη/η̄η̄ 扇区非零要求奇自对偶
  站形（Z2/Z2Z2，`exponentials.jl:748-752`）——U1/U1U1 的 BCS 已被
  `bcs_symmetry_guard`（D7，仅 :z2）在 kernel 构造端拒绝，PolarIF 无需重复
  守卫；normal 核无同味扇区，U1/U1U1 normal 走纯 ηη̄ Pairing 路线合法
  （if_gates ② 已有 U1 互换算例可复用）。
- **F5【包内 PolarExp 消费先例】** ExactTTI 装配链已用
  `GM.PolarExp(rtol=0.0/alg.cross_rtol, alg_mult=SVDCompression(trunc))`
  双用法（`exacttti.jl` `_exacttti_state`）——PolarIF 的参数形状与之一致，
  仅 rtol 语义通道不同（PolarIF 把 rtol 暴露为算法字段）。

## 2. 设计（源码改动 5 处，约 35 行）

### 2.1 算法对象——`src/influencefunctional/influencefunctional.jl`

算法 struct 与 `PartialIF`（:28）/`ExactTTI`（:48）同层，插在 `PartialIF`
之后：

```julia
"""
    PolarIF(; rtol=0.0, trunc=Defaults.alg_trunc_if())

PolarIF：fullfactor 的 SVD 秩一精确指数装配（expmult!×PolarExp）——同一
IF 算符 e^{Σ ā Δ a} 的第二种精确因子化（Pairing = 对角剥离 + 余阵 SVD；
AnomalousPairing = 味扇区折叠逐扇区 SVD）。适用范围与 PartialIF 一致
（任意轮廓 / kernel / 对称性；BCS 反常腿味在核键上，同一装配函数）。
`rtol` = PolarExp SVD 秩截断（相对 Frobenius 尾部² 口径；0 = 全数值秩，
机器精度精确；AnomalousPairing 跨扇区全局归一）；
`trunc` = 乘法压缩截断策略（Defaults.alg_trunc_if() = truncdimcutoff 缺省）。
"""
@kwdef struct PolarIF{Tr<:TruncationStrategy} <: InfluenceFunctionalAlgorithm
    rtol::Float64 = 0.0
    trunc::Tr = Defaults.alg_trunc_if()
end
```

要点：rtol 用具体 `Float64` 字段（同 ExactTTI.cross_rtol 风格），不进类型
参数——截秩是精度旋钮非结构信息；`trunc` 字段名与 PartialIF/ExactTTI 一致，
`hybriddynamics!` 入口的 `diagnostics[:truncation] = alg.trunc`
（influencefunctional.jl:100）自动成立。

### 2.2 装配函数——`src/influencefunctional/partialif.jl`

紧随 `_partialif_state`（:102-107）：

```julia
"PolarIF 组装（原地）：expmult!(gmps, fullfactor, PolarExp)——对角剥离/味扇区
折叠后逐秩一分量 mult!（精确因子化），乘法压缩截断 = alg.trunc；
返回 (state, [truncerr])（与 _partialif_state 同形）"
function _polarif_state(gmps, lattice::AbstractGrassmannLattice, kernel::HybridizationKernel,
                        alg::PolarIF)
    exp = GM.PolarExp(rtol=alg.rtol, alg_mult=SVDCompression(trunc=alg.trunc))
    state, err = GM.expmult!(gmps, _fullfactor(kernel, lattice), exp)
    return state, [err]
end
```

### 2.3 分发——`src/influencefunctional/influencefunctional.jl` `hybriddynamics!`

在 `alg isa PartialIF` 分支（:101-102）后加：

```julia
    elseif alg isa PolarIF
        state, truncerrs = _polarif_state(gmps, lattice, kernel, alg)
```

`else throw(ArgumentError(...))` 兜底不动。docstring（:80-94）算法分支描述
补一句：「PolarIF = 恒就地单链（同 PartialIF 契约，SVD 秩一因子化乘入）」。

### 2.4 BCS 便捷入口守卫——`src/influencefunctional/partialif.jl` `bcs_hybriddynamics`

:121 守卫放宽 + 文案更新（A2）：

```julia
    alg isa Union{PartialIF,PolarIF} || throw(ArgumentError(
        "BCS IF supports only PartialIF/PolarIF on the GM stack (D7), got $(typeof(alg))"))
```

docstring 同步（"装配走 fullfactor 轨道路径"叙述不变，算法枚举补 PolarIF）。

### 2.5 导出——`src/GTEMPO.jl:115`

```julia
export InfluenceFunctionalAlgorithm, PartialIF, PolarIF, ExactTTI
```

### 2.6 文件头注释

`influencefunctional.jl` 头注释（:1-16）partialif.jl 条目补 PolarIF 一行；
`changes.md` 落地时追加条目（见 §5）。

## 3. 明确不改动的面

- **标量型**：无需 `_complex_mps`——那是 ExactTTI 特有（拟合层恒 ComplexF64）；
  PolarIF 与 PartialIF 同路径，GM `expmult!` 内部
  `T = promote_type(scalartype(gmps), scalartype(factor))`，实链 × 实谱虚时核
  全程 Float64。
- **单链/多带契约**：fullfactor 单因子覆盖 `kernel.bands` 全部变量 ⟹ 多带
  带对角核与跨带核同恒单链就地装配（PartialIF 契约继承，无逐带元组）。
- **ExactTTI 守卫**：`_is_normal_conj` / 单带限制与 PolarIF 无交集，不动。
- **空核/全零块**：GM PolarExp 空因子表 → 恒等 no-op（truncerr=0，
  `exponentials.jl:649`）——无需 GTEMPO 侧特判。
- **docs/tutorials 不动**：教程头「顶点位移装配」陈旧文案是独立事项
  （normal 实时 W2 起已 raw），不混入本计划。
- **GM 侧零改动**（PolarExp/Pairing/AnomalousPairing 全现成）⟹ GM 套件不重跑。

## 4. 测试计划（验收闸门）

### T1 `test/api/if_gates.jl`：PolarIF 互换闸门（新 testset "IF 装配闸门：PolarIF"）

复用现有算例配置（`gd` = todense 逐元素最大差，金标准纪律）：

| # | 算例 | 配置（复用现有 ①②③） | 断言 |
|---|---|---|---|
| 1 | 虚时单模 μ=0/μ≠0 | `discretefermionicbath`/`DiracDelta`，N=4, δτ=0.25 | `gd(PolarIF, PartialIF) < 1e-12` |
| 2 | 混合连续谱 | Nt=3, Nτ=4 半圆谱 | `< 1e-12` |
| 3 | U1 normal 虚时 | `U1Symmetry()` 单模 | `< 1e-12`（F4 纯 ηη̄ 路线） |
| 4 | BCS 虚时反常 | ③ 的 latb/bathb（Δ=0.4, bands=2, Z2） | `< 1e-12`（AnomalousPairing 味扇区折叠路径覆盖） |
| 5 | 实时 DiracDelta 单模 | N=3, δt=0.15（③ 机器精度算例） | `< 1e-12`（实时 raw 直读——F1 无拟合档） |
| 6 | 跨带随机核 dense oracle | 「泛化 kernel」testset 的 `rand_kernel`，alg 换 `PolarIF(trunc=...D=512)` | vs `DT.exp(todense(fullfactor))` `< 1e-10`（虚时 4 算例 + 实时 1 算例 + 恒单链契约 2 断言） |

### T2 `test/api/if_gates.jl`（或 T1 内）：rtol 截秩 smoke

同算例 1 配 `PolarIF(rtol=0.1)`：结果有限、与 `rtol=0` 版有可观测差、
`diagnostics[:truncerrs]` 有限非负——验证截秩通道活着（数值锚实测回填，
不预设紧界）。

### T3 `test/bcs/runtests.jl`：BCS 实时 ED（新 testset）

复刻 `docs/tutorials/firstorderif/singlemodebath/bcs_realtime.jl` 模型缩
测试级格点：Anderson U=1、Δ=0.6、β=1.0、T=0.25、Nt=5（changes.md 记载该配置
ED 偏差 7.2e-5/1.1e-4 纯离散化档）。PartialIF 与 PolarIF 各跑 gt/lt vs ED
（容差 1e-3 档），并 PolarIF≡PartialIF 互拍（容差实测回填，预期 ≪ ED 容差）。

### T4 `test/normalbath/fewmodesbath/realtime.jl`：离散化阶语义复验

现有 `for order in (1,2)` ED 循环（容差 1e-4）加 PolarIF 一组——确认乘法序
改变不影响一/二阶离散化阶语义。需给 `test/normalbath/util.jl` `fermionic_setup`
（:45）加 `alg::InfluenceFunctionalAlgorithm=PartialIF()` kwarg 透传
（测试辅助面小改，缺省行为不变）。

### T5 全套回归（运维锚）

- `julia --project=. test/runtests.jl` 直接跑（Pkg.test 沙箱炸教训）；
  串行、跑前查用户并行 julia 进程、后台勿套 nohup&；
- 套件终点 = 末尾 Part4 testset 全绿即完整（以 WIP 后实际 testset 数为准）；
- GM 套件不重跑（零改动）。

## 5. 实施顺序与提交切分

1. **feat(if): PolarIF 装配算法（fullfactor × PolarExp）**——§2.1-2.6 全部
   源码 + T1/T2 闸门测试（测试是包的验收闸门，随功能同笔）；
2. **test: PolarIF BCS 实时 ED + normalbath 阶语义复验**——T3/T4；
3. `changes.md` 追加条目（存储约定/算法清单/测试增量），§7 回填实施状态。

## 6. 风险与边界

- **R1【rtol 跨扇区归一】**：AnomalousPairing 截秩阈值全局归一（尾部² 对
  三扇区 Σ|σ|² 总和比较，GM `exponentials.jl:735-737`）——BCS 配 rtol>0 时
  弱扇区可能整扇区截掉；`rtol=0`（推荐缺省）无此问题，docstring 已注明。
- **R2【性能不预设】**：PolarExp 逐秩一分量 mult!（步数 = 有效秩 + 对角剥离
  数）vs PartialExp 幂零级数逐行——步数与单步键维分布不同，孰快孰慢交给
  后续 performance/ 基准，不在本计划验收内。
- **R3【并行 WIP 基线】**：工作区有未提交修改（partialif.jl/influencefunctional.jl/
  correlationfunction.jl + 测试），与本计划同文件——实施前与用户确认
  rebase/提交基线，避免半路冲突。
- **R4【接口面最小】**：不加 Defaults 缺省、不改 `hybriddynamics` 缺省算法
  （仍 PartialIF）、不动测量层——PolarIF 纯增量算法对象，消费端
  `hybriddynamics(lat, kernel, PolarIF())` 即用。

## 7. 实施状态

**✅ 已落地（2026-11-19，两提交：`54c01a7a` feat + test 批）**。

- §2.1–2.6 全部落地（A1 = rtol 0.0 缺省、A2 = 守卫放宽，按计划推荐）；R3 的
  并行 WIP 基线 = 已随第七~九批提交（HEAD 18ce8f65），无冲突。
- T1/T2：PolarIF 对拍**并入 W3/W4 闸门既有算例**（复用 PartialIF 装配态，免
  重复装配）+ 独立 testset 只留 W3/W4 未覆盖面（BCS 反常味扇区 / dense
  oracle / 恒单链 / rtol smoke）。容差全部实测回填：缩格后虚时/U1/实时/混合
  互换均达机器精度档（1.2e-15–5.3e-15，压缩无损），BCS N=3 1.2e-14，
  oracle ≤ 3.2e-12。
- T3：**发现已知缺陷**（见 changes.md「PolarIF 装配算法批次」）——BCS 实时
  核块顶点移位口径与标准网格错位，ED 锚以 @test_broken 记录待裁决；互拍锚
  （1.5e-8）有效。ED 参照口径已验证（数守恒退耦 H0 = 杂质热态 ⊗ 浴热态，
  初态投影法：浴占据基底 +50 能量移位）。
- T4：fermionic_setup 加 alg kwarg；order=2 双算法复验（一阶 PolarIF 锚由
  T3 承接，免重复锚点）。
- 减时裁决（用户指示「不缩覆盖去重复」）：闸门对拍复用装配态 + 缩格 +
  oracle 裁剪，两闸门 testset 5.5 min → 0.5 min。
- 偏差记录：T1 断言按计划预设 1e-12 不可达（乘法序压缩路径差），改逐例
  实测锚（机器精度档～1e-8 视链宽）；imag 格 N=3 为热配对审计已知坏位型，
  用 N=4。
