# Td 模型 order 2 半步锚定修复计划（半区间中点冻结）

> **计划成文（2026-11-20，待实施）。** 缺陷：Weyl 2 阶的跨顶点参数化 S·W·S 对时变杂质
> 模型（`AbstractTdImpurityHamiltonian`）每胞两半步共用**胞左端点锚定**的同一矩阵
> ——第二半步的物理区间错位、内部合并步非连续区间，全局精度**名二阶实一阶**（O(δt)）。
> 修复在现有接口下 ~10 行可完成；**静态模型零回归由构造性保证**。

- 日期：2026-11-20（sysdynamics 半步实现审查中发现，用户提出修复计划）
- 基线：GTEMPO HEAD `a6672a84` **+ 合并计划（`fockmpo_dress_selfread_plan.md`）实施后的
  fockmpo.jl**（穿戴五步 + δ 自读已落地——本修复与其正交：只动传播子装配
  `_fock_propagator_at`/新增 `_fock_straddle_propagators` + `fock_sysdynamics` 的 mats
  构造两行，不碰穿戴/自读逻辑；实施前 rebase 到当期 HEAD）
- 波及面 = `src/sysdynamics/fockmpo.jl` 单文件（修复 + 测试），**无接口变更、无 src 新增
  导出符号**
- 关联：`fockmpo_dress_selfread_plan.md`（同文件叠加，正交）；`exactif_port_plan.md`
  （金标准 G1/G2 均静态模型，不受影响——边界见不做清单）

## 1. 缺陷定谳

**结构**：order 2 每格点 `S·W·S`（`fockmpo.jl:123`），S = `_fock_propagator_at`
（`fockmpo.jl:88-93`）的 `dt/2` 传播子；时不变模型相邻格点 S·S 合并成全步、两端剩首/末
半步 = Strang 序列（精确）。

**Td 缺陷**：`_fock_propagator_at` 对 order 2 只传胞物理左端点 t（+ 支 t = (j−1)·δt、
− 支 t = (Nt−j)·δt），`fock_sysdynamics` 把同一个 S 用了两次 ⟹

1. 第二半步应覆盖物理区间 [m_j, t_j]（m_j = 胞中点），实际用了 [t_l, m_j] 的矩阵——
   冻结 H 的锚点错（应锚 m_j 而用了 t_l）；
2. 内部"合并全步"= S(t_l)·S(t_l) 型同矩阵自乘——不是任何连续区间的演化。

**精度阶**：Td 传播子本身 = 冻结 H 近似（`fockmatrix.jl:93`：`exp(−i·h(t)·δt)`，t = 冻结
点）。锚定规则决定精度阶：区间**左端点**冻结 = 每门指数误差 O(δt²) ⟹ 2N 门累积
O(δt)（一阶）；区间**中点**冻结 = 每门 O(δt³) ⟹ 全局 O(δt²)（二阶）。现行缺陷态 =
左端点锚 ⟹ **Td order 2 实为一阶**（静态模型 h 与 t 无关，两半同矩阵，无缺陷）。

**Design A 否决记录**（先想到的"区间起点锚"方案不充分）：把两半步分别锚到 t_l 和
m_j（各自区间左端点）——修好了区间覆盖错位，但仍是左端点冻结 ⟹ 全局仍 O(δt)。
**必须用半区间中点（1/4 点）冻结**。

**Design C 否决记录（2026-11-20 用户问："只改 fock_propagator 内部为中点冻结、其他
不动"）**：**数值上修不好**——bug 的根是"一个矩阵当两个区间用"，S·W·S 两半调用
prop(δt/2, t_l) 传同一个 t；内部改中点冻结（exp(−ih(t+dt/2)·δt)）后两半都冻结在
t_l + δt/4：第一半恰好对、第二半仍错（[m_j, t_j] 应冻结 t_l + 3δt/4）——右半每胞
O(δt²) 误差 ⟹ 全局仍 O(δt)。**只要 mats 构造只传一个锚，内部怎么改都修不好；第二个
锚点必须在调用处出现**。变体"内部中点冻结 + 调用处传双锚"数值上 ≡ Design B，但代价
更高：fock_propagator 是 `AbstractTdImpurityHamiltonian` 的文档化契约（"t the left
endpoint of the physical time interval of the step"，`impuritymodels.jl:28-30`）——内部
中点化 = 契约语义变更（自定义 Td 模型照文档实现 exp(−ih(t)δt) 会静默不一致）、order 1
Td 数值漂移（斜率仍 ≈1——Lie 劈裂主导——但系数变、模型层单元测试破）、ED 参照等
其他消费者静默继承新语义。故**契约不动、锚点显式在调用处**（Design B）。

## 2. 修复设计（现有接口，~10 行）

**语义定谳（修复目标）**：order 2 的杂质侧 = 2N 个连续半区间门（全局 [t_0, t_N] 按
δt/2 切分，+ 支物理正向、− 支物理反向），每门冻结于**该半区间的中点**（中点规则，
每门 O(δt³)）⟹ 全局 O(δt²)；IF 源点保持物理胞中点 m_j；胞结构不变（S·W·S 跨顶点）。

```julia
# fockmpo.jl：order 2 专用（order 1 的 _fock_propagator_at 原样保留）
# 半步对：左半 [t_l, m] 冻结于 t_l + δt/4、右半 [m, t_r] 冻结于 t_l + 3δt/4
# （各半区间的中点——中点规则）。Const 分派忽略 t ⟹ 两半同矩阵 ⟹ 与现行 S·W·S
# 逐位一致（静态零回归构造性保证）。− 支轮廓物理反向 ⟹ 乘序镜像（与一阶
# S·W / W·S 的不对称同源）。
function _fock_straddle_propagators(lattice, model, b::Symbol, j::Int)
    dt = b === :τ ? lattice.δτ : lattice.δt
    t = b === :+ ? (j - 1) * dt : (b === :- ? (lattice.Nt - j) * dt : nothing)
    tm = t === nothing ? nothing : t + dt / 2
    S_lower = Matrix(_model_propagator(model, b, dt / 2, t === nothing ? nothing : t + dt / 4))
    S_upper = Matrix(_model_propagator(model, b, dt / 2, t === nothing ? nothing : t + 3dt / 4))
    return S_lower, S_upper
end

# fock_sysdynamics order 2 分支替换：
S_lower, S_upper = _fock_straddle_propagators(lattice, model, b, j)
mats = b === :- ? [S_lower * W * S_upper for W in weights] :
                 [S_upper * W * S_lower for W in weights]
```

**伴随的 doc-only 修订**：`fock_propagator(::Td, branch, dt, t)` 的接口 docstring
"t the left endpoint of the physical time interval of the step"（`impuritymodels.jl:29`）
放宽为"**t = H 的冻结时刻**（实现 = exp(−i·h(t)·δt)；锚定规则由调用方按离散化阶决定
——区间左端 = 一阶 Lie、半区间中点 = 二阶中点规则）"。纯措辞修订行为不变（现有全部
调用传的就是区间左端点，语义兼容）；Design B 传 1/4 点不是任何步的左端点，旧措辞下
构成契约误用，故必须同步修订。

**乘序镜像推导**（防实施走样，钉死）：矩阵积 = 右因子先作用。+ 支轮廓物理正向：
先下后半 [t_l, m]（S_lower 右位）、源在 m、再上后半 [m, t_r]（S_upper 左位）⟹
`S_upper·W·S_lower`；− 支轮廓物理反向：先上后半（S_upper 右位）、源在 m、再下后半
（S_lower 左位）⟹ `S_lower·W·S_upper`。Const 下两者都退化为现行 `S·W·S`。

**关键性质**：

- **静态零回归是构造性的**：`_model_propagator(::ConstImpurityHamiltonian, …, t)` 忽略
  t（`impuritymodels.jl:44-45`）⟹ S_lower ≡ S_upper ⟹ mats 与现行逐位一致——203 闸门/
  observables golden/九组二阶斜率 CSV 无需数值检查即保持（仍跑一遍确认）；
- **内部合并步恢复正确**：相邻胞交界 = S_upper_j·S_lower_{j+1}（时间正向），两半区间
  [m_j, t_j]、[t_j, m_{j+1}] 相邻接、各自中点冻结 ⟹ 中-中全步 [m_j, m_{j+1}] 的二阶
  近似（每半 O(δt³)）；
- **τ 支不变**：t = nothing ⟹ 两半同锚（Td 在 imag 无语义，hτ 静态）；
- **测量层自动继承**：inserted-K 重建走 sysdynamics ⟹ fix 在 fock_sysdynamics 内一次生效；
- **order 1 不动**：Td 的 order 1 = 左端点冻结全步 + 边界源 = 标准 Lie 一阶，本来就该
  一阶，无缺陷。

## 3. 验收

- **T0 静态零回归**：Const order 2 全套——203 接线闸门 + observables golden（U1GTEMPO2
  机器精度）+ 九组二阶斜率 CSV **逐字节**；order 1 全套不动。
- **T1 Td 斜率判别**（诊断 + 修复证明）：`TdImpurityHamiltonian` 线性驱动（h(t) = h₀ +
  t·h₁）× 单模 bath，n̄(t) 对 δt 收敛斜率：**修复前 order 2 ≈ 1**（bug 的数值定谳——
  实施前先跑，钉住缺陷）；**修复后 order 2 ≈ 2**；order 1 ≈ 1 前后不变。
- **T2 Td 机器精度锚**：同一 Td 模型 + 单模 bath——**ED 直接 Strang 门序列**（同 1/4 点
  锚定的半步对 + 分支镜像）vs `sysdynamics(order=2)` + IF 装配 n̄，相同 (N, δt) 调度下
  机器精度对齐（参照机 = 测试侧 helper，src 零新增；IFD `test/exact/_anchor_utils.jl`
  门序列机同款思路）。
- 边界回归：混合轮廓（τ 胞不变 + ± 胞修复后）、imag 轮廓（应逐位不变）。

## 4. 不做清单

- **exactif 严格核机器的 Td 面**：严格核 kicks 的 gap 相位取 εd(t_g) 与 Td 联合构造
  ——金标准 G1/G2 均静态模型，本计划不碰；未来 Td 严格核对齐需求出现时归
  `exactif_port_plan.md` 契约核对。
- order 1 锚定不动（标准 Lie 一阶，无缺陷）。
- 不新增接口/导出符号（修复 = fockmpo.jl 文件内部 ~10 行 + 测试）。
- Td 传播子本身的精度升级（冻结 H → 仿射精确指数）不在本计划——中点冻结已达
  二阶；若未来要更高阶再立项。

## 5. 工作量与顺序

~半天：缺陷前置诊断（T1 修复前斜率钉住 ≈1）→ 修复（~10 行）→ T0 静态零回归 →
T1 斜率 ≈2 + T2 机器精度锚 → 收尾。顺序硬约束：**先跑 T1 修复前斜率钉住缺陷，
再动源码**；T0 不过 = 回退。
