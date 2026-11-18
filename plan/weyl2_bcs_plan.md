# Weyl 二阶 × BCS bath 扩展计划

- 日期：2026-11-17
- 基线：GTEMPO HEAD `f3ee05c4` 系工作区
- 关联：`sysdynamics_fock_mpo_plan.md`（K 侧重构——本计划的 Weyl 标量结论须同步其 §3.1）；`if_expmult_backend_plan.md`（IF 装配，独立）；D7 裁决（BCS 仅 :z2）继续有效

## 0. 结论与裁决

**可行性：可行，且大部分机器已就位**（结构性论证 §1 F1-F3 + 冒烟实证 F4）。本计划 = 实验先行的验证与形式化，不是新机制开发——Weyl-2 的 K 侧机器对 BCS 零改动复用；唯一开放推导点 = Weyl 标量的 BCS 广义化（Phase 1 实验裁决）。

## 1. 可行性定谳

**F1【K 链 bath 无关】** `sysdynamics`/`sysdynamics_weyl2order` 的 K = f(lattice, model)——BCS bath 只进 IF（`bcs_partial_if`）。Weyl-2 的全部 K 侧机器（S·M·S 对称半步、q=1/2 measure、W₄）对 BCS **零新张量代码**。

**F2【BCS IF 与 normal 同机器、同变量字母表】** `bcs_partial_if`（influence.jl:652-681）= `partial_row` + `_site_charges`——与 normal PartialIF 完全同一机器；反常结构纯粹体现在 **BCSCorrelationFunction 的系数支撑位**（`_bcs_legs`：行 = (band1,ā)/(band2,ā)/(band2,a)，列 = (band1,a)/(band2,a)/(band2,ā)——āā/aa 同味对由 corr 的 cc 块供值）。:z2 parity 守卫内建。**aa/āā 配置在变量字母表内天然存在**（:z2 下 parity 偶、合法）。

**F3【一阶二阶同一套指标】** 离散化锚点（weyl-secondorder-discretization 记忆/代码实证）：Weyl 只改系数不改指标集 ⇒ BCS IF（一阶已验证 ED）与 Weyl-2 K **无新指标协调问题**。

**F4【冒烟实证（2026-11-17 实跑）】** bands=2、N=4、δt=0.1、real、AndersonIM(U=1, ϵ_d=−0.5)、`bcsbath(fermionicbath(DiracDelta(1.0, 0.5), β=1), Δ=0.3)`：
```
order1  站数=16  Z(K,IF) = 0.9460810747
order2  站数=16  Z(K,IF) = 0.9441024643     # 差 ~2e-3 = O(δt) 档，符合阶差预期
```
**现行代码无任何守卫阻止该组合**——Weyl-2 K + BCS IF 今天就能跑，量级合理。

**F5【唯一开放推导点 = Weyl 标量的 BCS 广义化】** 现行 `lscale = coupling_norm/4·Σ step²` 只含 normal α 和。BCS 的反常耦合（Δ）与反常道二阶累积是否对标量有 O(δt) 贡献 = **Phase 1 实验裁决**：naive 组合已 2 阶收敛 ⇒ 无贡献（成文）；若实测 1 阶 ⇒ 需要 Nambu 双高斯测度对称离散化残余的 Δ 修正项推导（子工作项，Plan §2.1）。
  - **✅ 已裁决（2026-11-17，bcsweyl/）**：反常道对 Weyl 标量**零贡献**——`coupling_norm(::AbstractBCSBath) = coupling_norm(bath.spectrum)`（每带裸谱强度，Δ/μ/β 不进；Bogoliubov 正则性 u²+v²=1，反常接触为异槽对、κ² 级 O(λ³)/格 污染不到 O(λ) 标量）。实证：平凡系统接触缺口提取 cn*(Δ)/Σα ≡ 1（≤4e-4，Δ∈[0,2]，二次拟合 N≤64）；generic H_S 收敛斜率 1.97–1.98（cn=Σα）vs 0.98（cn=0 对照）。推导与数据档案 `bcsweyl/README.md`；方法已落地 correlationfunction.jl + test/bcs 语义 testset。

**F6【现状缺口】** `test/bcs/` 目录为空（套件未建）；无 BCS order=2 测试与 goldens；ED 锚点函数现成（`bcs_ed` test/util.jl:219，含 U/Δ 全参）；旧栈 BCS ED 管线（docs/tutorials/.../bcs_realtime.jl）可移植。

## 2. 阶段与验收

- **Phase 1（实验先行，本计划科学核心）**：移植旧 bcs ED 管线到新包（fewmodes 单模 + `bcs_ed` 参照），测 greater/lesser 全时序相对误差的 log-log 斜率，N ∈ {4, 8, 16, 32}，order=1 vs order=2 对照（K 侧唯一变量；IF 同为 bcs_partial_if）。
  - 判据：order=2 斜率 ≈ 2.00 ⇒ naive 组合正确，Weyl 标量无 BCS 修正（推导结论成文进 sysdynamics 注释）；斜率 ≈ 1 ⇒ 反常标量修正存在，进入推导子项（Nambu 双高斯残余），推导后回归实验。
  - 同时验证 :z2 guard 与 U≠0（Anderson 相互作用）下的稳定性。
- **Phase 2（形式化）**：① BCS order=2 放开的形式化——`fock_sysdynamics`/`sysdynamics_weyl2order` 无 BCS 特例（维持现状）+ 文档声明 `coupling_norm` 语义 = **normal 通道 α 和**（Δ 不进 Weyl 标量，按 Phase 1 结论）；② `bcs_symmetry_guard` :z2 保留；③ `test/bcs/` 套件落地：fewmodes ED 对拍（order1/order2 双档）+ 冒烟配置回归。
- **Phase 3（goldens + 收尾）**：BCS goldens 建档（此前无 BCS golden——新档案）；changes.md；与 `sysdynamics_fock_mpo_plan.md` 的 Weyl 标量段同步。

## 3. 不做清单

- :u1/:u1u1 的 BCS（荷不守恒，物理禁止——guard 保留）；
- IF 侧 BCS 二阶化（IF 精度独立于 K 离散化阶——bcs_partial_if 维持）；
- 多带 BCS 的 fullfactor/tti_* 细节（归 IF 后端计划 W2 的 ContourBlockIndex 统一）；
- 缓存/效率（同 V2 口径）。

## 4. 工作量与顺序

Phase 1 ≈ 2 天（管线移植 + 收敛实验；若触发标量推导 +1~2 天）；Phase 2 ≈ 1 天；Phase 3 ≈ 半天。总 ≈ 3.5 天（无标量修正线）。顺序 1→2→3。
