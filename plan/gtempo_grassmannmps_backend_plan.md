# GTEMPO 后端移植计划：Z2Tensors 栈 → GrassmannMPS 包

- 日期：2026-11-13
- 基线：GTEMPO HEAD `99455a0`（77 组测试 0 fail），全量备份已置于 `backup/GTEMPO/`（404M，含 .git）
- 参考实现：**U1GTEMPO2** —— 已完成的同型移植（零依赖 GTEMPO，纯 GrassmannMPS 后端，929 tests 绿 + golden 六算例机器精度 + 与本包 extensions/secondorderif 端到端 ≤1.4e-11 对齐，含 Weyl 二阶）

## 0. 目标

把 GTEMPO 的张量引擎从「Z2Tensors + `src/grassmanntensor/` + `src/grassmannmps/`」三层旧栈整体替换为独立包 GrassmannMPS（TensorKit sector 栈），GTEMPO 保留业务 API 皮肤（模型、lattice、IF、observables 公开函数名不变），同时完成三项功能调整：

1. **删除 GrassmannOrdering 排布家族**，只支持 `extensions/secondorderif` 中的 lattice 排布（ChronoKeldysh / ChronoMixed / 虚时默认序）；
2. **暂停 electronphonon** 全部功能；
3. **新增 `sysdynamics_weyl2order`**：Weyl 二阶支持，仅改杂质（K 链）动力学构造，IF 与收缩零改动；
4. **对称性用户手选**（2026-11-13 裁决）：bridge/扇区层按对称性参数化，GrassmannMPS 四形（Z2 / Z2Z2 / U1 / ChargeCharge）全部支持，经选择器由用户指定；**BCS bath 仅支持最基本的 Z2**（详见 §3.4 与 D1/D7）。

## 1. 探查定谳（计划赖以成立的关键事实）

**F1【Weyl 二阶只需改 K 构造——用户设想成立】**
`extensions/secondorderif/INTEGRATION.md` §2 明证：「GTEMPO 的 `integrate` 是裸 Berezin 配对机器，不携带任何一阶/二阶信息。Weyl 二阶的全部修改都折叠在 K 链的局部张量里；IF 与收缩完全复用」。代码印证：IF 构造（`influencefunctional/`）不读 `sysdynamics`，两者只在 `integrate(lattice, K, IF...)` 变参收缩会合；secondorderif 的 `influence.jl` 仅 58 行转发、`contraction.jl` 仅 28 行转发。

**F2【U1GTEMPO2 = 现成参考实现】**
零依赖 GTEMPO（deps 只有 GrassmannMPS + TensorKit 栈）。可直接搬运的部件：
- `bridge.jl`：站形构造（ChargeCharge 交替 up/down + 键恒等双荷）、`pairize`、`LegacyTruncation`（逐位复刻旧 `select_dimensions` 截断语义）、`compress!`（QR 前扫 + SVD 回扫 + logscale 漂移记账）、`multiply`、`from_core`/`split_cell`、`monomial`；
- `system.jl` + `prepared_system.jl`：cell automaton K 构造（`local_weight` q=1/2 Weyl / q=1 一阶、`WeylInsertion`、`cell_maps`、`source_schedule`、插入算符与 ρ₀ 折叠）——order=1/2 同一台机器；
- `influence.jl`：`fit_exponentials`（矩阵铅笔）→ `memory_cell` → `split_cell` → `pairize` → `multiply` 的 ExactTTI 装配，及 PartialIF 参考路径；
- golden 对拍方法论：`dump_golden.jl`/`check_golden.jl` + 独立 oracle + 对齐扫描三层。

**F3【secondorderif 的 lattice 排布是什么】**
- 实时：`ChronoKeldyshOrdering`（`contours.jl:58-84`）——块序 `[(+,1)…(+,N),(−,N)…(−,1)]`、块内 (a,ā) 相邻对齐、AdjacentConjugation + TimeLocalLayout、K 链按时序本征序直铺**零 permute!**（旧交错序 N=4 双带 >9min → 1.3s，约 400×）；
- 混合：`ChronoMixedOrdering`（`contours.jl:133-174`）——块序 [τ 格 1..Nτ][rt 边界][+ 格][− 格]，块数恰等于 MixedGrassmannLattice 长度；
- 虚时：沿用主包默认 `A1Ā1B1B̄1`（`contours.jl:99`），单元 j 放顶点 j+1（顶点 1 = t=0 平凡边界）。
- 三者 subclass 主包 `GrassmannOrdering` 类型树，但 `_index` 在 extensions 私有；单元式排布（每时间步一个全宽源单元，无 QuAPI 额外浴传播端点）。

**F4【GrassmannOrdering 删除影响面】**
12 个具体序 + trait 三件（`ConjugationStyle`/`LayoutStyle`/`TimeOrderingStyle`）+ `changeordering`/`toadjacentordering` + TTIIF/ExactTTIIF/TDVPIF 的 ordering gate（`_AllowedRealGrassmannOrdering`，ttiif/realtime.jl:1 等）+ `sysdynamics_fast` 的 `_fastordering`（sysdynamics.jl:305-307）+ electronphonon 的 `FockOrdering` 映射（`similargrassmannordering`，随暂停一起走）+ test/runtests.jl:13-16 的私有序导入。

**F5【electronphonon 耦合面小】**
非 electronphonon 代码对它的引用只有：export 块（GTEMPO.jl:81-90，18 个符号）+ `retardedinteractdynamics_naive(_!)` export（:51）+ include（:194）+ 两处注释。测试面：`test/electronphonon/`、`test/api/electronphonon.jl`、`docs/tutorials/holstein/`。

**F6【Z2 单键宇并是已知结构缺陷，正好被本次换轨修复】**
τ 奇插入两点 ~1.7e-3 结构尾部（INTEGRATION.md §9），根因 = GTEMPO 单 Z2 键宇称合并每带电荷；修复 = grading 泛化，参考实现即 U1GTEMPO2 的 ChargeCharge 双荷键。Weyl 二阶 mixed 轮廓正是命中此缺陷的场景。

**F7【GrassmannMPS API 齐备】**
`mult!/compress!`（DMRG2/SVDCompression）、`swap!/reorder!`、`expmult!`（PartialExp/XTRG/PolarExp/TDVP2）、`SchurGMPSTensor/SchurGMPS`/`timeevo`(WI/WII)、引擎 `environments`/`ExpectationCache`/`expectation(op::SparseGMPS, env)`/`Zvalue`（berezin.jl，K≤5 链）、`coefficient_scaled`/`todense`、`truncrelerr/truncdimcutoff`、`vacuumstate`、因子族（Monomial/Pairing/ExpDecay 族）。

**F8【测量架构先例】**
U1GTEMPO2 最终采用「带插入 K 全量重建 + fresh environments」（solver.jl:151-179）；GM 引擎的 `expectation(op::SparseGMPS, env)` 支撑扫掠路径已在 GM engine.jl 验证（9/9）但未被 U1GTEMPO2 采用。两条路线需在 golden 对拍下裁决（D5）。

**F9【GrassmannMPS 对称性支持面——对称性可选的可行性基础】**
规范形注册表四形（abelian 扩展已交付、全套绿）：**Z2**（合并宇并，两朝向）、**Z2Z2**（每带宇并，自对偶两朝向）、**U1Irrep spinless**（单整数荷、异扇区无朝向）、**ChargeCharge**（U1×U1 双荷，:up/:down）；六站形 digit 语义统一（vac0/pair1/eta2/etabar3）、跨对称性同算符 dense 逐位相等有测试保证（GM test/abelian.jl 352 项）。物理约束：BCS 反常关联破坏粒子数守恒（每带荷与总荷皆破），宇称（±2 变化）恒守恒 → **BCS 唯一合法对称性 = Z2**——与用户裁决一致。

## 2. 目标架构

```
用户 API 皮肤（不变）：GrassmannLattice / AndersonIM… / sysdynamics / hybriddynamics /
                      integrate / gf / occupation / currents …
        │
业务层（保留+重接）：lattices(Chrono 三序) · sysdynamics 模型层 · correlationfunction(QuAPI)
        │
构造层（自 U1GTEMPO2 移植 + 对称性泛型化 §3.4）：bridge(站形/截断/pairize，四形可选)
                            · K cell automaton(order 1/2) · IF 装配(PartialIF/ExactTTI/BCS-Z2)
        │
引擎层（GrassmannMPS 包）：mult!/compress! · reorder! · environments/Zvalue/expectation
                          · SchurGMPS/timeevo(WI/WII) · expmult!
```

### 保留（零改动或只改接线）
| 部件 | 处置 |
|---|---|
| `sysdynamics/predefined/{siam,irlm,skmodel}.jl`、`generalimpurity/{fockhamiltonian,fockmatrix}.jl` | 零改动（`fock_propagator`/`fockmatrix` 纯数值，无张量依赖） |
| `correlationfunction.jl` + QuAPI/ImpurityModelBase 依赖 | 零改动 |
| `auxiliary/linalg.jl` | 零改动 |
| observables 物理公式层（gf/occupation/currents/nn/NEGF 定义） | 张量操作重接（Phase 5） |
| `sysinitstate.jl` 语义 | 重接：ρ₀/热态折进 K 构造单元 1 |

### 移植（U1GTEMPO2 → GTEMPO）
| 部件 | 落点 |
|---|---|
| bridge.jl 全套 | 新 `src/bridge.jl` |
| K cell automaton（system/prepared_system） | `src/sysdynamics/`（挂 sysdynamics 命名空间） |
| ExactTTI/PartialIF 装配 | `src/influencefunctional/` 重写 |
| Chrono 两 ordering + `_index` | `src/lattices/` |

### 删除
| 部件 | 理由 |
|---|---|
| Z2Tensors、FiniteMPSAlgorithms 依赖 | 后端退役（乘法/MPO/正交全走 GM） |
| `src/grassmanntensor/`（3 文件） | GM.GrassmannBackend + GM.@grassmann 替代 |
| `src/grassmannmps/`（11 文件旧栈） | GM 包替代 |
| `src/lattices/grassmannordering.jl` + `changeordering.jl` | 功能调整 #1（12 序 + trait + changeordering 全删，仅留 Chrono 三序 `_index`） |
| `src/mpo/` 大部 | SchurMPO/WII/稠密 MPO → GM.SchurGMPS/GM.timeevo；DMRGAlgorithm → GM.DMRG2（随 D3 裁决可整删） |
| `src/gvconnections/` | 边界闭合折进 K 构造单元（secondorderif 方式）；视 D6 留 bulkconnection 单函数 |
| `src/electronphonon/` | 功能调整 #2（暂停） |
| `ext/` 整目录（GTEMPOCUDAExt + ext/test，11 个 .jl，~76K） | **直接删除**（2026-11-13 用户裁决）：依赖 Z2TensorsCUDAExt 与 GTEMPO 20 个内部符号，随旧栈退役，不保留恢复桩 |
| sysdynamics 旧装配路径 | `_propagator_sparsegmps` 窗口 mult!、`sysdynamics_fast`/`_fastordering`、`baresysdynamics*`、8 个 `*_deprecated`、`termdynamics.jl`（已注释）、`fermionoperator.jl`（已注释）——K 构造统一走 cell automaton（D2） |
| `algorithms.jl` 的 DMRG1/DMRG2/SVDCompression/Orthogonalize | 改为 re-export GM 对应物 |

## 3. 三项功能调整落法

### 3.1 删除 GrassmannOrdering，只支持 secondorderif 排布
1. `ChronoKeldyshOrdering`/`ChronoMixedOrdering` 从 `extensions/secondorderif/src/contours.jl:58-174` 移入 `src/lattices/`（连 `_index` 重载与相邻共轭语义；若 trait 体系删除则 Chrono 类型直接携带固定行为，不再走 trait 分派）；
2. 虚时排布内联：`imaginarytime.jl:53` 的 `A1Ā1B1B̄1` `_index` 公式保留，序类型收编（建议命名 `ChronoImagOrdering`）；
3. 三个 lattice 类型保留（`ImagGrassmannLattice`/`RealGrassmannLattice`/`MixedGrassmannLattice`），构造器删 `ordering` kwarg（固定 Chrono）；`GrassmannLattice(; contour=:imag/:Keldysh/:mixed)` 工厂签名不变；
4. 删 `grassmannordering.jl` 全文 + `changeordering.jl` + 各 lattice 文件的旧 `_index` 方法 + IF 系 gate 收窄为单一序；
5. 连带删测试中 ordering 扫描（test/runtests.jl:13-16 私有序导入、test/api/lattice.jl ordering 组、multiorbital 的 `for ordering` 循环）。

### 3.2 暂停 electronphonon
1. `src/GTEMPO.jl` 删 include（:194）+ export 块（:81-90）+ `retardedinteractdynamics_naive(_!)`（:51）；
2. 删 `test/electronphonon/`、`test/api/electronphonon.jl`、runtests include 行；
3. `docs/tutorials/holstein/` 移出或删除；
4. 处置方式建议 **git rm**（backup/GTEMPO 已有全量副本，恢复路径写进 changes.md）；
5. 注意：`ImpurityModelBase`/`QuAPI` 依赖保留（bath 与 `GrassmannContourIndex` 仍需要）；`FiniteMPSAlgorithms` 依赖是否残留视 mpo/ 删除范围定。

### 3.3 sysdynamics_weyl2order
- 移植 U1GTEMPO2 `system.jl`/`prepared_system.jl` 的 cell automaton 机器 → `src/sysdynamics/`；
- API 双入口（同一台机器，只差 `local_weight` 的 q）：
  - `sysdynamics(lattice, model; trunc)` —— order=1 相干（q=1），默认路径；
  - `sysdynamics_weyl2order(lattice, model; trunc)` —— order=2 Weyl（q=1/2、Weyl 标量 `exp(±h²Σα/4)`、F 路由与跨带 XOR 补偿、`variant=:dag`）；
- ρ₀/初态折进单元 1、迹闭合折进单元 N → `boundarycondition!`/`bulkconnection` 对 K 链不再需要；
- IF 与收缩零改动（F1）；
- 验收：Weyl 收敛阶 2.00（`extensions/secondorderif/toulouse` 与 `singlemodebath` 基准搬入 validation）+ U=0 因子化恒等式（Z₂=Z₁²、K₂=K₁∧K₁）+ ED 对拍 + τ 奇插入两点尾检查（R3）。

### 3.4 对称性用户选择器（2026-11-13 裁决落地：Symbol 输入）

API：`GrassmannLattice(...; symmetry::Symbol=:z2)`——四个 Symbol（扇区类型映射到 GM 规范形；带→荷槽映射是 GTEMPO 层策略）：

| `symmetry` | GM 规范形 | 带→荷槽映射 | 适用 bands | 备注 |
|---|---|---|---|---|
| `:z2` | Z2Irrep | 全带合并宇并（= 旧 GTEMPO 语义） | 任意 | **默认值**；BCS 唯一合法选项 |
| `:z2z2` | Z2Z2 | band1/2 各占一个宇并槽 | ≤2 | 每带宇并分辨 |
| `:u1` | U1Irrep spinless | 全带合并总粒子数 | 任意 | 数守恒模型；异扇区无朝向 |
| `:u1u1` | ChargeCharge | band1→:up 槽、band2→:down 槽 | ≤2 | 每带荷分辨（τ 尾缺陷修复通道，U1GTEMPO2 同款交替站形） |

- **传递方式**：`symmetry` 随 lattice 携带，`sysdynamics`/`hybriddynamics`/`integrate`/观测构造统一从 lattice 读取。ordering 与 symmetry 正交：前者管链排布（站的位置），后者管扇区标记（键的荷结构），互不干涉。
- **Symbol 校验**：未知 Symbol → `ArgumentError`（合法集 `{ :z2, :z2z2, :u1, :u1u1 }`，不做静默回退）。
- **显式拒绝（禁静默降级）**：`:z2z2`/`:u1u1` × bands>2 → `ArgumentError`；`:u1` × 非数守恒问题（如 BCS）→ `ArgumentError`（φ 荷合并静默降级对角 U1 已有否决先例）。
- **站形推导**：bridge 按 symmetry 推导站形与荷槽——`:u1u1` 走 U1GTEMPO2 同款交替站形（奇站 :up/偶站 :down）；`:z2`/`:u1` 全站同形（单荷槽，无需交替）；`:z2z2` 两朝向交替。
- **验收新增「跨对称性一致电池」**：同一物理问题在所有适用对称性下各跑一遍，观测互相一致（截断容差内）且合 ED oracle；**Z2 通道另锚 golden**（旧后端唯一扇区同构通道，稠密系数直接可比）。

## 4. 设计裁决点

### D1 对称性方案 —— ✅ 已裁决（2026-11-13 用户）：对称性用户手选，四形全支持；BCS 仅 Z2
- **裁决内容**：bridge/扇区层按对称性参数化，不硬编码 ChargeCharge；用户经选择器（§3.4）在四形中任选；BCS bath 路径仅接受 Z2（构造器级检查）。
- **实现要点**：站形/扇区代数/`LegacyTruncation` 荷排序（(value, charge, order) 跨荷预算）对荷类型泛型——GM abelian 扩展已保证四形统一 digit 语义与跨对称性 dense 相等（F9）；带→荷槽映射是 GTEMPO 层策略（§3.4 表）。
- **约束**：`:z2z2`/`:u1u1` 每带分辨要求 bands≤2，超限显式拒绝（禁静默降级）；BCS 反常关联破坏粒子数守恒、宇称恒守恒 → 仅 Z2。
- **Kanamori norb≥2（4 带）**：荷分辨形不可用，但用户选 `:z2`/`:u1` 仍可跑（合并标记 = 旧语义/总荷）；恢复每带荷分辨 = GM 注册表扩 U1^M（未来独立小任务，不在本计划范围）。
- **覆盖风险**：U1GTEMPO2 先例只验证了 ChargeCharge 通道；其余三通道靠 GM abelian 测试电池 + 本包跨对称性一致电池（Phase 5 验收）+ ED oracle 兜底（R9）。

### D2 K 构造统一走 cell automaton【建议采纳】
旧窗口 SparseGMPS 路径（`_propagator_sparsegmps` + `sysdynamics_fast` + bare 族）删除；order=1/2 同一机器（§3.3）。收益：删 ~900 行旧路径、K 构造直接产出 Chrono 本征序（零 permute）。

### D3 IF 算法保留集【待裁决】
- 建议：**第一批** PartialIF（naive/hybrid）+ ExactTTI（U1GTEMPO2 `influence.jl` 路线）；**第二批** TTIIF 经 GM.SchurGMPS + timeevo(WII) 移植；**XTRGIF** → GM.`expmult!(XTRG)` 现成可接（GM 已有完整 XTRG 实现）；**TDVPIF 暂停**。
- 影响面：TTIIF/XTRGIF 保留则 `src/mpo/` 需保留薄壳（或直接消费 GM 的 Schur/timeevo）；TDVPIF 暂停则删 `tdvpif/`。

### D4 integrate 引擎【建议采纳】
`integrate(alg, lattice, x0, x1...)` 内部走 GM.`environments` + `Zvalue`（K + IF1 + IF2 三链，K≤5 支持充足）；`ExactIntegrate` zipup 退役；**BMPSIntegrate 删除**（Chrono 排布恒为 AdjacentConjugation，AC 语义由 GM 引擎覆盖）。`IntegrationAlgorithm` export 收敛。

### D5 测量架构【待裁决】
- 主路线建议：GM 引擎 `expectation(op::SparseGMPS, env)`（支撑扫掠；GTEMPO 的 environments 第一链 = 算符目标链约定与 GM 引擎同型）；
- fallback：inserted-K 重建（U1GTEMPO2 先例，solver.jl:151-179）；
- Phase 5a 两条都对拍 golden，按精度/性能定案后删另一条。

### D6 gvconnections 处置【待裁决】
建议：`boundarycondition` 随 K 构造闭合退役；`bulkconnection`（观测算符跨步重叠）若 observables 仍需要则保留单函数桥。

### D7 BCS —— ✅ 已裁决（2026-11-13 用户）：保留，仅支持 Z2
- `bcsinfluencefunctional/` 随 Phase 3 重接（`hybriddynamics!` 对 `BCSCorrelationFunction` 的 partial-IF 乘路径）；对称性检查：BCS 入口断言 `symmetry === :z2`，否则 `ArgumentError`（物理依据：反常关联破坏粒子数守恒，仅宇称守恒，见 F9）。
- GM.`AnomalousPairing` 因子可用（Z2 守卫已自对偶化，Z2 通道原生支持）。
- 验收：一组反常通道 golden/ED 对拍。

## 5. 分阶段实施

> 依赖关系：Phase 0 完全独立先行；1→2→3→4→5 串行，每阶段以 golden 对拍为闸门；Phase 6 收尾。
> 建议在 git 分支 `grassmannmps-backend` 上进行。

### Phase 0 —— 基线冻结与对拍基建（不动 src）
- 当前套件全跑留档（77 组绿基线）；
- **golden fixtures**：`test/fixtures/dump_golden.jl` 在旧后端一次性 dump 算例矩阵（imag/real/mixed × Anderson/Toulouse × bands 1/2 × ExactTTIIF/PartialIF × 一阶；6~9 例），每组含：lattice 位置映射、K 链 `todense` 系数 + logscale、IF 链系数、`Zvalue`、occupation/gf 全端点；
  ⚠️ **对称性换轨后键维序列必然不同**——golden 对拍对象是稠密系数/物理量（无截断算例 1e-10 档，参照 U1GTEMPO2 port_golden 容差），键维只做新实现内部记账，不做跨栈对拍；
- golden 锚定 **Z2 通道**：新实现选 `symmetry=:z2` 时扇区结构与旧后端同构（全带合并宇并），稠密系数直接可比；其余对称性通道以 GM 跨对称性 dense 相等 + ED oracle 锚定（Phase 5 验收电池）；
- alignment sweep 基建：vs `backup/GTEMPO` 扫描脚本（参照 U1GTEMPO2 `align_sweep`/`compare_sweep` 模式）；
- **验收**：golden JSON 入 `test/fixtures/port_golden/`，`check_golden.jl` 对旧实现自洽通过。

### Phase 1 —— 底座（bridge + lattice 简化 + 依赖换轨）
- `Project.toml`：删 Z2Tensors、FiniteMPSAlgorithms（视 mpo 处置），加 GrassmannMPS（dev path）；
- 新 `src/bridge.jl`（U1GTEMPO2 bridge.jl 移植：站形/pairize/LegacyTruncation/compress!/multiply/from_core/split_cell/monomial）+ **对称性泛型化**（原版硬编码 ChargeCharge 交替站形；泛化为按 `symmetry::Symbol` 推导站形与荷槽映射，§3.4 四形）+ Symbol 校验；
- export 皮肤对表：`const GrassmannMPS = GM.GMPS`、`const SparseGMPS = GM.SparseGMPS`、`truncdimcutoff`/`DMRG2`/`leftorth!`/`rightorth!`/`mult!`/`canonicalize!` 等逐一映射到 GM（GTEMPO.jl export 名单不变，下游教程/测试零改动）；
- Chrono ordering 移入 + GrassmannOrdering 家族删除（§3.1）；
- 删 `grassmanntensor/`、`grassmannmps/` 旧栈；`@grassmann` re-export GM 版；
- **验收**：包可加载、真空态/randomgmps 烟测（四形各一）、对称性校验拒绝路径单测（bands>2 × 荷分辨形、BCS × 非-Z2）、lattice `_index` 单元测试（位置映射 vs Phase 0 golden）、`vacuumstate`/`coefficient_scaled` 往返。

### Phase 2 —— K 构造（sysdynamics + sysdynamics_weyl2order）
- 模型层保留接线；cell automaton 机器移植（§3.3）；删旧装配路径（D2）；
- `sysinitialstate`/`systhermalstate` 折入单元 1；
- **验收**：一阶 K `todense` 对拍 golden；`sysdynamics_weyl2order` vs secondorderif 参考数值（`extensions/secondorderif/validation/` 的 rt_fastpath/u0_factor 三件套直接引用，不搬运）；迹归一化 ∏(1+c̄c) 恒等式。

### Phase 3 —— IF 构造
- PartialIF/ExactTTI 装配路径移植；bath 核对接：`KernelCorrView`（secondorderif kernels.jl:122-149，QuAPI 接口桥）或 `correlationfunction` 直连（INTEGRATION §10.3 评估项，取改动小者）；
- BCS（`bcsinfluencefunctional/`）重接：partial-IF 乘路径 + **Z2-only 守卫**（D7 裁决）；
- **验收**：golden IF 系数对拍 + `if_band_invariance`（IF 带不变性）+ ExactTTI vs PartialIF 互检 + BCS 反常通道 golden/ED 对拍一组。

### Phase 4 —— integrate 引擎
- `integrate` → GM.`environments`/`Zvalue` 接线（D4）；
- **验收**：golden Z 对拍 + U=0 因子化三件套（Z₂=Z₁²、K₂=K₁∧K₁）+ 一阶 vs 二阶共享 IF 检验 + `extensions/secondorderif` 的连续浴对齐脚本复跑（≤1e-10 档）。

### Phase 5 —— observables（分三批）
- **5a**：`gf`/`greater`/`lesser`/`occupation` + 归一化（D5 两路线对拍定案）；
- **5b**：`cached_*` 快路径（GM `ExpectationCache` 映射；multintegrateband 语义由引擎支撑扫掠替代）；
- **5c**：currents/heat/nn/NEGF；
- **验收**：golden 物理量对拍 + Toulouse 解析收敛基准 + `test/api/quapi_crosscheck.jl` 保留组 + **跨对称性一致电池**（同一物理问题 × 四形各跑，观测互一致且合 ED，§3.4）；允许 5c 降级「先精确版、快路径后补」。

### Phase 6 —— 功能收尾
- electronphonon 摘除（§3.2）、**CUDA ext 直接删除**（`git rm -r ext/`；清 Project.toml 残留 weakdeps/extras 条目；grep 全包确认无 `Z2TensorsCUDAExt`/CUDA 条件路径残留）、docs/tutorials 改写（ordering kwarg 已删、新 K API）、`changes.md` 记录（含 electronphonon/TDVPIF 恢复路径与 CUDA 移除说明——未来若需 GPU 走 GrassmannMPS 侧方案，另立项）；
- 测试套件重组：api/normalbath/bcs(视 D7)/multiorbital(≤2 带保留)；
- **总验收**：全套绿 + vs `backup/GTEMPO` alignment sweep 报告 + Weyl 收敛阶 2.00 报告 + 性能基线（参照 U1GTEMPO2 longchain CSV 格式）。

## 6. 风险表

| # | 风险 | 缓解 |
|---|---|---|
| R1 | 截断语义换轨（`select_dimensions` 跨荷预算 → MAK TruncationStrategy） | `LegacyTruncation` 已逐位复刻（U1GTEMPO2 D1/D2 先例）；对拍键 = 无截断算例 + CompressionReport 记账 |
| R2 | 费米符号约定重标定 | secondorderif T1–T7 标定陷阱是 Z2 栈约定（graded⊗ OFF、oddlater=:none 等），GM/TensorKit 栈约定不同；**以 U1GTEMPO2 代码为准，禁止照抄 secondorderif 符号位**（其端到端 ≤1.4e-11 对齐已证明两套各自自洽） |
| R3 | τ 奇插入两点尾 | Z2 下 ~1.7e-3 结构尾（F6）；ChargeCharge 下预期消失——列入 Phase 2 验收（复跑 `diag_tail_investigation.jl`） |
| R4 | Chrono 排布下 IF 构造 | `hybriddynamics` 与 K 必须同一 lattice 实例（(a,ā) 相邻是 AC 积分前提）；Phase 3 用同一 `weyl_lattice` 构造 |
| R5 | observables 广度回归成本 | cached_*/NEGF 机器庞大；分三批 + golden 兜底；5c 允许先精确版 |
| R6 | 性能回归 | U1×U1 块维小于 Z2 合并块、键空间分辨更多；以 U1GTEMPO2 longchain 基线（`validation/records/port_*_longchain.csv`）为参照，Phase 4/5 各跑一次 benchmark |
| R7 | 命名冲突 | GTEMPO export 类型名 `GrassmannMPS` vs 包名：`import GrassmannMPS as GM` + `const` 别名；GM 对 norm/dot 的重载关系按 GM 模块 docstring 处理（选择性 import，勿整包 using） |
| R8 | 私有符号依赖 | GM 引擎私有名（`_boundary_pair`/`_bond_op_left`/`_absorb_*`/`_rescale_sites!`/`_left_step!`）按 U1GTEMPO2 纪律「只 import 不扩方法」；升级 GM 时回归 |
| R9 | 对称性泛型化的四形覆盖风险 | U1GTEMPO2 先例只验证 ChargeCharge 通道；其余三通道依赖 GM abelian 测试电池（352 项 + 跨对称性 dense 相等）+ 本包跨对称性一致电池 + ED oracle 三层兜底；BCS 强制 Z2 缩小暴露面 |

## 7. 待用户裁决清单（开工前需确认）

1. ~~**D1**~~：✅ 已裁决（2026-11-13）——对称性用户手选四形全支持，BCS 仅 Z2（§3.4/D1）；
2. ~~**D7**~~：✅ 已裁决（2026-11-13）——BCS 保留，仅 Z2；
3. **D3**：IF 保留集（建议：PartialIF + ExactTTI 首批；TTIIF/XTRGIF 二批；TDVPIF 暂停）；
4. **D5**：测量架构两路线对拍后定案 —— 是否同意先双轨再收敛；
5. **D6**：gvconnections 是否随 K 闭合退役；
6. 分支策略：`grassmannmps-backend` 分支 or 直接 main（backup 已在）；
7. ~~选择器命名~~：✅ 已裁决（2026-11-13）——Symbol 输入 `:z2`/`:z2z2`/`:u1`/`:u1u1`，默认 `:z2`（§3.4）。

## 8. 规模粗估

| 阶段 | 相对规模 | 主要不确定性 |
|---|---|---|
| Phase 0 | 小（脚本为主） | 算例矩阵选取 |
| Phase 1 | 中偏大（bridge ~400 行移植 + **对称性泛型化** + lattice 删改 + 依赖换轨） | export 对表遗漏面、四形站形推导 |
| Phase 2 | 中（cell automaton ~850 行移植改造） | ρ₀/闭合折入与 golden 容差 |
| Phase 3 | 中偏大（ExactTTI ~350 行 + PartialIF 重接 + BCS Z2-only） | KernelCorrView vs correlationfunction 直连 |
| Phase 4 | 小（引擎接线） | BMPSIntegrate 用户存量 |
| Phase 5 | **大**（observables 13 文件 + cached 机器） | 5b/5c 快路径是否逐个保留 |
| Phase 6 | 小 | — |
