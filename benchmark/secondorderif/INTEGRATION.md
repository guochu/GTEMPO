# SecondOrderWeyl 二阶 Weyl 算法：成果总结与 GTEMPO 集成指南

> 2026-11-07 定稿。本文档面向**将来把二阶 Weyl 算法集成进 GTEMPO 主库**的
> 实施者：完整记录当前成果状态、架构复用边界、二阶形式的全部标定规格
> （含"看起来错其实对"的约定清单）、验证电池、金标准对照状态与已知限制。
> 历史判决链见 `HISTORY.md`（十二轮），本文只做集成视角的终态总结。

---

## 1. 定位与当前状态总览

`benchmark/secondorderif/src/SecondOrderWeyl.jl`（1605 行）是 **Weyl 二阶
展开（U₂(λ) = exp(λ²c/4)·e^{−λH_S/2}·Q_W[e^{−λH_Bᴳ}]·e^{−λH_S/2}）在
GTEMPO Z2 分级栈上的完整实现**，覆盖三种 contour（实时 / 混合 / 虚时），
与 GTEMPO 一阶（anti-Wick）共用 IF、收缩引擎与 lattice 基础设施。

| 能力 | 状态 | 精度 |
|---|---|---|
| rt 两点（greater/lesser 全网格，含 (+,0) 与 − 支） | ✅ ED 闭合 | 1.54e-3 @Nt=8（O(h²) 干净） |
| rt 两点 vs 金标准 U1GTEMPO2 | ✅ 对齐 | Nt=8 最差误差逐位相同（1.359e-3）；首格 O(h²) 细节差 |
| 混合两点 vs ED / vs U1 | ✅ 闭合 / ✅ **机器精度逐位** | 6.8e-4→1.7e-4；\|mine−U1\| = 3.7e-14 / 7.5e-12 |
| τ 两点全网格（含 j==0 平移路由） | ✅ ED 闭合 | 1.71e-3 @N=16（已知最优 Z2 值，见 §9 限制） |
| occupation（三轮廓 × 逐带） | ✅ ED + U1 对齐 | mixed 逐位 ≤2e-8；连续谱机器精度 |
| 收敛阶 scaling（Toulouse/Anderson × rt/τ/mixed） | ✅ | 2.001/1.996、2.002/2.004、2.005/2.004 |
| 连续谱浴（cos² 有界谱） | ✅ 对齐 + 收敛阶闭合 | U1 对齐 ≤1.4e-11；τ order-2 逐对精确 2.000 |
| U=0 因子化恒等式族 | ✅ | Z₂=Z₁²、K₂=K₁∧K₁（0/400）、**M2≡M1 逐点**（机器锐利判据） |

---

## 2. 架构与复用边界（核心结论）

**GTEMPO 的 `integrate` 是裸 Berezin 配对机器，不携带任何一阶/二阶信息。**
Weyl 二阶的全部修改都折叠在 K 链的局部张量里；IF 与收缩完全复用。本模块
本身就是这个解耦的证明：

| 层 | 本模块实现 | 复用 GTEMPO |
|---|---|---|
| IF 构造 | `influence.jl`（58 行转发） | `hybriddynamics`/`partialif_hybrid`（一、二阶共用同一 IF） |
| 收缩 | `contraction.jl`（28 行转发） | `integrate`（ExactIntegrate zipup，分级 Berezin 配对） |
| lattice | 三个 Weyl contour **就是** GTEMPO lattice | `ImagGrassmannLattice`/`RealGrassmannLattice`/`MixedGrassmannLattice` |
| 单元核 | `realtime_kernel_quapi`/`mixed_kernel_quapi`/τ `Δτ` | QuAPI 包 `Δt`/`Δm`/`Δτ`（离散谱精确极点装配） |
| 观测量机制 | `measurement.jl`（140 行，插入约定见 §5） | 插入 → K 重建 → 同一 `integrate`（与 GTEMPO `gf` 机制同构） |
| **系统因子 K** | **`system.jl`（692 行）+ `tensors.jl`（166 行）** | —— **唯一的新内容 = "sysdynamics_weyl"** |

其它复用：`ImpurityModelBase`（bath/DiracDelta/DiscreteSpectrum、
ToulouseIM/AndersonIM、`fock_propagator`）、QuAPI 包。

---

## 3. 文件清单（src/，共 1605 行）

| 文件 | 行数 | 内容 | 集成处置 |
|---|---|---|---|
| `system.jl` | 692 | **二阶核心**：`local_weight`/`_slot_operators`/`_fock_annihilators`、`WeylInsertion`、`system_factor`（_rt/_mixed/_tau 三装配） | **搬运（sysdynamics_weyl 实体）** |
| `tensors.jl` | 166 | `CellCore`/`split_core`（稠密格核 → Z2 site 张量的精确 SVD 拆分）/`tile_cells!` | **搬运**（或对齐 GTEMPO sysdynamics 的核拆分设施） |
| `measurement.jl` | 140 | `occupation`/`raw_correlation`/`greater`/`lesser`/`matsubara`（含 j==0 平移路由） | **搬运**（Weyl 测量 API + 约定） |
| `contours.jl` | 197 | 三 contour 薄壳 + `weyl_lattice` + **`ChronoKeldyshOrdering`/`ChronoMixedOrdering`（自定义，GTEMPO 无对应物）** | **搬运 ordering**；contour 薄壳可映射到既有 API |
| `kernels.jl` | 157 | `HybridizationKernel`、`α_weight`、rt/mixed QuAPI 核、`KernelCorrView` | 大半可省：Weyl lattice 是真 GTEMPO lattice，合并版可试 `GTEMPO.correlationfunction(bath, lattice)` 直连 `hybriddynamics`（需核对带 gate 与 getindex 语义，§6-T7） |
| `solver.jl` | 59 | `Problem`/`Solution`/`solve` | 与 GTEMPO 既有 Problem/Solution 形状对齐即可 |
| `influence.jl` | 58 | IF 逐带装配 wrapper | 复用（或并入 GTEMPO 逐带 IF 工具） |
| `contraction.jl` | 28 | `partition_function`/`expectation` 薄壳 | 删除（直接用 GTEMPO） |
| `diagnostics.jl` | 35 | `fit_order`/`convergence_report` | 可选 |

**搬运清单合计 ≈ 1160 行 + 验证电池（§7）**。

---

## 4. 二阶 K 构造完整规格（system.jl）

### 4.1 Weyl 局部权重（`local_weight`）

每带 4 槽基 `[q·I − n, a†, a, W₄]`，按带序乘积展开 4^M 个 2^M×2^M 矩阵：

- **q = 1/2（Weyl）/ 1（一阶 anti-Wick）** —— 常数项的序差；
- **W₄ = I/2 + n（Weyl）/ n（一阶）** —— 高斯核的 Weyl 序移动
  （`mode=:gaussian`，ED 标定 = `validation/calibrate.jl`）；
- **a 槽 ↦ +d†（`variant=:dag`）**——吸收 GTEMPO 分级 trace 配对
  （`g_trace_phy`）在奇宇称扇区多出的 −1。**:std（−d†）仅供对照，勿用**。
- ⚠️ **不能用 2-bit 槽码直接索引 `local_weight` 的 4^M 列表**（M≥2 会取错
  矩阵并重复计入 qI−n）——必须经 `_slot_operators` 的逐带乘积装配。

JW 约定（`_fock_annihilators`）：带 b 的算符穿越所有更低索引带（band 1
无符号；band 2 带 (−1)^{n₁}）——与 U1GTEMPO2 `fock_operators` 对齐。

### 4.2 单元装配

- 二阶：`mats[m] = sc · S · M(m) · S`（对称半步，O(h³)/格）；
  一阶：`S·M`（+支）/ `M·S`（−支）（Lie，O(h²)/格）。
- **Weyl 标量**：τ 格 `sc = exp(+h²·Σ_b α_b/4)`、rt 格 `sc = exp(−h²·Σ_b α_b/4)`
  （全局 exp(λ²c/4) 的单元化；λ = ±ih）。注意符号随轮廓段翻转。
- ρ₀：τ = I（热迹，boundary-0 插入折入 `R = O·R`）；rt = 因子化初态对角
  （用户经 `Problem(…; rho0)` 传入，(+,0) 切位插入左乘）。
- 穿线语义：`pairs = [(s,t) for t in 1:dim for s in 1:dim]`（t 大股 s 快，
  U1 约定）；内部格 `(tL==tR) ? A[sL,sR] : 0`（t 线恒等穿线）；fold
  `Σ_{s0} A[sL,s0]·R[s0,tL]`；close `A[tR,sR]`；N=1 退化格 `Tr[A·R]`。

### 4.3 F 路由与跨带补偿（十二轮标定的核心，勿动）

- **宇称路由**：`inpar = pairpar(l) ⊕ F_j`、`outpar = pairpar(l) ⊕ F_{j−1}`
  （第一轴朝晚端）；`F_j` = 切位奇插入累积宇称（含 boundary-0 的 F₀）。
  每带版本 `Fb_τ[b]` 同步累积（供补偿用）。
- **rt/mixed 逐元素补偿**（`compexp`，每格每键元素）：
  `exp = ⊕_b c_b·(p_b(l)⊕F_b) ⊕ ⊕_{b≠b'} c_b·p_{b'}(l)`
  （c_b = 槽码宇称、p_b(l) = 键态带 b 的 s⊕t 宇称）。
- **τ 补偿（TXCOMP，F-gated）**：`acc = ⊕_b c_b·Fb_b ⊕ 交叉项(出键)`，
  仅在 `F_before ≠ 0` 的格施加；交叉项 = `⊕_{b≠b'} c_b·p_{b'}(出键)`。
- **τ fold R 线交叉项（第十一轮 V1）**：fold 输出键真空 ⟹ 键公式交叉项
  恒 0；交叉项按 **R 线 (s0→tL) 逐 s0** 拆开 fold 求和施加（与内部格公式
  同构）。fold 门控 = `F₀ ≠ 0`。N=1 退化格天然因子化，**无补偿**。
- **oddlater = :none（三路径皆然）**：分级收缩已供给全部奇算符交叉符号；
  显式再加会重复计入（τ 上表现为 δt 无关的线性漂移）。
- ⚠️ **graded⊗ 符号（U1 `local_weight` 里有 `(−1)^{popcount(a)·popcount(b)}`）
  在本 Z2 管线中必须为 OFF**——分级收缩 + `:dag` 变体已供给，再加即重复
  计入（第六轮组合扫描实证；M=1 不敏感，M≥2 致命）。U1 侧为 ON 是其
  Gaussian 测度 + U(1) 管线的自有约定，不可移植。

### 4.4 插入机械（`WeylInsertion`）

- `WeylInsertion(branch, boundary, matrix; side=:left)`；cut 位置：τ j → j、
  + j → j（mixed: Nτ+j）、− j → 2N−j+1（mixed: Nτ+2Nt−j，(−,0) = 末格）。
- `cutop(k)`：该切位插入矩阵**左乘**（O 在左 = 轮廓更晚）——同切位多算符
  的 0⁺ 语义由 **push 序**承载（先 push = 轮廓更早）。
- **oddlater 不施**；插入后 f_k 进 F 路由与补偿门控。

### 4.5 存储（零 permute!）

- rt：`ChronoKeldyshOrdering`（自定义 ordering；块序 [(+,1)…(+,N),(−,N)…(−,1)]，
  块内 ysv 序）——K = 时序链本征序，MPS 免 permute!（N=8 双带 1.3s，
  旧 permute! 路线 >9min，~400×）。
- mixed：`ChronoMixedOrdering`（块序 [τ 边界][τ 格 1..Nτ][rt 端点][+ 格][− 格]，
  块数恰 = MixedGrassmannLattice 长度）。
- 装配循环的 ysv↔site 带映射：`ysv[v] = (m >> (2·(M−((v−1)÷2+1)) + (v−1)%2)) & 1`
  ——slot code 带序 = band 1 高位；直接 `(m>>(v−1))&1` 在 M≥2 **互换两带**
  （带对称观测量不可见的历史陷阱，第十轮补遗定位）。

---

## 5. 观测量 API 与约定（measurement.jl）

全部经「带插入重建 K → 同一 `integrate` → /Z」实现：

- `occupation(sol, j; branch, band)`：j≥1 折 (+,j) 边界 :left；j=0 折 ρ₀。
- `raw_correlation(sol, i, j; kind)`（rt/mixed，API 承诺 i≥j）：
  - greater = ⟨d(+,tᵢ) d†(+,tⱼ)⟩：**push 序 ddag 先、d 后**（同切位 0⁺ =
    ⟨dd†⟩ = 1−n；旧序给 ⟨n⟩）；
  - lesser = ⟨d†(tⱼ) d(tᵢ)⟩：d@(+,i)、**ddag@(−,j) side=:right**；
  - i<j 返回轮廓时序积类量（不在承诺语义内）。
- `matsubara(sol, i, j=0)`（τ）：两插入均 :left；**j==0 平移路由（M≥2）**：
  i<N → (i+1, 1)；i=N → (1,1) 等时（G(β,0)=+⟨d†d⟩ 环等时同号）。平移等价
  G(τᵢ,0) ≡ G(τᵢ₊₁,1) 在链上**机器精度**成立（M=1 direct vs shifted
  ≤5.6e-16）。τ 链 G(β,0) = +⟨d†d⟩（无时间序反号）。
- rt 等时核 (i,i) = **n(tᵢ)（演化后的占据）**——参考实现见
  `diag_u1_mixed_twopoint.jl` 的 `edn`（因子化 ρ₀ 与 H 不对易，两时关联
  必须显式演化，不能套平衡态平移不变公式）。

---

## 6. 标定陷阱清单（"看起来错其实对"——集成时禁止"清理"）

| # | 约定 | 表面现象 | 为什么是对的 |
|---|---|---|---|
| T1 | `_slot_operators` 无 graded⊗ 符号 | U1 的同函数有 | Z2 分级收缩 + `:dag` 变体已供给；再加 = 重复计入（第六轮组合扫描） |
| T2 | `oddlater_mode = :none` | U1 施加显式 oddlater | 分级收缩已供给；`:u1` 模式带 δt 无关 τ 线性漂移 |
| T3 | τ 补偿 F-gated 且无 `c_b·p_b` 项 | 与 rt/mixed compexp 形状不同 | τ 两轴路由的楔积约定与 rt/mixed 不同（等时端到端 = 基线真值）；改动前必跑等时基线 |
| T4 | fold 格 N=1 无补偿、fold R 线有逐 s0 交叉项 | 不对称 | N=1 的 Tr[A·R] 天然因子化（U=0 M2≡M1 实证）；R 线是 fold 的"键" |
| T5 | `variant=:dag`（a ↦ +d†） | U1 用 :std（−d†） | 吸收 `g_trace_phy` 奇扇区 −1；ED 标定 |
| T6 | matsubara j==0 走平移路由而非直接折 R | 多一步间接 | 热环循环迹平移等价链上机器精度；绕开 Z2 尾部（§9） |
| T7 | `KernelCorrView` getindex：a-ā → +Δ、ā-a → −Δ、同 conj → 0 | 看似多余的一层 | 镜像 QuAPI `RealCorrelationFunction` 语义（跨支负号内嵌块值，读数处补回）；合并版若改用 `GTEMPO.correlationfunction` 直连需逐位核对 |

**总原则**：任何符号/约定的改动必须过 §7 的验证电池；跨 MPS 系数比较只信
**比值判据**（`coeffs()` 经典展平对 F≠0 链不保 Z2 荷守恒，插入链禁用系数
表对账）。

---

## 7. 验证电池（回归保护，集成后必须全绿）

全部位于 `validation/`（ed.jl = ED 参考依赖）：

| 脚本 | 内容 | 当前基线 |
|---|---|---|
| `twopoint_m2_check.jl` | rt/τ M=1/M=2 全网格 vs ED | rt ✅、τ M=1 3.67e-5、M=2 1.71e-3 |
| `mixed_twopoint_check.jl` | 混合 M=1 全网格 + M=2 等时/跨格 | 6.82e-4 ✅ |
| `twopoint_band_check.jl` | 逐带 occupation + 实时等时两点 | ③ 斜率 2.002、④ 1.06e-5 |
| `rt_fastpath_check.jl` | rt 双带端到端 + 计时 | 1.34e-6@Nt=8 |
| `kernel_paths_check.jl` | QuAPI 统一路径端到端（δ + 双 δ） | ✅ |
| `u0_factor_mixed.jl` / `u0_factor_rt.jl` / `k_factor_mixed.jl` | U=0 因子化三件套 | Z₂/Z₁²=1.0；楔积 0/400 |
| `diag_j0_factor.jl` | **M2≡M1 机器锐利判据**（τ 列） | i=1/N=2 列 ~1e-16 |
| `diag_u1_mixed_twopoint.jl` | 混合/rt 三框架金标准对照 | 混合 3.7e-14/7.5e-12 |
| `diag_u1_tau_twopoint.jl` | τ 三框架对照 | U1 干净、本侧 1.7e-3（§9） |
| `continuous_bath_scaling.jl` / `continuous_bath_{u1,mine}.jl` | 连续谱收敛阶 / 对齐 | 2.000 / ≤1.4e-11 |
| `if_band_invariance.jl` / `k_factor_mixed.jl` | IF/K 带不变性探针 | 比值判据 spread≈0 / 0/400 |
| `toulouse|singlemodebath/convergence*.jl` + `plot_scaling.py` | scaling 基线 | 2.00x 全绿 |
| `calibrate.jl` / `unit_weyl_step.jl` | 单元级约定标定档案 | — |
| `diag_tail_investigation.jl` / `diag_j0_gauge3.jl` / `twopoint_m2_frontier.jl` | 尾部判决链档案 | — |

判决链档案：`HISTORY.md`（十二轮全文）、`mixed_u1_alignment.md`、
`continuous_bath_alignment.md`、`u1_tau_twopoint_alignment.md`。

---

## 8. 金标准对照状态（U1GTEMPO2，~/Documents/Missile/U1GTEMPO2）

U1GTEMPO2 是**独立栈**（直接 TensorKit，自有 GMPS/`berezin_vertex` 配对；
键标签 = `U1Irrep⊠U1Irrep⊕FermionParity` 每带电荷），实现同一 Weyl 二阶。

- **occupation**：mixed 逐位一致（≤2e-8）；连续谱 ≤1.4e-11。
- **混合两点**：机器精度逐位（3.7e-14 / 7.5e-12）。
- **rt 两点**：离散化精度内一致（Nt=8 最差误差逐位相同）。
- **τ 两点**：本侧 1.7e-3 vs U1 干净（1e-4 级）——唯一差距，见 §9。
- 对照方法论：**只做端到端观测量对比**（U1 管线约定自洽，表级对比是历史
  陷阱）；U1 lesser 默认分支 = d@(+,i), d†@(−,j)，勿覆盖为 +/+（会结构性
  变成 greater）。

---

## 9. 已知限制：τ 奇插入两点的 ~1e-3 结构尾部（Z2 grading 上限）

**现象**：τ（虚时）奇插入两点链相对 ED/U1 有 ~1e-3 结构性残差（固定 τ 不随
N 收缩，非离散化）；等时/occupation/rt/混合**不受影响**（F≡0 链 num/Z 同
几何；rt/混合的链结构不触发）。

**根因**（第十二轮定性，五路线穷尽 + IF₂ 解耦决定性实验）：(交错
Weyl-K ⊗ IF₁ ⊗ IF₂) 分级收缩的**值**与 M=1 乘积之差，跨带 bilinear
（K₁-slots × IF₂-selections）。**K 侧符号在数学上不可修**（比例失真
Σφ_{IF₂}(−1)^h 非 ±1）；收缩系统序、格点置换、IF 构造变体均逐位不变
（收缩值对表示变换不变）。

**根子**：GTEMPO 栈的单 Z2 键宇称把每带电荷 (1,0)/(0,1) 合并进同一扇区；
U1GTEMPO2 的每带电荷标签（`U1Irrep⊠U1Irrep⊕FermionParity`）天然分辨。
GTEMPO 核心 Z2Irrep 硬编码 36 处（integration/observables/grassmannmps/
lattices/sysdynamics）。

**修复规格**（集成后可选立项）：GTEMPO 核心的 grading 泛化
（Z2Irrep → 抽象 sector），完整参考实现 = U1GTEMPO2 的 `tensors.jl`
（`split_cell`/`from_core` 电荷标签拆分、`local_weight` 的 graded⊗ 符号、
显式 oddlater——注意这些在 U1 栈内自洽，移植时须连带其收缩约定）。
升级落地前，1.71e-3 为已知最优 Z2 值，验收线（2e-3）内。

---

## 10. 集成步骤建议

1. **前置**：跑通 §7 验证电池（基线数字以本文档为准）。
2. **搬运**（建议顺序）：contours.jl 的两个 ordering → tensors.jl →
   system.jl（`system_factor` 改名 `sysdynamics_weyl` 挂入 GTEMPO
   sysdynamics 命名空间）→ measurement.jl（作为 Weyl 测量 API）→
   solver 薄壳合并。
3. **kernel 层评估**：先试 `GTEMPO.correlationfunction(bath, weyl_lattice)`
   直连 `hybriddynamics`（绕开 `KernelCorrView`），与 §7 的
   kernel_paths_check 逐位对照后再决定去留。
4. **回归**：§7 全电池 + scaling 四条逐位对照本文档基线。
5. **不搬运**：contraction.jl（直接用 GTEMPO）、diag_* 判决链档案（历史
   价值，留 benchmark 或归档）。
6. **风险点**：§6 陷阱清单里的七条约定是高频"顺手清理"事故区——任何
   改动先过 M2≡M1 机器判据 + 等时基线 + U=0 三件套。
