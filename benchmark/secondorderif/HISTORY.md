# singlemodebath — 单模 Anderson 模型的 Weyl 二阶 GTEMPO 基准

基于 GTEMPO 框架（Z2 GMPS）实现 `benchmark/weyl-secondorder.html` 中的
**Weyl 二阶算法**（`weyl_tti` 方案），结构与 `~/Documents/Missile/U1GTEMPO2`
包完全对齐，但全部构建在 GTEMPO 已有机制之上。

## 算法

模型：H = H_S + H_B，H_B = H_bath + Σ_α(d†_α F_α + F†_α d_α)。单元（源）
变量 (a_j, ā_j)，每单元恰好一个全宽源单元（无额外浴端点）。

- **二阶（Weyl）**：`U₂(λ) = exp(λ²c/4)·e^{−λH_S/2}·Q_W[e^{−λH_Bᴳ}]·e^{−λH_S/2}`，
  局部权重 `M_W = qI−n + a·(系数) + ā·(系数) + aā·W₄`，q = 1/2；
- **一阶**：anti-Wick 权重（q = 1）+ Lie 全步；
- 两者**共用同一 IF**（物理浴、相同 N 个单元）；
- K 为 automaton 链：虚拟键携带 (第一成员, 第二成员) Fock 对（实时每顶点
  折叠正反支，"X → A₊XA₋"），单元 1 折叠 ρ₀、单元 N 闭合迹；
- K 侧无需 U1GTEMPO2 的 Gaussian 测度适配项之外的处理——经 ED 校准，
  GTEMPO 的收缩配对约定与 U1GTEMPO2 的 Gaussian 测度适配等价
  （aā 槽 = I/2+n），但耦合槽 a ↦ **+d†**（U1 的 −d† 差一个由
  `g_trace_phy` 分级 trace 的奇宇称 −1 产生的符号，见 `src/system.jl`）。

## 模块结构（对齐 U1GTEMPO2 的 src/）

| U1GTEMPO2        | 本目录                | 说明 |
|------------------|----------------------|------|
| contours.jl      | src/contours.jl      | N 全宽单元轮廓 → GTEMPO 格点映射 |
| kernels.jl       | src/kernels.jl       | Weyl N×N 浴核（复用 QuAPI Δt/Δτ）+ coupling_norm |
| exponentials.jl  | （并入 kernels）     | 单模浴单元积分闭式 |
| tensors.jl       | src/tensors.jl       | 单元核 → Z2 site tensors（逐变量分级精确 SVD） |
| system.jl        | src/system.jl        | `local_weight` / `system_factor`（Weyl/一阶 K） |
| influence.jl     | src/influence.jl     | IF：**复用 GTEMPO hybriddynamics**（核视图） |
| contraction.jl   | src/contraction.jl   | **复用 GTEMPO integrate** |
| solver.jl        | src/solver.jl        | Problem / Solution / solve |
| measurement.jl   | src/measurement.jl   | occupation / greater / lesser / matsubara（插入重建 K） |
| prepared_system  | （并入 system）      | 基准规模无需环境缓存 |
| scaling/persistence/hierarchical | —   | 基准范围外 |

复用的 GTEMPO 功能：`GrassmannLattice`、`vacuumstate`、`hybriddynamics`
（IF 全量复用，经 `KernelCorrView` 喂单元级核）、`integrate`、
`swapband`（多带复制）、`fock_propagator`（精确 Fock 传播子）、
`TruncateDimCutoff`。

## 运行

```bash
# 在 GTEMPO 仓库根目录
julia --project=. benchmark/secondorderif/singlemodebath/validation/convergence.jl  # 收敛阶基准
julia --project=. benchmark/secondorderif/singlemodebath/validation/calibrate.jl    # 单元校准
```

## 已验证结果

**① 虚时平衡 ⟨n⟩（β=1，单带 Toulouse + Dirac-δ 单模浴），ED 参考**
（`results/convergence.csv`）：

| order | N=5      | N=10     | N=20     | N=40     | 斜率  |
|-------|----------|----------|----------|----------|-------|
| 1     | 6.0e-4   | 2.1e-4   | 8.0e-5   | 3.4e-5   | 1.38  |
| **2 (Weyl)** | **3.8e-4** | **9.4e-5** | **2.3e-5** | **5.9e-6** | **2.00** |

Weyl 二阶收敛阶精确为 **2**（步长减半误差 ÷4）。两者共用同一 IF——
精度提升全部来自 Weyl 局部权重 + 对称半步 + 标量补偿，与推导一致。

**虚时 Matsubara G(τ,0)（N=10，d†(0) 经 boundary-0/ρ₀ 折叠）**：
G(0.2/0.4/0.6/0.8) = 0.5526/0.4800/0.4251/0.3832 vs ED
0.5526/0.4801/0.4252/0.3832（≤1e-4，全点一致）。

**实时 greater/lesser（dt=0.1，N=4，全部边界 vs ED）**：

| t   | G(t,0) mine | G(t,0) ED | L(t,0) mine | L(t,0) ED |
|-----|-------------|-----------|-------------|-----------|
| 0.1 | −0.0471−0.6642i | −0.0466−0.6649i | +0.0234+0.3298i | +0.0231+0.3302i |
| 0.2 | −0.0928−0.6585i | −0.0922−0.6551i | +0.0461+0.3270i | +0.0458+0.3253i |
| 0.3 | −0.1382−0.6513i | −0.1357−0.6390i | +0.0686+0.3234i | +0.0674+0.3173i |
| 0.4 | −0.1830−0.6424i | −0.1763−0.6170i | +0.0909+0.3190i | +0.0876+0.3064i |

**两带总 ⟨n₁+n₂⟩（U=0.8，τ，N=16）**：0.6094 vs ED 0.6118（0.4%）。

## 本轮修复的四个 bug

1. **实时核约定**：QuAPI `Δt` 的原始值是顶点约定（`_vs_rt_acc` 节点
   映射），不是 Weyl 全宽单元积分值。改为从**精确浴迹多项式**直接
   提取单元级核（单模 Dirac-δ 浴解析精确；单元演化算符
   `E = e^{lam·H_b}(I+lam·Q'+(lam·Q')²/2)`，Q' 带海森堡共轭，
   O(N²) 直取线性系数）后，实时 greater/lesser 全边界对齐。
2. **两点宇称轴路由**：内部单元的第一轴（朝晚端 bond j↔j+1）= F_j、
   第二轴（朝早端 bond j−1↔j）= F_{j−1}——此前装反；修正后奇插入
   分子按 |ys| = f_j（本格奇插入奇偶）正确筛选。
3. **虚时折叠单元的 R**：τ 折叠单元此前硬编码 R = I——boundary-0
   插入（d†(0)）与 ρ₀ 从未进入虚时链。
4. **多带 IF 复制**：`influence_functional` 用 `mpsI = swapband(...)`
   覆盖了带 1 的 IF（第二带耦合丢失）。改为 IF 列表 + 变参收缩。

## 开放项（未对齐，诚实记录）

1. **实时单带：已完全闭合（2026-11-05）**。二阶收敛阶 = **1.996 ≈ 2.00**
   （err 1.9e-6@δt=0.05）；occ/gt/lt 全时刻 1e-6/1e-4 级 = U1 端到端同级。
2. **③ 两带 τ U≠0：已闭合（2026-11-05 第六轮，见下方修复记录 24–27）**。
   二阶收敛阶 = **2.004**（err 7.45e-5@N=16 vs ED 0.6117976）；order 1 = 1.61。
   U=0 因子化检验：Z₂/Z₁² = 1.0000000000（order 1 偏差 8e-15）、两带 K₂ 与
   K₁∧K₁ 逐支撑一致（order 1 摊散 3.3e-15、order 2 8.7e-14，含 Weyl 标量）。
3. **实时等时两点 0⁺ 约定 + 逐带 ⟨n_b⟩ 插入：已闭合（2026-11-05
   第七轮，见修复记录 31–33）**：等时核 = ED 0⁺ 值 ≤7.6e-6（t=0 边界
   精确 1e-15）；τ 逐带 ⟨n₁⟩/⟨n₂⟩ = ED 0.3058988、收敛斜率 **2.002**；
   实时逐带 + 等时两点 ≤2.65e-4（Nt=2）。
4. **实时两带**：正确性已确认（第六轮）：N=2 err 3.8e-4、N=3 err 2.3e-4
   （t=0.4，order 2）。**性能开放项**：`_system_factor_rt` 的 `permute!`
   （TEMP→格点 swap 门链）在 M≥2、N≥4 显著变慢（预存问题，与符号修复
   无关；τ 路径不用 permute!）——见 `validation/realtime_tb_sanity.jl` 头注。
5. **τ N=1 退化闭合已修**（修复记录 10）；单带 τ N≥2 基线不受影响。

## 配对层判决史（防走回头路）

- **pairing_G ≡ pairing_u1（2026-11-05 实证，推翻 11-04 的相反判决）**：
  用合成单项式 MPS 逐 (SK,SI) 对直接测量双方 integrate 的符号表
  （`gtable_measure.jl` / `gtable_u1.jl`）——313/313 个非零对**逐位相同**
  （最大差 0.0）。此前"配对差"判决是被 IF 跨支符号差污染的间接推断。
- 部件级对齐清单（全部实测）：核 = U1 逐位（±1 映射）；IF ≡ U1 逐位
  （Nt=2 表 66/66，修复 17 后）；K（无插入）≡ U1（66/66，reorder_sign
  补偿后）；K（插入 n̂）≡ U1（41+41 全 +1）；配对 ≡ U1（313/313）。
- 交叉核对工具：`gtable_measure.jl`（GTEMPO 侧 G 表 + 单项式合成 MPS），
  `gtable_u1.jl`（U1 侧同表），`kcoeff_ins.jl`（插入 K 表，SIDE 环境变量
  两侧复用），`kpermute_wedge.jl`（置换楔检验——置换无罪）。

## 2026-11-04 会话修复记录

1. **τ Matsubara 恒零**：折叠单元键标签把 boundary-0 奇宇称 F₀ 重复计入
   （f_j 并入后又经 F_before⊕f_j 抵消）。修复：f_j 保持纯净。
2. **oddlater 双重计入（τ）**：分级收缩已内部供给，U1 式显式符号再加
   = 双重计入（退耦标定 6 位一致）。已删。**实时路径同样标定：:none**
   （退耦 gt/lt/occ 全部 1e-14/15；:incl/:strict 产生 (−1)^k 交替）。
3. **实时核浴演化指数**：Ebath = e^{−λH_b}（原 + 反）；F' 共轭同步翻转。
4. **实时核热浴权重**：轮廓迹以 ρ_b = e^{−βh_b} 起始、tr(e^{−βh_b}) 归一。
5. **单元多项式精确化**：线性 −λ∫F(s)ds = −F(1−e^{−λω})/λω；二次有序三角。
6. **p-h 镜像浴混淆（重大澄清）**：此前 U1 对比误用 DiscreteBath(−ω)
   （p-h 共轭浴，⟨n_b⟩ = 0.731）。正确参照 (+ω)：
   **U1GTEMPO2 端到端 = ED 至 1e-5（全时间点）= 观测量级金标准**。
7. **QuAPI Δt 再定性**：其值为单元级全宽积分（_gⱼₖ 矩形/_gⱼⱼₖₖ 三角）；
   且 GTEMPO 读取经 `_vs_rt_acc` 顶点位移（首行零 = 边界顶点 ā 腿无
   配对）。GTEMPO 一阶 = 顶点式架构；U1 = 单元式架构。
8. **B′ 存储约定（现行）**：U1/QuAPI 存储符号经本管线补偿，早期时间
   精度 2.4e-3 → 2.1e-5。

## 2026-11-04 第二轮（IF/K 系数交叉核对）修复记录

9. **实时接触项（有序三角）物理 bug（重大）**：`realtime_kernel_trace`/
   `_cell_E` 的同单元 aā 二次项此前用宗量 `λω·dt` + f/g 通道配对相反
   （"经验标定配置 A"）——**IF 单项式系数交叉核对定位**：核 diff 显示
   非对角（矩形）条目与 U1 逐位一致（±1 存储映射），唯一数值差异 =
   对角接触（虚部差 dt 倍且反号：ours −8.3e-6 vs U1 +8.3e-5）。修复：
   宗量 = `λω`（= U1 的 e·li；∫₀¹du(1−u)e^{±λωu}）、f 通道 tri(+λω)、
   g 通道 tri(−λω)。修复后核与 U1 mode_block **逐位一致**（64 条目
   全部 ±1 映射，偏差 ≤ 1.1e-14）。旧标定之所以"早期 1e-5 看起来对"，
   是错接触 × 补偿翻转恰好抵消的巧合。
10. **τ N=1 退化闭合**：`_system_factor_tau` 在 `j==1==N` 时两侧都真空
    闭合，只读 A[1,1]，热迹丢失（n(β) 直接错）。改为全 Fock 迹
    Tr[A·R]。N≥2 不受影响（闭合格的 outsize=pairdim 本就是全迹）。
11. **存储开关**：`realtime_kernel_trace(; storage=)`：`:B′`（默认，
    K 约定所配）与 `:u1`（:B′ 复合自逆映射 S：同支 i<j 翻、跨支 i>j
    翻 = U1 mode_block 原始存储）。**端到端标定实验**：`:u1` 直配产生
    均匀 +1e-2 偏移——K 侧 oddlater/槽适配是围绕 B′ 标定的，不可单改
    IF 存储（补偿耦合）。
12. **交叉核对基础设施**（validation/）：`if_coefficients_mine.jl` /
    `if_coefficients_u1.jl`（双侧全单项式系数导出；GTEMPO 侧为键扇区
    行走批量 DP，振幅 = 真空归一后的链收缩——canonical form 的 scaling
    按 `n^(1/L)` 几何均值吸收范数，静态推导不可靠）、
    `if_coeff_rowprobe.jl`（单行 partial-IF 探针：提取器与读数规则
    0.0e+00 吻合）、`if_coefficients_compare.jl`（变量双射：v =
    4(Nt−j)+(支偏移)+(槽偏移)，两侧逐位同序）、`kcoeff_*.jl`（K 侧
    同流程，全局尺度拟合比较）。
    **关键结论**：IF/K 两个多项式都与 U1 精确 ±1 对齐（无复数散布）；
    GTEMPO canonical MPS 的"振幅 = scaling × 原始链收缩"不成立，
    真空归一是唯一可靠途径。

## 2026-11-04 第三轮（约定穷举 + 配对差实证）记录

13. **IF 存储切 :u1（默认）**：B′ 三角翻号破坏单模浴核的秩 1/Toeplitz
    结构——IF 四阶 Pfaffian 扇区结构性偏离（多出 1429 个 U1 精确为零的
    实质单项式）；:u1 保持物理核的符号结构。KernelCorrView/partialif
    读取路径不变（跨支 Keldysh 补偿与 QuAPI 存储约定相容）。
14. **K 格间重排补偿（reorder_sign，默认 true）**：K 系数表实证本构造
    相对 U1 多出 (−1)^{C(m,2)}（m = 奇选取格数）——已按
    Σ_q [n_q odd]·[inpar] 在核元素循环逐格补偿；补偿后 K ≡ U1（Nt=2
    系数表逐位 +1）。`system_factor(; reorder_sign=false)` 恢复旧行为。
15. **pairing_G ≠ pairing_u1 实证**：U1 多项式原样（IF=:u1 + K≡u1）喂
    GTEMPO 管线 → δt 完全平坦偏移、第一格边界即错；而 U1 自己的配对 +
    同一组多项式 = ED。至此剩余差异**唯一**位于配对层符号语义。
16. **permute! 楔检验（kpermute_wedge.jl）**：TEMP→格点置换 = 楔正确
    重展开 × 全局常数（91/91 单项式均匀）——置换层洗清。注意钩子必须
    `deepcopy`（permute! 原位变更，引用型钩子会得到同一对象）。

## 2026-11-05 第四轮（实时闭合）修复记录

17. **IF 跨支 Keldysh 预翻符号（最后一缺口）**：`KernelCorrView` 读数
    预翻跨支符号，抵消 `partialif_hybrid` 消费端的 QuAPI 式跨支 (−1)
    补偿——净读数 = U1 核原始值。Nt=2 IF 表实证：修复前 32/68 反号
    （恰好全部跨支二次项），修复后 **66/66 全 +1**。修复来源链：G 表
    逐对测量证明配对逐位相同 ⟹ 排除配对层 ⟹ 插入 K 表证明 ≡ U1 ⟹
    唯一未直比过的部件 = :u1 下的 IF 多项式 ⟹ 直接对比暴露跨支符号。
    修复后：管线 Z = 0.99999（此前 0.98325），实时二阶阶 = 2.00。
18. **验证状态（最终）**：实时 occ err [3.0e-5, 7.6e-6, 1.9e-6]
    （δt=0.2/0.1/0.05，阶 2.00）；N=4 序列 occ/gt/lt 全时刻 1e-6/1e-4 级；
    τ 二阶 2.002 不变；两带 τ 基线不变。order=1 实时斜率 0.82
    （N=2 已受二阶项污染，合理）。

## 2026-11-05 第五轮（双带部件级对齐进行中）记录

19. **`_slot_operators` M≥2 槽矩阵 bug（重大）**：旧实现用 2-bit 槽码索引
    `local_weight` 的 4^M 组合码列表——单带恰好重合、M≥2 取错矩阵且重复
    计入 qI−n 因子。重写为从每带 4 槽基矩阵直接构造（结构 = U1
    local_weight 的 locals×graded-sign）。签名改为
    `_slot_operators(M, order; mode)`。
20. **JW 约定平移**：`_fock_annihilators` 符号掩码改为「带 b 穿越更低索引
    带 b′<b」（band1 无符号、band2 带 (−1)^{n₁}）= U1 fock_operators；
    旧实现穿越更高索引带（恰好反转）。单带不变。
21. **τ 尾格楔补偿**：`_system_factor_tau` 闭合格（j=N）乘
    (−1)^{popcount(槽码)}（M≥2）→ 双带 K ≡ U1（36/36 全 +1）。
22. **`KernelCorrView` τ 跨格预翻（M≥2）** → 双带 IF₁ ≡ U1（6/6）。
    注意：单带 τ 无此差——预翻会破坏单带基线（实证 3.6e-2），条件必须
    限定 M≥2。
23. 工具：`tb2_compare.jl`（双带 τ 三表，SIDE 环境变量两侧复用）、
    `gtable_measure.jl`/`gtable_u1.jl`（配对符号表直测）、
    `kcoeff_ins.jl`（插入 K 表）、`u1w_decomp.jl`、`basis_check.jl`、
    `realtime_dtconvergence.jl`。

## 2026-11-05 第六轮（两带 τ 闭合）修复记录

**方法论转折**：第五轮的"最后缺口（partialif M≥2 多格乘积楔）"被证明是
**比较器伪影**——跨框架系数表对比时，mine 侧 dump 把 site 位重排到 mu 约定
（cell 升序）但 site 序是 cell 降序，每跨 2 格支撑差 (−1)^{C(k,2)} 宇称因子；
第五轮为迎合它引入的补偿（IF τ 跨格预翻、K τ 格间楔补偿）全部是假象补偿。
本轮全部改用**约定无关的本征检验**：site 基解析指数引擎（20-bit 外代数）、
U=0 因子化恒等式、配对结合律、逐支撑 Python 参照链。

24. **IF 构造实证精确（wedge_probe.jl）**：每行 partialif_hybrid、逐级
    mult! 累积、完整 per-band IF 全部 = e^{+Δ}（0/66 失配；单带对照同样
    全对）。U1 表实证 = e^{−Δᵀ}（mu 约定 66/66）——两框架互为**取向反转**
    管线，各自 internally 一致；单带 τ 喂 −Δᵀ 平坦偏差 1.5e-2
    （tau_orientation_verdict.jl），证明 e^{+Δ} 是本管线唯一正确取向。
    ⟹ **KernelCorrView 的 M≥2 τ 跨格预翻移除**（实时跨支 Keldysh 翻保留）。
25. **K τ 格间楔补偿移除**：(−1)^{Σ_{j≥2}Nⱼ}（M≥2）与 IF 预翻同批引入、
    同源。移除后 U=0 因子化恢复（0.9925→见 27 后 1.0）。
26. **Weyl 标量改为总耦合**：sc = exp(λ²·Σ_b α_b/4)（×bands）——两带
    K₂/K₁∧K₁ 曾差 exp(−Nτh²α/4) = 0.96923（逐支撑实证）。
27. **_slot_operators 移除带间 graded⊗ 符号 (−1)^{|a||b|}（根因）**：
    U1 的该符号属于其 Gaussian 测度配对 + 取向反转管线（其 IF = e^{−Δᵀ}），
    照搬进本仓库裸 Berezin 管线 = 符号重复计入。JW 穿越
    （_fock_annihilators，带 b 穿越更低带）保留。正确组合由
    k_factor_probe.jl 的 4 组合扫描实证：graded=off + jw=on ⟹ 联合链
    vs K₁∧K₁ 逐支撑 0 失配。M=1 该符号本就不生效，单带基线无损。
28. **配对层结合律实证（pair_assoc_probe.jl）**：integrate(K,I₁,I₂) =
    integrate(K∧I₁, I₂)（20-site 交错支撑，差 0.0）。
29. **验证终态**：③ 两带 τ U≠0 二阶 = **2.004**（err 7.45e-5@N=16）、
    order 1 = 1.61；U=0 同阶 2.004；单带 τ 2.002、实时 1.996 回归无损；
    Z₂/Z₁² = 1.0（order 1 偏差 8e-15）。新工具：`wedge_probe.jl`（IF 逐级
    site 基探针）、`z_factor_probe.jl`（U=0 因子化）、`k_factor_probe.jl`
    （K₂ vs K₁∧K₁ 逐支撑）、`tau_flip_verdict.jl`/`tau_orientation_verdict.jl`/
    `tau_u0_verdict.jl`（ED 裁决）、`pair_assoc_probe.jl`（配对结合律）、
    `results/delta_tau4.csv`（τ 核导出）。
30. **跨框架对比教训（重要）**：U1GTEMPO2 与本仓库是取向反转的两套自洽
    管线（IF e^{−Δᵀ} vs e^{+Δ}；其表比较必须做置换宇称校正 + 约定显式化）；
    部件级"≡ U1 全 ±1"式判决若比较器不带宇称校正即不可靠。本征锚 =
    解析指数引擎 + 因子化恒等式 + ED 端到端。

## 2026-11-05 第七轮（等时两点 0⁺ + 逐带占据）修复记录

31. **实时等时两点 0⁺ 约定（greater 同切位 push 序）**：`raw_correlation`
    的 greater 原先 push [d, ddag]——`cutop` 后 push 的矩阵乘在左 = 轮廓
    更晚 ⟹ 同切位（i==j，含 j=0 折 ρ₀）给出 d†·d = ⟨n⟩；0⁺ greater 应取
    d(t⁺) 更晚 ⟹ 改为 push [ddag, d]（等时核 = ⟨d d†⟩ = 1−⟨n⟩）。lesser
    两腿异支（(+,i)/(−,j)）不同切位，不受影响。实证：修正前等时 greater
    核 = n(t)（0.33056@t=0.2），修正后 = 1−n(t)（0.66944，ED ≤7.6e-6）；
    t=0 边界精确到 1e-15；不等时 gt(t,0)/lt(t,0) 回归无损（≤1.55e-4）。
32. **逐带 occupation**：`occupation(sol, j; band=b)` 插入 n̂_b = d†_b d_b
    （缺省仍为总占据）。τ 两带 ⟨n₁⟩ = ⟨n₂⟩ = ED 0.3058988
    （1.49e-4@N=8 → 3.73e-5@N=16，斜率 **2.002**）；实时逐带 ≤2.65e-4
    （Nt=2）。第三轮已对齐的插入 K 表（n̂₁@τ ≡ U1 41/41）在此兑现为
    ED 级观测量。
33. **occupation 的 j=0 边界接线**：旧实现 j=0 返回空插入（恒 1，占位）；
    改为 n̂ 折进 ρ₀（tau：R = n̂·R；实时：cutop(0)·ρ₀，机器早已在位）。
    实证：τ @j=0 = @j=8 = 0.3058616（循环一致）；实时 @t=0 = 0.2888640
    vs 因子化初态精确值 0.2886863（Nt=2 离散 1.8e-4）。
    验证工具：`validation/twopoint_band_check.jl`（①②③④ 四节全绿）。

## 2026-11-06 第十轮（M≥2 混合/实时工件闭合：每带 XOR 补偿）修复记录

**判决链（全程本征锚，无跨框架比较）**：

34. **IF 判决（无罪）**：if_band_invariance.jl 复跑 + 升级为比值判据——
    双带 IF₁（band-1 嵌入）vs 单带 IF 的 20 个可比单项比值 spread = 3.3e-15
    （机器精度）⟹ IF₁ ≡ IF_single × 全局标量，标量在 num/Z 的 vacuum 归一
    收缩里相消。第九轮"~2e-4 结构性偏移"是探针误读：绝对容差 1e-9 撞上
    coeffs() 的 scaling(mps) 几何均值归一（两个 MPS 的标量本不可比）。
    （教训：canonical MPS 的全局标量不进物理；跨 MPS 比较必须用比值判据。）
35. **U=0 因子化判决（K 有罪）**：diag（现 u0_factor_mixed.jl）——U=0 双带
    mixed 的 Z₂/Z₁² ≈ 0.9960（K 级差 ~3.5e-3 且不随 N 衰减）、⟨n₁+n₂⟩ vs
    2×单带 ~1.8e-3 不收敛；逐带对称 |n₁−n₂| ~ 1e-16（排除带不对称）。
36. **逐支撑楔积判决（定位到补偿形式）**：k_factor_mixed.jl——K₂ vs K₁∧K₁
    （graded 楔积，inv = 跨格 band1→band2 反转数）：96/400 单项**纯反号**，
    且全部（仅）deg₁≥2∧deg₂≥2。代数解释：总宇称补偿指数 (c₁⊕c₂)(p₁⊕p₂)
    相对每带 XOR c₁p₁⊕c₂p₂ 恰差交叉项 **c₁p₂⊕c₂p₁**（只在两带都选中时
    非零）——M=1 标定的 paritybits(m)·inpar[l] 在 M≥2 的系统性缺口。
37. **修复（system.jl，mixed + rt 两路径）**：格间重排补偿改每带 XOR
    exponent = ⊕_b c_b·(p_b(l_b) ⊕ F_b)；c_b = 带 b 槽位奇偶、p_b = 带 b 线
    对宇称、F_b = 每带 prefix 奇插入宇称（新 _isoddband；⊕_b F_b = 总 F₀，
    M=1 精确退化为旧形式 ⟹ 单带基线逐位不变）。occupation 的 n̂ 类偶插入
    F_b ≡ 0；奇插入（greater/lesser）的每带路由随修复自动成立。
38. **验证终态（全部绿）**：
    - k_factor_mixed.jl：**0/400**；
    - u0_factor_mixed.jl：Z₂/Z₁² = 1.0、vs 2×单带机器精度、order 2 斜率
      **2.020**、|n₁−n₂| ≤ 5e-15；
    - u0_factor_rt.jl：rt 同款全绿 + ED 绝对值（杂质真空⊗热浴初态）二阶
      斜率 1.993（ED 初态的 kron 因子序 = [bath2,bath1,band2,band1]，ed.jl
      杂质在低位）；
    - Anderson mixed U≠0：order 2 斜率 0.641 → **2.004**（err 6.4e-5@N=16），
      与 U1GTEMPO2 的 4 点（order 1/2 × N=4/8）全部逐位一致；
    - rt 两带 U≠0：斜率 1.927 → **2.001**（err 8.3e-8@Nt=32；rt_fastpath
      N=2 3.8e-4 → 2.1e-5）；
    - 连续谱浴 U1 对齐复测：实时 4.8e-5 → **≤1.4e-11**（12 点全机器精度；
      第九轮"指数拟合容差"实为 rt 工件）；
    - 单带回归逐位无损：Toulouse τ 2.002 / rt 1.996 / mixed 2.005；
      twopoint_band_check ①②③④ 全绿（④ 2.65e-4 → 1.06e-5）；
      kernel_paths_check 全绿。
    - 工具：validation/{k_factor_mixed,u0_factor_mixed,u0_factor_rt}.jl。

## 2026-11-06 第十轮补遗（两点函数全面闭合：ysv 带交换 + 双向跨带交叉补偿）

39. **混合 − 支插入切位映射 off-by-one（修复）**：`(−, j)` 的 cut 由
    `Nτ+2Nt−boundary+1` 修正为 `Nτ+2Nt−boundary`（与 rt 已验证公式逐位对应；
    旧式把 (−,0) 落到 L+1 被静默丢弃 ⟹ lesser(t,0) 全错）。转折点约定：
    τ 边界 Nτ ≡ (+,0)（cut Nτ）；− 边界 Nt ≡ (+,Nt)（cut Nτ+Nt）。
40. **重大根因（rt/mixed 装配）：ysv↔site 带映射互换（修复）**。装配循环
    `ysv = (m >> (v−1)) & 1` 直接把 slot code 的低位当 site 低带——而 slot
    code 带序 = band 1 高位（bits 2(M−b), 2(M−b)+1 ↔ band b）⟹ **M≥2 两条
    带在装配中被互换**。τ 路径走 `_code_from_ys`（正确换算）无此病；M=1
    退化不变。**带对称性掩蔽**：U=0 因子化 0/400（wedge 对带交换不变）、
    occupation/逐带/连续谱对齐（对称观测量）全部照常通过——bug 潜伏三轮，
    直到带插入 d₁⊗I（两点函数）打破带对称才暴露。
41. **补偿缺失双向跨带交叉项（修复）**：ysv 修复后楔积对账暴露剩余 ±1 符号
    （96/200，GF(2) 拟合证实非逐 pair 结构）。最终正确补偿（rt/mixed 同修）：
    exponent = ⊕_b c_b·(p_b(l_b)⊕F_b) **⊕ ⊕_{b≠b'} c_b·p_{b'}**
    （双向跨带交叉项 = 楔积 inv 的逐元素实现：band b 变量 × band b' 线宇称；
    之前「每带 XOR 修复得 0/400」实为带交换 bug 与补偿的巧合抵消——修复
    ysv 后无插入对账重新失配，交叉项补上才真正归零）。M=1 退化不变。
42. **判决链与三次自我纠错（重要教训）**：
    - 「coeffs() 对 F 偏移 MPS 失效」**撤回**：单插入提取为空 = 奇算符迹的
      费米宇称恒零（Tr[ddag·ρ]=0），是正确行为；双插入提取可靠；
    - 「无插入 0/400 = 因子化成立」**撤回**：0/400 是带交换 bug 与补偿的
      巧合抵消（wedge 对带交换对称）——修复 ysv 后对账重新失配；
    - 「rt α→0 常数路径病」**撤回**：ED 参考初态用错（rt 物理初态 = rho0
      杂质真空⊗热浴，非耦合热态；mixed 有 τ 热化段才用耦合热态）——修正
      参考后 rt α→0 机器精度。
    - 「α→0 精确 ⟹ 装配正确」**撤回**：α→0 值只测 mask-0 常数项（slot
      补偿不进入），不构成对补偿机制的检验。
    - 「Z2 表示原理性冲突、需 Z2×Z2」**撤回**：现有 Z2 框架内的符号簿记
      （双向交叉项）即为完整修复。
    - τ 窄口判决：τ 无插入楔积 96/400 反号 = τ 两轴路由的楔积约定与
      rt/mixed 不同（非 bug——等时端到端证明「缺 rt/mixed-inv 版」才是 τ
      基线真值）⟹ τ 修复 = **F 区门控**跨带补偿（交叉项仅施加于奇插入
      后 F≠0 的格，保持第六轮无插入/等时基线逐位不变）；fold-R（j==0）
      的逐 s0 修正执行但端到端无效（fold R 线符号流未解）——唯一窄口。
43. **验证终态（全部 ED 端到端）**：
    - mixed：M=2 (i≥j) 全网格 greater/lesser（含 (+,0) 转折点、− 支插入）
      6.8e-4（二阶离散水平，N=8 等时 1.5e-4 收敛）——**完全闭合**；
    - rt：α→0 机器精度 + 全网格（因子化热 ρ0_imp2 初态）Nt=4 2.98e-3 →
      Nt=8 1.54e-3（二阶收敛）——**完全闭合**；occupation 四象限
      （{因子化热,真空}×{U=0,U≠0}）全部离散水平；
    - τ：F 区门控跨带补偿——等时 3.73e-5（恢复第六轮基线）、j≥1 两点
      大部 ≤1e-3（(8,4)/(8,2)/(6,2) 全绿；N=16 (15,1) 1.58e-3——大 Δτ
      尾部残留小残差，非干净二阶，如实记录）；
    - τ 窄口：**j==0（ddag 折入 fold 的 R）**——基于 p₂ 的 fold 逐 s0
      修正**恒 0 无效**（fold 非零 R 项只有 s0=3，其 band-2 位 = 0——
      修正不在 fold 元素的 p₂ 结构上），matsubara 对 j==0 显式抛错
      （force 绕过）；M=1 全网格 + M=2 等时/j≥1 可用（N=16：(15,1)
      1.58e-3，j≥1 域收敛确认）；
    - 全回归逐位/收敛无损：Anderson mixed 2.004（逐位）、Toulouse mixed
      2.005、τ 2.004/2.002、rt 两带 2.001；twopoint_band_check ①–④ 全绿
      （④ 1.06e-5）；rt_fastpath 全绿；楔积 plt 304/96（τ F 门控基线形态，
      等时端到端验证正确）；
    - API 语义：raw_correlation/matsubara 承诺区间 i≥j（i<j 返回轮廓时序
      积类量）；matsubara 拦截 j==0。
    - 工具：validation/{twopoint_m2_check, mixed_twopoint_check,
      twopoint_m2_frontier}.jl。
    - 工具：validation/{twopoint_m2_check, mixed_twopoint_check,
      twopoint_m2_frontier}.jl。

## 2026-11-06 第十一轮（j==0 窄口闭合：fold R 线修复 + 热环平移路由）记录

1. **判决链（U=0 因子化判据全程锚定）**：① 复现：j==0 列 M=2 err 2.6e-2，
   **U=0 同量级 ⟹ 与传播子无关，纯补偿/F 路由结构病**；band-2 与 band-1
   逐位同差（带对称）。② **M2−M1 = 全部误差**（diag_j0_factor.jl：M=2(U=0)
   j==0 列减 M=1 列 = 逐位等于 M=2 的 ED 残差）⟹ U=0 因子化在奇插入 j==0
   拓扑被破坏（等时/occupation 的偶插入折入不受影响——既有检验盲区）。
   ③ N=1 退化格 M2≡M1 逐位相等 ⟹ 病不在 fold 矩阵代数（Tr[A·R] 天然因子化）。
2. **coeffs() 经典展平对 F≠0 链不忠（工具教训）**：融合树展平求和不约束
   Z2 荷一致，F 移位链的中间键可有不一致路径混入——插入链的系数表级对账
   不可信（diag_j0_gauge3.jl 的「结构失配」实为提取器工件）。gauge 测量
   （无插入/等时）仍可靠：χ₀ = c₁@k·c₂@k′（k<k′，96/400 复现第十轮记录）。
3. **fold R 线修复（V1，system.jl）**：fold 输出键真空 ⟹ 键公式交叉项退化
   恒 0；「键」在 fold 处的对应物 = R 线 (s0→tL)——交叉项
   ⊕_{b≠b'}c_b·p_{b'}(s0,tL) 拆开 fold 求和逐 s0 施加（F 项 slot-only 仍由
   公共块施加）。判决：i=1（d 也在 fold 格）拓扑 U=0 M2−M1 达 **1e-15/16
   机器精度**（N=2/4/8）；j≥1/等时/M=1 基线逐位无损（twopoint_m2_check
   基线 1.58e-3 复现）。
4. **热环平移路由（measurement.jl）**：Matsubara 环的循环迹平移等价
   G(τᵢ,0) ≡ G(τᵢ₊₁,1)——**链上机器精度实证**（M=1：direct vs shifted
   d−s ≤ 5.6e-16 全网格，含 i=N ↔ 等时点 (1,1)）。M≥2 的 j==0 路由到
   (i+1,1)（i<N）/ (1,1)（i=N，G(β,0)=⟨d†d⟩ 环等时同号），全程落在
   F₀=0 已验证扇区；fail-loud 防守移除。
5. **闭合终态**：twopoint_m2_check **M=2 全网格（含 j==0）= 1.71e-3
   （(15,0)）✅**——与第十轮 j≥1 扇区的已知尾部水平（(15,1) 1.58e-3）
   同层。全回归逐位无损：U=0 三件套（Z₂/Z₁²=1.0、楔积 0/400）、mixed/rt
   两点、twopoint_band ③ 2.002/④ 1.06e-5、rt_fastpath、四条 scaling
   （τ 2.002/rt 1.996+2.001/mixed 2.005+2.004）。
6. **新发现（下轮入口）：内部格 ~1e-3 结构性尾部**。U=0 机器锐利判据
   （M2≡M1）暴露：所有含被门控内部格的 τ 两点拓扑（(·,1) 路由扇区、
   F₀-起源直接扇区、第十轮已验证的 j≥2 扇区的 (15,1) 大 Δτ 尾部）共享
   ~1e-3 级缺口，固定 τ 下不随 N 收缩（非离散化）。14 特征 GF(2) 邻域
   穷举（diag_j0_phiscan.jl，Hamming ≤2 from cur=385）**无命中** ⟹
   缺陷不在逐元素符号层——候选：内部格 Z2 扇区路由 / Fb 更新时序 /
   fold-cutop 格（j==1 拓扑）的补偿缺位。修复后 M=2 τ 两点应全面达
   M=1 的 5e-5 水平。变体扫描方法论教训：neg 变体（全交叉反号）≡ 每元素
   常数 −1 ⟹ 在 num/Z 相消 = 无效变体；扫描设计前先做不变性分析。
7. 工具：validation/diag_j0_{repro,factor,wedge,gauge,gauge3,phiscan}.jl
   （factor = U=0 M2≡M1 判据主工具；其余为判决链档案）。
8. **连续谱浴收敛阶基准（遗留项②，validation/continuous_bath_scaling.jl）**：
   模型同 continuous_bath_{mine,u1}.jl（cos² 谱 ∫A=0.5、Anderson U=0.8）；
   τ（N≤64）+ rt（N≤32）× orders 1/2，Richardson 自洽参考。终态：
   τ order-2 逐对斜率 **2.000 精确**、order-1 → 1.000；rt order-2 → 2.000、
   order-1 → 1.000；两轮廓 Richardson 参考互洽（~4e-7）；全部数值与既有
   continuous_bath_mine.csv 逐位一致。连续谱基准至此**闭合 ✅**。

## 2026-11-07 第十二轮（内部格结构尾部：机制定性 + K 侧不可修性证明）记录

1. **判决总表**（U=0 机器锐利判据 M2≡M1，diag_j0_factor.jl）：
   - K 侧符号三连败：① 14 特征 GF(2) Hamming≤2 邻域穷举（106 候选）
     无命中；② 32 结构变体（c₁q₂/c₂q₁/c₁ρ₂/c₂ρ₁ × 门控）全部聚在
     1.9~3.6e-3 平台；③ 解析 χ₀ 形式（bare-joint = 楔积×(−1)^χ₀ 的逐格
     键化，χ₀ = ⊕_{k<k′}c₁@k·c₂@k′）**过矫正**（残差反号 ~+2e-3）
     ——「F≠0 时 bare-joint 与楔积之差非 slot-only」是 χ₀ 推导失效根因。
   - 表示变换三连败：④ 收缩系统序 (K,IF₁,IF₂)↔(K,IF₂,IF₁) 逐位相同；
     ⑤ 格点置换 changeordering(BandMajorTime)（GTEMPO 精确 swap 网络）
     **逐位相同**（收缩值对顺序不变——范畴论必然，配对顺序不是自由
     约定）；⑥ IF₂ 构造 swapband（上游模式）vs 逐带直构——逐位相同
     （2026-11-05「swapband 零真空振幅」弃用判决与此互证：两条 IF₂ 等）。
   - **决定性实验：IF₂ 解耦**——IF₂→真空 MPS ⟹ M2−M1 = **8.9e-16 机器
     精确**（且与补偿无关：IF₂ 真空使全部 c₂ 补偿项死掉，纯 F 项恰好
     正确）。失真 = 张量网络**值**的属性，绑定 IF₂ 物理内容。
2. **判定**：尾部 = (交错 Weyl-K ⊗ IF₁ ⊗ IF₂) 分级收缩的值与 M=1 乘积之
   差，跨带 bilinear（K₁-slots × IF₂-selections）结构。**K-only 符号在
   数学上不可能精确修复**（比例失真 = Σ_{mIF₂}φIF₂(−1)^h(mK,mIF₂) 非
   ±1 函数）。也非 K-local：等时/occupation（F≡0，num/Z 同几何）天然
   精确——病只在奇插入链的 num/Z 失配几何上显形。
3. **上游对照**：GTEMPO 一阶形式论（sysdynamics K）同物理量
   （μ=0.7,U=1.0,α=0.5 双带 τ 两点，δτ=0.1）relerr = **2.9e-4**（纯
   离散化，干净）——但其 K 是另一形式论（anti-Wick、非 Weyl 联合
   automaton），不可移植；上游自身多轨道 τ 两点验收线 = relerr<1e-2。
   本基准 1.7e-3（N=16 全网格）为已知最优二阶 Weyl 值。
4. **修复路径（未来，均为大工程）**：① GTEMPO 分级融合约定（update_pair/
   compensate_twists/g_trace_phy 层）的多系统 surgery；② 二阶 K 的带顺序
   重构（U 交互单元 Trotter 化为带间链 + 显式 swap 网络）。判据已永久
   化：diag_j0_factor.jl（M2≡M1）+ diag_tail_investigation.jl（本判决
   链档案，5 条路线的可复现脚本）。
5. **方法论教训**：① 「收缩顺序/表示变换修复符号」类假设先用**逐位
   不变性**（permutation invariance）检验——值不变则整条路线秒杀；②
   单变量解耦实验（IF₂→真空）的判决力远大于盲扫——先设计判别实验再
   扫描；③ χ₀ 解析推导务必先验证前提（F≡0 的测量外推到 F≠0 是本轮
   过矫正根因）。

## 2026-11-07 第十二轮补遗（金标准对照：U1GTEMPO2 无此病；根因最终定位）记录

1. **判决（diag_u1_tau_twopoint.jl，三框架同模型对比）**：U1GTEMPO2 的
   τ 两点**无结构尾部**——U=0 双带列与单带离散化地板逐位相同
   （4.805e-5…1.468e-4），U=0.8 全列 2.0e-5~1.65e-4 干净；本基准同列
   3.0e-4~2.9e-3（尾部）。**金标准自身不受影响**。
2. **根因最终定位（修正 7782c0b 的「收缩引擎架构固有」表述）**：两框架
   **不共用收缩引擎**——U1GTEMPO2 是独立栈（直接 TensorKit，自有
   GMPS/berezin_vertex 配对），键标签 = `U1Irrep⊠U1Irrep⊕FermionParity`
   （**每带独立粒子数电荷**）；本基准的 GTEMPO 栈 = 单 Z2 宇称，
   (1,0)/(0,1) 合并同扇区 ⟹ 跨带融合符号错标 = 尾部。K 侧符号不可修
   证明不变，但正确修复 = **grading 升级**（非符号）。
3. **升级规格已知**：GTEMPO 核心 Z2Irrep 硬编码 36 处
   （integration/observables/grassmannmps/lattices/sysdynamics），泛化
   = 核心重构；完整参考实现 = U1GTEMPO2 的 SpinSector 方案
   （tensors.jl `split_cell`/`from_core` 电荷标签拆分、`local_weight`
   的 graded⊗ 符号 `(−1)^{popcount(a)·popcount(b)}`、显式 oddlater）。
   在升级落地前，本基准 1.71e-3 的 τ 两点为已知最优 Z2 值，与金标准
   的差距（~1e-3，奇插入链）记录为已知限制。
4. 工具：`validation/diag_u1_tau_twopoint.jl` +
   `validation/u1_tau_twopoint_alignment.md`（判决档案）。

## 文件

```
src/SingleModeWeyl.jl   模块入口
src/contours.jl         WeylRealTime / WeylImagTime / weyl_lattice
src/kernels.jl          HybridizationKernel + KernelCorrView（hybriddynamics 接口）
src/tensors.jl          CellCore / split_core / tile_cells!（Z2 分级精确拆分）
src/system.jl           local_weight / system_factor / WeylInsertion
src/influence.jl        influence_functional（复用 hybriddynamics）
src/contraction.jl      partition_function / expectation（复用 integrate）
src/solver.jl           Problem / Solution / solve
src/measurement.jl      occupation / greater / lesser / matsubara
src/diagnostics.jl      fit_order / convergence_report / maxbond
validation/ed.jl                        单模 Anderson ED 参考（bit 基 JW）
validation/calibrate.jl                 单元级校准（约定核对）
validation/convergence.jl               收敛阶基准（①②③）
validation/realtime_dtconvergence.jl    实时 δt 收敛探针（开放项 1 判据）
validation/if_coefficients_mine.jl      IF 全单项式系数导出（本侧）
validation/if_coefficients_u1.jl        IF 全单项式系数导出（U1 侧）
validation/if_coeff_rowprobe.jl         单行 partial-IF 探针（提取器自检）
validation/if_coefficients_compare.jl   IF 系数比较器（双射 + 失配分类）
validation/kcoeff_mine.jl / kcoeff_u1.jl / kcoeff_compare.jl
                                        K 侧同流程交叉核对
results/convergence.csv 基准数据
results/ifcoeff_*、kcoeff_*             交叉核对数据与报告
```
