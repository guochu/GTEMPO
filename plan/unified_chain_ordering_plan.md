# 统一链存储序计划（v3）：热段降序 + 实时段 ± 成对升序——纯置换路线落地

- 日期：2026-11-20（v3：Phase 0 探针推翻 Route B 必要性，纯置换路线实施完成）
- 基线：GTEMPO `df8a9af1`（if-expmult-backend；插入符号化 + 穿戴内联 +
  荷账自读已落地，全套件绿）
- **状态：已实施（Phase 0/1 完成，Route A 纯置换；v2 的合并超算子设计 §3 撤销备用）**
- 关联：`lattices.jl` 头注释「排布固定 Chrono」；backup/GTEMPO 旧栈 ordering
  家族（± 逐时刻相邻布局先例）；plan/fockmpo_dress_selfread_plan.md

## 0. 设计（用户裁决，2026-11-20）

统一存储模板（三 contour 的 lattice 寻址与 GMPS 链存储）：

```
[ τ Nτ, τ Nτ−1, …, τ 2, τ 1 ]  ‖  [ (+,1), (−,1), (+,2), (−,2), …, (+,Nt), (−,Nt) ]
        热段（τ 降序）                    实时段（± 成对、时间升序）
```

物理读法：混合链 = 「热态制备段（β→0）+ ρ₀ 接缝 + 实时演化段」（τ 段/实段
分界键 = ρ₀ 切位，⟨U†ρ₀U⟩ 同构）；实时段 ± 成对 = 时刻局域 + 右端追加增长
语义（backup makestep TimeAscending 同款，为 InfiniteGMPS 留结构）。
imag 零变更；real/mixed 真变更；混合 τ 段降序为模板本意（已确认）。

## 1. 路线定谳（Phase 0 探针实证，2026-11-20）

**Route A（纯置换）实证成立**——K 链与 IF 链按同一 perm 经 `GM.reorder!`
（= 相邻 swap! 序列）重排后，Z/插穿观测值三轮廓全部机器精度不变
（real 零插入/greater 双插/lesser 含 (−,0)/n̂+初态、imag 反序对照、mixed τ 奇
双插+rt 偶插，Δ ≤ 5e-15）。

**机理**：K 链编码的是标量 Grassmann 多项式（Fock 迹闭合已烙进系数），链站
序只是多项式因子化细节；swap! 按费米号做变量重标号，多项式值不变——IF 链
同 perm 重排后逐站配对一致。v2 的「算符乘积序」担忧不成立（那是对裸 FockMPO
的直觉，不是对已转换 GMPS 多项式的）；**v2 §3 合并超算子格设计撤销**（对值
保持而言不需要；「杂质含 Lindblad 时合并格是双线耦合的天然容器」的判断
仍记录在案，留作未来 Liouvillian 扩展的设计参考）。

**IF 侧等价性实证**：直接按新序装配（cell_block 抽象驱动）≡ 旧序装配 +
reorder!（todense 机器精度一致）——IF 装配无需任何改动。

**约束不变量**（全部保持）：`chronological_cells`/切位语义/K 装配公式/
(a,ā) 相邻 AC 前提/imag contour 一切（含 ExactTTI）。

## 2. 实施（已完成）

| 项 | 内容 |
|---|---|
| `src/lattices.jl` | real `_index`（块 = 2i−1/2i）、mixed `_index`（τ → Nτ−i+1、± → Nτ+2i−1/2i）、`cell_block` 同步、头注释/docstring |
| `src/sysdynamics/sysdynamics.jl` | `_storage_cells`：real/mixed 改 ± 成对交错 |
| K converter（fockmpo.jl） | **零改动**——perm 由 `chrono`×`_storage_cells` 自动重派生 |
| `test/api/port_smoke.jl` | `_index` 位置断言按新序更新 |
| 其余测试 | 零改动全绿（值类测试天然序无关） |

## 3. 踩坑记录

- `_storage_cells` 首版误写成「+ 整块再 − 整块升序」（分支分段而非 ± 交错）
  ——它产生的 perm 是自洽的恒等式变形，闸门/δ 锚全数通过，**只有 golden
  值锚抓住**（Z = 0.9588 vs 0.99985）。教训：链序类改动的验收必须包含
  绝对物理值锚（golden），相对一致性测试（对拍/自洽）对此类 bug 不敏感；
- `GM.reorder!` 的 perm 语义 = 新链第 k 站 = 原链第 perm[k] 站（旧用法恒等/
  反序均为自逆置换，方向不可见——交错序是首个方向敏感用例，本次语义无恙）。

## 4. 验收

Part 1 API 全绿（golden U1GTEMPO2 机器精度锚 + ED 电池 + 闸门 300 项 +
δ≡0 锚）；Part 2–4 全套件见套件记录。物理值逐点不变。

## 5. 不做清单

- 不恢复 backup 的 GrassmannOrdering 类型层；不动 `chronological_cells`/
  切位语义/imag 装配；不实施 Lindblad（v2 §3.4 设计红利记录保留——合并格
  作为 Liouvillian 容器的判断不受本路线影响，另行立项）；不做增长式 API。

## 6. 实施状态

**✅ 已实施（2026-11-20，本会话）**。工作树待提交。

---

## 7. imag converter `reorder!` 不可去性定谳（2026-11-20 晚，探针实证）

折叠装配（foldedmpo.jl）做到 real/mixed 装配内序直出后，imag 支仍保留
「时序核发射 + 末端 `GM.reorder!`」。用户再问可否装配内直发；结论：**在
v3 存储模板（imag τ 降序，"imag 零变更"）下不可去**，已数值定谳：

- **恒等式**（四配置：b=1 Nτ=3/4、order 1/2、b=2 Nτ=2 Kanamori；机器精度
  ≤1.6e-15）：
  `todense(A)[m⃗] = (−1)^{C(k(m⃗),2)} · todense(C)[变量块反排(m⃗)]`
  其中 A = deprecated converter 输出（存储布局），C = 同一 imag 核发射
  **去掉 reorder!**（时序布局），k = 该赋量下奇占据格点数（格 = 2·bands
  连续变量）。folded imag 支与 deprecated todense 逐位 Δ = 0.0。
- **无符号直发残差 O(1)（1.8–2.3）**，且带 −1 号的系数占非零系数多数
  （6/8、14/18、20/28、13/29）——既非全局相位（同链 ±1 并存），也非逐格
  可分解（翻单格奇偶只在 k ≡ 2,3 (mod 4) 时翻号 ⟹ 依赖其它格点的占据）
  ⟹ 不可能写进任何逐格核值；只能经键承载的 graded 辫（`GM.reorder!` =
  分级相邻 `swap!` 序列，逐融合块补费米交叉号）执行。等价去处 = bond
  携带奇格点计数 mod 4（键膨胀 ×4 + 自制机器），贵于 reorder!，不取。
- **imag 为何是唯一例外**：real 直发靠 ± 折叠恒等式（𝓛_j = A⁻_j(·)A⁺_j
  把两支并成单格，装配序 = 存储序，无反排）；mixed τ 段反序被 vec/bra
  双线恒等式（τ₁ vec 头 + τⱼᵀ bra 侧追加）**代数吸收**。imag 单线环两类
  装置都没有——ρ₀ 邻接结构又把 FockMPO 数据序物理地钉死在时序（演化方向
  不可反），于是「时序装配 → 降序存储」的 graded 反排无代数出口。
- **去 reorder! 的唯一干净出路 = 虚轴存储序翻转**（imag `_index`/`cell_block`/
  `_storage_cells` 三处 k−i+1 → i ⟹ 存储 = 时序，converter 直发，perm 恒等）。
  代价：违 v3「imag 零变更」用户裁决（imag 与 mixed 热段 τ 降序「链头 = β 端」
  的模板对称性破裂）；ExactTTI 虚轴环族按新方向重推（「链站序与核胞序相反」
  烙在 ch1/ch2/ch3 通道分配，7 配置×gt/lt dense 对拍须全部重走）；golden 布局
  锚重锚（值类锚不变，按 §3 踩坑须配绝对值锚）；partialif/measurements 的
  cell_block 驱动面审计。**属待裁接口语义变更，未裁决不实施。**

探针逻辑（/tmp 易失，故录此）：deprecated imag 发射逐行镜像、去掉 reorder!，
`pairize` 后 `GM.Dense.todense`（轴 = 变量级，每轴 2 维、占用 = index−1），
变量块反排逐条目比对 + k 逐条目计数。

### §7.1 用户两条解除路线的探针裁决（2026-11-20 深夜；探针 2–6）

用户质疑：①阻碍是否来自边界条件烙在 FockMPO（若然 mixed 也该做得到）；②
是否因限定 FockMPO 直积态——放松后先转 GMPS 再补边界/或合格装配。两条均
被实验否决，且 §7 的论证由此**补强为一般性结论**：

- **①边界条件无关**。反序号是体现象：两格开放链（无迹闭合/ρ₀）即出现
  `x₂x₁ = −x₁x₂` 的偶-偶交叉号；(s,t) 读法的头/尾形态在穷举中作为自由
  取向覆盖，无组合幸存。mixed 之所以"看似做得到"也**不是**因为边界——见
  ③，其 τ 段直发撞的是同一堵墙。
- **②放松直积买不到**。三层实验：
  a. **vec/bra 路由**（mixed τ 段同款装置：Tr 循环对称 + `Σ_s X_{ss}` 闭合格
     + δ_{s,v}·A[w,t] bra 中间格 + vec 头，降序直发、pair 荷键沿存储序走）
     ——产物是**共轭多项式**（η↔ā 全站镜像，支持集错位），非带号表；
  b. **64 组合穷举**（中间格 4 种 tie × 4 种矩阵元取向 × 闭/头格各 2 种转置）
     ——零幸存；
  c. **一般性证明（二格双线性论证）**：正确表 d(m₂,m₁) = (−1)^{f₂f₁}·cᵀ，
     若站点张量 = 逐格读法（值 = 单个传播子元，可乘逐条目符号 γ₁(l,m)γ₂(l,m)），
     则 γ 之积须与键态 l 无关 ⟹ 分解为每站槽函数 ⟹ XOR≠AND，不可能。槽级
     oddlater 型符号同理（线性 vs 二次型）。
- **键维测量（重要否定结果）**：正确降序表的标准键表示**存在且不胀**——
  b=1 bonddims(A)=[4,4,4,4,4] = C；b=2 A=[4,16,16,16,16,16,4] < C（split_cell
  无截断满数值秩 64/37/7，reorder 的 1e-14 截断反而压缩）；U1U1 同（A ≤ C）。
  即"×4 寄存器"猜想错误——信息装得下标准键空间，但**只能以 SVD 重划分形态**
  装入：A 的站点张量是 Schmidt 规范混合值（相邻格内容互相渗透，Nτ=2 全块
  数值实证），非任何逐格函数。**这个重划分正是 swap 网络计算的东西**——
  在 converter 内"直接装配"它 = 把 reorder! 搬进装配内部，无净收益。
- **③mixed τ 段连带发现（给 folded WIP 会话）**：folded mixed 的 τ 段降序
  直发与 imag 同构受阻（vec/bra 装置不吸收 graded 反序号——探针 2 即该装置
  在纯虚时的实现，失败形态=共轭镜像）。mixed 未过门的原因清单应加上这一
  项；τ 段若要无 reorder!，与 imag 共享同一裁决（虚轴翻转或接受重划分）。
- **结论不变**：现行"τ 降序"模板下，imag（与 mixed τ 段）converter 的
  graded 辫不可去掉、只可换位置；结构性解除 = 虚轴存储序翻转（§7 末，待裁）。

### §7.2 ED 裁判：reorder 前后系数不等价（2026-11-20 深夜；探针 7–8）

用户质疑"reorder 只是置换顺序，构造时对齐排布即可；升序/降序结果应一致，
reorder 前后系数应完全等价"。裁决分两层：

- **事实层：`GM.reorder!` 不是纯置换**。最小两站演示（`monomial` 链 x₁…x₃）：
  升序写 +1，降序带号写 −1，`GM.swap!` 后 −1（与带号降序逐位相等），
  朴素纯对调站点张量 → 0.0（键荷失配）。`_swap_gate` 对每个融合树块显式补
  `(−1)^{奇(p₁)·奇(r₂)}`（2026-11-13 修复史同源）。
- **物理层：ED 三方对拍**（单模浴 Nτ=8 β=0.8 PartialIF，同一管线/IF 因子/
  插入约定，仅换 K 链 = 生产带号降序 vs 朴素无号降序）：
  - ⟨n⟩：生产 0.428897 vs ED 0.429038 ✓；朴素 0.421671 ✗；
  - G(τβ)：生产 −0.428897 vs ED −0.429038 ✓；朴素 −0.548954 ✗（28% 偏，
    且破坏平移结构）；
  - Z：生产 1.8328 vs 朴素 1.6069（12% 偏）；
  - 系数表自检：生产 = 带号块反排（3e-15）、≠ 无号块反排（1.92）。
- **概念澄清**："升序/降序结果一致"在物理层成立（Phase 0 值不变实证），
  这正是虚轴翻转可行的依据；但"系数完全等价"不成立——系数表是多项式在
  **给定生成元写序下的坐标**，写序翻转 = 外代数的 graded 基变换
  （转移矩阵 = (−1)^{C(k,2)}，非恒等）。逐站读坐标的收缩机器必须配对
  正确的那张表；无号表配存储约定消费端 = ED 可证地错。

---

## 8. v7 虚轴翻转实施（2026-11-20 深夜，用户裁决后落地）

**裁决**：用户拍板"虚轴定义按正序排列"——纯虚时 lattice 存储 = 时序升序，
converter 装配内序天然一致，imag `reorder!` 删除。

| 项 | 内容 |
|---|---|
| `lattices.jl` | imag `_index`（块 = i）、`cell_block`（imag = j）、头注释；**mixed τ 段不动**（仍 Nτ−i+1，v3 模板） |
| `sysdynamics.jl` | `_storage_cells` imag = [τ 1…Nτ] |
| `foldedmpo.jl` | folded `fockmpo2gmps` imag 支删 perm/reorder 块（时序直出）；docstring 同步 |
| `fockmpo.jl` | deprecated converter perm 由 `_storage_cells` 自适应 = 恒等（reorder! 退化为无操作，跨布局兼容保留）；docstring 注记 |
| `exacttti.jl` | **虚轴环族按新方向重推**（§8.1）；零核列 λ=0 通道不变（带序在翻转下不变） |
| `test/api/port_smoke.jl` | imag `_index` 断言按升序更新（格 1 = 块 1） |

自动适配面（零改动）：partialif/measurements/bridge/influencefunctional
（全部经 `cell_block`/`_index` 间接寻址）、golden U1GTEMPO2 锚（值类）。

### §8.1 ExactTTI 虚轴环族（升序重推版）

取向读号规则（由旧降序家族的正确性反推钉死）：原取向（conja=ca, conjb=cb,
asite=ba, bsite=bb）dense 读 = +系数；翻转取向（cb,ca / bb→ba）读 = −系数
（旧 ch2 的上三角覆盖需求 ⟹ 读 = −系数 实证）。

- ch1 前向族（原取向）：α₁ = **−w·q^{n−1}**、节点 **1/q**——覆盖 {jr < jc}
  （链前向 = 上三角），距离 d 系数 upper[d] = −lower[n−d]；
- ch2 后向族（翻转取向）：α₂ = **−w/q**、节点 **q**——覆盖 {jr > jc}（下三角），
  翻转读 −系数 ⟹ 贡献 +lower[m]；
- ch3 λ=0（同胞折叠取向随带序，带序在翻转下不变 ⟹ 物化条件同旧）：
  ba < bb ⟹ Σα₃ = M[1,1] − T₁；ba > bb ⟹ Σα₃ = −M[1,1] − T₂；
  ba == bb ⟹ Σα₃ = **M[1,1] − T₁ + T₂**（对角读 (T,F)→−T₁、(F,T)→+T₂）。
- **符号定谳路径**：手推 ba==bb 首版写反（+T₁−T₂）→ 4×4×7 符号暴力搜索
  （oracle = ET vs PartialIF todense 逐位 + Z）唯一幸存 (s1,s2,c3) =
  (−,−,M−T1+T2)，Δ表 = 1.6e-15 —— 与"对角读"独立推导一致。
- 零核列 λ=0 单发路径与热配对审计（upper = −reverse(lower)，核侧恒等式）
  不变。

### §8.2 验收

Part 1 API 全绿（golden U1GTEMPO2 机器精度锚 + ED 电池 + FockMPO 闸门 300 +
ET≡PartialIF 虚时 Z rtol 1e-8）；Parts 2–4 见套件记录。物理值逐点不变
（imag ⟨n⟩ ED 误差保持离散化常数档 1.5e-4~7e-4）。

### §8.3 遗留

- mixed τ 段降序保留（v3 模板）——与 imag 方向不再一致（裁决接受）；
  mixed 折叠 τ 段直发若立项仍需按 §7.1-③ 处理（reorder 或重推）。
- golden fixture 的"虚时 compact 通道 [格 Nτ…1] 天然同序"备忘注已过时
  （compact 通道代码不在 fixture 中，仅历史注记）。

### §8.4 折叠 converter 删 trunc 输入（2026-11-21 用户裁决）

折叠 `fockmpo2gmps` 删除 `trunc`/`truncerrs` kwargs 与两处 `GM.truncate!`
后处理（imag + real/mixed 支）——converter 成为纯翻译（签名 =
`fockmpo2gmps(f, lattice)`）。套件从不经此管道传 trunc（grep 全测试面零调用；
IF 算法的 trunc = SVDCompression 算法对象字段，另一机制不动）。冒烟：折叠
imag vs deprecated todense Δ = 3.3e-16。

**未动项与原因**：deprecated `fockmpo2gmps_deprecated` 的 `reorder!`/
`trunc` 保留——其 reorder 是 real/mixed 布局的承重件（chrono→storage graded
重排；v7 只把 imag 段 perm 变恒等），删除 = real/mixed 链布局错位 = 套件红；
它随折叠管线过门接线后整体退役。`sysdynamics` 公开 API 的 trunc kwarg 管道
（转发 deprecated）同保留——公开接口变更另行裁决。


---

## 附录（v2 历史存档：合并超算子格设计——本路线下撤销，留作 Liouvillian 扩展参考）

> v2 曾按「± 相邻 ⟹ per-支独立装配不便 ⟹ 合并超算子格」的推演把 Route B
> （超算子重构）立为主线。Phase 0 探针推翻其前提：K 链的 GMPS 多项式值与
> 站序解耦（swap! 值保持），纯置换即可。以下设计要点存档：
>
> - 合并格 𝓛_j（ket/bra 双线同时演化，单项式腿 4^{2B}，FockMPO 站数
>   real 2Nt→Nt）：对**密度矩阵型**演化对象（ADM/Liouvillian 折叠
>   ρ_j = A⁻_j ρ_{j−1} A⁺_j）仍是最自然的容器；
> - 折叠恒等式：Tr[A⁺₁⋯A⁺_N A⁻_N⋯A⁻₁] = Tr[A⁻_N⋯A⁻₁·𝟙·A⁺₁⋯A⁺_N]
>   （轮廓积的循环平移 = 逐步超算子演化）——杂质含 **Lindblad 耗散**时
>   D[L] 的双线耦合块（ket L·L† + bra 共轭 + 跨线交叉）以合并格为天然容器，
>   per-支独立装配结构上无法表达跨支耦合——**Liouvillian 扩展立项时的
>   装配路线即此**；
> - 定谳点 D1（键荷差标签双线保守性）/D4（插入半支寻址）/D5（δ 双线账）
>   的分析在 Liouvillian 立项时直接复用。
