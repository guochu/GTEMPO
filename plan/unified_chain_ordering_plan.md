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
