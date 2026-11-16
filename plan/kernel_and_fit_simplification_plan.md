# HybridizationKernel 轻量化 + fit_exponentials 与 ExpExp 去重裁定计划

- 日期：2026-11-16
- 基线：GTEMPO HEAD `c01f446c` + 工作区
- 关联：`measurement_cache_refactor_plan.md`（缓存键引用 kernel.blocks 的语义——Phase B 实施时以本计划 Part 1 的 accessor 为准）；其余计划无文本冲突

## 0. 已裁决

- Part 1（kernel.blocks 物化退役 → QuAPI corr 惰性视图）**立项实施**（用户 2026-11-16 提案，本计划论证可行）。
- Part 2（fit_exponentials vs ExpExp）**裁定为"保留现状 + 文档注明"**（推荐案 A；rebase 案 B 的代价与路径一并列出，供用户翻案）。

## 1. Part 1：HybridizationKernel 轻量化

### 1.1 探查定谳

**F1【QuAPI 侧已是平移不变存储——GTEMPO 物化是纯损失】**
- `fermionic_Δt` 返回 `RealCorrelationFunction(CorrelationMatrix(.-c,.-b), CorrelationMatrix(.-a,.-d))`（QuAPI realtime/fermionic.jl:153）：只存交叉对，四个 branch 块共享四个 Toeplitz 向量 O(N)，构造时一次算完（注释明言此前每块四重复制翻倍 quadgk 成本）；
- `RealCorrelationFunction{CorrelationMatrix}` 统一约定（realcorr.jl:14-130）：交叉对携带全部信息，`branch(x; b1, b2)` 从 Toeplitz 向量**符号组合派生**任意块（:176-197），构造器还校验 Keldysh 恒等式（:97-104）；
- `Δτ` → `ImagCorrelationFunction`（单 `CorrelationMatrix`，ηⱼₖ/ηₖⱼ 按 |i−j| 查表，corrmat.jl:37-54）、`Δm` → `MixedCorrelationFunction`（Gττ Toeplitz + 实时交叉 Toeplitz，mixedcorr.jl:12-40）。
- GTEMPO `hybridization`（influence.jl:61-83）把每个 branch 视图物化成稠密 `Matrix{ComplexF64}`：O(N²) 内存 + 整块复制 + **丢弃 Toeplitz/统一约定结构**——三项全是损失，零收益（物化副本不被二次修改）。

**F2【branch 视图逐位精确】** `CorrelationMatrix` 的 `getindex(i,j)` 是 |i−j| 查表（corrmat.jl:54），branch 派生是纯符号组合——惰性视图与现稠密块**逐位一致**，golden 不重锚。

**F3【消费方清单（全部只读索引）】**
- `partial_if`（influence.jl:147）：逐系数 `blocks[(b, cell[1])][j, cell[2]]`；
- ExactTTI（influence.jl:336-339+）：整列读取 `blocks[(:+,:)][i+1, 1]` + `validate_tti(blocks[...])`；
- currents（measurement.jl:354-355/388-389）：`ηpp = kernel.blocks[(:+,:+)]` 后逐元素 `[k,j]` 循环；
- BCS 路线**不经过 HybridizationKernel**（`bcs_partial_if(lattice, corr; …)` 直通 corr）——不在范围内。

**F4【validate_tti 可简化】** `validate_tti(blocks[(a,b)])` 校验输入块的时间平移不变一致性；换 corr 存储后 Toeplitz 是**结构性保证**（QuAPI 构造器已校验统一约定/Keldysh 恒等式）——该检查退化为恒真，可删（保留一处注释指向 QuAPI 构造器校验）。ExactTTI 的 `fit_exponentials` 只消费整列，与存储形式无关。

**F5【窗口语义】** QuAPI 对象带 N+1 边界向量；现 kernel 取 `[1:Nt, 1:Nt]` / `[1:Nτ, 1:Nτ]` 窗口（influence.jl:66/71/80）——惰性视图须保留窗口（`view(branch_view, 1:N, 1:N)` 的 SubArray 保持 O(1) 查表索引）。

### 1.2 设计

```julia
struct HybridizationKernel{L,C}
    lattice::L
    corr::C                                   # QuAPI Imag/Real/MixedCorrelationFunction
    coupling_norm::Union{Nothing,Float64}
    diagnostics::Dict{Symbol,Any}
end
```

- **兼容 accessor（零消费者 churn 的关键）**：`Base.getproperty(k, :blocks)` 惰性构建 `Dict{Tuple{Symbol,Symbol},AbstractMatrix}`，值 = 窗口化 branch 视图（F5）。现有 `kernel.blocks[(b1,b2)][i,j]` 语法与数值**逐位不变**；视图按 (b1,b2) 缓存于该惰性 Dict（每对一次构建）。
- 构造端：`hybridization(lattice, bath)` 删三处物化（influence.jl:66/71/80），直接存 `corr`；形状/有限性校验移到惰性视图首建或保留在构造端（对 branch 视图 `size` 抽查一次）。
- `validate_tti` 删除（F4），ExactTTI 入口留注释指 QuAPI 构造器校验。
- currents 两处（measurement.jl:354/388）改局部物化 `ηpp = Matrix(kernel.blocks[(:+,:+)])`——函数内一次性、O(N²) 瞬态与现状相同，循环性能不变。
- `correlationfunction.jl` 的跨支负号约定注记（:24-34）**原样有效**：blocks 视图物化的就是 QuAPI branch 的值，读点负号处理不动。

### 1.3 收益

- 内存：实轴 4×O(N²) 稠密 → O(N) 向量；虚轴 1×O(N²) → O(N)；混合同比例；
- 构造端删全部物化复制（长轮廓 N=数万时 N² 复制是秒级开销）；
- 平移不变结构保留（未来长时拟合/解析通道可直接吃 Δk）；
- 跨块符号/Keldysh 约定单一来源（QuAPI），GTEMPO 不再维护物化副本。

### 1.4 验收

1. 全测试绿（含 quapi_crosscheck、fewmodesbath ED 对拍、golden 全套——值逐位不变是硬闸门）；
2. 新增惰性视图逐位测试：`Matrix(kernel.blocks[(b1,b2)])` 与旧物化路径对拍（Imag/Real/Mixed × 全 branch 对）；
3. 内存基准：performance/ 加 N scaling（如 N=512/2048 的 kernel 驻留内存）before/after。

## 2. Part 2：fit_exponentials 与 ExpExp——裁定保留，不 rebase

### 2.1 探查定谳

**F6【同族但非重复：GTEMPO 版 = U1GTEMPO2 逐字移植 + 三处领域扩展】**
- 核心同为 Hua–Sarkar 矩阵铅笔（Hankel SVD → 移位算符特征值 → LS 系数）；
- GTEMPO/U1GTEMPO2 版（influence.jl:173，diff 确认逐字同源）有三处 ExpExp **没有**的扩展：
  1. **双端锚定**（`rightanchors`：|q|>1 的极点折算为 :right 锚 `q^{L−d}`，保证全部传播节点模长 ≤ 1——IF 因子稳定性的关键，ExpExp 无此概念）；
  2. **行子采样**（`pencil_rows=128`，长块 Hankel 行数上限）；
  3. **rtol 秩扫描 + 退化窗回退**（r=1..maxterms 早退；L ≤ maxterms 时单位根精确插值回退——注释明言"精确有界展开，非回退到其它 IF 算法"）。
- ExpExp 的迭代/剪枝在 `exponential_expansion` 包装层（stepsize 级联 + cut），其 `matrix_pencil(s, p)` 是单秩、无锚定的标准型——两边移位算符的数值路径不同（eigvals(U1\U2) vs U′H₁V′D⁻¹），数学等价但**浮点位级漂移**。

**F7【依赖面】** TEMPO 依赖 ExpExp；GTEMPO/U1GTEMPO2 不依赖——fit_exponentials 的谱系是 U1GTEMPO2 专门化铅笔，非 ExpExp vendor 拷贝。

### 2.2 裁定

**案 A（推荐，已选）：保留 fit_exponentials 现状**，文件头补两行文档（与 ExpExp 的同族关系 + 三处领域扩展清单），不引入 ExpExp 依赖。理由：改动收益 ≈ 核心 ~15 行的去重；代价 = 因子位级漂移 → ExactTTI 全链重锚 golden + 双栈对拍重跑——性价比倒挂。

**案 B（备选，供翻案）**：GTEMPO 增 ExpExp 依赖，`fit_exponentials` 重写为 ExpExp 核心 + 本地锚定/扫掠/回退包装。安全路径：先做 ExpExp 侧"anchored MatrixPencil"算法扩展（锚定逻辑进 ExpExp，通用受益）→ GTEMPO 包装替换 → golden 重锚（ExactTTI 六基准 + 双栈对拍新基线）。工作量 ≈ 2–3 天，收益 = 消除算法双维护。若裁决走 B，建议作为独立计划与 golden 重锚批次一起排期。
**【2026-11-16 翻案】案 B 已采纳并细化：见 `ExpExp/plan/expexp_anchored_pencil_plan.md`（含默认参数对齐、锚定三元组独立函数、共享循环 vs 自带扫的性能裁定、MAK/LinearAlgebra SVD 位精确前置闸门）。**

## 3. 实施顺序与验收

- Part 1 单独一提交：实现 → §1.4 三闸门（全绿 + 惰性视图逐位对拍 + 内存基准）。
- Part 2 文档两行随 Part 1 提交或单独提交；案 B 不排期（除非翻案）。
- 与 `measurement_cache_refactor_plan.md` 的顺序：建议本计划 Part 1 先行（缓存键在 accessor 语义上写，避免两改一层）。

## 4. 不做清单

- BCS 浴 kernel 化（corr 直通已是最轻形态）；
- ExpExp 依赖引入（除非案 B 翻案）；
- `hybridization(lattice, corr)` accessor 路线复活（T7 已裁退役， throwing 现状保持）。
