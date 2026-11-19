# BCS/ExactTTI 深化计划：is_tti 门控 + 任意多带推广 + 去顶点移位 + Weyl2×BCS

> **状态：立项待实施（v3，2026-11-19）**。演化线：v1 = ExactTTI 支持 BCS
> 轨道核（直读面 → 拟合面）；v2 = 用户裁决深化——TTI 假设改为
> `is_tti` 显式块级检验，band 结构守卫退役，适用面 = 任意多带核；
> **v3 = 用户裁决新增两点**：①**去掉 BCS 实轴顶点移位**（A2 就此裁决，
> 比 O2-lite 更彻底——移位整体移除而非搬到读点，二阶装配吃未移位 IF）；
> ②**Weyl 二阶 × BCS 并入本计划**（吸收 `weyl2_bcs_plan.md`——其唯一
> 开放推导点 F5 已由 bcsweyl 裁决解决，剩余 = 端到端收敛验证与建测试）。
> 用户指定 R/W 两阶段**先行**。

- 日期：2026-11-19（v1/v2/v3 同日成文）
- 基线：GTEMPO 分支 `if-expmult-backend`（第九批 `18ce8f65` 后）+ 工作区
  并行 WIP——实施前确认基线
- GM 前置：同 ExactTTI 现行前置集（`expdecayexp-same-orbital-fix` 已推送）；
  **GM 侧零改动**
- 关联计划：`weyl2_bcs_plan.md`（**被本计划 Phase W 吸收/取代**）；`polarif_plan.md`
  （测试算例共享）；`exactif_port_plan.md`（严格 IF 消费侧，无冲突）；
  `ifd_backend_switch_plan.md`（is_tti 的容器类型经 reexport 在域内，顺序无关）

## 0. 裁决面

**已裁决**：

- v2（用户 2026-11-19）：is_tti 门控 + band 守卫退役 + 任意多带适用面；
- **v3-①（用户 2026-11-19）：BCS-RT 去顶点移位**——blocks 改存 raw 味块
  branch 视图，移位（`_vs_rt_acc` 重索引）在 GTEMPO 消费面**整体移除**
  （不搬读点——PartialIF/ExactTTI/Weyl2 全消费端统一吃 raw）；构造期
  跨支 Keldysh 号折叠保留（CM 整块取负，保 Toeplitz 性）；
- **v3-②（用户 2026-11-19）：Weyl2×BCS 并入**（Phase W，见 §5）；原
  A2（raw 存储仲裁）就此关闭。

**待裁**：

- **A1【阶段范围与顺序】**：推荐 R（去移位）→ W（Weyl2×BCS）先行（用户
  已指"先"），A′（ExactTTI is_tti+模式+守卫退役+direct 泛化）、B′（fit
  泛化+同味探针）随后；四阶段相互独立可插队。
- **A4【is_tti 容差缺省】**：推荐 `rtol=1e-12`（相对块幅值）——物理核
  同源装配逐位相同不触发；备选精确 `==`。

## 1. 核心设计：is_tti + 通道模式选择（v2 承继）

### 1.1 `is_tti`（exacttti.jl，或 correlationfunction.jl 公共面）

```julia
"is_tti(M)：块的时间平移不变性（Toeplitz 性）。CorrelationMatrix 结构性
成立恒 true；一般矩阵逐对角线检验常数性（相对块幅值的 rtol 容差，
零块恒 true；矩形块按 τ↔rt 交叉块形状支持）。"
is_tti(::CorrelationMatrix) = true
function is_tti(M::AbstractMatrix; rtol::Real=1e-12)
    n1, n2 = size(M); amax = maximum(abs, M)
    iszero(amax) && return true
    for d in -(n2 - 1):(n1 - 1)
        i0 = max(d + 1, 1); j0 = i0 - d
        v0 = M[i0, j0]; i = i0 + 1; j = j0 + 1
        while i <= n1 && j <= n2
            isapprox(M[i, j], v0; rtol=rtol, atol=rtol * amax) || return false
            i += 1; j += 1
        end
    end
    return true
end
```

### 1.2 通道模式选择（is_tti 的消费方式）

`tti_imag_factors`/`tti_real_factors` 入口逐块检验，模式写进 audit
（`tti_mode::Symbol` ∈ `(:fit, :direct)`）：

| 通道 | fit 模式（TTI 专用机器） | direct 模式（回退，精确对任意块） |
|---|---|---|
| **imag** | ESS 极点级联（热配对审计 + 环族拟合；§3.2） | 全矩阵 Pairing 直读（混合格 τ 段恒 direct——容器限制非 TTI 缺席） |
| **real** | (g,l) 列拟合重构（Keldysh 审计；§3.3） | 原始块 Pairing 直读（混合格 rt 块既有先例） |
| **cross** | ——（无 TTI 假设） | 恒 key 直读 + SVD 秩一 |

- 两模式对各自输入块**精确**；is_tti 失败不报错、降级 direct + audit 钉住；
- **同味环族就绪度门**：imag fit 对 (T,T)/(F,F) 味对额外要求几何已验证
  （§3.2 探针闸门）——未验证阶段同味对自动 direct（可分味对混合）；
- real fit 仅纯实时格。

### 1.3 守卫重构

`_tti_band_guard` 三重限制（单带/拒反常/拒跨带）退役 → `_tti_guard`
（band 范围 + 任意键味对；TTI 前提下沉到两因子构造器入口的 is_tti）。
装配路由：**band 对角 normal 核 → 逐带元组**（现行不动，golden 锚）；
**含跨带/反常键 → 单链**（PartialIF 第九批先例）。D7 的 IF 算法面叙述
退役；对称性守卫（BCS 仅 :z2，构造端）不变。

## 2. Phase R：BCS-RT 去顶点移位（先行，~0.5–1 天）

### 2.1 改动

`correlationfunction` BCS-RT 分支（`correlationfunction.jl`
`_bcs_dense_block`，:222-229）退役，替换为 imag/mixed 同款 raw branch
视图 + 构造期跨支折叠：

```julia
W = _window_block(branch(_bcs_flavor(corr, dba, dbb); b1=b1, b2=b2), n1, n2, T)
(keldysh && b1 !== b2) && (W = -W)     # 跨支 Keldysh 号保留（CM 整块取负）
```

- **跨支号保留在构造期**（整块取负 = `CorrelationMatrix` 向量取负，
  Toeplitz 性不变 ⟹ is_tti ✓）——与现行指数系数约定仅差移位本身，
  符号簿记零重推导；
- `_bcs_dense_block`/GrassmannContourIndex 访问器消费删（GTEMPO 从此
  **不再消费任何移位读取**——QuAPI 包内 `_vs_rt_acc` 语义不动）；
- 存储约定注释（correlationfunction.jl 头 + docstring）与 changes.md
  更新：「BCS 实时 = 访问器物化（顶点移位内嵌，bcsweyl 锚定口径）」
  标注被本裁决取代。

### 2.2 效益（本轮讨论的动因闭环）

- **Weyl2×BCS 的对齐前提**：二阶格点的接触门/系数结构按 raw 口径锚定
  （IF_constant≡Weyl 接触门锚点）——未移位 IF 使二阶装配免于接触行
  重分配的交互风险（上一轮定性：移位与二阶正交，但接触门交互是
  BCS-未验证组合的唯一可疑面；去移位后结构性消除）；
- **ExactTTI-BCS 实时 fit 模式解锁**：raw CM 通过 is_tti ⟹ §3.3 逐味
  Keldysh 审计 + 拟合重构可用（此前物化块只能 direct）；
- **BCS 电流 η 张力闭环**：`_kernel_eta`/currents 读到 raw = η 禁移位
  禁令（`inserted_k.jl:246` ED 锚）自动满足；
- normal/BCS 存储约定统一（BCS 三轮廓全 raw CM，与 normal 同构）。

### 2.3 代价与接受度

BCS 一阶实轴误差常数变差（bcs_ed 标定：vertex 比 raw 好 4–16×）——
**阶不变**（raw 仍一阶收敛），精度恢复走 Phase W 的 order=2（与 normal
同款决策路径）。imag/mixed 轮廓零影响（本就 raw）。受影响锚点 =
`docs/tutorials/firstorderif/*/bcs_realtime.jl` 常数档（顺带修正其
过时的「顶点位移装配」头注释）与未来 realtime ED 容差设定。

### 2.4 验收闸门（Phase R）

- 结构：BCS-RT blocks 全为 `CorrelationMatrix`（is_tti 恒 true）；
  同支 raw 值 ≡ 未移位访问器值（去重索引对照，一次性探针）；
- 物理：BCS 一阶实时 ED 斜率仍 ~1（常数变差预期内，容差按标定放宽）；
  imag/mixed BCS 测试零变化；
- 全套回归 + golden 按值（normal 不涉、BCS 虚时不涉—— golden 面不变）。

## 3. ExactTTI 深化（v2 承继）：band-pair 泛化因子构造

### 3.1 键遍历与站表

有序 band 对 = kernel.blocks 键去重投影 `(ba, bb, ca, cb)`（第七批 conj
键自完备）；`_kernel_block(k,b1,b2,band)` 推广 `_kernel_flavor_block`
（normal 同带 = 特例）；站表 `_station(lattice, cell, band)` 已逐带。

### 3.2 imag ESS：跨带与同味环族

- **(T,F) 对（任意 band 对）**：环几何与现行逐字相同（`ExpDecaySum` 
  asite=ba/bsite=bb 前向族 + 环绕族 + 接触修正）；热配对审计逐对；
- **(T,T)/(F,F) 对——推导提案（dense probe 定谳）**：同味下三角经反交换
  改写为前向通道——**环绕族 = `ExpDecaySum(T,T; asite=bb, bsite=ba)`**
  （站反转），系数 = −（转置下三角列）；两通道有序积集不重叠、无同站
  碰撞；接触修正公式同构。探针不过 ⟹ 同味对滞留 direct（就绪度门）。

### 3.3 real 拟合重构（fit 模式）

逐味对 (g,l) 列拟合 + 四块重构：现行 Keldysh 公式逐味复用——QuAPI BCS
每味块 = 标准 `fermionic_Δt` 产物（谱权重 u²/uv/vu/v² + 能隙色散），
内部恒等式构造端逐味校验 ⟹ 无新轮廓代数；**Phase R 后 BCS-RT raw 块
使本路线全解锁**。Keldysh 残差审计逐对化。

### 3.4 direct 模式 / 标量型

逐键 Pairing 直读（混合格 rt/τ 现行 direct 语义不变；手造非 TTI 核同
路径）；物化标量 = kernel T（第八批动态标量型；BCS-RT/Mixed 核结构性
ComplexF64 ⟹ 真空态已 promote，无需 `_complex_mps`）。

## 4. Phase A′ / B′（ExactTTI 深化的两段交付）

- **Phase A′（~1.5 天）**：is_tti 函数 + 模式选择 + `_tti_guard` 退役 +
  band-pair 泛化 direct 面 + 装配路由（单链/元组）+ 测试 T0/T1/T2。
- **Phase B′（~2–3 天）**：fit 面泛化——跨带 (T,F) ESS/拟合 + 同味环族
  推导探针（§3.2）+ real 逐味拟合重构与审计（§3.3）+ 测试 T3/T4。

## 5. Phase W：Weyl 二阶 × BCS（吸收 `weyl2_bcs_plan.md`，~2.5–3 天）

**前置全部就绪**：

- **F1【K 链 bath 无关】**：`sysdynamics(...; order=2)` 的 K = f(lattice,
  model)——BCS 只进 IF，Weyl-2 K 侧机器（S·M·S 半步、q=1/2 measure、W₄）
  零新张量代码；
- **F5【Weyl 标量 BCS 已裁决】**（bcsweyl/，2026-11-17）：反常道对 Weyl
  标量零贡献——`coupling_norm(::AbstractBCSBath)` = 每带裸谱强度
  （Δ/μ/β 不进；Bogoliubov 正则性；实证 cn\*(Δ)/Σα ≡ 1 ≤ 4e-4 +
  generic H_S 斜率 1.97–1.98 vs cn=0 对照 0.98）——原计划唯一开放
  推导点关闭；
- **F2【BCS IF 与 normal 同机器】**：fullfactor + conj 键（第七~九批）
  承接原 `bcs_partial_if` 语义，aa/āā 支撑位 :z2 下天然合法；
- **F3【一阶二阶同一套指标】**：Weyl 只改系数不改指标集；
- **Phase R【raw 基础】**：二阶装配吃未移位 IF（§2.2 对齐前提）。

**工作项**：

1. **BCS 实时 ED 管线移植 + 阶数实验（科学核心）**：`bcs_ed` 锚
   （test/util.jl 现成）+ fewmodes 单模配置，greater/lesser 全时序
   相对误差 log-log 斜率，N ∈ {4, 8, 16, 32}，order=1 vs order=2 对照
   （IF 同为 BCS PartialIF）。**判据：order=2 斜率 ≈ 2.00**（Weyl 标量
   无 BCS 修正的端到端确证；斜率 ~1 ⟹ 进入 Nambu 残余推导子项——
   F5 已裁决零贡献，此为保险闸门非预期路径）。同时验证 :z2 guard 与
   U≠0（Anderson）稳定性；
2. **接触门交互实证**：raw 基础下二阶斜率判据本身即接触门配合的验证
   （§2.2 风险闭环的实验面）；
3. **形式化**：`fock_sysdynamics` 无 BCS 特例（维持现状）+ 文档声明
   `coupling_norm` 语义 = normal 通道 α 和（Δ 不进标量）+ `bcs_symmetry_guard`
   :z2 保留；
4. **test/bcs 实时套件 + goldens**：order1/order2 双档 ED 对拍 + 去移位
   后的一阶常数档案（与 Phase R §2.4 容差设定合一）；BCS goldens 新建档
   （此前无 BCS golden）；
5. **电流 η raw 口径 BCS 验证**（Phase R 红利的实验面）。

**不做清单**（承原计划）：:u1/:u1u1 BCS（物理禁止，guard 保留）；IF 侧
二阶化（IF 精度独立于 K 离散化阶）；缓存/效率。

## 6. 测试计划（验收闸门总表）

- **T0** is_tti 单元（CM=true；Toeplitz 稠密=true 含噪声容差；randn=false；
  矩形两例；零块）；
- **T1** 模式选择（normal real=:fit；手造 randn 稠密=:direct；现行
  "反常核响抛"语义退役的回归面清点）；
- **T2** 互换闸门（跨带 TTI 随机核 2/3/4 带全耦合、多轨道 BCS——ExactTTI ≡
  PartialIF 机器精度；BCS-RT raw 后 real 走 fit 模式的逐味互换）；
- **T3** 同味环族 dense probe（§3.2 提案定谳；不过则就绪度门生效且 T2 仍绿）；
- **T4** ED 锚（BCS 虚时 ED 参数化 ExactTTI；跨带 TTI 核 U≠0 toy ED）；
- **T5** golden 按值（exactttiif 6/6、weyl 8/8——normal 单带现行路径零变）；
- **T6** Phase R 闸门（§2.4：raw 结构探针 + 一阶斜率 + imag/mixed 零变化）；
- **T7** Phase W 闸门（§5.1：order1 ~1 / order2 ~2.00 斜率 + :z2/U≠0 +
  BCS goldens 新档）；
- **T8** 全套回归（`julia --project=. --heap-size-hint=8G test/runtests.jl`
  直跑；GM 套件不重跑）。

## 7. 风险

- **R1【GM 同味通道物化覆盖】**：同味 ExpDecaySum 站反转 + 同味 Pairing
  的 GTEMPO 侧首次消费——T3 probe 与 T2 oracle 先行；
- **R2【同味环族符号簿记】**：§3.2 反交换负号若与 GM ESS 单项式序约定
  不合——probe 定谳，失败退 direct（不阻塞）；
- **R3【守卫退役的回归面】**：现行 ExactTTI 对反常/跨带核响抛——退役后
  改走 direct/fit，依赖抛错语义的现行测试需清点（T1）；
- **R4【fit/direct 混合模式】**：同一调用可分味对混合——audit 逐对钉住
  `tti_mode`，两模式精度语义文档钉死（fit = 拟合残差档、direct = 逐位
  原值档）；
- **R5【去移位的一阶常数回归】**：BCS 一阶实时常数变差 4–16×（标定
  支撑、阶不变）——容差设定与文档须同步，避免误判为 bug；
- **R6【跨支折叠符号】**：Phase R 保留构造期取负——若 raw branch 派生
  号与现行折叠号在某支对上不一致，T6 的 ED 斜率与探针暴露（一次性
  核对，预期零差异——访问器内部同用 branch 派生值）；
- **R7【Weyl2 斜率闸门失败】**：order2 斜率 ~1 ⟹ Nambu 残余推导子项
  （F5 已裁决零贡献，此为保险路径，触发则 +1–2 天）；
- **R8【WIP 基线】**：实施前确认基线（correlationfunction.jl 在列）。

## 8. 实施状态

**✅ Phase R / Phase W 已落地（2026-11-19，提交 `e1f7f291` R + Phase W 批；
基线 = 十批后 HEAD 02faecb5，R8 无冲突）**。A1 = 按用户指示 R→W 先行。

**✅✅ Phase A′ / B′ 已落地（2026-11-19，同日第五批；混合轴面按用户裁决退役
——比计划 §1.2 更彻底：混合格响抛而非 direct 承接，`tti_cross_factors` 与
`ExactTTI.cross_rtol` 删除）**：

- `is_tti` 公共谓词（CM 结构性 true / 稠密 Toeplitz 检验，rtol=1e-12 = A4
  采纳推荐值）+ 味对分组 `(ba, ca, bb, cb)` fit/direct 门控 + `_tti_band_guard`
  守卫退役 + 装配路由（band 对角 normal 逐带元组 / 跨带反常恒单链）；
- **虚时环族广义化**（§3.2 (T,F) 跨带 + 同味环族一并落地——推导以「通道 =
  链前向对 + 异味等时折叠」+ 链站序反转（`cell_block(:τ,j) = k−j+1`）钉死，
  dense 算符对拍：单对 bitwise + 完整 BCS 四键 3.5e-18，即 T3 探针通过）；
  **追加发现**：GM ExpDecayExp（WII 顶点序）只物化胞内前向对——asite>bsite
  通道的同胞折叠被结构性排除 ⟹ ch3 等时接触 λ=0 通道恒取带序前向取向收编
  全部同胞账目（ba==bb 退化为原式 `M[1,1] − Σα₁ + Σα₂`）；
- 实时逐味 (g,l) 拟合重构 + 逐味 Keldysh 审计（§3.3；BCS 四味对全部 :fit
  模式，残差 0.0）；audit 失败（非 TTI/热配对违反/Keldysh 违反）降级 direct
  不报错（R3/R4 处置：audit 钉 `mode` + 逐对残差）；
- `fit_exponentials` NaN 健壮性修复（铅笔退化列 [a,0,−a] 的 `err > tol` 轮盘）；
- 测试：if_gates 新 24 项门控（跨带/同味/BCS 虚实时互换、direct 降级、
  is_tti 单元、混合拒绝、元组/单链契约）；全套件 53 testset ~809 项绿 +
  exactttiif golden 6/6 绿（normal 单带路径零变）。

**Phase R（T6/T3）**：
- raw 直存落地；**R6 探针否决计划预设**——「跨支构造期取负」删除：bcs_Δt
  味块支视图自带与 normal 相同的指数系数 Keldysh 号（Δ=0 四支对逐位 0.0），
  取反即错。结构与物理闸门全绿：块全 CM、Δ=0 逐位 ≡ normal（maxd==0.0）、
  T3 ED 锚翻绿（1.4e-4/6.5e-5 = 一阶标定档 ✓ δt 减半减半）；polarif 批
  @test_broken 的缺陷即移位口径错位，就此闭环。R5（一阶常数回归）按标定
  预期，无额外回归。

**Phase W（T7 + goldens）**：
- 套件闸门：order=2 较 order=1 偏差收缩 ≥10×（实测 ~40×，Nt=8
  1.2e-4 → 1.2e-5）——Weyl-2 × BCS 端到端成立，:z2/U≠0（Anderson）稳定；
- `validation/bcs_realtime_scaling.jl`：order=1 斜率 1.07–1.09 ✓；
  **order=2 渐近斜率 ≈ 2.00 未达**——Nt≥8 撞 ~1e-5 δt 无关地板（紧截断
  不变、随 Nt 线性缓增 = 每格常数残差特征），疑二阶接触门 × BCS 交互
  （R7 触发，待排查；normal bath 九组基准无此地板）；
- **标定方法两陷阱（已修）**：ED 参照须按时刻对齐切片（`gt_ed[(1:Nt).*stride]`
  ——首版半时刻错位伪造 order 无关大误差）；斜率口径 = 绝对偏差（relerr
  分母随窗口变化伪造斜率）；
- BCS goldens 新档 12/12 绿（raw 口径四组合，1e-6 档）；
- 排除项：K 链截断（紧截断逐位不变）、IF 截断（D512 不变）、ED 初态口径
  （三口径同偏差）、反常通道（Δ=0 同偏差）——地板在 order=2 装配的接触/
  测度层，非 BCS 反常专属（normal 二阶同框架待对照定位）。
