严格 IF 实时核理论：Keldysh 轮廓上的逐步推导
===

**benchmark/exactif/theory/keldysh.md**（2026-11-07，详导版）
配套实现：`src/NativeExact.jl`；验证入口：`validation/verify_compact.jl`、
`verify_gtlt_align.jl`、`verify_align_u15.jl`、`scaling_compact_models.jl`、`zidentity.jl`。

**本文约定**（学习自 `benchmark/weyl-secondorder.html` 的讲解标准）：

- 每一个新变量、新算符在首次出现处定义；推导不跳跃——每个等号给理由；
- 先在最简单的情形（单费米模、单 bath 模、单接口）把每条规则**显式算出**，
  再陈述多模推广及推广成立的理由；
- 明确区分两种陈述：**精确恒等式**（无任何近似，数值上机器精度成立）与
  **近似/约定**（Trotter 分割、截断、以及绝对符号这类由证书锁定的簿记）；
- 绝对符号（Berezin 积分取向、整体正负）是纯约定；实现里的绝对符号由
  Z 恒等式证书与 Fock 链证书**锁定**。本文的推导只涉及相对结构与相对符号，
  并在每处标注"由哪个证书锁定"。

---

## 0. 模型与目标

### 0.1 模型

杂质模型（单带推导，多带在 §4.7 推广）：

\[
H = H_{\mathrm{imp}} + H_q ,
\]

逐项定义：

- **H_imp**：杂质哈密顿量，只含杂质算符。单带时
  `H_imp = ε_d·n`（`n = d†d`，`d` 为杂质湮灭算符）；双带 Anderson：
  `H_imp = ε_d(n₁+n₂) + U·n₁n₂`。**关键性质：H_imp 一般不是二次型**
  （U 项是四费米算符）——这是杂质部分必须走精确门、不能并入高斯积分的
  根本原因（§2.1）。
- **H_q ≡ H_bath + H_hyb**：二次型部分。逐项：
  `H_bath = Σ_α ε_α c_α†c_α`（bath 模式 α 的能级 `ε_α`、湮灭算符 `c_α`）；
  `H_hyb = Σ_α (d†V_α c_α + c_α†V_α d)`（耦合幅 `V_α`，复数）。

`H_q` 是单粒子二次型：`H_q = Σ_{ij} h_{ij} γ_i† γ_j`，其中 `γ_i` 是全部
湮灭算符（`γ_0 = d`、`γ_{α} = c_α`），单粒子矩阵

\[
h = \begin{pmatrix} 0 & V^\dagger \\ V & H_{\mathrm{bath}} \end{pmatrix} ,
\]

**杂质对角元为 0**——杂质的自身演化（含 ε_d 与 U）全部归入 `H_imp`。
这个分配不是唯一可能，但它使 `h` 的 bath 侧只含 bath 自由度，是 §2 高斯
积分成立的条件；杂质自身演化由 §4 的精确门承担。

### 0.2 初态与目标观测量

初态可分离：

- bath 热态：`ρ_b = e^{−βH_bath}/Z_b`（β 逆温、`Z_b = tr e^{−βH_bath}`）；
- 杂质初态：`ρ_0 = e^{−βH_imp}/tr(e^{−βH_imp})`（热态；任意初态同样支持）。

目标：Keldysh 轮廓 `0 → T`（前向 + 支）`→ 0`（后向 − 支）上的实时观测量：

- 末态占据数 `n(T)`：`ρ(T) = Tr_bath[U ρ_tot(0) U†]` 的 `ρ[1,1]` 类对角元；
- 双时 Green 函数（裸关联形式；物理 gt/lt 差一个 `∓i` 因子）：
  `greater(k) = ⟨c(T) c†(τ_k)⟩`、`lesser(k) = ⟨c†(τ_k) c(T)⟩`，
  采样时刻 `τ_k` 见 §4.6；
- 全部观测量都是**比值**：分子 = 带插入的收缩，分母 `Z` = 无插入收缩——
  归一化常数（§2.6 的 norm）在比值中消去。

### 0.3 总路线

| 步骤 | 积掉的对象 | 方法 | 输出 | 实现位置 |
|---|---|---|---|---|
| （1） | 物理 bath 的相干场 | 逐 cell 精确高斯积分 + Schur | 接口变量上的 Gaussian `G` | `contour_G` / `compact_delta` |
| （2） | IF（接口）Grassmann 变量 | 带源高斯积分（Pfaffian 恒等式） | 杂质轨迹腿上的消元核 `D` | `delta_eff` 恒等式 |
| （3） | ——（组装） | 调度门链 + 消元核 + 测度/闭合 | 收缩 → 读出 | `sysdynamics_exact!` + `hybriddynamics` + `gtlt_compact` |

（1）（2）是**精确恒等式**（数值验证：（2）逐系数 ≤2.2e-16；（1）的 G 与
IFGTEMPO 逐位一致）。（3）携带全部 Trotter 误差：order 1/2 全局 O(δt²)、
order 4 全局 O(δt⁴)（§4 给出推导与数值）。

记号约定：每条 Grassmann 腿一个生成元；`x̄` 型腿记 creation 类、`x` 型腿记
annihilation 类；接口（IF）变量对记 `a = (ā, a)`，杂质轨迹变量对记
`x = (x̄, x)`；多带时带 flavor 下标。

---

## 1. 预备：费米相干态与积分规则

本节全部在**单费米模**上手推（Fock 空间只有 `|0⟩`、`|1⟩ = d†|0⟩` 两个态）；
多模推广在 §1.4 给出理由。本节的每条规则都在后续章节反复使用。

### 1.1 相干态定义与基本恒等式

**定义**（本文采用的正号约定）：

\[
|x⟩ = e^{x d†}|0⟩ = |0⟩ + x|1⟩ ,
\qquad
⟨\bar x| = ⟨0| e^{-\bar x d} = ⟨0| - \bar x⟨1| .
\]

`x, x̄` 是**外部 Grassmann 生成元**：相互反对易（`x x̄ = −x̄ x`）、与物理
算符 `d, d†` 反对易、平方为零（`x² = 0`）。指数展开只有两项：
`e^{x d†} = 1 + x d†`（`(x d†)² = 0`）。

**恒等式 （1）（ket 本征关系）**：`d|x⟩ = x·|x⟩`。

推导：`d|x⟩ = d(|0⟩ + x|1⟩) = d|0⟩ + x·(d|1⟩)`。`d|0⟩ = 0`；
`d|1⟩ = d d†|0⟩ = (1 − d†d)|0⟩ = |0⟩`。故 `d|x⟩ = x|0⟩`；而
`x|x⟩ = x|0⟩ + x²|1⟩ = x|0⟩`。（证毕）

**恒等式 （2）（bra 本征关系）**：`⟨\bar x| d^\dagger = -\bar x·\langle\bar x|`
对 `|0\rangle` 分量成立，即 `\langle\bar x|d^\dagger|0\rangle = -\bar x·\langle\bar x|0\rangle = -\bar x`。

推导：`\langle\bar x|d^\dagger|0\rangle = (\langle0| - \bar x\langle1|)d^\dagger|0\rangle = \langle0|d^\dagger|0\rangle - \bar x\langle1|d^\dagger|0\rangle`。
`\langle0|d^\dagger|0\rangle = \langle0|1\rangle = 0`（`d^\dagger|0\rangle = |1\rangle`）；`\langle1|d^\dagger|0\rangle = \langle1|1\rangle = 1`。
故 `= 0 - \bar x·1 = -\bar x`。（证毕）

（（2）的负号与（1）的正号是同一约定的两面；它们将在 §1.3 的迹规则里相互配合。）

**恒等式 （3）（重叠）**：`\langle\bar x|x\rangle = 1 - \bar x x`（= `e^{-\bar x x}`）。

推导：`\langle\bar x|x\rangle = (\langle0| - \bar x\langle1|)(|0\rangle + x|1\rangle)`
`= \langle0|0\rangle + \langle0|1\ranglex - \bar x\langle1|0\rangle - \bar x\langle1|1\ranglex`
`= 0 + 0 - 0 - \bar x·1·x`。末项 `\langle1|x\rangle = \langle1|d^\dagger|0\rangle = \langle1|1\rangle = 1`，
故 `= -\bar x x`。（证毕）

### 1.2 分辨率与逐对积分规则

**逐对积分规则**（声明为约定，紧接着用分辨率验证其一致性）：

| 单项式 | 积分值（`\intd\bar x\, dx` 次序） |
|---|---|
| `1` | 0 |
| `x̄` | 0 |
| `x` | 0 |
| `x̄·x` | **-1** |

**分辨率**：

\[
\int d\bar x\, dx\; e^{-\bar x x}\, |x\rangle\langle\bar x| = 单位算符 .
\]

验证（展开 `e^{-\bar x x} = 1 - \bar x x`，作用到两组基矢）：

- 作用于 `|0\rangle`：被积式 `= e^{-\bar x x}(|0\rangle + x|1\rangle)(\langle0| - \bar x\langle1|)`。
  保留积分值非零的项（须含恰好一个 `x̄` 与一个 `x`）：`|0\rangle\langle0|` 项的系数
  `-\bar x x` 积分 `= -(-1) = +1` \to 贡献 `|0\rangle\langle0|` （成立）；
  `|0\rangle\langle1|` 项：含 `x` 但需配 `x̄`——`-x\bar x\langle1|` 的 `x\bar x = -\bar x x`，
  积分 `= +\bar x x\langle1| = -\langle1|` \to 贡献 `-|0\rangle\langle1|·(-1) = +|0\rangle\langle1|` （成立）；
  `|1\rangle\langle0|` 项：`+x|1\rangle\langle0|(1-\bar x x)`：`x` 项积分 `\intx = 0`、`-x\bar x x`
  中 `x\bar x x = 0`（两个 `x`）\to 0 （成立）；
  `|1\rangle\langle1|` 项：`-x|1\rangle\langle1|\bar x(1 - \bar x x)`：`x\bar x = -\bar x x` \to
  `+\bar x x·\bar x = \bar x·\bar x·x = 0`……直接展开
  `-x\bar x(1 - \bar x x) = -x\bar x = +\bar x x` \to 积分 `= -1` \to 贡献
  `-|1\rangle\langle1|·(-1) = +|1\rangle\langle1|` （成立）。
- 作用于 `|1\rangle`：同法逐项（`\langle\bar x|1\rangle = -\bar x` 进入展开）\to 合计 `|1\rangle\langle1|` （成立）。

两条基矢都还原，分辨率成立。**`\int\bar x x = -1` 的负号正是使分辨率为 +1 的
约定**；改用 `dx\,d\bar x` 次序则所有符号翻转、公式组同样自洽——绝对符号是
约定，相对结构不是。

### 1.3 迹规则与反周期边界

**迹规则**：

\[
\operatorname{Tr}[O] = \int d\bar x\, dx\; e^{-\bar x x}\; \langle-\bar x|\,O\,|x\rangle ,
\]

bra 端边界取 `⟨−\bar x|`（**反周期边界**）：把 `x̄ \to -x̄` 代入 bra 定义，
`\langle-\bar x| = \langle0|e^{+\bar x d} = \langle0| + \bar x\langle1|`（展开中 `(+)` 号与 §1.1
同理由）。

**验证**（单模、任意 `\alpha`，含 `Tr[单位算符] = 2` 的检查）：
`e^{\alphan}|x\rangle = |0\rangle + e^{\alpha}x|1\rangle`；
`\langle-\bar x|0\rangle = 1`；`\langle-\bar x|1\rangle = +\bar x`（`\langle-\bar x|d^\dagger|0\rangle = +\bar x`，
由（2）代入 `x̄ \to -x̄`：bra 定义里 `(+)` 号给出）。
`\langle-\bar x|e^{\alphan}|x\rangle = 1 + e^{\alpha}x\bar x = 1 - e^{\alpha}\bar x x`。
乘测度因子 `e^{-\bar x x} = 1 - \bar x x`：

\[
(1 - \bar x x)(1 - e^{\alpha}\bar x x) = 1 - (1 + e^{\alpha})\bar x x + (\bar x x)^{2}(= 0) .
\]

积分：`0 − (1 + e^{α})·(−1) = 1 + e^{α}`。Fock 直算 `Tr[e^{αn}] = e^{α·0} + e^{α·1} = 1 + e^{α}`
（成立）。特别地 `α = 0` 给 `Tr[单位算符] = 2` （成立）。

**为什么边界必须是反周期的**：迹对单模求和 `Σ_{n=0,1}⟨n|O|n⟩` 在相干态表示
中把 bra 端的 `x̄` 与 ket 端的 `x` 配对；配对规则（§1.2）要求两个变量以
`(x̄, x)` 顺序相邻，而 bra 定义自然给出 `⟨+x̄|`——为落到 §1.2 的积分规则须
把 `x̄` 换成 `−x̄`（Grassmann 变量取负 = 翻转其线性项符号）。这就是"反周期"
的全部内容；多模时逐模反周期。

### 1.4 多模推广

多模（1 杂质 + m bath）时每模一对 (x̄, x)，外部生成元全体相互反对易；
分辨率与迹规则是逐模乘积。**分次乘积注意**：多模算符串重排的费米子符号 =
置换符号；本文涉及多模置换绝对符号处均由 §5 的证书锁定，推导只使用
"每生成元在单项式中至多出现一次"的幂零性。

### 1.5 高斯映射引理（cell 传播子的相干矩阵元）

**引理**：对单粒子矩阵 `h` 与 `U = e^{-iτh}`（分块：杂质行/列 = 1，bath = 其余），
二次型算符 `B(τ) = e^{-iτh}` 的相干态矩阵元是外部变量的高斯型指数，
其系数恰为 `U` 的块：

\[
\langle\bar x', \bar a'\,|\,B(\tau)\,|\,x, a\rangle \propto
\exp\big[\, \bar x' U_{11} x \;+\; \bar x' U_{12} a \;+\; \bar a' U_{21} x \;+\; \bar a' U_{22} a \,\big] ,
\]

`∝` 号吸收的常数前因子最终进入 `Pf(G)`（由数值 Pfaffian 计算）。

**单模（单 bath 模 + 杂质）验证**：取 `h` 对角（`V = 0`）：
`U = diag(1, e^{−iτ ε})`，即 `U_{11} = 1`、`U_{22} = e^{−iτ ε}`、交叉块 0。
引理给出矩阵元 `∝ exp[\bar x'·1·x + \bar a' e^{-i\tau\epsilon} a]`。直接展开：
`e^{-i\tauh} = e^{-i\tau\epsilon c^\daggerc} = 1 + e^{-i\tau\epsilon} c^\daggerc`（生成元幂零截断），矩阵元
`= \langle\bar x', \bar a'|x, a\rangle + e^{-i\tau\epsilon}\langle\bar x', \bar a'|c^\daggerc|x, a\rangle`
`= e^{-\bar x' x - \bar a' a} + e^{-i\tau\epsilon}·(\bar x' x)(\bar a' a)·(费米符号)`。

两条幂零展开在"最多一个 `x̄`、一个 `x`、一个 `\bar a`、一个 `a`"的单项式上
与引理的指数展开一致（指数展开 `= 1 + \bar x'x + e^{-i\tau\epsilon}\bar a'a + 高阶(0)`），
常数因子 `e^{-\bar x'x - \bar a' a}` 正是 §1.2 型的边界分辨率——在链的合成中
成为相邻 pair 的测度（§2.2）。含 `V` 耦合的单模情形由同一幂零展开逐系数
验证；其整体正确性由 §5 的两条机器精度证书锁定。

**本文使用规则**：§2.3 的 Q-kernel 缝合处按引理写入 `U` 块系数；常数前因子
进 `Pf` 通道。**方法论**：结构由推导给出，绝对符号由证书锁定。

### 1.6 高斯积分规则（带源）

**规则**（标准费米高斯积分；n=2 展开验证见 §3.2）：

\[
\int \mathcal{D}a\; e^{\frac12 a^{\mathsf T} G a + j^{\mathsf T} a} = \Pf(G)\; e^{-\frac12 j^{\mathsf T} G^{-1} j} ,
\]

`a` = Grassmann 生成元列向量、`G` = 反对称系数矩阵（`G^ᵀ = −G`）、`j` = 外部
源列向量、`Pf(G)` = `G` 的 Pfaffian（`Pf(G)² = det G`）。多模情形即多对
(x̄, x) 的乘积形式，§2.4 由它推导 Schur 公式。

---

## 2. 第一步：积掉 bath

### 2.1 链的搭建：pair 变量与逐对测度

把实时演化按调度切成 cell：`U 链 = A(ℓ_{N_c+1}) B(τ_{N_c}) A(ℓ_{N_c}) … B(τ_1) A(ℓ_1)`。
每个 `B(τ_c)` 是 bath 自由度的高斯映射；相邻 `B` 的 bath 相干场通过**边界
分辨率**（§1.2 的多模版）缝合。由 §1.3 的迹规则（多模版），整条链的 bath 迹：

\[
\operatorname{Tr}_{\mathrm b}\big[ B(τ_{N_c}) … B(τ_1) \big] =
\int \mathcal{D}(\bar \gamma, \gamma)\;
e^{\frac12 (\bar \gamma, \gamma)^{\mathsf T} \mathcal{B} (\bar \gamma, \gamma)}
\times \prod_{\text{pair}} (-\bar \gamma \gamma)\ \big|_{\text{热边/迹边闭合}} ,
\]

逐符号定义：

- `(γ̄_α(j), γ_α(j))`：第 `j` 个 cell **边界**上的 bath pair 变量（前向支
  `+`、后向支 `−` 各一套；branch 下标略）；
- `B 块`：由各 `B(τ)` 的 `U_{22}` 块（bath 内部传播）组装的 bath 块矩阵
  （组装算法见 §2.4）；
- `∏(−γ̄γ)`：逐对测度。

**为什么每 pair 必须带测度 `−γ̄γ`**：相邻两次分辨率闭合（§1.2）恰好把链
缝合——单模验证：`∫ e^{-\bar x x} \langle-\bar x|x\rangle · ([-\bar x|` 链 `|x\rangle])` 型
双分辨率闭合正是 §1.2 分辨率的成对形式；缺一个测度相当于缺一次分辨率，
实现上表现为 `Z \approx \tfrac{1}{2}·(\text{正确值})`（每个缺失边界丢一半相干态重叠）——
与"\tau=0 恒等门不可跳"（§4.5）同源：恒等门与测度都是分辨率的载体。

### 2.2 热边与迹闭合

热边：bath 热态的相干表示（§1.5 引理单模验证同型）：

\[
\rho_b = Z_b^{-1} e^{-\betaH_{\mathrm{bath}}} \;\Rightarrow\; ⟨\bar \gamma(1,-)| \rho_b |\gamma(1,+)\rangle \propto e^{-\beta \sum_\alpha \epsilon_\alpha \bar \gamma_\alpha(1) \gamma_\alpha(1)} ,
\]

（`Z_b^{-1}` 进 `ν`，见 §2.6。）迹闭合：整条 bath 链的两端由 §1.3 的反周期
边界规则闭合——compact 布线中实现为**共享 trace pair**（远端一对 bath pair
由 forward 末核与 backward 末核共用）+ 末个 forward 核输出行的 `⟨−ψ̄|`
翻转；其绝对符号由 Z 证书锁定。

### 2.3 杂质接口作为"源"：Q-kernel 的四个块

每个 cell `B(τ_c)` 通过 `H_hyb` 把**杂质接口腿**耦合进高斯（§1.5 引理的
四个块；接口 `c` 的杂质腿对 = `(\bar x_c, x_c)`，`x_c` 即该接口的轨迹腿）：

- `U[1,1]`：接口 pair `(\bar x_c, x_c)` 之间的系数。**`U[1,1] = 1`**：
  由 `h` 的杂质对角元为 0（§0.1）——即接口分辨率恒等；
- `U[1,1+\alpha] = \bar x_c \times` bath 输入 eta：**入射幅**（`d^\daggerV_\alpha` 项：杂质产生、
  bath 湮灭）；
- `U[1+\alpha,1] =` bath 输出 bar `\times x_c`：**出射幅**（`V_\alpha d` 项）；
- `U[1+\alpha,1+\beta]`：bath 内部传播（对角 `\epsilon_\alpha`；多模 bath 时含 `H_bath` 块）。

**为什么杂质对角元可以是 0**：`B(\tau_c)` 只承担 `H_q` 部分的演化；杂质自身
演化由骨架门 `A(\ell) = e^{-i\ellH_imp}` 承担（§4.1）。

### 2.4 组装算法与 Schur 消元的推导

实现（`compact_delta`）把**接口变量与 bath 变量放进同一个大反对称矩阵**
`A`（接口块 `E`、交叉块 `C`、bath 块 `B`；`put_G` 逐项：上三角 `+= v`、
下三角 `-= v`；§2.1 的 `-\gammā\gamma` 测度与 §2.2 的热/迹边都作为线性项进入 `A`），
然后**一步 Schur 消掉全部 bath 变量**：

\[
\int \mathcal{D}(\text{bath 变量})\; e^{\frac12 (\text{接口},\text{bath})^{\mathsf T} A (\text{接口},\text{bath})} = \Pf(B)\; e^{\frac12 (\text{接口})^{\mathsf T} \big(E + C B^{-1} C^{\mathsf T}\big) (\text{接口})} ,
\]

**推导**（分块配方法，从 §1.6 规则直接读出）：把指数按块写出

\[
\tfrac12 \eta_{\mathrm{if}}^{\mathsf T} E \eta_{\mathrm{if}}
+ \eta_{\mathrm{if}}^{\mathsf T} C\, b_{\mathrm{bath}} + \tfrac12 b^{\mathsf T} B b ,
\]

对 bath 变量 `b` 用 §1.6 规则（源 = `C^{\mathsf T}\eta_{\mathrm{if}}`，线性项
系数矩阵 = `B`）：积分给 `Pf(B)·exp(+\tfrac{1}{2}\eta_{\mathrm{if}}^{\mathsf T} C B^{-1} C^{\mathsf T} \eta_{\mathrm{if}})`。
合并接口块即得 Schur 公式。
（反对称块约定使 `C`/`C^{\mathsf T}` 的摆法有一个记法差；数值上实现为
`effective = E + C(B \setminus C^{\mathsf T})`，正确性由 §5 证书锁定。）

随后 `G = antisym(effective)`（反对称化：Gaussian 指数只需反对称部分）。
**对照**：G 与 IFGTEMPO `influence_matrix_G` 逐位一致（order 2 4\times4、
order 4 含负 cell 12\times12，max|ΔG| = 0）。

（等价观点：也可以**逐 cell 顺序 Schur**——每乘一个 cell 把该 cell 的 bath
变量积掉、接口集合向前扩展一格。两种次序结果相同——高斯积分次序无关；
实现取"一次大 Schur"。）

### 2.5 U(1) 结构与同 kind 零

`H_q` 数量守恒 \Rightarrow 装配的每一项都连接一个 creation 类槽位与一个
annihilation 类槽位（§2.3 四块全是跨类对）\Rightarrow 大 `A` 与其 Schur
结果保持结构 \Rightarrow

\[
G_{(\bar,\bar)} = G_{(a,a)} = 0 \quad\text{对所有接口对（全 order 数值验证 0）}。
\]

这条性质是 §3 行扫表示可行性的关键。

### 2.6 归一化 ν

Schur 抛掉的常数（`Pf(B)` 等）连同 bath 热态归一化构成

\[
\nu = Z_N/Z \;=\; \frac{\operatorname{tr}[\,\text{bath 链的二次量子化迹}\,]}
{\operatorname{tr}\big[e^{-\betaH_{\mathrm{bath}}}\big]} ,
\]

（`bath_norm`；分子 = 链上逐 cell `U_{22}` 块乘积的二次量子化迹——即 bath
配分函数在链传播下的值。）闭合配分函数 `Z = integrate · norm`；
比值型观测量中 `ν·Pf` 消去。

**本步小结**：bath 积掉后，剩余多项式只含杂质轨迹腿：

\[
P(x) = Z^{-1}\; e^{\frac12 x^{\mathsf T} G x}\; K[x]\cdot(\text{测度/闭合已并入}),
\]

——在 compact 实现里这是最终形态（接口腿即轨迹腿）；在两扇区记法中它还带
接口交叉测度，其消去是 §3 的内容（两种记法的等价由 §5 证书链锁定）。

---

## 3. 第二步：积掉 IF Grassmann 变量

### 3.1 接口交叉测度作为源（逐分量）

两扇区记法（`a` = 接口 IF 变量、`x` = 杂质轨迹变量，二者由接口交叉测度
连接）：

\[
P(a,x) = \nu\cdot\Pf(B)\; e^{\frac12 a^{\mathsf T} G a}
\cdot \prod_c e^{-\bar a_c\, x^{\eta}_c - \bar x^{\bar}_c\, a_c}
\cdot K[x]\cdot M_0 .
\]

逐符号：`c` 遍历接口（= cell 边界，`c = 1..N_c`）；`x^{η}_c` = 接口 `c` 上
杂质轨迹的 annihilation 腿（前向支），`x^{\bar}_c` = creation 腿（后向支）；
交叉测度是**两扇区构造的接口恒等分辨率**（每个 U_imp–U_bath 接口一个，逐
接口独立因子——它的非独立性正是 §2 的 Schur 要一次做完的原因）。

把交叉测度读成 Gaussian 的**线性源**：

\[
j[\bar a_c] = -x^{\eta}_c ,\qquad j[a_c] = -\bar x^{\bar}_c ,
\]

即**源 = partner 杂质腿**（bar↔eta 交叉；"partner 槽位"映射的定义）。

### 3.2 带源高斯积分（Pfaffian 恒等式的 n=2 展开）

**规则**（§1.6；n=2 展开验证如下）：对单接口，`a = (ā, a)` 两个生成元，
`G = [[0, g],[−g, 0]]`（反对称 2×2，`g = G_{\bar a, a}` 即该接口的高斯幅），
源 `j^{\mathsf T} a = j_{\bar}\bar a + j_a a`，其中 `j_{\bar} = -\theta(\thetā 的源 = impurity-eta 腿 \theta)`、
`j_a = -\thetā(\thetā 腿)`——把杂质腿
也记为生成元 \theta（它们是**同一 Grassmann 代数**的成员，这是路径积分的机制：
所有腿都是生成元）。

逐项展开（全部只到二阶；更高阶因生成元重复而为零）：

- `exp(\tfrac{1}{2}a^{\mathsf T}Ga)`：`\tfrac{1}{2}a^{\mathsf T}Ga = \tfrac{1}{2}(g\bar a·a + a·(-g)\bar a) = \tfrac{1}{2}(g\bar a a + g\bar a a)`
  = `g\,\bar a\, a`（用了 `a\bar a = -\bar a a`），故 `= 1 + g\,\bar a\, a`；
- `exp(j^{\mathsf T}a) = 1 + j_{\bar}\bar a + j_a a + j_{\bar}j_a\,\bar a a·(符号)`。
  `(j^{\mathsf T}a)^{2}` 项：`j_{\bar}\bar a·j_a a` 重排符号分析——两个源的积
  `(-\theta)(-\thetā)\bar a a = \theta\thetā\bar a a`：把 `\thetā` 移过 `\bar a`（反号）再移过
  `a`（反号）……净符号 `(+1)`，得 `+\theta\thetā\bar a a`；
- **乘积**：`(1 + g\bar a a)(1 + \theta\thetā\bar a a)` 中含 `\bar a a` 的项：
  `g\bar a a` 与 `\theta\thetā\bar a a`；二者之和 = `(g + \theta\thetā)\bar a a`。
  注意 4 生成元项 `(g\bar a a)(\theta\thetā\bar a a) = 0`（`\bar a` 重复）。
- **积分**（按 §1.2 型逐对规则 `\int\bar a\, a = +1`，符号即 `elim_sign` 的
  约定通道）：`\int(g + g\theta\thetā)\bar a a = g·(1 + \theta\thetā)`。

**结果**：单接口消元给出 `g·(1 + \theta\thetā)`——即杂质腿对上的核
`g·\theta\thetā` 加常数 `g`。**一般 n**：`exp(j^{\mathsf T}a)` 展开到二阶
（更高阶为零），其二阶项 `\tfrac{1}{2}(j^{\mathsf T}a)^{2}` 与
`\tfrac{1}{2}a^{\mathsf T}Ga` 的交叉把 `\tfrac{1}{2}j^{\mathsf T}G^{-1}j` 吸收进指数
（标准费米高斯积分），得规则式；杂质腿对的核由 `G^{-1}` 的 partner 读出给出
（§3.3）。

**本步是精确恒等式**：实现层面由 Step B 暴力 Berezin 积分逐系数对照锁定
（\leq2.2e-16，orders 1/2/4 全验）。

### 3.3 消元核 D 的读出（partner 映射）

一般 `N_c` 接口、4\timesN_c 个杂质腿（每接口 `x̄⁺, x⁺, x̄⁻, x⁻` 四条）：

\[
D[\lambda,\lambda'] = -\tfrac12\big(G^{-1}_{ik} - G^{-1}_{ki}\big),\qquad
i = \text{slot(}\bar\lambda\text{)},\ k = \text{slot(}\lambda'\text{)} ,
\]

逐符号：`λ, λ'` = 两条杂质轨迹腿；`slot(λ̄)` = `λ` 的 **partner IF 槽位**
（§3.1 的源映射：bar↔eta 交叉——`λ = x̄⁺(c)` 的 partner 是 IF 槽位
`(c,+,eta)`，`λ' = x⁺(c′)` 的 partner 是 `(c′,+,bar)`）；`−½(·−·)` 是把
`G⁻¹` 反对称化后按 (bar→eta) 有向读出的**反对称化因子**；绝对符号
`elim_sign = −1` 由 Step B 证书锁定（纯约定）。

### 3.4 同 kind 零与行分解（Q1 判决）

§2.5 的 U(1) 结构经 Schur 与求逆保持 ⇒ `G⁻¹` 同 kind 元为零 ⇒ `D` 同 kind
元为零（全 order 数值验证 0 个）⇒

\[
e^{\sum D} = \prod_{\text{bar } r}\Big(1 + \theta_{\bar r}\sum_c D[r,c]\,\theta_c\Big)
\]

按行分解——每行一个 partialmpo、键维 2，`hybriddynamics` 行扫完整承载
（Q1 判决）。**消元精确、无近似**；代价是 `D` 非局域（IF 键维 4^{N_c}，
生产上键维截断）。

### 3.5 两种记法的关系（一步 vs 两步）

compact 实现把 §2（Schur）与 §3（Pfaffian 消元）**融合为一次大 Schur**：
接口槽位直接取为轨迹腿（交叉测度成为 Gaussian 的对角恒等对 `U[1,1] = 1`，
其与轨迹腿的连接由骨架链承担）。两种记法的最终多项式逐位等价——等价性的
证书链：（1）G vs IFGTEMPO 逐位；（2）消元核 vs 暴力积分逐系数；（3）管线
vs Fock 离散链机器级；（4）管线 vs IFGTEMPO 机器级。§3 的两步记法用于推导
与验证，§2 的融合记法用于实现。

---

## 4. 杂质部分与 1/2/4 阶组装

### 4.1 杂质核 K 链

杂质扇区**不走高斯积分**（`H_imp` 含 `U` 项，非二次型）——走精确 Fock 门：

\[
K[x] = \langle\bar x_0^{+}|\rho_0|x_0^{-}\rangle
\prod_{t=1}^{N_c+1} \langle\bar x_t^{+}|A(l_t)|x_{t-1}^{+}\rangle
\prod_{t=N_c+1}^{1} \langle\bar x_{t-1}^{-}|A(l_t)^\dagger|x_t^{-}\rangle ,
\qquad A(\ell) = e^{-i\ell H_{\mathrm{imp}}} .
\]

`A(ℓ)` 对多带含 `U` 时是 `2^{bands}` 维 Fock 矩阵的精确指数——**门内无
近似**；近似只在门与 cell 的分割方式（调度）里。

### 4.2 调度表

| order | cells（进 G） | links（进 K） | 全局阶 |
|---|---|---|---|
| 1 | `fill(δt, N)` | `[δt×N, 0]`（末位 0：开放端在末时刻） | O(δt²)（机理 §4.4） |
| 2 | `fill(δt, N)` | `[δt/2, δt×(N−1), δt/2]` | O(δt²)（对称，§4.3 推导） |
| 4 | `[a,b,a]δt × N` | `[a/2, c, c, (a,c,c)×(N−1), a/2]` | O(δt⁴)（§4.3 推导） |

`links` 比 `cells` 多一个（链两端各一门）。order 1 links 末位 0 的恒等门
不可跳（§4.5）。

### 4.3 回文与阶条件（order 2 与 4 的推导）

记 `X = −iτ₁A`、`Y = −iτ₂B`、`Z = −iτ₃A`（τᵢ = 各 A/B 因子时长）。对
`e^{X}e^{Y}e^{Z}`，BCH 展开的指数二阶项：

\[
\tfrac12[X,Y] + \tfrac12[Y,Z] + \tfrac14\big([X,Z] + [Z,X]\big) .
\]

order 2 的链 `X = −iδtA/2`、`Y = −iδtB`、`Z = −iδtA/2`：`[X,Z] = 0`（同为
`A` 向），且 `[Y,Z] = [Y,X]·(τ 相等)`，于是

\[
\tfrac12[X,Y] + \tfrac12[Y,X] = 0 :
\]

**δt² 项恰好相消**（"对称/回文"的代数含义：首尾两个 A-半步相等使两个
commutator 贡献互抵）。残留从 δt³ 起 ⇒ 单步 O(δt³)、N 步全局 O(δt²)。

**order 4**：回文 Suzuki `S4 = S2(aδt)·S2(bδt)·S2(aδt)`，合并链（相邻 A 因
子相加：`a/2 + b/2 = c` 等）后逐宏步

\[
A(a/2)\, B(a)\, A(c)\, B(b)\, A(c)\, B(a)\, A(a/2) ,
\qquad c = \frac{a+b}{2} .
\]

权重条件：

\[
2a + b = 1 \quad(\text{一阶：总时长})；
\qquad 2a^3 + b^3 = 0 \quad(\text{三阶：BCH 第三项的正比系数为零})，
\]

解出 `a = 1/(2−∛2) ≈ 1.2252`、`b = 1−2a ≈ −1.4505 < 0`。**回文性的角色**：
组合中 S2 因子序正读倒读一致 ⇒ 所有奇幂 BCH 余项成对相消（§4.3 的二阶
推导是其二阶实例）；**非回文序 `S2(a)²·S2(b)`** 的 δt³ 项不再满足该条件
——实测比值序列 12.7→9.0（滑向 8，即 p→3）是其失效的数值表现；修正为
回文序后比值恢复 16（N=16 定点 15.97/16.12）。**负 cell 的必然性**：
`2a³ + b³ = 0` 要求 `b < 0`（两个正立方之和不能为零）——三阶相消必须有
反向步，接口格点因此非均匀且时间倒卷（macro 内边界
`t + a·δt → t + (a+b)·δt = t − 0.225·δt`）。

### 4.4 order 1 全局二阶的机理

order 1 的裸分裂 `A(δt)·B(δt)` 单步只有 O(δt)——但严格核管线的实测全局阶
是 **p = 2**。机理：分裂的领头误差 `[A,B]δt²/2` 是"杂质门与耦合演化的交错
缺失"；严格 cell `B(δt) = e^{-iδt h}` 的耦合顶点幅对 `V` 的任意幂精确，它
把这项交错**重新生成到所有阶**（D 的二阶交叉项正是 `[A,B]δt²` 的内容），
残留从更高阶起 ⇒ 全局 O(δt²)。冻结核（顶点冻结在 cell 内）没有这个重求和
⇒ 全局 O(δt)（实测 p ≈ 1）。最终裁判是数值：order 1 比值 4.00、
p = 2.01（Toulouse）/1.98（Anderson）。

### 4.5 组装约定逐条（每条的理由）

| 约定 | 内容 | 理由 |
|---|---|---|
| 每 gap 一门 | τ=0 恒等门不可跳 | 恒等传播子张量携带该 pair 的相干态分辨率 `∏(1+c̄c)`，是测度的一部分；跳过丢测度（Z 偏 ~½） |
| 接口反测度 `cc = −1` | 接口 pair 乘 `(1 − x̄x)` | 隐式 pair 分辨率被消元重复计入一次；幂零代数 `(1−x̄x)(1+x̄x) = 1` 精确抵消 |
| 归一化 `normdelta` | `ν·Pf(G)` 进 scaling | 闭合 `Z = integrate·norm`；比值型观测量自动消去 |
| 消元符号 `sgn = −1` | `M = −Ā` 的 partner 读 | Step B 暴力积分证书锁定（绝对符号纯约定） |
| 无 siteflip | M 不预翻 | 行扫 machinery（partialmpo）内部处理位点序符号；预翻双重 |
| 门符号 | 前向 `e^{−iℓH}`、后向 `e^{+iℓH}` | Keldysh 两支的物理演化方向 |

### 4.6 读出约定（闭环的出口；每条的理由）

- **n(T)**：折返点自由腿对 `(a⁺_k, ā⁻_k)` 的 2 腿 GTerm 投影。理由：轮廓
  折返把前向末腿与后向末腿配对，插入 `θ_{a⁺}θ_{ā⁻}` 与折返双线性共享腿、
  按幂零性自动相消——与稠密 ρ(T) 投影同一约定。`occupation` 的 insert_n
  带 O(δt) 有效时刻偏移，勿用；
- **⟨n₁n₂⟩**：band-1 + band-2 两条折返腿 GTerm 的乘积投影（幂零性保证乘积
  = 4 腿单项式）；
- **gt/lt**：`c(t_f)` @ 开放端 backward-eta 腿；
  `greater[k+1] = ⟨c(t_f)c†(k·δt)⟩` 用 GTerm(fwd\_bar@层 k+1, csite)；
  `lesser[k+1] = ⟨c†(k·δt)c(t_f)⟩` 用 GTerm(bwd\_bar@层 k, csite)；
  裸关联、Z 归一。**有效半步时刻**：门链把算符时刻推到相邻门的中点——
  `G^>` 的插入腿在层 k+1（时间 (k+½)δt）、`G^<` 在层 k（(k−½)δt）；与
  IFGTEMPO 同约定读出逐位对齐；
- **ρ(T)**：开放腿提取（FormB 路线）或折返跨越投影。

### 4.7 多带推广

bands > 1：IF 因子逐带复制（每带独立 bath 与 D）；杂质门用多带 Fock 传播子
（`2^{bands}` 维含 U 项）——**U 项把两带的轨迹腿耦合进同一门张量**（非每带
独立）；折返逐 band 施加。已验证：双带 Anderson U=1.5 的 n/⟨n₁n₂⟩/gt/lt 与
IFGTEMPO 直接同模型对齐（≤8.6e-14）。

---

## 5. 闭环验证阶梯（每个数字的出处）

| 阶梯 | 判据 | 结果 | 出处 |
|---|---|---|---|
| 恒等式级 | ε_d=0：Z=1（任意 order/N；测度/簿记证书） | 最差 4.6e-13 | zidentity |
| 代数级 | 暴力 Berezin ∫Da == ν·Pf·exp(D) | ≤2.2e-16 | formb_dense Step B |
| 结构级 | 消元核同 kind 项 / bar→eta 充分性 | 0 个 / PASS | samekind_check / restrictd_probe |
| 离散链级 | n/⟨n₁n₂⟩/ρ/gt,lt vs 同调度 Fock 直算 | ≤3.8e-15 / 3.8e-14 | verify_compact / formb_dense Step C |
| 对齐级 | vs IFGTEMPO 直接同模型（U=1.5 双带） | ≤2.2e-14（n）、≤8.6e-14（gt/lt） | verify_align_u15 / verify_gtlt_align |
| 连续极限级 | vs ED（有效时刻），order 1/2/4 | p = 2/2/4；比值→4/4/16 | scaling_compact_models / realtime |

各 order 一览（n(T) 标度）：Toulouse p = 2.01/2.01/3.94（渐近段
2.00/2.00/3.99）；Anderson p = 1.98/1.98/3.94（渐近段 1.997/1.997/3.990）；
order 4 比值序列 14.4→15.6→15.9→16.0（N=16 定点 15.97/16.12）。

两条已知的非误差项（预期行为）：①order 4 粗端（δt=0.4）误差仅略优于
order 1——Suzuki 负 cell 的大误差常数的预渐近区；②quapi 路线 order 1 N=1
的 n₀ 偏差——调度表"末位补零"约定在 N=1 清掉唯一 cell 的边缘情况。

## 6. 边界与延拓

- **quapi 路线**（点基 Δ）：仅 orders 1/2——Suzuki 的杂质门插在 bath cell
  之间，点基 Δ 的点对点直连假设不成立；
- **连续谱 bath**：原生管线接受有限模 (H_bath, V)（δ 谱 = 精确）；真连续谱
  需先谱离散化（求积节点，误差随模式数收敛）——IFGTEMPO 的连续谱路线
  （Volterra + 谱求积，无模式化）在均匀 δt 网格上工作，order 4 的非均匀
  cells 目前显式拒绝；
- **IFGTEMPO 的角色**：全部严格性验证以它为独立参照（n/⟨n₁n₂⟩/ρ/gt/lt 的
  机器级对齐）；src/ 与活跃验证链本身 IFGTEMPO-free。
