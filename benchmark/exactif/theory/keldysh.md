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

$$H = H_{\mathrm{imp}} + H_q ,$$

逐项定义：

- **H_imp**：杂质哈密顿量，只含杂质算符。单带时
  $H_{\mathrm{imp}} = \epsilon_d\, n$（$n = d^\dagger d$，$d$ 为杂质湮灭算符）；
  双带 Anderson：$H_{\mathrm{imp}} = \epsilon_d(n_1+n_2) + U\,n_1 n_2$。
  **关键性质：H_imp 一般不是二次型**（U 项是四费米算符）——这是杂质
  部分必须走精确门、不能并入高斯积分的根本原因（§2.1）。
- **H_q ≡ H_bath + H_hyb**：二次型部分。逐项：
  $H_{\mathrm{bath}} = \sum_\alpha \epsilon_\alpha c_\alpha^\dagger c_\alpha$
  （bath 模式 α 的能级 $\epsilon_\alpha$、湮灭算符 $c_\alpha$）；
  $H_{\mathrm{hyb}} = \sum_\alpha (d^\dagger V_\alpha c_\alpha + c_\alpha^\dagger V_\alpha d)$
  （耦合幅 $V_\alpha$，复数）。

$H_q$ 是单粒子二次型：$H_q = \sum_{ij} h_{ij} \gamma_i^\dagger \gamma_j$，
其中 $\gamma_i$ 是全部湮灭算符（$\gamma_0 = d$、$\gamma_{\alpha} = c_\alpha$），
单粒子矩阵

$$h = \begin{pmatrix} 0 & V^\dagger \\ V & H_{\mathrm{bath}} \end{pmatrix} ,$$

**杂质对角元为 0**——杂质的自身演化（含 $\epsilon_d$ 与 $U$）全部归入
$H_{\mathrm{imp}}$。这个分配不是唯一可能，但它使 $h$ 的 bath 侧只含 bath
自由度，是 §2 高斯积分成立的条件；杂质自身演化由 §4 的精确门承担。

### 0.2 初态与目标观测量

初态可分离：

- bath 热态：$\rho_b = e^{-\beta H_{\mathrm{bath}}}/Z_b$（β 逆温、
  $Z_b = \operatorname{tr} e^{-\beta H_{\mathrm{bath}}}$）；
- 杂质初态：$\rho_0 = e^{-\beta H_{\mathrm{imp}}}/\operatorname{tr}(e^{-\beta H_{\mathrm{imp}}})$
  （热态；任意初态同样支持）。

目标：Keldysh 轮廓 $0 \to T$（前向 + 支）$\to 0$（后向 − 支）上的实时观测量：

- 末态占据数 $n(T)$：$\rho(T) = \operatorname{Tr}_{\mathrm{bath}}[U \rho_{\mathrm{tot}}(0) U^\dagger]$ 的对角元；
- 双时 Green 函数（裸关联形式；物理 gt/lt 差一个 ∓i 因子）：
  $\mathrm{greater}(k) = \langle c(T) c^\dagger(\tau_k)\rangle$、
  $\mathrm{lesser}(k) = \langle c^\dagger(\tau_k) c(T)\rangle$，
  采样时刻 $\tau_k$ 见 §4.6；
- 全部观测量都是**比值**：分子 = 带插入的收缩，分母 $Z$ = 无插入收缩——
  归一化常数（§2.6 的 $\nu \cdot \operatorname{Pf}$）在比值中消去。

### 0.3 总路线

| 步骤 | 积掉的对象 | 方法 | 输出 | 实现位置 |
|---|---|---|---|---|
| （1） | 物理 bath 的相干场 | 逐 cell 精确高斯积分 + Schur | 接口变量上的 Gaussian $G$ | `contour_G` / `compact_delta` |
| （2） | IF（接口）Grassmann 变量 | 带源高斯积分（Pfaffian 恒等式） | 杂质轨迹腿上的消元核 $D$ | `delta_eff` 恒等式 |
| （3） | ——（组装） | 调度门链 + 消元核 + 测度/闭合 | 收缩 → 读出 | `sysdynamics_exact!` + `hybriddynamics` + `gtlt_compact` |

（1）（2）是**精确恒等式**（数值验证：（2）逐系数 ≤2.2e-16；（1）的 G 与
IFGTEMPO 逐位一致）。（3）携带全部 Trotter 误差：order 1/2 全局 $O(\delta t^2)$、
order 4 全局 $O(\delta t^4)$（§4 给出推导与数值）。

记号约定：每条 Grassmann 腿一个生成元；$\bar x$ 型腿记 creation 类、$x$ 型腿记
annihilation 类；接口（IF）变量对记 $a = (\bar a, a)$，杂质轨迹变量对记
$x = (\bar x, x)$；多带时带 flavor 下标。

---

## 1. 预备：费米相干态与积分规则

本节全部在**单费米模**上手推（Fock 空间只有 $|0\rangle$、$|1\rangle = d^\dagger|0\rangle$
两个态）；多模推广在 §1.4 给出理由。本节的每条规则都在后续章节反复使用。

### 1.1 相干态定义与基本恒等式

**定义**（本文采用的正号约定）：

$$|x\rangle = e^{x d^\dagger}|0\rangle = |0\rangle + x|1\rangle , \qquad ⟨\bar x| = ⟨0| e^{-\bar x d} = ⟨0| - \bar x⟨1| .$$

`x, x̄` 是**外部 Grassmann 生成元**：相互反对易（`x x̄ = −x̄ x`）、与物理
算符 `d, d†` 反对易、平方为零（`x² = 0`）。指数展开只有两项：
`e^{x d†} = 1 + x d†`（`(x d†)² = 0`）。

**恒等式 （1）（ket 本征关系）**：$d|x\rangle = x\cdot|x\rangle$。

推导：$d|x\rangle = d(|0\rangle + x|1\rangle) = d|0\rangle + x\cdot(d|1\rangle)$。
$d|0\rangle = 0$；$d|1\rangle = d\,d^\dagger|0\rangle = (1 - d^\dagger d)|0\rangle = |0\rangle$。
故 $d|x\rangle = x|0\rangle$；而 $x|x\rangle = x|0\rangle + x^2|1\rangle = x|0\rangle$。（证毕）

**恒等式 （2）（bra 本征关系）**：$\langle\bar x| d^\dagger = -\bar x\cdot\langle\bar x|$
对 $|0\rangle$ 分量成立，即 $\langle\bar x|d^\dagger|0\rangle = -\bar x\cdot\langle\bar x|0\rangle = -\bar x$。

推导：$\langle\bar x|d^\dagger|0\rangle = (\langle0| - \bar x\langle1|)d^\dagger|0\rangle = \langle0|d^\dagger|0\rangle - \bar x\langle1|d^\dagger|0\rangle$。
$\langle0|d^\dagger|0\rangle = \langle0|1\rangle = 0$（$d^\dagger|0\rangle = |1\rangle$）；
$\langle1|d^\dagger|0\rangle = \langle1|1\rangle = 1$。
故 $= 0 - \bar x\cdot 1 = -\bar x$。（证毕）

（（2）的负号与（1）的正号是同一约定的两面；它们将在 §1.3 的迹规则里相互配合。）

**恒等式 （3）（重叠）**：$\langle\bar x|x\rangle = 1 - \bar x x$（= $e^{-\bar x x}$）。

推导：$\langle\bar x|x\rangle = (\langle0| - \bar x\langle1|)(|0\rangle + x|1\rangle)$
$= \langle0|0\rangle + \langle0|1\rangle x - \bar x\langle1|0\rangle - \bar x\langle1|1\rangle x$
$= 0 + 0 - 0 - \bar x\cdot 1\cdot x$。末项 $\langle1|x\rangle = \langle1|d^\dagger|0\rangle = \langle1|1\rangle = 1$，
故 $= -\bar x x$。（证毕）

### 1.2 分辨率与逐对积分规则

**逐对积分规则**（声明为约定，紧接着用分辨率验证其一致性）：

| 单项式 | 积分值（$\int d\bar x\, dx$ 次序） |
|---|---|
| 1 | 0 |
| $\bar x$ | 0 |
| $x$ | 0 |
| $\bar x \cdot x$ | **-1** |

**分辨率**：

$$\int d\bar x\, dx\; e^{-\bar x x}\, |x\rangle\langle\bar x| = \text{identity} .$$

验证（展开 $e^{-\bar x x} = 1 - \bar x x$，作用到两组基矢）：

- 作用于 $|0\rangle$：被积式 $= e^{-\bar x x}(|0\rangle + x|1\rangle)(\langle0| - \bar x\langle1|)$。
  保留积分值非零的项（须含恰好一个 $\bar x$ 与一个 $x$）：$|0\rangle\langle0|$ 项的系数
  $-\bar x x$ 积分 $= -(-1) = +1$，贡献 $|0\rangle\langle0|$（成立）；
  $|0\rangle\langle1|$ 项：含 $x$ 但需配 $\bar x$——$-x\bar x\langle1|$ 的 $x\bar x = -\bar x x$，
  积分 $= +\bar x x\langle1| = -\langle1|$，贡献 $-|0\rangle\langle1|\cdot(-1) = +|0\rangle\langle1|$（成立）；
  $|1\rangle\langle0|$ 项：$+x|1\rangle\langle0|(1-\bar x x)$：$x$ 项积分为零、
  $-x\bar x\, x$ 中 $x\bar x\, x = 0$（两个 $x$）→ 0（成立）；
  $|1\rangle\langle1|$ 项：$-x|1\rangle\langle1|\bar x(1 - \bar x x)$：$x\bar x = -\bar x x$，
  $+\bar x\, x\cdot\bar x = \bar x\cdot\bar x\cdot x = 0$……直接展开
  $-x\bar x(1 - \bar x x) = -x\bar x = +\bar x x$，积分 $= -1$，贡献
  $-|1\rangle\langle1|\cdot(-1) = +|1\rangle\langle1|$（成立）。
- 作用于 $|1\rangle$：同法逐项（$\langle\bar x|1\rangle = -\bar x$ 进入展开），合计 $|1\rangle\langle1|$（成立）。

两条基矢都还原，分辨率成立。**$\int\bar x\, x = -1$ 的负号正是使分辨率为 +1 的
约定**；改用 $dx\,d\bar x$ 次序则所有符号翻转、公式组同样自洽——绝对符号是
约定，相对结构不是。

### 1.3 迹规则与反周期边界

**迹规则**：

$$\operatorname{Tr}[O] = \int d\bar x\, dx\; e^{-\bar x x}\; \langle-\bar x|\,O\,|x\rangle ,$$

bra 端边界取 $\langle-\bar x|$（**反周期边界**）：把 $\bar x \to -\bar x$ 代入 bra 定义，
$\langle-\bar x| = \langle0|e^{+\bar x d} = \langle0| + \bar x\langle1|$（展开中 (+) 号与 §1.1
同理由）。

**验证**（单模、任意 $\alpha$，含 $\operatorname{Tr}[\text{identity}] = 2$ 的检查）：
$e^{\alpha n}|x\rangle = |0\rangle + e^{\alpha}x|1\rangle$；
$\langle-\bar x|0\rangle = 1$；$\langle-\bar x|1\rangle = +\bar x$
（由（2）代入 $\bar x \to -\bar x$：bra 定义里 (+) 号给出）。
$\langle-\bar x|e^{\alpha n}|x\rangle = 1 + e^{\alpha}x\bar x = 1 - e^{\alpha}\bar x x$。
乘测度因子 $e^{-\bar x x} = 1 - \bar x x$：

$$(1 - \bar x x)(1 - e^{\alpha}\bar x x) = 1 - (1 + e^{\alpha})\bar x x + (\bar x x)^2(= 0) .$$

积分：$0 - (1 + e^{\alpha})\cdot(-1) = 1 + e^{\alpha}$。Fock 直算
$\operatorname{Tr}[e^{\alpha n}] = e^{\alpha\cdot 0} + e^{\alpha\cdot 1} = 1 + e^{\alpha}$
（成立）。特别地 $\alpha = 0$ 给 $\operatorname{Tr}[\text{identity}] = 2$（成立）。

**为什么边界必须是反周期的**：迹对单模求和 $\sum_{n=0,1}\langle n|O|n\rangle$ 在相干态表示
中把 bra 端的 $\bar x$ 与 ket 端的 $x$ 配对；配对规则（§1.2）要求两个变量以
$(\bar x, x)$ 顺序相邻，而 bra 定义自然给出 $\langle+\bar x|$——为落到 §1.2 的积分规则须
把 $\bar x$ 换成 $-\bar x$（Grassmann 变量取负 = 翻转其线性项符号）。这就是"反周期"
的全部内容；多模时逐模反周期。

### 1.4 多模推广

多模（1 杂质 + m bath）时每模一对 $(\bar x, x)$，外部生成元全体相互反对易；
分辨率与迹规则是逐模乘积。**分次乘积注意**：多模算符串重排的费米子符号 =
置换符号；本文涉及多模置换绝对符号处均由 §5 的证书锁定，推导只使用
"每生成元在单项式中至多出现一次"的幂零性。

### 1.5 高斯映射引理（cell 传播子的相干矩阵元）

**引理**：对单粒子矩阵 $h$ 与 $U = e^{-i\tau h}$（分块：杂质行/列 = 1，bath = 其余），
二次型算符 $B(\tau) = e^{-i\tau h}$ 的相干态矩阵元是外部变量的高斯型指数，
其系数恰为 $U$ 的块：

$$\langle\bar x', \bar a'\,|\,B(\tau)\,|\,x, a\rangle \propto
\exp\big[\, \bar x' U_{11} x + \bar x' U_{12} a + \bar a' U_{21} x + \bar a' U_{22} a \,\big] ,$$

∝ 号吸收的常数前因子最终进入 $\operatorname{Pf}(G)$（由数值 Pfaffian 计算）。

**单模（单 bath 模 + 杂质）验证**：取 $h$ 对角（$V = 0$）：
$U = \operatorname{diag}(1, e^{-i\tau\epsilon})$，即 $U_{11} = 1$、$U_{22} = e^{-i\tau\epsilon}$、交叉块 0。
引理给出矩阵元 $\propto \exp[\bar x'\cdot 1\cdot x + \bar a' e^{-i\tau\epsilon} a]$。直接展开：
$e^{-i\tau h} = e^{-i\tau\epsilon c^\dagger c} = 1 + e^{-i\tau\epsilon} c^\dagger c$（生成元幂零截断），矩阵元
$= \langle\bar x', \bar a'|x, a\rangle + e^{-i\tau\epsilon}\langle\bar x', \bar a'|c^\dagger c|x, a\rangle$
$= e^{-\bar x' x - \bar a' a} + e^{-i\tau\epsilon}\cdot(\bar x' x)(\bar a' a)\cdot(\text{费米符号})$。

两条幂零展开在"最多一个 $\bar x'$、一个 $x$、一个 $\bar a'$、一个 $a$"的单项式上
与引理的指数展开一致（指数展开 $= 1 + \bar x'x + e^{-i\tau\epsilon}\bar a'a +$ 高阶(0)），
常数因子 $e^{-\bar x'x - \bar a' a}$ 正是 §1.2 型的边界分辨率——在链的合成中
成为相邻 pair 的测度（§2.2）。含 $V$ 耦合的单模情形由同一幂零展开逐系数
验证；其整体正确性由 §5 的两条机器精度证书锁定。

**本文使用规则**：§2.3 的 Q-kernel 缝合处按引理写入 $U$ 块系数；常数前因子
进 $\operatorname{Pf}$ 通道。**方法论**：结构由推导给出，绝对符号由证书锁定。

### 1.6 高斯积分规则（带源）

**规则**（标准费米高斯积分；n=2 展开验证见 §3.2）：

$$\int \mathcal{D}a\; e^{\frac12 a^{\mathsf T} G a + j^{\mathsf T} a} = \Pf(G)\; e^{-\frac12 j^{\mathsf T} G^{-1} j} ,$$

$a$ = Grassmann 生成元列向量、$G$ = 反对称系数矩阵（$G^{\mathsf T} = -G$）、$j$ = 外部
源列向量、$\operatorname{Pf}(G)$ = $G$ 的 Pfaffian（$\operatorname{Pf}(G)^2 = \det G$）。多模情形即多对
$(\bar x, x)$ 的乘积形式，§2.4 由它推导 Schur 公式。

---

## 2. 第一步：积掉 bath（compact_delta 逐步导览）

本节对照实现代码逐步走一遍"积掉 bath"：输入是模型参数（$H_{\mathrm{bath}}$、$V$、$\beta$）
与调度表（`sched.cells` = $N_c$ 个 cell 时长、`sched.links` = $N_c{+}1$ 个门时长），
输出是接口变量上的 Gaussian 核 $G$ 与归一化。实现 = `NativeExact.compact_delta`。

### 2.0 代码总览：一个大反对称矩阵 A

`compact_delta` 的全部工作就是组装一个 $L \times L$ 的反对称矩阵 $A$，

$$L = \underbrace{4 N_c}_{\text{接口腿}} + \underbrace{2m(2N_c + 1)}_{\text{bath 腿}} ,$$

（$m$ = bath 模式数；$2m$ = 每 pair 的变量数），然后一次 Schur 消元把 bath
块积掉。代码骨架（`NativeExact.jl`）：

```julia
nif = 4 * nl                          # 接口腿数（4 条/接口）
total = nif + 2m * (2nl + 1)          # 接口 + bath 腿
A = zeros(ComplexF64, total, total)
h = _coupled_h(V, Hbath)
Us = [exp(-1im * cell .* h) for cell in cells]   # 逐 cell 单粒子传播子
# …… 逐项填充 A（下文逐条）……
E = A[1:nif, 1:nif]; C = A[1:nif, nif+1:end]; B = A[nif+1:end, nif+1:end]
effective = E + C * (B \ transpose(C))           # Schur：积掉 bath
G = 0.5 .* (effective .- transpose(effective))   # 反对称化
```

### 2.1 变量布局（每个槽位的定义）

**接口腿**（前 $4N_c$ 个槽位）：每个接口 $c = 1..N_c$ 有 4 条腿。代码：

```julia
_ifslot(c, br, stat) = 4*(c-1) + (br === :+ ? 1 : 3) + (stat === :bar ? 0 : 1)
```

即接口 $c$ 的槽位 $4(c{-}1){+}1 \mathrel{..} 4(c{-}1){+}4$ 依次是
$\bar x^{+}(c), x^{+}(c), \bar x^{-}(c), x^{-}(c)$
（$\bar x^{+}(c)$：前向支接口 $c$ 的 creation 腿，其余类推）。这些腿**同时**是
（1）杂质轨迹变量（§4.1 的门链挂在其上）与（2）高斯核 $G$ 的变量——两扇区
记法中它们是两组变量由交叉测度绑定，紧凑记法直接合一（等价性见 §3.6）。

**bath 腿**（后面的 $2m(2N_c{+}1)$ 个槽位）：组织成 pair。每个 pair 是一组
bath 相干场的 $(\bar\gamma, \gamma)$ 对——`pair.bars[k]` = bath 模式 $k$ 的 $\bar\gamma$ 槽位、
`pair.etas[k]` = $\gamma$ 槽位。共 $2N_c + 1$ 对（代码 `pair_data`）：

| pair | 个数 | 角色 |
|---|---|---|
| `forward_pairs[c]`，$c = 1..N_c$ | $N_c$ | forward 链的边界：$B(\tau_c)$ 的输入 etas 与 $B(\tau_{c+1})$ 的输出 bars |
| `backward_pairs[c]` | $N_c$ | backward 链同构 |
| `trace_pair` | 1 | 链末端的共享 trace pair（§2.3） |

三个块 $E$（接口×接口）、$C$（接口×bath）、$\mathcal{B}$（bath×bath）由下面
§2.2–§2.3 的逐项填充给出。

### 2.2 逐 cell 的 Q-kernel 缝合（四个 U 块）

每个 forward cell $c$ 在代码中是一次 `_q_add!` 调用：

```julia
_q_add!(A, _ifslot(c, +, :bar), _ifslot(c, +, :eta),
        forward_pairs[c+1].bars,   # 输出 bath bars
        forward_pairs[c].etas,     # 输入 bath etas
        Us[c])                     # U = exp(-i τ_c h)
```

`_q_add!` 放置四个块（每行一个块；`put_G!(A, i, j, v)` = $A[i,j]\mathrel{+}=v$ 且
$A[j,i]\mathrel{-}=v$，即写入反对称元）：

| 代码行 | 块 | 系数 | 物理含义 |
|---|---|---|---|
| put_G!(bar_imp, eta_imp, U[1,1]) | 接口 pair | $U_{11}(\tau_c)$ | 单粒子**透射幅**：杂质模式穿过 cell $c$（含 bath 虚交换的重整化路径） |
| put_G!(bar_imp, bath_etas[j], U[1,1+j]) | 接口 bar × bath eta | $U_{1,1+j}$ | **入射幅**：$d^\dagger V_\alpha$ 项——bath 粒子被杂质吸收 |
| put_G!(bath_bars[i], eta_imp, U[1+i,1]) | bath bar × 接口 eta | $U_{1+i,1}$ | **出射幅**：$V_\alpha d$ 项——粒子从杂质发射回 bath |
| put_G!(bath_bars[i], bath_etas[j], U[1+i,1+j]) | bath×bath | $U_{22}$ 块 | bath **内部传播**（含 $H_{\mathrm{bath}}$ 的对角与跃迁） |

**这四个块为什么恰好是这些系数**：它们是单粒子传播子 $U = e^{-i\tau_c h}$ 的
四个象限——cell 的二次型演化 $B(\tau_c)$ 作用在（杂质 + bath）单粒子振幅上，
输出振幅 = $U \times$ 输入振幅；相干态表示中这个线性映射写成外部变量的高斯
指数，指数系数就是 $U$ 的四个块（§1.5 引理；单模验证见该节）。

**透射幅的诚实陈述**：$U_{11}(\tau_c) = [e^{-i\tau_c h}]_{11}$ 一般 $\neq 1$（它含
bath 虚交换的重整化路径，单模时 $= \cos(\Omega\tau) - i(\epsilon/2\Omega)\sin(\Omega\tau)$，
$\Omega^2 = \epsilon^2/4 + |V|^2$）；$\tau\to 0$ 极限才趋于 1。旧版文档"接口分辨率恒等"的说法
指的正是这个 $\tau\to 0$ 极限。

### 2.3 backward 链、热边、迹 pair、逐对测度

**backward 链**（代码第二段循环）：后向支是前向的厄米共轭演化——同样的
pair 结构，核系数取 $\operatorname{adjoint}(U)$，接口腿用 $\bar x^{-}/x^{-}$：

```julia
_q_add!(A, _ifslot(c, -, :bar), _ifslot(c, -, :eta),
        backward_pairs[...].bars, backward_pairs[...].etas, adjoint(Us[c]))
```

**热边**（bath 初态）：前向链与后向链在 $t{=}0$ 由 bath 热态闭合（§1.3 迹规则
的多模应用）：

```julia
for i in 1:m, j in 1:m
    put_G!(A, forward_pairs[1].bars[i], backward_pairs[1].etas[j], thermal[i, j])
end   # thermal = exp(-β H_bath)（§1.3 单模验证的多模版）
```

**迹 pair 与 ⟨−ψ̄| 翻转**：最后一个 forward cell 的输出不是新 pair 而是
`trace_pair`（代码 `output_pair = step==nl-1 ? trace_pair : ...`），且该核的
输出行系数**翻号**（代码 `trace_bra` → `_q_add!` 内对输出行取负）。这正是
§1.3 迹规则的**反周期边界** $\langle-\bar\psi|$：bra 端生成元整体取负 = 输出行翻号。
共享 trace pair 的含义：forward 链与 backward 链在 $T$ 端共用同一对边界变量，
把两条链闭合成一个迹。

**逐对测度 $-\bar\gamma\gamma$**（代码 `put_G!(A, pair.bars[mode], pair.etas[mode], -1)`）：
每个 bath pair 上的边界分辨率（§1.2 分辨率的多模版）；绝对符号 −1 由 Z
证书锁定。

### 2.4 Schur 消元：积掉 bath 的那一步

至此 $A$ 是（接口 | bath）分块的反对称矩阵：

$$
A = \begin{pmatrix} E & C \\ -C^{\mathsf T} & B \end{pmatrix},
\qquad
\exp = \tfrac12 \eta_{\mathrm{if}}^{\mathsf T} E \eta_{\mathrm{if}}
+ \eta_{\mathrm{if}}^{\mathsf T} C\, b + \tfrac12 b^{\mathsf T} B b .
$$

对 bath 变量 $b$ 用 §1.6 的带源高斯积分规则（源 = $C^{\mathsf T}\eta_{\mathrm{if}}$，
二次型系数 = $B$）：

$$
\int \mathcal{D}b\; (\cdots) = \Pf(B)\;
e^{\frac12 \eta_{\mathrm{if}}^{\mathsf T}\big(E + C B^{-1} C^{\mathsf T}\big)\eta_{\mathrm{if}}} ,
$$

（推导 = 分块配方：对 $b$ 平移 $b \to b - B^{-1}C^{\mathsf T}\eta_{\mathrm{if}}$
使交叉项消失，常数部分给 $\Pf(B)$。）**这就是"积掉 bath"**：全部 bath
Grassmann 变量被精确积掉，剩余多项式只含接口变量，其高斯核为

$$
G = \operatorname{antisym}\big(E + C B^{-1} C^{\mathsf T}\big)
= \tfrac12\big(\text{effective} - \text{effective}^{\mathsf T}\big) ,
\qquad \text{effective} = E + C B^{-1} C^{\mathsf T} .
$$

代码即 `effective = E + C * (B \ transpose(C))` + 反对称化两行。

**为什么一次 Schur 就把 bath 积干净**：bath 变量在 $A$ 中只通过 $B$ 块相互
耦合、通过 $C$ 块与接口耦合——配平方后 bath 部分成为独立高斯，其积分把
接口块重整化为 $E + CB^{-1}C^{\mathsf T}$。多 cell 的链结构（每个 bath pair 出现在相邻
两个核中）全部体现在 $C$/$B$ 的具体元里，Schur 一次性处理。

**证书**：所得 $G$ 与 IFGTEMPO influence_matrix_G **逐位一致**（order 2
的 4×4、order 4 含负 cell 的 12×12，max|ΔG| = 0）——即本步的簿记（pair
布局、热边、迹翻转、测度符号）与 IFGTEMPO 完全相同。

### 2.5 U(1) 结构、归一化 ν 与本步小结

- **U(1) 同 kind 零**：$H_q$ 数量守恒 ⇒ 上述四个块全是 (bar, eta) 跨类对
  ⇒ $A$ 与 Schur 后的 $G$ 都没有同 kind 元（全 order 数值验证 0）。这条
  性质是 §3.5 行扫表示可行性的关键。
- **归一化 ν**：Schur 抛掉的 $\Pf(B)$ 与逐对测度常数、bath 热态归一化合成

$$
\nu = Z_N/Z = \frac{\operatorname{tr}[\text{bath 链的二次量子化迹}]}{\operatorname{tr}[e^{-\beta H_{\mathrm{bath}}}]} ,
$$

（`bath_norm`；分子 = 逐 cell $U_{22}$ 块乘积的二次量子化迹。）闭合
$Z = (\text{收缩})\cdot\nu\cdot\operatorname{Pf}(G)$；比值型观测量中自动消去。

- **本步小结**：$A$ 装配 + 一次 Schur 后，剩余多项式只含接口腿：

$$
P(x) = (\text{常数})\; e^{\frac12 x^{\mathsf T} G x}\; K[x]\cdot(\text{测度/闭合已并入}) ,
$$

——在两扇区记法中它还带接口交叉测度与独立 IF 变量，其消去是 §3 的内容
（两种记法的等价由 §3.6/§5 证书链锁定）。

---

## 3. 第二步：积掉 IF Grassmann 变量（消元核 D）

### 3.1 本步做什么（两种记法）

第一步得到的 $G$ 定义在**接口变量**上。两扇区记法中接口变量是独立的
IF Grassmann 变量 $a$，与杂质轨迹变量 $x$ 由交叉测度绑定：

$$
P(a,x) = \nu\Pf(\mathcal{B})\;
e^{\frac12 a^{\mathsf T} G a}
\prod_c e^{-\bar a_c\, x^{\eta}_c - \bar x^{\bar}_c\, a_c}
\, K[x]\, M_0 .
$$

逐符号：$c$ 遍历接口；$x^{\eta}_c$ = 接口 $c$ 的 annihilation 腿（前向支）、
$\bar x^{\bar}_c$ = creation 腿（后向支）；$M_0$ = $t{=}0$ 的分辨率测度。交叉测度是
**两扇区构造的接口恒等分辨率**（每个 U_imp–U_bath 接口一个，逐接口独立因子
——它们的非独立性正是 §2 的 Schur 要一次做完的原因）。

**本步**：把 IF 变量 $a$ 精确积掉——用 §1.6 的带源高斯积分规则，源由交叉
测度提供（§3.2）。积掉后 $a$ 消失、$D$ 落在轨迹腿上。紧凑实现把这一步与
§2 融合为一次大 Schur（§3.6 说明等价）。

### 3.2 交叉测度作为源（partner 映射表）

把交叉测度读成 Gaussian 的线性源——逐分量：

| IF 槽位 $s$ | 源 $j_s$（= partner 杂质腿） |
|---|---|
| $\bar a_c^{+}$（IF bar，前向） | $-x^{\eta}(c,+)$（杂质 annihilation 腿） |
| $a_c^{+}$（IF eta，前向） | $-\bar x(c,+)$（杂质 creation 腿） |
| $\bar a_c^{-}$、$a_c^{-}$ | 同构（后向支） |

即**源 = partner 杂质腿**（bar↔eta 交叉；"partner 槽位"映射的定义，代码
`partner_var(t, br, stat) = imp(t, br, 对(stat))`）。

### 3.3 带源高斯积分的 n=2 展开（规则验证）

**规则**：$\int \mathcal{D}a\; e^{\frac12 a^{\mathsf T} G a + j^{\mathsf T} a}
= \Pf(G)\; e^{-\frac12 j^{\mathsf T} G^{-1} j}$。

**n=2 展开**：单接口，$a = (\bar a, a)$，$G = \begin{pmatrix} 0 & g \\ -g & 0 \end{pmatrix}$
（反对称 2×2，$g = G_{\bar a, a}$ 即该接口的高斯幅），源
$j^{\mathsf T} a = j_{\bar a}\,\bar a + j_a\, a$，其中 $j_{\bar a} = -x^{\eta}$、
$j_a = -\bar x$（$x, \bar x, \bar\theta, \theta$ 是杂质腿生成元——**它们与 $a$ 同属一个
Grassmann 代数**，这是路径积分的机制：所有腿都是生成元）。逐项展开
（全部只到二阶；更高阶因生成元重复而为零）：

- $\exp(\tfrac12 a^{\mathsf T}Ga)$：$\tfrac12 a^{\mathsf T}Ga = \tfrac12(g\bar a\,a + a\cdot(-g)\,\bar a) = g\,\bar a\,a$
  （第二项用 $a\bar a = -\bar a\,a$），故 $= 1 + g\,\bar a\,a$；
- $\exp(j^{\mathsf T}a) = 1 + j_{\bar a}\,\bar a + j_a\,a + j_{\bar a}\,j_a\,\bar a\,a\times(\text{符号})$。
  $(j^{\mathsf T}a)^2$ 项：$j_{\bar a}\,\bar a\cdot j_a\,a$ 重排符号分析——两个源的积
  $(-x)(-\bar x)\,\bar a\,a = x\bar x\,\bar a\,a$：把 $\bar x$ 移过 $\bar a$
  （反号）再移过 $a$（反号）……净符号 $(+1)$，得 $+x\bar x\,\bar a\,a$；
- **乘积**：$(1 + g\bar a\,a)(1 + x\bar x\,\bar a\,a)$ 中含 $\bar a\,a$ 的项：
  $g\,\bar a\,a$ 与 $x\bar x\,\bar a\,a$；二者之和 = $(g + x\bar x)\,\bar a\,a$。
  注意 4 生成元项 $(g\bar a\,a)(x\bar x\,\bar a\,a) = 0$（$\bar a$ 重复）。
- **积分**（按逐对规则 $\int\bar a\,a = +1$，符号即 elim\_sign 的约定通道）：
  $\int(g + x\bar x)\,\bar a\,a = g\,(1 + x\bar x/g)$。

**结果**：单接口消元给出 $g\,(1 + x\bar x/g)$——杂质腿对上的核 $g\,\bar x\,x$
（即消元核 $D = g\,\bar x\,x$）加常数 $g$；绝对符号由证书锁定。
**一般 n**：$\exp(j^{\mathsf T}a)$ 展开到二阶（更高阶为零），其二阶项
$\tfrac12(j^{\mathsf T}a)^2$ 与 $\tfrac12 a^{\mathsf T}Ga$ 的交叉把 $\tfrac12 j^{\mathsf T}G^{-1}j$ 吸收进
指数（标准费米高斯积分），得规则式；杂质腿对的核由 $G^{-1}$ 的 partner
读出给出（§3.4）。

**本步是精确恒等式**：实现层面由 Step B 暴力 Berezin 积分逐系数对照锁定
（≤2.2e-16，orders 1/2/4 全验）。

### 3.4 消元核 D 的读出（partner 映射）

一般 $N_c$ 接口、$4N_c$ 个杂质腿（每接口 $\bar x^{+}, x^{+}, \bar x^{-}, x^{-}$ 四条）：

$$
D[\lambda,\lambda'] = -\tfrac12\big(G^{-1}_{ik} - G^{-1}_{ki}\big),\qquad
i = \operatorname{slot}(\bar\lambda),\ k = \operatorname{slot}(\lambda') ,
$$

逐符号：$\lambda,\lambda'$ = 两条杂质轨迹腿；$\operatorname{slot}(\bar\lambda)$ =
$\lambda$ 的 **partner IF 槽位**（源映射：bar↔eta 交叉——
$\lambda = \bar x^{+}(c)$ 的 partner 是 IF 槽位 $(c,+,\eta)$，
$\lambda' = x^{+}(c')$ 的 partner 是 $(c',+,\operatorname{bar})$）；$-\frac12(\cdot-\cdot)$
是把 $G^{-1}$ 反对称化后按 (bar$\to$eta) 有向读出的**反对称化因子**；
绝对符号 $\mathrm{elim\_sign} = -1$ 由 Step B 证书锁定（纯约定）。

### 3.5 同 kind 零与行分解（Q1 判决）

§2.5 的 U(1) 结构经 Schur 与求逆保持 $\Rightarrow$ $G^{-1}$ 同 kind 元为零
$\Rightarrow$ $D$ 同 kind 元为零（全 order 数值验证 0 个）$\Rightarrow$

$$
e^{\sum D} = \prod_{\text{bar } r}\Big(1 + \theta_{\bar r}\sum_c D[r,c]\,\theta_c\Big)
$$

按行分解——每行一个 partialmpo、键维 2，hybriddynamics 行扫完整承载
（Q1 判决）。**消元精确、无近似**；代价是 $D$ 非局域（IF 键维按
$4^{N_c}$ 增长，生产上键维截断）。

### 3.6 融合记法与两步记法的等价（证书链）

compact 实现把 §2（Schur）与 §3（Pfaffian 消元）**融合为一次大 Schur**：
接口槽位直接取为轨迹腿（交叉测度成为 Gaussian 的对角恒等对，其与轨迹腿的
连接由骨架链承担）。两种记法的最终多项式逐位等价——等价性的证书链：
（1）$G$ vs IFGTEMPO 逐位；（2）消元核 vs 暴力积分逐系数；（3）管线
vs Fock 离散链机器级；（4）管线 vs IFGTEMPO 机器级。§3 的两步记法用于
推导与验证，§2 的融合记法用于实现。

---

## 4. 杂质部分与 1/2/4 阶组装

### 4.1 杂质核 K 链

杂质扇区**不走高斯积分**（$H_{\mathrm{imp}}$ 含 $U$ 项，非二次型）——走精确 Fock 门：

$$
K[x] = \langle \bar x_0^{+}|\rho_0|x_0^{-}\rangle
\prod_{t=1}^{N_c+1}\langle \bar x_t^{+}|A(\ell_t)|x_{t-1}^{+}\rangle
\prod_{t=N_c+1}^{1}\langle \bar x_{t-1}^{-}|A(\ell_t)^\dagger|x_t^{-}\rangle,
\qquad A(\ell) = e^{-i\ell H_{\mathrm{imp}}} .
$$

$A(\ell)$ 对多带含 $U$ 时是 $2^{\mathrm{bands}}$ 维 Fock 矩阵的精确指数
——**门内无近似**；近似只在门与 cell 的分割方式（调度）里。

### 4.2 Trotter 调度表

| order | cells（进 $G$） | links（进 $K$） | 全局阶 |
|---|---|---|---|
| 1 | $\operatorname{fill}(\delta t, N)$ | $[\delta t\times N,\,0]$ | $O(\delta t^2)$（机理 §4.4） |
| 2 | $\operatorname{fill}(\delta t, N)$ | $[\delta t/2,\,\delta t\times(N{-}1),\,\delta t/2]$ | $O(\delta t^2)$（对称） |
| 4 | $[a,b,a]\delta t\times N$ | $[a/2,\,c,\,c,\,(a,c,c){\times}(N{-}1),\,a/2]$ | $O(\delta t^4)$ |

links 比 cells 多一个（链两端各一门）。order 1 links 末位 0 的恒等门
不可跳（§4.5）。

### 4.3 回文与阶条件（order 2 与 4 的推导）

记 $X = -i\tau_1 A$、$Y = -i\tau_2 B$、$Z = -i\tau_3 A$。对
$e^{X}e^{Y}e^{Z}$，BCH 展开的指数二阶项：

$$
\tfrac12[X,Y] + \tfrac12[Y,Z] + \tfrac14\big([X,Z] + [Z,X]\big) .
$$

order 2 的链 $X = -i\delta t A/2$、$Y = -i\delta t B$、$Z = -i\delta t A/2$：
$[X,Z] = 0$（同为 $A$ 向），且 $[Y,Z] = [Y,X]\cdot(\tau\text{ 相等})$，于是

$$
\tfrac12[X,Y] + \tfrac12[Y,X] = 0 :
$$

$\delta t^2$ 项恰好相消（"对称/回文"的代数含义：首尾两个 A-半步相等使两个
commutator 贡献互抵）。残留从 $\delta t^3$ 起，单步 $O(\delta t^3)$、
$N$ 步全局 $O(\delta t^2)$。

order 4：回文 Suzuki $S_4 = S_2(a\delta t)\,S_2(b\delta t)\,S_2(a\delta t)$，
合并链后逐宏步

$$
A(a/2)\,B(a)\,A(c)\,B(b)\,A(c)\,B(a)\,A(a/2),
\qquad c = \frac{a+b}{2} .
$$

权重条件：

$$
2a + b = 1\ (\text{一阶：总时长}),\qquad 2a^3 + b^3 = 0\ (\text{三阶相消}),
$$

解出 $a = 1/(2-\sqrt[3]{2}) \approx 1.2252$、$b = 1-2a \approx -1.4505 < 0$。
回文性的角色：组合中 S2 因子序正读倒读一致，所有奇幂 BCH 余项成对相消；
非回文序 $S_2(a)^2 S_2(b)$ 的 $\delta t^3$ 项不满足该条件——实测比值序列
$12.7\to9.0$（滑向 8，即 $p\to3$）是其失效的数值表现；修正为回文序后
比值恢复 16（$N{=}16$ 定点 15.97/16.12）。负 cell 的必然性：
$2a^3 + b^3 = 0$ 要求 $b < 0$——三阶相消必须有反向步，接口格点因此
非均匀且时间倒卷。

### 4.4 order 1 全局二阶的机理

order 1 的裸分裂 $A(\delta t)\,B(\delta t)$ 单步只有 $O(\delta t)$——但严格
核管线的实测全局阶是 p = 2。机理：分裂的领头误差 $[A,B]\delta t^2/2$ 是
"杂质门与耦合演化的交错缺失"；严格 cell $B(\delta t) = e^{-i\delta t\,h}$ 的
耦合顶点幅对 $V$ 的任意幂精确，它把这项交错重新生成到所有阶（$D$ 的二阶
交叉项正是 $[A,B]\delta t^2$ 的内容），残留从更高阶起，全局 $O(\delta t^2)$。
冻结核没有这个重求和，全局 $O(\delta t)$（实测 $p\approx1$）。

### 4.5 组装约定逐条

| 约定 | 内容 | 理由 |
|---|---|---|
| 每 gap 一门 | $\tau = 0$ 恒等门不可跳 | 恒等传播子携带分辨率 $\prod(1+\bar cc)$，是测度的一部分 |
| 接口反测度 | $\mathrm{cc} = -1$ | 幂零代数 $(1-\bar xx)(1+\bar xx) = 1$ 精确抵消 |
| 归一化 | $\nu\operatorname{Pf}(G)$ 进 scaling | 闭合 $Z$；比值型观测量自动消去 |
| 消元符号 | $\mathrm{sgn} = -1$ | Step B 证书锁定 |
| 无 siteflip | $M$ 不预翻 | partialmpo 内部处理位点序 |
| 门符号 | 前向/后向共轭 | Keldysh 两支方向 |

### 4.6 读出约定

- n(T)：折返点自由腿对 GTerm 投影
- 双占据：双折返腿 GTerm 乘积
- gt/lt：c(t_f) @ backward-eta + GTerm 对（有效半步时刻）

### 4.7 多带推广

IF 逐带复制；门用多带 Fock 传播子。

---

## 5. 闭环验证阶梯（每个数字的出处）

| 阶梯 | 判据 | 结果 | 出处 |
|---|---|---|---|
| 恒等式级 | $\epsilon_d=0$: $Z=1$ | $4.6\times10^{-13}$ | zidentity |
| 代数级 | $\int\mathcal{D}a = \nu\Pf\cdot e^{D}$ | $\le 2.2\times10^{-16}$ | Step B |
| 离散链级 | vs 同调度 Fock 直算 | $\le 3.8\times10^{-15}$ | verify_compact |
| 对齐级 | vs IFGTEMPO（U=1.5 双带） | $\le 8.6\times10^{-14}$ | verify_align_u15 |
| 连续极限级 | vs ED, order 1/2/4 | $p = 2/2/4$ | scaling_compact_models |

---

## 6. 边界与延拓

- quapi 路线：仅 orders 1/2
- 连续谱 bath：order 4 非均匀 cells 显式拒绝
- IFGTEMPO：验证参照
