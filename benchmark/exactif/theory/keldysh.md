# 严格 IF 实时核：Keldysh 轮廓上的全流程闭环理论

**benchmark/exactif/theory/keldysh.md**（2026-11-06）
配套实现：`src/NativeExact.jl`；验证入口：`validation/verify_compact.jl`、
`verify_gtlt_align.jl`、`verify_align_u15.jl`、`scaling_compact_models.jl`、`zidentity.jl`。
姊妹文档：[MATHEMATICS.md](MATHEMATICS.md)（形态 B 的消元原理与稠密验证，本文是其
实时 Keldysh 全流程版）；实现对照 IFGTEMPO（对齐参照，非依赖）。

---

## 0. 模型、目标与总览

杂质模型（单带或多带）：

$$
H = H_{\mathrm{imp}}[d] + \underbrace{H_{\mathrm{bath}}[c] + \textstyle\sum_\alpha (V_\alpha d^\dagger c_\alpha + h.c.)}_{\equiv H_q}
$$

bath 处于热平衡 `ρ_b ∝ e^{−βH_bath}`，杂质初态 `ρ₀`（热态，与 bath 可分离）。
目标：Keldysh 轮廓 `0 → T (+)`、`T → 0 (−)` 上的实时观测量——`n(T)`、`ρ(T)`、
双时 Green 函数 `G^>(t,τ) = ⟨c(t)c†(τ)⟩`、`G^<(t,τ) = ⟨c†(τ)c(t)⟩`。

**全流程 = 两次精确高斯积分 + 一次 Trotter 组装**，除 Trotter 分割本身（与
可选的键维截断）外无任何近似：

| 步骤 | 积掉的对象 | 输出 | 实现位置 |
|---|---|---|---|
| ① 积掉 bath | 物理 bath 模式的相干场 | 接口变量上的两扇区 Gaussian `G` + 归一化 `ν·Pf(G)` | `contour_G` / `compact_delta` 的 G 装配 |
| ② 积掉 IF Grassmann 变量 | 接口（influence-functional）变量 | 杂质轨迹腿上的消元核 `D` | `delta_eff` 恒等式（`delta_eff_factors`） |
| ③ 组装 | —— | 杂质核 `K`（order 1/2/4 调度）× `D` × 测度/闭合 → 收缩 → 读出 | `sysdynamics_exact!` + `hybriddynamics` + `gtlt_compact` |

第①②步是**精确恒等式**（数值上分别验证到 Pfaffian/系数级 ≤1e-16）；第③步
携带全部 Trotter 误差，其阶数由调度表决定（1/2 阶全局 O(δt²)、4 阶全局 O(δt⁴)）。

记号：每条 Grassmann 腿一个生成元 θ；`x̄`/`θ̄` 型腿记 creation 类，`x`/`θ` 记
annihilation 类；接口（IF）变量记 `a = (ā, a)`，杂质轨迹变量记 `x = (x̄, x)`。

---

## 1. 第一步：积掉 bath——两扇区 Gaussian 与 Schur 消元

### 1.1 耦合二次型 cell 与单粒子传播子

把 `H_q`（bath + 耦合；**杂质能级移入 H_imp，故 `h` 的杂质对角元为 0**）写成
单粒子形式 `H_q = Σ_{ij} h_{ij} γ_i† γ_j`（mode 0 = 杂质，mode α = bath α）：

$$
h = \begin{pmatrix} 0 & V^\dagger \\ V & H_{\mathrm{bath}} \end{pmatrix},
\qquad B(\tau) = e^{-i\tau h},
$$

每个 Trotter cell 的二次型演化 `B(τ_c)` 是相干场上的**精确高斯映射**：
`exp(½a'ᵀW(τ)a + …)`，其中 `W(τ)` 由 `B(τ)` 的矩阵元给出。任意时长（含
Suzuki 负 cell）都能算——这是 order 4 原生支持的来源。

### 1.2 两扇区 Q-kernel 缝合与 G 装配（compact 布线）

bath 相干场配对成 per-(branch, node) 的 (γ̄, γ) pair；forward 链沿
`(+,1)→(+,2)→…→(+,N_c+1)`、backward 链反向；热边连接 `(+,1)` bars × `(−,1)`
etas（`exp(−βH_bath)`）；远端共享 trace pair（末个 forward 核输出行带
`⟨−ψ̄|` 翻转）；每 pair 测度 `−γ̄γ`。每个 forward cell 把 `B(τ_c)` 的矩阵元
缝合进大反对称矩阵 `A`（`put_G`：上三角 `+= v`、下三角 `−= v`）：

$$
A[\bar x_c,\, x_c] \mathrel{+}= U[1,1],\quad
A[\bar x_c,\, a'_{\alpha}] \mathrel{+}= U[1,1+\alpha],\quad
A[\bar a'_{\alpha},\, x_c] \mathrel{+}= U[1+\alpha,1],\quad
A[\bar a'_{\alpha},\, a_\beta] \mathrel{+}= U[1+\alpha,1+\beta]
$$

（`U[1,1] = 1`：杂质能级在门里，`h` 的杂质对角元为 0。）

### 1.3 Schur 消元与 G 的性质

`A` 分块为（接口 | bath）：`A = [[E, C], [−Cᵀ, B]]`。**积掉 bath** = 对 bath
块 `B` 作 Schur 消元：

$$
G \;=\; \operatorname{antisym}\!\big(E + C\,B^{-1}C^{\mathsf T}\big),
\qquad \operatorname{antisym}(M) = \tfrac12(M - M^{\mathsf T}),
$$

得到**接口变量上的两扇区 Gaussian** `G`（compact 布线 = IFGTEMPO
`_analytic_contour_G`，order 2 与含负 cell 的 order 4 已逐位验证一致）。

性质（全部数值验证）：
- **U(1) 结构**：`H_q` 数量守恒 ⇒ `G` 的非零元只在 (bar, eta) 跨类槽位对上，
  同 kind（bar,bar）/(eta,eta) 块恒为零 ⇒ `G⁻¹` 同样 ⇒ 消元核 `D` 无同 kind
  项（全 order 数值验证 0.000e+00）；
- 归一化 `ν = Z_N/Z`：bath 链传播的迹比（`bath_norm`）；闭合配分函数
  `Z = (\text{收缩}) · ν · Pf(G)`；
- `Pf(G)` 用 Parlett–Reid O(n³) 数值 Pfaffian（勿用 2ⁿ 展开）。

---

## 2. 第二步：积掉 IF Grassmann 变量——Pfaffian 恒等式与消元核 D

### 2.1 接口交叉测度作为线性源

积掉 bath 后的两扇区多项式（`a` = 接口 IF 变量，`x` = 杂质轨迹变量）：

$$
P(a, x) \;=\; \nu\, e^{\frac12 a^{\mathsf T} G a}\;
\prod_{t} e^{-\bar a_t x_t - \bar x_t a_t}\; K[x]\, M_0 ,
$$

接口交叉测度 `exp(−ā_t x_t − x̄_t a_t)`（每接口每 branch 独立因子）正是
Gaussian 的**线性源**：slot `(t, br, bar)` 上的源为 `−x(t, br, eta)`，slot
`(t, br, eta)` 上的源为 `−x̄(t, br, bar)`——即**源 = partner 杂质腿**
（bar↔eta 交叉）。

### 2.2 Pfaffian 恒等式与消元核

Grassmann Gaussian 恒等式（n=2 手算核验）：

$$
\int \mathcal{D}a\; e^{\frac12 a^{\mathsf T} G a + j^{\mathsf T} a}
\;=\; \Pf(G)\; e^{-\frac12 j^{\mathsf T} G^{-1} j} .
$$

代入源 `j`（源变量成对出现，`j_i j_k` 给出杂质腿对的单项式），积掉全部 IF
变量后：

$$
\int \mathcal{D}a\; P(a,x)
\;=\; \nu\,\Pf(G)\; e^{\sum_{\lambda<\lambda'} D[\lambda,\lambda']\,
\theta_\lambda \theta_{\lambda'}}\; K[x]\, M_0 ,
\qquad
D[\lambda,\lambda'] = -\tfrac12\big(G^{-1}_{ik} - G^{-1}_{ki}\big),
$$

其中 `i, k` 是 `λ, λ'` 的 partner IF 槽位（即 `D = −Ā` 在 partner 映射下的
读出，`Ā = (G⁻¹ − G⁻¹ᵀ)/2`）。**这一步是精确恒等式**：暴力 Berezin 积分
与解析式逐系数一致（≤2.2e-16，orders 1/2/4 全验）。

### 2.3 D 的结构：bar→eta 有向块

`G⁻¹`（同 `G`）只有 (bar, eta) 跨类元 ⇒ `D` 只含 **bar→eta 有向对**，同 kind
单项式（bar–bar / eta–eta）恒为零（全 order 验证 0 个）。因此

$$
e^{\sum D} \;=\; \prod_{\bar r} \Big(1 + \theta_{\bar r}
\textstyle\sum_c D[r,c]\,\theta_c\Big)
$$

按行分解（每行一个 partialmpo，键维 2）——**partialif 行扫完整承载消元核**
（Q1 判决：全 order 0 个同 kind 单项式）。这是 GTEMPO `hybriddynamics` 能直接
消费 D 的机制基础。IF 变量不再出现；消元不损失任何信息，代价是 `D` 的非局域性
（IF 键维按 4^{N_c} 增长，生产上用键维截断控制）。

---

## 3. 杂质部分与 1/2/4 阶完全组装

### 3.1 杂质扇区：K 因子

杂质核是**精确 Fock 传播子**的链（任意 `ε_d`、`U`，逐门精确）：

$$
K[x] \;=\; \langle \bar x_0^+|\,\rho_0\,|x_0^-\rangle
\prod_{t=1}^{N_c+1} \langle \bar x_t^+|\,A(\ell_t)\,|x_{t-1}^+\rangle
\prod_{t=N_c+1}^{1} \langle \bar x_{t-1}^-|\,A(\ell_t)^\dagger\,|x_t^-\rangle,
\qquad A(\ell) = e^{-i\ell H_{\mathrm{imp}}},
$$

`N_c` = cell 数（接口数），`ℓ_t` = 门时长（`N_c + 1` 个）。前向链止于末层
creation 腿、后向链止于末层 annihilation 腿——开放端读出挂在这两条腿上。

### 3.2 Trotter 调度表（order 1/2/4）

`A`/`B` 分割 `e^{-i(A+B)\delta t}`：`A`-因子（杂质门，进 K）与 `B`-cell（bath
+ 耦合，进 G）交错合并：

| order | cells（进 G） | links（进 K） | 全局阶 |
|---|---|---|---|
| 1 | `fill(δt, N)` | `[δt,…,δt, 0]`（末位 0：开放端在末时刻） | **O(δt²)**：严格核精确顶点幅已重求和 `[A,B]δt²` 首项 |
| 2 | `fill(δt, N)` | `[δt/2, δt,…,δt, δt/2]`（Strang，对称） | **O(δt²)** |
| 4 | `[a, b, a]δt × N` | `[a/2, c, c, (a,c,c)×(N−1), a/2]`，`c = (a+b)/2` | **O(δt⁴)** |

（order 1 与 order 2 对角 `H_imp` + 对角初态下逐位同链——已知恒等式。）

### 3.3 order 4：回文 Suzuki 的必要性

标准 Suzuki 四阶 = **回文组合** `S4 = S2(aδt)·S2(bδt)·S2(aδt)`：

$$
2a + b = 1,\qquad 2a^3 + b^3 = 0
\;\Longrightarrow\;
a = \frac{1}{2-\sqrt[3]{2}} \approx 1.2252,\quad
b = 1 - 2a \approx -1.4505 < 0,\quad
c = \frac{a+b}{2} < 0 .
$$

三个条件缺一不可：
1. `2a + b = 1`：一阶正确性；
2. **回文**（时间反演对称）：组合中 S2 因子序正读倒读一致——奇幂（O(δt³)）
   BCH 余项恒零的充要条件；**非回文序 `S2(a)²S2(b)` 携带 O(δt³)**，渐近比值
   滑向 8（实测 12.7→9.0 序列）——曾实现的错误序；
3. `2a³ + b³ = 0`：三阶相消，与回文序联立解出上权重。

负 cell（`b < 0`）是三阶相消的必然代价：接口格点非均匀且时间倒卷（macro 内
边界 `t + a·δt → t + (a+b)·δt`）。回文链（合并相邻 A 后）

$$
A(a/2)\,B(a)\,A(c)\,B(b)\,A(c)\,B(a)\,A(a/2)
$$

逐步构造（cell 边界逐一 Schur）即得 order-4 的 `G` 与 `K`。实测：比值序列
14.4→15.6→15.9→**16.0**（N=16 定点 15.97/16.12，渐近拟合 p = 3.98/3.99）。

### 3.4 原生格上的完全组装与锁定约定

GTEMPO 原生实时格（`GrassmannLattice(N = N_c+1, δt, :real, bands)`：node 1 =
t=0、node 2..N_c+1 = 接口、node k = 开放端）上的三个因子：

1. **IF 因子** `mpsIs`：`compact_delta` 的 `M` 矩阵（D 的 bar→eta 有向读）经
   `hybriddynamics` 行扫；
2. **骨架** `mpsK`：按调度的杂质门（**每 gap 一门，τ=0 恒等门不可跳**——
   恒等传播子携带相干态分辨率 `∏(1+c̄c)`，是测度的一部分）+ 接口反测度
   `(1 − x̄x)`（补偿隐式 pair 分辨率的重复施加；**锁定 `cc = −1`**）+
   Keldysh 折返 `boundarycondition!`（`exp(−x̄⁻x⁺)`）+ `systhermalstate!`
   （`ρ₀ = e^{−βH_imp}/tr`）；
3. **归一化**：`normdelta = ν·Pf(G)` 进 scaling；`Z = integrate·norm`。

消元核与骨架的约定**成对锁定**（Z=1 与 Fock 双证书标定）：
`counter_c = −1`、`sgn = −1`、无 siteflip、无 node-1 修正。

### 3.5 读出约定（闭环的出口）

- **n(T)**：折返点自由腿对 `(a⁺_k, ā⁻_k)` 的 2 腿 GTerm 投影（与稠密 ρ(T)
  投影同约定；`occupation` 的 insert_n 带 O(δt) 有效时刻偏移，勿用）；
- **⟨n₁n₂⟩**：band-1 + band-2 两条折返腿 GTerm 的乘积投影；
- **gt/lt**：`c(t_f)` @ 开放端 backward-eta 腿；
  `greater[k+1] = ⟨c(t_f)c†(k·δt)⟩` 用 GTerm(fwd_bar@层 k+1, csite)，
  `lesser[k+1] = ⟨c†(k·δt)c(t_f)⟩` 用 GTerm(bwd_bar@层 k, csite)；
  裸关联、Z 归一；**有效半步时刻** `G^> ≈ (k+½)δt`、`G^< ≈ (k−½)δt`；
- **ρ(T)**：开放腿提取（FormB 路线）或折返跨越投影。

---

## 4. 闭环验证阶梯与数值汇总

| 阶梯 | 判据 | 结果 |
|---|---|---|
| 恒等式级 | ε_d = 0：`Z = 1`（任意 order/N；测度/簿记的最锐利证书） | 最差 **4.6e-13** |
| 代数级 | 暴力 Berezin ∫Da == ν·Pf·exp(D)（系数级） | ≤2.2e-16 |
| 离散链级 | n / ⟨n₁n₂⟩ / ρ / gt/lt vs 同调度 Fock 直算 | ≤**3.8e-15 / 3.8e-14** |
| 对齐级 | vs IFGTEMPO 直接同模型（U=1.5 双带，热态 4×4） | n/⟨n₁n₂⟩ ≤**2.2e-14**、gt/lt ≤**8.6e-14** |
| 连续极限级 | vs ED（有效时刻），order 1/2/4 | **p = 2 / 2 / 4**；比值 →4/4/16（N=16：15.97/16.12） |

各 order 一览（Toulouse/单模 bath 与 Anderson U=1.5/双带，n(T) 标度）：
Toulouse p = 2.01/2.01/3.94（渐近段 2.00/2.00/3.99）、
Anderson p = 1.98/1.98/3.94（渐近段 1.997/1.997/3.990）。

两点已知的非误差项（预期行为，勿误判）：
- order 4 粗端（δt=0.4）误差仅略优于 order 1——Suzuki 负 cell 的大误差常数
  所致的预渐近区，斜率 4 在细端兑现；
- quapi 路线 order 1 N=1 的 n₀ 偏差——调度表"末位补零"约定在 N=1 时清掉唯一
  cell 的边缘情况（该路线 order 4 本就不支持，见下）。

## 5. 边界与延拓

- **quapi 路线**（`exactif_native_quapi`，点基 Δ）：仅 orders 1/2——Suzuki 的
  杂质门插在 bath cell 之间，点基 Δ 的点对点直连假设不成立；
- **连续谱 bath**：原生管线接受有限模 `(H_bath, V)`（δ 谱 = 精确）；真连续谱
  需先谱离散化（求积节点，误差随模式数收敛）——IFGTEMPO 的连续谱路线
  （Volterra + 谱求积，无模式化）在均匀 δt 网格上工作，order 4 的非均匀
  cells 目前显式拒绝；
- **IFGTEMPO 的角色**：全部严格性验证以它为独立参照（n/⟨n₁n₂⟩/ρ/gt/lt 的
  机器级对齐）；`src/` 与活跃验证链本身 IFGTEMPO-free。
