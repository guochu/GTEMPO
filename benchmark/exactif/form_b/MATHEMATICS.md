# 形态 B（消元式单扇区）的数学原理

本文档给出 `form_b/formb_dense.jl`（稠密系数级验证）与 `src/FormB.jl`
（GTEMPO GrassmannMPS 移植）所实现表示的完整数学推导：从 IFGTEMPO 的
两扇区 Grassmann 路径积分出发，把 bath 影响 IF 的生成元**精确**积掉，
得到只剩一套杂质轨迹变量的等价表示，并说明它与 GTEMPO pair 收缩机制
的逐条对接。所有公式均被代码验证（验证入口见文末表格）。

约定：单 flavor 记号；`±` 为 Keldysh 前向/后向支；bar 为产生侧标签
（不是数值共轭）；`θ_λ` 按 lattice site 递增序排 monomial。

---

## 1. 设定：离散 Trotter 链与两扇区变量

### 1.1 算符链

总 Hamiltonian 分解为杂质踢与二次胞：

$$H = H_{\rm imp} + H_Q,\qquad A(\ell)=e^{-i\ell H_{\rm imp}},\qquad Q(c)=e^{-icH_Q},$$

其中 $H_Q$ 的杂质对角为零（杂质能级归踢算符携带），$H_{\rm imp}$ 含
$\epsilon$ 与 $U$。Trotter 调度（order 1/2/4，`trotter_schedule`）给出
胞长 $c_1..c_{n_c}$ 与链节长 $\ell_1..\ell_{n_c+1}$（order 2 的 Strang
合并半步、order 4 的 Suzuki $S_2(a\delta t)^2S_2(b\delta t)$、$2a+b=1$，
合并内层 $\ell=(a+b)/2$ 可为负）。离散链为（右端先作用）：

$$U \;=\; A(\ell_{n_c+1})\,Q(c_{n_c})\,A(\ell_{n_c})\cdots Q(c_1)\,A(\ell_1),$$

$$\rho^{\rm disc}_{\rm imp}=\operatorname{Tr}_B\!\big[U\,(\rho_0\otimes\rho_B)\,U^\dagger\big].$$

实时 Keldysh 围道：前向支携带 $U$、后向支携带 $U^\dagger$。
**迹保持**：所有因子幺正 ⟹ $Z=\operatorname{tr}\rho^{\rm disc}_{\rm imp} =\operatorname{tr}\rho_0=1$——这是区分结构正确性的最锐利证书
（对任意参数、任意 order 成立）。

### 1.2 变量与 lattice

每个 $U_{\rm imp}$–$U_{\rm bath}$ 接口插入一个 Grassmann 恒等 resolution
$1=\int d\bar x\,dx\,e^{-\bar xx}\,|x\rangle\langle\bar x|$。接口 $t$ 的
生成元分两扇区（每 flavor 8 个 site，`timelocallattice` 序）：

$$\big(\bar a_t^{+},\,x_t^{+},\,\bar x_t^{+},\,a_t^{+},\;\;
     \bar a_t^{-},\,x_t^{-},\,\bar x_t^{-},\,a_t^{-}\big)$$

- $a,\bar a$：**IF 变量**（bath 响应），G 矩阵槽位序为
  $\mathrm{slot}\,0,1,2,3=(\bar a_t^{+},a_t^{+},\bar a_t^{-},a_t^{-})$，即
  $\mathrm{flat}=4(t{-}1)+\mathrm{slot}{+}1$（1-based）；
- $x,\bar x$：**IMP 轨迹变量**（ket/bra 腿），$x_t^{\pm}$ 同时是
  $\ell_{t+1}$ 的输入 ket、$\bar x_t^{\pm}$ 是 $\ell_t$ 的输出 bra；
- layer 0：$(\bar x_0^{+},x_0^{+},\bar x_0^{-},x_0^{-})$，初态核
  $\langle\bar x_0^{+}|\rho_0|x_0^{-}\rangle$；
- layer $N{+}1$（开放腿）：$\bar x_{N+1}^{+},\,x_{N+1}^{-}$。

相邻 site 恰为一个恒等 resolution 的 (bar, eta) 对——这一布局使每接口
的交叉测度天然落在 pair 上。

---

## 2. 两扇区多项式

$$P(a,x)\;=\;\nu\,e^{\frac12 a^{T}Ga}\;
   \prod_{t=1}^{n_c} e^{-\bar a_t^{+}x_t^{+}-\bar x_t^{+}a_t^{+}
                            -\bar a_t^{-}x_t^{-}-\bar x_t^{-}a_t^{-}}\;
   e^{-\bar x_0^{+}x_0^{+}-\bar x_0^{-}x_0^{-}}\;K[x]$$

各因子：

1. **IF Gaussian** $e^{\frac12 a^TGa}=e^{\sum_{i<j}G_{ij}a_ia_j}$：
   $G$ 为反对称 IF 系数矩阵（`influence_matrix_G`）——bath 围道大矩阵
   $A=\begin{pmatrix}E&C\\-C^T&B\end{pmatrix}$（bath 生成元 + 热核 +
   逐 pair 恒等测度 + 反周期 trace 边界）对 bath 块的 Schur 补
   $G=\mathrm{antisym}(E+CB^{-1}C^{T})$；其元素即论文六类相关函数族。
   $\nu=Z_N/Z$：传播后热 bath 迹对未传播迹之比。
2. **接口交叉测度** $M_t=e^{-\bar a_t x_t-\bar x_t a_t}$：每个恒等
   resolution 的 ket 变量归左核、bra 变量归右核，因此测度把 IF 与 IMP
   扇区**交叉**配对——这是"两扇区"的来源（IFGTEMPO 推导文档 §3、9）。
3. **layer 0 测度** $M_0=e^{-\bar x_0^{+}x_0^{+}-\bar x_0^{-}x_0^{-}}$。
4. **$K[x]$**：踢核 $\langle\bar x_t^{\pm}|A(\ell_t)|x_{t-1}^{\pm}\rangle$
   （前向 $e^{-i\ell\epsilon}$ 对角、后向 $A^\dagger$）+ 初态核。核的
   单项式符号配方：$\mathrm{coeff}(\mathrm{bra},\mathrm{ket}) =(-1)^{C(\mathrm{ket},2)}K[\mathrm{bra}{+}1,\mathrm{ket}{+}1]$
   （$C(n,2)=n(n{-}1)/2$，双占据符号），local 序 (bar…, eta…)。

物理量读出：

- **开放核** $\langle\bar x_{N+1}^{+}|\rho(T)|x_{N+1}^{-}\rangle$：
  积掉全部内部变量、保留开放腿，系数按上述配方的逆映射回矩阵；
- **迹** $Z$：反周期边界
  $\operatorname{Tr}\rho=\int d\bar x\,dx\,e^{-\bar xx}\langle-\bar x|\rho|x\rangle$
  ——把 $\langle\bar x|\to\langle-\bar x|$ 的翻转并入闭合符号后，即
  闭合因子取 $(1+\bar x_{N+1}^{+}x_{N+1}^{-})$（系数 **+1**；取 −1 是
  周期迹，$Z\ne1$）。

---

## 3. 消元：IF 变量的精确高斯积分（核心）

### 3.1 源嵌入

把全部含 $a$ 的指数项集中，测度成为 IF 生成元上的**线性源**：

$$\frac12 a^{T}Ga+\sum_t\big(-\bar a_t x_t-\bar x_t a_t\big)
 =\frac12 a^{T}Ga+j^{T}a,\qquad
 j_i=-\theta_{\lambda(i)},$$

即 IF bar slot（$\bar a_t^{\pm}$）的 partner 是 IMP eta 腿
$x_t^{\pm}$、IF eta slot（$a_t^{\pm}$）的 partner 是 IMP bar 腿
$\bar x_t^{\pm}$（bar↔eta 交叉，§1.2 的配对）。

### 3.2 高斯恒等式

$$\boxed{\;\int Da\,\exp\!\Big(\tfrac12 a^{T}Ga+j^{T}a\Big)
 =\operatorname{Pf}(G)\,\exp\!\Big(-\tfrac12 j^{T}G^{-1}j\Big)\;}$$

配方：$-\tfrac12\xi A\xi+\eta\xi=-\tfrac12(\xi-\eta A^{-1})A(\xi-\eta A^{-1}) +\tfrac12\eta A^{-1}\eta$，Grassmann 积分的平移不变性把平方项积出
$\operatorname{Pf}$。$n=2$ 直接核验：
$\int e^{A_{12}a_1a_2+j_1a_1+j_2a_2}=A_{12}+j_1j_2$，与公式一致。

代入 $j_i=-\theta_{\lambda(i)}$，展开

$$-\tfrac12 j^{T}G^{-1}j
 =-\tfrac12\sum_{i,k}G^{-1}[i,k]\,\theta_{\lambda(i)}\theta_{\lambda(k)},$$

取递增 monomial $\theta_\lambda\theta_{\lambda'}$（$\lambda<\lambda'$，
$i=\mathrm{slot}(\lambda)$、$k=\mathrm{slot}(\lambda')$）的系数：
$(i,k)$ 项贡献 $-\tfrac12G^{-1}[i,k]$，$(k,i)$ 项经
$\theta_{\lambda'}\theta_\lambda=-\theta_\lambda\theta_{\lambda'}$
贡献 $+\tfrac12G^{-1}[k,i]$，合计

$$\boxed{\;D[\lambda,\lambda']
 =-\tfrac12\big(G^{-1}[i,k]-G^{-1}[k,i]\big)\;}\qquad(\lambda<\lambda')$$

——只依赖 $G^{-1}$ 的**反对称部分**（对称部分因 $\theta$ 反对易自动
消失）。这就是 README `Δ_eff = −½PᵀG⁻¹P` 的反对称多项式形式。

### 3.3 消元多项式与性质

$$P_{\rm elim}(x)\;=\;\nu\operatorname{Pf}(G)\;
   e^{\sum_{\lambda<\lambda'}D[\lambda,\lambda']\,\theta_\lambda\theta_{\lambda'}}\;
   K[x]\;M_0$$

- **精确恒等式**：$\int Da\,P=P_{\rm elim}$，无近似、无信息损失；
  全部接口测度被吸收——消元表示中**内部接口不得再出现任何测度**。
- $D$ 只作用在接口轨迹腿上（layer 0 / 开放腿无源，对应行全零）；
- $D$ 稠密、非局域（$G^{-1}$）——与 QuAPI 冻结核 $\Delta t$ 的
  稀疏时序结构不同；这是同一离散链的另一种多项式实现，二者不可逐
  元素互换。骨架测度 $-\bar ax$ 的吸收表现为 $D$ 在同接口相邻
  slot 对上的"测度块"（$\approx-\tfrac12$）；
- **归一化** $\nu\operatorname{Pf}(G)$ 替代精确核的 prefactor。

### 3.4 收缩语义：投影积分

稠密表示下的全收缩 = 逐生成器提系数（$\int d\theta\,\theta f\to f$，
顺序任意）。因为每个 monomial 中每个生成元至多出现一次，这一投影
**与次序无关、无定向符号**——测度权重已全部显式在多项式因子里
（测度因子、闭合因子），与 GTEMPO pair sweep 的"逐 pair 带权规则"
$f_0-f_{be}$ 是同一运算的两种分布方式（见 §5）。

### 3.5 验证（`formb_dense.jl`）

| 步骤 | 内容 | 结果 |
|---|---|---|
| A | 稠密两扇区副本：$Z=1$ 且 $\rho(T)$ = Fock 离散链 | $2.6\times10^{-18}$ / $1.1\times10^{-16}$ |
| B | 暴力 $\int Da\,[e^{\frac12aGa}\times\text{测度}]$ 逐系数对照 $\nu\operatorname{Pf}\cdot e^{D\theta\theta}$ | $\le2.2\times10^{-16}$（含 order 4） |
| C | 消元式单扇区收缩 vs Fock（orders 1/2 × N≤4、order 4 × N=1） | $\le3.3\times10^{-16}$ |
| D | 消元式两自旋 vs Fock vs IFGTEMPO GMPS | 三方 $\le2.2\times10^{-14}$ |

---

## 4. 因子分解的精确性（工程规则，两处踩坑来源）

1. **逐恒等拆分**：$M_t=e^{-\bar ax-\bar xa}$ 是**四个独立恒等**的指数
   之积。不同恒等的双线性必须作为**独立因子**相乘——合并成一张单项式
   表会使同一因子内的单项式互相组合丢失（指数 $e^{\sum}$ 的组合结构
   被破坏），多项式直接塌为零。
2. **行因子化**：$e^{\sum_{i<j}G_{ij}a_ia_j} =\prod_i\Big(1+\theta_i\sum_{j>i}G_{ij}\theta_j\Big)$（行首唯一 +
   幂零性 ⟹ 精确，`get_if_gmps` 逐行 partialmpo 同款）。行的乘积中
   每个 perfect matching 恰出现一次、带正确符号——即
   $\operatorname{Pf}$ 的展开。Δ_eff 因子同法。
3. **数值 Pfaffian 用 Parlett–Reid**（$O(n^3)$，斜对称 Gauss 消元到
   三对角，$\operatorname{Pf}=\prod$ pivots）：pair 因子稠密展开既丢
   monomial 合并的排列符号（order 4 偏 $8.3\times10^{-6}$，被显式
   perfect matching 求和判定），又需 $2^{4n_c}$ 内存（大 N 即爆）。

---

## 5. 与 GTEMPO pair 机制的数学对接（GMPS 移植）

### 5.1 GTEMPO 的逐 pair 规则

`GrassmannTransferMatrix` 按**相邻 site pair** 扫描，`g_trace_phy` 对
每对施加**带权 graded 规则** $f_0-f_{be}$（即 $\int d\bar x\,dx\, e^{-\bar xx}f$ 的分辨率加权积分 + 内部 distributed graded 符号）。
两扇区机器的 pair 布局恰好逐对落在恒等 resolution 上，故该隐式权重
就是物理测度——`build_impurity_path` 不需要显式 M 因子。
（注意区分：GTEMPO **自身单扇区管线**的测度在 `sysdynamics` 里以
overlap-removal GTerms $\prod(1+\bar cc)$ 施加——将来 FormB 库内集成
时的补偿应接在该测度语义上；本 bench 移植走 IFGTEMPO 收缩路径。）

### 5.2 消元表示的对接

紧凑 IMP-only lattice（`timelocallattice` 滤掉 :IF，相邻 pair 结构
保持）上的 pair 语义逐对翻译：

| pair | GTEMPO 隐式规则 | 消元表示需要 | 处理 |
|---|---|---|---|
| layer 0 $(\bar x_0^{\pm},x_0^{\pm})$ | resolution | $M_0$ | 恰好相等，不放因子 |
| 开放 $(\bar x_{N+1}^{+},x_{N+1}^{-})$ | resolution | 反周期闭合 $(1+\bar x^{+}x^{-})$ | 恰好相等（+1 形式），不放因子 |
| 接口 $(x_t^{\pm},\bar x_t^{\pm})$ | resolution | **无测度**（已吸收进 $D$） | 乘幂零反测度 $(1-x_t\bar x_t)$（系数 **−1**），$(1-\bar xx)(1+x\bar x)=1$ 精确抵消 |

反测度折进 Δ_eff 因子所在的同一 GMPS。**因子排列**取
$(\Delta,\text{bwd},\rho_0{+}\text{前向})$（三个非平凡 GMPS，喂
GTEMPO 任意因子数的 pair transfer）；$Z$ = `l_LL × transfer × scalar`，
$\rho(T)$ 复用 `gtempo_open_kernel_direct`/`density_matrix_direct`
（单 flavor 自写 2 开放腿版——gtempo 版硬编码 4 条开放腿）。

### 5.3 graded 读出符号

多项式乘积与因子分组无关，但 pair sweep 的 graded 符号**依赖因子
排列**——消元格与 timelocallattice 的扇区结构不同，开放系数提取的
某些扇区差一个费米子符号。按 $Z=1$ / $\rho$=Fock 双证书标定（在全部
order/N/模型上一致，是固定读出约定而非拟合）：

- 单 flavor：直接映射；
- 两 flavor：乘 $\mathrm{diag}(1,1,1,-1)$——仅 $|\uparrow\downarrow\rangle$
  对角扇区差一个符号（与两自旋提取中"中间交换"符号
  $(-1)^{\mathrm{bit}_1\cdot\mathrm{bit}_2}$ 同性质）。
- 排列与读出符号**联合标定**：换因子排列必须重标
  （单 flavor 在 $(\Delta,\rho_0{+}\text{fwd},\text{bwd})$ 排列下需
  $(-1)^{\mathrm{bra}\cdot\mathrm{ket}}$，在当前排列下恒等）。

### 5.4 归一化与核

- $\nu\operatorname{Pf}(G)$ 折进 GMPS scaling（GTEMPO 按
  $\mathrm{scaling}^{n_{\rm sites}}$ 计）；
- 单 flavor 对角核（踢/对角初态）是纯双线性
  $(1-\,e^{\mp i\ell\epsilon}\,\theta_\eta\theta_\beta)$，
  partialmpo 直构（含 local→global 排列符号）；两 flavor 用
  `impurity_kernel_gmps` 的 Z2 分级构造（电荷零扇区；**奇宇称核单项式
  不在其表示内** ⟹ 初态密度须对角——IFGTEMPO 全部测试同类）。

---

## 6. 精确性边界与规模

- **精确档**（notrunc）：$D$ 非局域 ⟹ bond 按 $4^N$ 增长——逐位验证
  只到小 N；
- **截断档**（`truncdimcutoff(D=512, ϵ=1e-12)`，同生产用法）：截断误差
  $\sim10^{-12}$ 量级，远低于最小 Trotter 误差，标度不受影响；
- order 4 的 lattice 层数 $n_c=3N$（每宏步 3 个 subcell），凡按
  $4n_c$ 计槽位的函数必须传 $n_c$ 而非宏观 $N$（否则静默只取前 4 槽
  ——症状：$\rho$ 形状全对但整体差一个归一化因子；判据：$Z\ne1$）；
- order 4 的固定-$T$ 平台行为**依赖初态**：$|0\rangle\langle0|$ 初态
  相邻比值 $\sim10$–$12$（preasymptotic），占据对角初态恢复 $\sim16$
  ——与已发表基准对比标度时必须同 $T$ 同初态。

---

## 7. 验证入口

| 声明 | 验证位置 |
|---|---|
| §1 算符链 / trace 保持 | `formb_dense.jl: fock_chain_*`，电池 Z=1 |
| §2 因子结构 | `formb_dense.jl: Step A`（两扇区副本 == Fock） |
| §3.2 高斯恒等式 | $n{=}2$ 直算 + `Step B`（暴力 vs 解析，≤2.2e-16） |
| §3.3 消元收缩 | `Step C`（≤3.3e-16）、`Step D`（三方 ≤2.2e-14） |
| §3.4/迹闭合 | `Step A` 旋钮网格（closure_k=+1 唯一存活） |
| §5.2 反测度 | `src/verify_formb.jl: calibrate`（−1 唯一存活，Z≈5e-15） |
| §5.3 读出符号 | `calibrate` + 全电池（各 order/N/模型一致） |
| §4.3 Pfaffian | `pf_matchings`（完美匹配显式）+ `pf_mine`（行因子 Berezin）三方 ≤2.2e-16 |
| Trotter 阶数 | `src/scaling_formb.jl` → `results/`（p = 2.00/2.00/3.99/4.16） |
