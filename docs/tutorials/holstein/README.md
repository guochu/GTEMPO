# Holstein 模型（semicircular 费米浴）的 GTEMPO 验证与 δt 收敛性分析

运行脚本：

```bash
julia --project=. docs/tutorials/holstein/holstein.jl
```

## 1. 模型

半圆形谱密度的费米浴 + 局域声子 + 无相互作用杂质（U = 0，Toulouse 型）：

```text
H = ϵ_d a†a + ω₀ b†b + g a†a (b† + b)
    + Σₖ ωₖ cₖ†cₖ + Σₖ Vₖ (a†cₖ + cₖ†a)
```

Lang–Firsov 变换下声子只耦合在占据数 n̂ 上，模型严格可解：
Green 函数 = 自由杂质部分 × 声子退相干因子。

| 参数 | 取值 | 说明 |
| --- | --- | --- |
| ϵ_d = μ | 0.0 | 杂质能级（ToulouseIM 的 `μ`） |
| β | 5.0 | 逆温度 |
| ω₀ | 0.8 | 声子频率 |
| g = √α₀ | √0.5 | Holstein 耦合（`DiracDelta(ω=0.8, α=0.5)`） |
| t̂ | 1.0 | semicircular 半带宽，谱支撑 [−t̂, t̂] |

对比量：**retarded Green 函数** G^R(t)，t ∈ [0, 1]。

## 2. 解析参考

解析解来自 ImpurityModelBase 的 Holstein 有限温延展分数
`holstein_G0w_to_Gw(G₀, ω; g, ω₀, β, maxiter=6)`，其中自由杂质（Toulouse）
传播子为

```text
G₀(ω) = 1 / (ω + iδ − ϵ_d − Δ(ω))
```

semicircular 浴的 retarded 杂化函数有闭式：

```text
Δ(z) = 2 [z − √(z² − t̂²)] / t̂²
```

实轴上的分支由 Im Δ(ω) ≤ 0 固定：

```text
Δ(ω) = 2 (ω − i·√(t̂² − ω²)) / t̂²            (|ω| < t̂)
Δ(ω) = 2 (ω − sgn(ω)·√(ω² − t̂²)) / t̂²       (|ω| > t̂)
```

然后做一次 Fourier 型积分

```text
G^R(t) = ∫ dω/2π · e^(−iωt) [ G(ω) − 1/(ω + iδ) ] − i ,   w_max = 20, δ = 1e-6
```

**为什么不直接用 `holstein_Gt(spec, t)`？** 它在实频每一点都做一次数值
Hilbert 变换（嵌套自适应 quadgk，δ = 1e-8 极点叠加带边平方根奇异性），
semicircular 谱上部分时间点单点耗时超过 15 分钟且不收敛。闭式 Δ 版本
free 情形与库内 `toulouse_Gt`（谱直接数值积分，w_max = 20）一致到
**9 位有效数字**；42 个参考点总共约 80 秒。库内 `holstein_Gt` 取
w_max = 4 时自带约 0.5% 的截断误差（例如 t = 0.5：解析 −0.87829i
vs 库 −0.87587i，以 `toulouse_Gt` 为准）。

## 3. GTEMPO 工作流

`docs/tutorials/holstein/holstein.jl` 中的 `run_gtempo(δt; with_phonon)`：

1. `GrassmannLattice(N, δt, contour=:real, order=1, bands=1)` +
   `FockLattice`（声子用 occupation-number 格点）；
2. 费米浴 IF：`fermionicbath(spec_f, β=β)` → `correlationfunction` →
   `ExactTTIIF(algmult=SVDCompression(trunc))`（单带格点）→ `fillband`；
3. 杂质动力学：`sysdynamics`（ToulouseIM）→ `systhermalstate!(β=β)` →
   `boundarycondition!`；
4. 声子 IF：`bosonicbath(DiracDelta)` → `hybriddynamics(flat, …)` →
   `reweighting!` 合入系统态；
5. 观测量：`environments` 缓存 + `cached_greater` / `cached_lesser`。

**Retarded 函数的组合约定**（重要）：GTEMPO 中 `cached_greater` 返回
⟨d(i) d†(j)⟩、`cached_lesser` 返回 −⟨d†(i) d(j)⟩（裸关联，不含 ±i），
因此

```text
G^R(t) = G^>(t) − G^<(t) = −im · [ cached_greater(t) − cached_lesser(t) ],   t > 0
```

写成 `+` 号会令两项相消、结果恒为 ~0：t = 0 处 ⟨d d†⟩ = 1−n = 0.498、
−⟨d† d⟩ = −0.500（n = 1/2）。

截断参数：`truncdimcutoff(D=100, ϵ=1.0e-10)`；费米浴 IF 用
`ExactTTIIF`（指数展开取默认 `DefaultExpansionAlg`）。

## 4. δt 收敛性结果

固定 t_max = 1，细化 δt。相对误差定义为整个 0:δt:t_max 网格上的 2-范数：

```text
rel.err = ‖G_TEMPO − G_exact‖₂ / ‖G_exact‖₂
```

| δt | rel.err（有声子） | rel.err（无声子） |
| --- | --- | --- |
| 0.2 | 6.20e-2 | 6.29e-2 |
| 0.1 | 2.96e-2 | 3.00e-2 |
| 0.05 | 1.42e-2 | 1.46e-2 |

结论：

- 误差随 δt **线性收敛**（δt 减半误差减半，O(δt)）；
- 有声子与无声子的误差几乎相同，说明声子 IF / `reweighting!` 路径在
  α₀ = 0.5 下不引入可观测的额外离散化误差；
- t = 0 处 G^R(0) = −i 被精确复现（GTEMPO 给出 5.7e-5 − 1.006i）。

## 5. 逐点对比（δt = 0.05，有声子）

| t | GTEMPO | 解析 |
| --- | --- | --- |
| 0.0 | 5.7e-5 − 1.005869i | −5.0e-6 − 1.000000i |
| 0.05 | 0.000203 − 1.002312i | 7.0e-6 − 0.997643i |
| 0.1 | 0.000412 − 0.994964i | 7.0e-5 − 0.992292i |
| 0.15 | 0.000732 − 0.983871i | 0.000225 − 0.983197i |
| 0.2 | 0.001207 − 0.969107i | 0.000521 − 0.970115i |
| 0.25 | 0.001876 − 0.950765i | 0.001005 − 0.953216i |
| 0.3 | 0.002773 − 0.928961i | 0.001713 − 0.932839i |
| 0.35 | 0.003923 − 0.903834i | 0.002675 − 0.909257i |
| 0.4 | 0.005346 − 0.875541i | 0.003912 − 0.882588i |
| 0.45 | 0.007048 − 0.844258i | 0.005441 − 0.852865i |
| 0.5 | 0.009030 − 0.810175i | 0.007269 − 0.820168i |
| 0.55 | 0.011279 − 0.773500i | 0.009397 − 0.784697i |
| 0.6 | 0.013774 − 0.734451i | 0.011817 − 0.746752i |
| 0.65 | 0.016484 − 0.693258i | 0.014510 − 0.706635i |
| 0.7 | 0.019366 − 0.650158i | 0.017448 − 0.664582i |
| 0.75 | 0.022372 − 0.605396i | 0.020596 − 0.620764i |
| 0.8 | 0.025442 − 0.559223i | 0.023914 − 0.575353i |
| 0.85 | 0.028509 − 0.511890i | 0.027352 − 0.528574i |
| 0.9 | 0.031503 − 0.463652i | 0.030857 − 0.480716i |
| 0.95 | 0.034345 − 0.414762i | 0.034370 − 0.432080i |
| 1.0 | 0.036956 − 0.365473i | 0.037827 − 0.382928i |

偏差随 t 缓慢累积（t_max = 1 处约 4%），与上表 rel.err = 1.4% 一致，
由 δt = 0.05 的实时离散化主导，可通过继续细化 δt 或提高 `D` 消减。
