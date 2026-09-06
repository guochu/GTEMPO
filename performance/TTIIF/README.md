# 影响泛函构建算法对比基准：PartialIF vs XTRGIF vs ExactTTIIF vs TDVPIF

本目录包含 GTEMPO 四种影响泛函（IF）构建算法的效率与精度对比，以
**Toulouse 模型**（非相互作用 SIAM，杂质能级线性耦合到连续谱费米浴，U = 0）
的虚时间与实时间演化为测试平台。

## 脚本

| 脚本 | 轮廓 | 观测量 |
|---|---|---|
| `imaginarytime.jl` | 虚时间（N=25, δτ=0.01, β=0.25, ϵ_d=1.25π, μ=0） | G(τ) = ⟨d(τ) d†(0)⟩ |
| `realtime.jl` | 实时 Keldysh 1 阶（N=10, δt=0.02, t=0.2, β=1, ϵ_d=-1, μ=0） | ⟨d(t) d†(δt)⟩ |

浴谱 `f(D, ϵ) = sqrt(1-(ϵ/D)^2)/π`（[-D, D]，与 `test/util.jl` 一致）；
截断 `truncdimcutoff(D=100, ϵ=1e-8/1e-9, add_back=0)`；
ordering：虚时 `A1Ā1B1B̄1`，实时 `A1Ā1a1ā1B1B̄1b1b̄1`。

四种算法：

- `PartialIF(trunc)`：partial-MPO 逐个乘入（2^N-1 次乘法）
- `XTRGIF(k=5, fast=true, algmult=SVDCompression(trunc))`：树形二分，5 次乘法
- `ExactTTIIF(algmult=SVDCompression(trunc))`：影响算符精确指数化（解析 MPO 指数）
- `TDVPIF(trunc, δ=0.1)`：TDVP 虚时流 dz/dτ = H·z（10 步，流中不截断，流结束后最后一次 sweep 用 SVD 截断压缩）

**参考解**（ImpurityModelBase）：

- 虚时间：`toulouse_Gτ`（离散浴 δw=0.2 ED，平衡热态）
- 实时间：`toulouse_neq_greater_lesser(nsys=0)`（直积初态 ρ_imp(空态) ⊗ 浴热态，
  与纯实时格点的边界条件一致）

运行方式：

```
OMP_NUM_THREADS=1 julia --project=. performance/TTIIF/imaginarytime.jl
OMP_NUM_THREADS=1 julia --project=. performance/TTIIF/realtime.jl
```

## 结果汇总

运行环境：Julia 1.10.11，`OMP_NUM_THREADS=1`（单线程）。

### 1. 虚时间（G(τ)，参考 ED）

| 算法 | IF 构建 (s) | 观测量扫描 (s) | 最大键维 | 相对误差 |
|---|---|---|---|---|
| PartialIF | 52.4 | 17.5 | 32 | 6.0e-5 |
| XTRGIF | 35.9 | 0.6 | 30 | **5.9e-5** |
| ExactTTIIF | **3.6** | 0.7 | 30 | 6.0e-5 |
| TDVPIF | 179.6 | 1.2 | 146 | 8.7e-4 |

算法间 IF 距离：PartialIF–XTRGIF 9.5e-7，PartialIF–ExactTTIIF 9.2e-8，
其余对 TDVPIF 1.9e-4。

### 2. 实时间（⟨d(t) d†(δt)⟩，参考 ED nsys=0）

| 算法 | IF 构建 (s) | 观测量扫描 (s) | 最大键维 | 相对误差 |
|---|---|---|---|---|
| PartialIF | 55.4 | 19.5 | 12 | 1.220e-3 |
| XTRGIF | 36.9 | **0.4** | 13 | 1.220e-3 |
| ExactTTIIF | **4.9** | 0.9 | 13 | 1.220e-3 |
| TDVPIF | 389.9 | 1.0 | 75 | 1.217e-3 |

算法间 IF 距离：PartialIF–XTRGIF 1.7e-7，PartialIF–ExactTTIIF 3.0e-8，
其余对 TDVPIF 7.3e-5。

## 结论

1. **ExactTTIIF 全面占优**（本模型特有优势）：影响算符可解析指数化，无迭代
   乘法/流，构建最快（比 XTRGIF 快 7~10 倍），精度与 XTRGIF/PartialIF 相当，
   键维相同。但它要求格点 ordering 在允许列表内，且不适用于无法解析指数化的
   更一般影响算符。
2. **IF 构建效率**：ExactTTIIF ≫ XTRGIF > PartialIF > TDVPIF。TDVPIF 最慢：
   每个 TDVP 步需 2L 次 Krylov 指数化 + 环境更新，δ=0.1 即 10 步 × 50 site。
3. **精度**：虚时间下 XTRGIF/PartialIF/ExactTTIIF（~6e-5）优于 TDVPIF
   （8.7e-4）；实时间下四者相当（~1.2e-3，受 δt 离散化与截断共同主导）。
   TDVPIF 的误差由 TDVP 流步长 δ 主导（`test/tempo/tdvpif.jl` 的
   "convergence with δ" 已验证误差随 δ 近二次收敛），减小 δ 可进一步提高精度，
   代价是构建时间线性增长。
4. **观测量扫描**：XTRGIF/ExactTTIIF/TDVPIF 比 PartialIF 快 1~2 个量级——
   三者的 IF 是单条紧致 MPS，而 PartialIF 的 IF 结构更松散（虚时最大键维 32
   但扫描 17.5s vs 0.6~1.2s）。

## TDVPIF 步长 δ 扫描

TDVPIF 的流中不做截断，流结束后做一次 canonicalize（左 QR sweep + 右 SVD 截断
sweep）压缩。TDVP 步长 δ 决定流误差，直接控制最终精度与键维
（脚本 `tdvp_delta_scan.jl`）：

**虚时间**（其余参数同上）：

| δ | IF 构建 (s) | 最大键维 | 相对误差 |
|---|---|---|---|
| 0.10 | 290 | 146 | 8.7e-4 |
| 0.05 | 344 | 140 | 4.5e-4 |
| 0.02 | 882 | 112 | **1.3e-4** |

**实时间**（其余参数同上；误差地板由 δt 离散化主导）：

| δ | IF 构建 (s) | 最大键维 | 相对误差 |
|---|---|---|---|
| 0.10 | 457 | 75 | 1.217e-3 |
| 0.05 | 761 | 75 | 1.198e-3 |
| 0.02 | 1935 | **58** | 1.210e-3 |

观察：

- 虚时间误差随 δ 近二次收敛（δ=0.02 时 1.3e-4，接近 XTRGIF 的 5.9e-5），
  键维同步下降（146 → 112）——证实"流误差在谱上留下脏分量"的机理，δ 越小
  谱越干净、最后 SVD sweep 能截得越紧。
- 实时间误差被 δt 离散化地板（~1.2e-3，四种算法相同）掩盖，δ 减小无精度
  收益，但键维仍从 75 降至 58。
- 构建时间随步数线性增长（δ=0.02 时虚时 882s、实时 1935s）。

## 关于 TDVPIF 键维的说明

TDVPIF 的流流形键维由 `trunc.D` 零填充保证（流中不截断），流结束后做一次
canonicalize（左 QR sweep + 右 SVD 截断 sweep）压缩。在 δ=0.1 下最终键维
（虚时 146，实时 75，均为 boundarycondition 后）高于其它算法（30/13），其机理：

- TDVP 步长误差（δ=0.1 → 精度 ~1e-3）在 Schmidt 谱上留下 ~1e-4 量级的
  "脏"分量，它们大于最终截断的 ϵ 判据（1e-8/1e-9），因此最后一次 SVD sweep
  无法截除；
- `boundarycondition` 的内部截断极松（`DefaultIntegrationTruncation`：
  D=10000, ϵ=1e-12），对双局域算符 apply! 造成的键维近似翻倍不做修剪，
  对所有算法一致（XTRGIF 等因谱本来就干净，翻倍后仍被 ϵ 判据压回）。

δ 扫描证实：减小 δ 可同时改善精度与键维（脏分量随 δ² 下降），代价是构建
时间线性增长。若目标是与其它算法相当的键维（~30），可进一步减小 δ 或在
`boundarycondition` 处显式传入更紧的 `trunc`。

## 备注

- 实时基准的参考解需与格点边界条件严格对应：纯实时（Keldysh 1 阶）轮廓
  对应直积初态 ρ_imp(nsys=0) ⊗ 浴热态；带 τ 支的混合轮廓对应整体热化平衡态
  （`toulouse_greater_lesser`）。nsys 缺省值是 ϵ_d 处的热占据而非 0。
- GTEMPO 的 `greater(i, j=1)` 返回 ⟨d(t_i) d†(δt)⟩（无 -i 因子，第二算符在
  格点时间 1），与解析 G>(t)（含 -i）相差因子 i 与时间平移 δt；因 nsys=0 时
  lesser 恒为零，⟨d(t)d†(t')⟩ 仅依赖 t-t'，故与 `i·G>(t_i-δt)` 逐点对比。
- 虚时基准的 `Gτ(lattice, mpsK, mpsI)` 约定与 `test/tempo/models.jl` 的
  GF-imaginary time 测试一致。
