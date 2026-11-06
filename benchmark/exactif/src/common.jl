# =============================================================================
# exactif/common.jl —— IF 严格核方法（exact_Δt）公共部分
#
# 参考：N. Ng, G. Park, A. J. Millis, G. K.-L. Chan, D. R. Reichman,
#       PRB 107, 125103 (2023)（MPS-IF），及 QuAPI 包 exactfermionic.jl。
#
# 方法：用严格核 exact_Δt（QuAPI 的 exact_fermionic_Δt：耦合二次型 cell 的
# 单粒子传播子 + bath Grassmann 高斯 Schur 消元，见 QuAPI/src/exactfermionic.jl）
# 替代 QuAPI 冻结核 Δt（cell 三角积分），在 GTEMPO 实时 Keldysh 管线中：
#   - IF：严格核的稠密高斯 Δ 经行扫 partialmpo 装配（partialif 同构，无 Prony）；
#   - K：杂质演化按 Trotter 调度放置（order 1/2/4），bath cell 由严格核精确处理。
#
# Trotter 调度（cells = 二次型 cell 长度 ×N_c，links = 杂质 link 时长 ×(N_c+1)）：
#   order 1: e^{-iHδt} ≈ A(δt) B(δt)，cells = fill(δt, N)
#   order 2: e^{-iHδt} ≈ A(δt/2) B(δt) A(δt/2)，cells = fill(δt, N)
#   order 4: S4(δt) = S2(aδt)² S2(bδt)，2a+b=1（b<0），cells = repeat([a,a,b]δt, N)
#
# 精度结构（已由 IFGTEMPO 对连续 Fock 精确演化的数值验证）：
#   冻结核 + 均匀 cell 链        → 全局 O(δt)
#   严格核 + 均匀 cell 链 (1/2)  → 全局 O(δt²)（精确顶点幅已重求和 [A,B]δt² 项）
#   严格核 + Suzuki 链 (4)       → 全局 O(δt⁴)
#
# 符号约定（已由 ε_d=0 空杂质 Z·prefactor == 1 的机器精度恒等式标定）：
# 适配器按读序节点取 v = Δ[ā(i 节点), a(j 节点)]，跨支块内嵌 Keldysh 负号
# （v·(-1)^[branch(i)≠branch(j)]），与 QuAPI 统一实时约定一致；GTEMPO 行扫
# 的跨支补偿读回物理块。
# =============================================================================
using LinearAlgebra, Printf
using GTEMPO
using Z2Tensors: TruncationScheme
using ImpurityModelBase
# test/util.jl 的 ED 参考解（singlemode_ed / greater_lesser_ed）；ordering 表引用
# 了两个未导出的 ordering 结构体
using GTEMPO: A1B1B̄1Ā1, A1B1ā1b̄1Ā1B̄1a1b1

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(GT_ROOT, "test", "util.jl"))
include(joinpath(GT_ROOT, "test", "normalbath", "util.jl"))

# ------------------------------------------------------------- Trotter 调度 --
"""
    trotter_factors(order, N, δt) -> TrotterFactors

`order ∈ (1,2,4)` 的 Trotter 因子调度。GTEMPO 单扇区格点上，每对相邻切片
（pair）承载一个 Trotter 因子：

- order 1（同 pair 放置）：每个 pair 同时承载杂质门 A(δt) 与二次型 cell B(δt)
  （A·B 同区间乘积，无需交错），`npairs = N`；
- order 2（交错放置）：合并链 A(δt/2)·[B(δt)·A(δt)]^{N−1}·B(δt)·A(δt/2) 的
  2N+1 个因子交替占据 pair（奇 = 杂质门，偶 = cell），`npairs = 2N+1`；
- order 4（交错放置）：Suzuki 链 S4(δt) = S2(aδt)²S2(bδt)（2a+b=1，b<0）合并
  后每宏步 7 因子 A B A B A B A，`npairs = 6N+1`。

字段：
- `gate_durations`：每 pair 的杂质门时长（0 = 无门）；
- `cell_durations`：每 pair 的二次型 cell 长度（0 = 恒等 cell；传给 exact_Δt
  时以 0 cell 补齐，使核节点与格点切片一一对应）；
- `node_times`：各切片的物理时刻标签（order 1/2 由因子时长平铺 [0, T]；
  order 4 取杂质门累加时刻——Suzuki 子步的形式时长不物理平铺）；
- `obs_layers`：取数切片（order 4 仅取门后节点，其时刻即 links 累加）。
"""
struct TrotterFactors
    npairs::Int
    gate_durations::Vector{Float64}
    cell_durations::Vector{Float64}
    node_times::Vector{Float64}
    obs_layers::Vector{Int}
end

function trotter_factors(order::Integer, N::Integer, δt::Real)
    N >= 1 || throw(ArgumentError("number of macro steps must be at least one"))
    dtf = Float64(δt)
    if order == 1
        npairs = N
        gates = fill(dtf, N)
        cells = fill(dtf, N)
        node_times = collect(0:dtf:(N * dtf))
        obs_layers = collect(1:(N + 1))
    elseif order == 2
        npairs = 2N + 1
        gates = zeros(npairs); cells = zeros(npairs)
        links = vcat(dtf / 2, fill(dtf, N - 1), dtf / 2)      # A 因子时长（N+1 个）
        gates[1:2:end] .= links
        cells[2:2:end] .= dtf                                  # B cells（N 个）
        # 合并链的节点物理平铺在半步网格上：pair 1 → [0, δt/2]（A₁），pair 2 →
        # [δt/2, δt]（B₁），…（与 IFGTEMPO 标定的层时刻 (m−½)δt 约定一致）；
        # 末端边界与 T 重合（截断）
        node_times = min.([(j - 1) * dtf / 2 for j in 1:(npairs + 1)], N * dtf)
        obs_layers = collect(1:(npairs + 1))
    elseif order == 4
        a = 1 / (2 - cbrt(2))
        b = 1 - 2a
        c = (a + b) / 2
        # 合并链的 A 因子时长（3N+1 个）与 B cell 时长（3N 个）
        links = vcat(a * dtf / 2, a * dtf, c * dtf)
        for _ in 2:N
            links = vcat(links, c * dtf, a * dtf, c * dtf)
        end
        links = vcat(links, b * dtf / 2)
        cellseq = repeat([a * dtf, a * dtf, b * dtf], N)
        npairs = 2 * length(cellseq) + 1                       # 6N+1，奇 A 偶 B
        gates = zeros(npairs); cells = zeros(npairs)
        gates[1:2:end] .= links
        cells[2:2:end] .= cellseq
        @assert length(gates[1:2:end]) == length(cellseq) + 1
        # Suzuki 子步的形式时长不物理平铺：节点时刻取杂质门累加时刻
        #（links 平铺 [0, T]；观测仅取门后节点），cell 后节点取相邻门时刻中点
        gatenodes = vcat(0.0, cumsum(links))                   # 3N+2 个门时刻
        node_times = [(iseven(j) ? gatenodes[j ÷ 2 + 1] :
                        (gatenodes[(j - 1) ÷ 2 + 1] + gatenodes[(j + 1) ÷ 2 + 1]) / 2)
                      for j in 1:(npairs + 1)]
        obs_layers = collect(2:2:npairs+1)                     # 门后节点
    else
        throw(ArgumentError("unsupported Trotter order $order (use 1, 2 or 4)"))
    end
    @assert length(node_times) == npairs + 1
    @assert length(obs_layers) <= npairs + 1
    return TrotterFactors(npairs, gates, cells, node_times, obs_layers)
end

# --------------------------------------------------- 严格核 exact_Δt 适配器 --
"""
    ExactIFKernel

严格核 exact_Δt 在 GTEMPO 格点腿空间上的稠密高斯 Δ（QuAPI `ExactGrassmannIF`
的 hmat）+ 整体 prefactor。作为 `QuAPI.AbstractRealCorrelationFunction` 的
子类型直接供 `hybriddynamics!` 的行扫消费。
"""
mutable struct ExactIFKernel <: QuAPI.AbstractRealCorrelationFunction
    pos::Dict{GrassmannContourIndex,Int}
    Δ::Matrix{ComplexF64}
    prefactor::ComplexF64
end
GTEMPO.TK.scalartype(::Type{<:ExactIFKernel}) = ComplexF64
branches(::ExactIFKernel) = (:+, :-)

function Base.getindex(x::ExactIFKernel, i::GrassmannContourIndex, j::GrassmannContourIndex)
    (i.conj == j.conj) && return zero(ComplexF64)
    # 读序节点规则（与冻结核访问器的 block[i.j, j.j] 模式同构）：
    # 系数 = Δ[ā(i 节点), a(j 节点)]，i/j 谁带 conj 均按读序取
    v = x.Δ[x.pos[GrassmannContourIndex(i.j, i.branch, i.band, true)],
            x.pos[GrassmannContourIndex(j.j, j.branch, j.band, false)]]
    # 统一实时约定：跨支块内嵌 Keldysh 负号（行扫读数处补偿回来）
    return (i.branch != j.branch) ? -v : v
end

"""
    exact_Δt(bath, cells::AbstractVector{<:Real}) -> ExactIFKernel

严格核 exact_Δt：bath + 耦合的二次型 cell 演化被精确处理（单粒子传播子
W(τ) = exp(-iτh)），杂质相干腿在节点插入，bath 腿经 Schur 精确积掉。
`cells[k]` 是第 k 个 cell 的长度（0 = 恒等 cell；order 4 的 Suzuki 内层 cell
为负）。收缩时必须把 `prefactor` 乘进 Z（见 QuAPI exactfermionic.jl）。
"""
function exact_Δt(bath, cells::AbstractVector{<:Real})
    ex = exact_fermionic_Δt_cells(bath.spectrum; β=bath.β, cells=collect(Float64.(cells)), μ=bath.μ)
    pos = Dict{GrassmannContourIndex,Int}(l => k for (k, l) in enumerate(ex.hmat.basis))
    return ExactIFKernel(pos, ex.hmat.Δ, ex.prefactor)
end

# ------------------------------------------------------ 按调度的系统演化 K --
"""
    sysdynamics_factors!(gmps, lat, model, gates; trunc) -> GrassmannMPS

按因子调度放置杂质传播门（`sysdynamics` 的调度版）：pair i 上
`gate_durations[i] ≠ 0` 时把 `fock_propagator(model, branch, τ)` 放在切片
(i → i+1) 之间（前向支沿时间方向，后向支 bra/ket 互换、时长取同一表）。
时不变模型仅按时长缓存 SparseGMPS。
"""
function sysdynamics_factors!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice,
                              model::AbstractImpurityHamiltonian, gates::TrotterFactors;
                              trunc::TruncationScheme=GTEMPO.Defaults.alg_trunc_sys())
    (lattice.Nt == gates.npairs) ||
        throw(DimensionMismatch("lattice must carry Nt == npairs slice pairs"))
    # 注意：τ = 0 的门**不能跳过**——恒等传播子张量携带该 pair 的相干态分辨率
    # 重叠 ∏(1 + c̄c)，是 Trotter 化相干态路径积分测度的一部分（每个 Trotter
    # 因子各有一套分辨率）；跳过会丢测度（Z 偏 ~1/2）。
    cache = Dict{Tuple,Any}()
    for branch in (:+, :-)
        for i in 1:gates.npairs
            fm = GTEMPO.fock_propagator(model, branch, gates.gate_durations[i])
            sparse = GTEMPO._propagator_sparsegmps(lattice, fm, i, branch, cache, nothing)
            gmps = mult!(gmps, sparse, trunc=trunc)[1]
        end
    end
    return gmps
end

# ----------------------------------------------------------------- 完整管线 --
"""
    exactif_pipeline(bath, model, N, δt; order, trunc, β) ->
        (lat = lat, mpsK = mpsK, mpsI = mpsI, kernel = kernel, factors = factors)

严格核管线的实时 Keldysh 演化：exact_Δt 行扫 IF（cell 调度零补齐到每个
pair）+ 按调度的杂质门 + Keldysh 折返闭合 + 杂质热态初态。`kernel` 是严格核
对象（收缩 Z 须再乘 `kernel.prefactor`）；取样时刻见 `factors.node_times` /
`factors.obs_layers`。
"""
function exactif_pipeline(bath, model::AbstractImpurityHamiltonian, N::Integer, δt::Real;
                          order::Integer=2, trunc::TruncationScheme=Defaults.alg_trunc_if(),
                          β::Real=bath.β, bands::Int=1)
    factors = trotter_factors(order, N, δt)
    lat = GrassmannLattice(N=factors.npairs, δt=δt, contour=:real, bands=bands)
    kernel = exact_Δt(bath, factors.cell_durations)
    mpsI1 = hybriddynamics(lat, kernel; trunc=trunc, band=1)
    mpsIs = GrassmannMPS[mpsI1]
    for band in 2:bands
        push!(mpsIs, swapband(mpsI1, lat, 1, band; trunc=trunc))
    end
    mpsK = sysdynamics_factors!(vacuumstate(ComplexF64, lat), lat, model, factors; trunc=trunc)
    for band in 1:bands
        mpsK = boundarycondition!(mpsK, lat, band=band, trunc=trunc)
    end
    mpsK = systhermalstate!(mpsK, lat, model; trunc=trunc, β=β)
    return (lat=lat, mpsK=mpsK, mpsIs=mpsIs, kernel=kernel, factors=factors)
end

"""
    frozen_pipeline(bath, model, N, δt; trunc, β) -> (lat, mpsK, mpsI)

冻结核（QuAPI Δt）参考管线：`correlationfunction` + `hybriddynamics` +
`sysdynamics`（均匀全步长杂质演化，order 1 放置）。
"""
function frozen_pipeline(bath, model::AbstractImpurityHamiltonian, N::Integer, δt::Real;
                         trunc::TruncationScheme=Defaults.alg_trunc_if(), β::Real=bath.β,
                         bands::Int=1)
    lat = GrassmannLattice(N=N, δt=δt, contour=:real, bands=bands)
    corr = correlationfunction(bath, lat)
    mpsI1 = hybriddynamics(lat, corr; trunc=trunc, band=1)
    mpsIs = GrassmannMPS[mpsI1]
    for band in 2:bands
        push!(mpsIs, swapband(mpsI1, lat, 1, band; trunc=trunc))
    end
    mpsK = sysdynamics(lat, model; trunc=trunc)
    for band in 1:bands
        mpsK = boundarycondition!(mpsK, lat, band=band, trunc=trunc)
    end
    mpsK = systhermalstate!(mpsK, lat, model; trunc=trunc, β=β)
    return lat, mpsK, mpsIs
end

"""
    gtlt_exactif(pipe; band=1) -> (gt, lt, ts)

严格核管线的 greater/lesser 系列取数：取样切片为 `factors.obs_layers`，对应
物理时刻 `factors.node_times[obs_layers]`；Z 含核 prefactor（归一化后消去）。
"""
function gtlt_exactif(pipe; band=1)
    (; lat, mpsK, mpsIs, kernel, factors) = pipe
    Z = integrate(lat, mpsK, mpsIs...) * kernel.prefactor
    gt = [-im * gf(lat, (GrassmannContourIndex(k, conj=false, branch=:+, band=band),
                         GrassmannContourIndex(1, conj=true, branch=:+, band=band)),
                    mpsK, mpsIs...; Z=Z) for k in factors.obs_layers]
    lt = [-im * gf(lat, (GrassmannContourIndex(1, conj=true, branch=:-, band=band),
                        GrassmannContourIndex(k, conj=false, branch=:+, band=band)),
                   mpsK, mpsIs...; Z=Z) for k in factors.obs_layers]
    ts = factors.node_times[factors.obs_layers]
    return gt, lt, ts
end

"""
    ntT(pipe) -> Float64

末时刻杂质占据数 ⟨n̂(T)⟩ 的**折返跨越投影**提取：在折返点自由腿对
`(a⁺_k, ā⁻_k)` 上施加 2 腿 `GTerm` 单项式（稠密参考 ρ(T) 投影的同一约定；
折返双线性与插入共享腿，其贡献按 Grassmann 幂零性自动相消）。相比
`occupation`（清零位点张量偶宇称块的插入约定，带 O(δt) 有效时刻偏移），
该提取对严格核链与 Fock 直算一致到截断噪声量级（N=1..8 逐点验证）。
"""
function ntT(pipe)
    (; lat, mpsK, mpsIs, kernel) = pipe
    Z = integrate(lat, mpsK, mpsIs...) * kernel.prefactor
    k = lat.k
    mpsN = copy(mpsK)
    pos_a = lat[GrassmannContourIndex(k, conj=false, branch=:+, band=1)]
    pos_ab = lat[GrassmannContourIndex(k, conj=true, branch=:-, band=1)]
    GTEMPO.apply!(GTEMPO.GTerm(pos_a, pos_ab, coeff=1.0), mpsN)
    return real(integrate(lat, mpsN, mpsIs...) * kernel.prefactor / Z)
end

"""
    gtlt_frozen(lat, mpsK, mpsIs) -> (gt, lt)

冻结核管线的 greater/lesser 系列（无 prefactor，标准 gtlt_series 公式）。
"""
function gtlt_frozen(lat, mpsK, mpsIs::Vector{<:GrassmannMPS})
    Z = integrate(lat, mpsK, mpsIs...)
    gt = [-im * gf(lat, (GrassmannContourIndex(k, conj=false, branch=:+, band=1),
                         GrassmannContourIndex(1, conj=true, branch=:+, band=1)),
                    mpsK, mpsIs...; Z=Z) for k in 1:lat.k]
    lt = [-im * gf(lat, (GrassmannContourIndex(1, conj=true, branch=:-, band=1),
                        GrassmannContourIndex(k, conj=false, branch=:+, band=1)),
                   mpsK, mpsIs...; Z=Z) for k in 1:lat.k]
    return gt, lt
end

# ------------------------------------------------------------------ 误差工具 --
relerr(a, b) = norm(a - b) / norm(a)
_m(v) = sum(v) / length(v)
# 误差取自然对数、横轴取 -log2(N)（∝ log δ）：斜率即「误差对 δ 的阶」
slope(hs, es) = (sum((hs .- _m(hs)) .* log.(es)) / sum(abs2, hs .- _m(hs))) / log(2)

# 与具体基准无关的 matplotlib 画图脚本：argv = csv, png, title, xlabel, ylabel。
# csv 列：series,dt,relerr；每个 series 画一条 log-log 线并标注最小二乘拟合阶，
# 另画锚在首条曲线首点的 O(δt)/O(δt²)/O(δt⁴) 参考线。
const PY_PLOTTER = raw"""
import sys, csv, math
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

csv_path, png_path, title, xlabel, ylabel = sys.argv[1:6]
series = defaultdict(list)
with open(csv_path) as f:
    for row in csv.DictReader(f):
        series[row["series"]].append((float(row["dt"]), float(row["relerr"])))

fig, ax = plt.subplots(figsize=(7.6, 5.2))
colors = ["#7f7f7f", "#d62728", "#1f77b4", "#2ca02c", "#9467bd"]
dts = None
for i, ((name, pts), c) in enumerate(zip(sorted(series.items()), colors)):
    pts.sort()
    dt = [p[0] for p in pts]
    er = [p[1] for p in pts]
    dts = dt
    ax.loglog(dt, er, "o-", color=c, label=name)
    x = [math.log(d) for d in dt]
    y = [math.log(e) for e in er]
    n = len(x)
    xm, ym = sum(x) / n, sum(y) / n
    p = sum((xi - xm) * (yi - ym) for xi, yi in zip(x, y)) / \
        sum((xi - xm) ** 2 for xi in x)
    ax.annotate(f"order = {p:.2f}", xy=(dt[-1], er[-1]),
                xytext=(4, -13 if i % 2 == 0 else 9),
                textcoords="offset points", fontsize=9, color=c)

# O(dt^p) 参考线：锚在第一条曲线的首/尾点
_, pts0 = sorted(series.items())[0]
d1, e1 = min(pts0)
d2 = max(dts)
for p, ls in ((1, "--"), (2, ":"), (4, "-.")):
    ax.loglog([d1, d2], [e1, e1 * (d2 / d1) ** p], ls, color="gray", lw=1,
              label=rf"$O(\delta t^{p})$")

ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_title(title, fontsize=11)
ax.grid(True, which="both", alpha=0.3)
ax.legend(loc="best", fontsize=9)
fig.tight_layout()
fig.savefig(png_path, dpi=150)
print("图 ->", png_path)
"""

function plot_png(csv, png, title, xlabel, ylabel)
    py = tempname()
    open(py, "w") do io
        write(io, PY_PLOTTER)
    end
    try
        run(`python3 $py $csv $png $title $xlabel $ylabel`)
    finally
        rm(py; force=true)
    end
end

"""
    scaling_main(runner; Ns, T, csvfile, pngfile, title, D, ϵ)

通用标度驱动。`runner(N; trunc)` 跑一条管线，返回
`Pair("序列名", 相对误差)` 的有序向量。扫描 Ns：打印误差表（含相邻误差比的
阶数）、写 CSV（series,dt,relerr）、拟合 log-log 斜率并画标度图。
"""
function scaling_main(runner; Ns, T, csvfile, pngfile, title, D=200, ϵ=1.0e-12)
    trunc = truncdimcutoff(D=D, ϵ=ϵ)
    hs = -log2.(collect(Ns))
    names = String[]
    errs = Vector{Vector{Float64}}()   # errs[j][i]：第 j 条序列在 Ns[i] 的误差
    open(csvfile, "w") do io
        println(io, "series,dt,relerr")
        for (i, N) in enumerate(Ns)
            t0 = time()
            res = runner(N; trunc=trunc)
            if isempty(names)
                append!(names, first.(res))
                errs = [Float64[] for _ in names]
                @printf("%-6s%-12s%s\n", "N", "δt",
                        join(rpad("$(names[j]) relerr", 16) for j in eachindex(names)))
            end
            es = last.(res)
            for (j, e) in enumerate(es)
                push!(errs[j], e)
                println(io, names[j], ",", T / N, ",", e)
            end
            ostr = join((i == 1 ? rpad("  --", 8) :
                         rpad(@sprintf("%5.2f", log2(errs[j][i-1] / es[j])), 8)
                         for j in eachindex(es)))
            @printf("%-6d%-12.6g%s%s  (%.1fs)\n", N, T / N,
                    join(rpad(@sprintf("%.6e", e), 16) for e in es), ostr, time() - t0)
        end
    end
    println("拟合收敛阶：",
            join((rpad("$(names[j]) p = $(@sprintf("%.2f", slope(hs, errs[j])))", 26)
                  for j in eachindex(names))))
    plot_png(csvfile, pngfile, title, "δt = T / N", "relative error of ĝ>/ĝ< series")
    println("数据 -> ", csvfile)
    return errs
end
