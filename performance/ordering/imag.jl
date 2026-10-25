# ============================================================================
# performance/ordering/imag.jl
#
# 单带费米 IF（虚时 / Matsubara contour）的 OSEE 随链上 site 的分布
#
# 两件事：
#  1. 小时间步数下，把 GTEMPO `hybriddynamics` 实际构造出的 GrassmannMPS 的
#     bond 熵（`entanglement_entropy`）与 QuAPI + GaussianStates 的协方差
#     快速 OSEE（`osee(FermionicGaussianState(fermionic_bdg_matrix(corr)))`）
#     逐 bond 对比，确认两者一致。
#  2. 用快速 OSEE（不必构造 MPS）在更大时间步数下扫出每种 Grassmann
#     ordering 的 OSEE-vs-site 分布，写 CSV + SVG 图，作为挑选 ordering 的依据。
#
# 关键映射（与 GTEMPO 内部 `_latticeindex` 一致）：correlation function 的
# 节点序由 `QuAPI.Dense.contour_ordering` 给出；GTEMPO lattice 的虚时分支比
# correlation function 多一个边界节点（kτ = Nτ + 1），故
#
#     correlation 节点 p  ↔  lattice 节点 i = p + 1     （:τ 分支）
#
# lattice 上 i = 1 的 τ 节点没有物理自由度（恒等张量）。每个格点 (i, conj)
# 映射到 Gaussian 模式：conj ? n+k : k（k 是该节点在 contour_ordering 中的位置）。
#
# 运行（GTEMPO 项目环境；GaussianStates 由脚本挂到 LOAD_PATH）：
#   JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. \
#       performance/ordering/imag.jl
# ============================================================================

using GTEMPO
using QuAPI
using QuAPI.Dense: contour_ordering
using LinearAlgebra
using Printf

# GaussianStates.jl 是独立包，不加入 GTEMPO 的依赖；挂到 LOAD_PATH 即可
const GS_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", "GaussianStates"))
GS_ROOT in LOAD_PATH || push!(LOAD_PATH, GS_ROOT)
using GaussianStates
using GaussianStates: FermionicGaussianState, fermionic_bdg_matrix, osee

const OUTDIR = joinpath(@__DIR__, "results")
isdir(OUTDIR) || mkpath(OUTDIR)

# 单带虚时 ordering（GTEMPO 只导出第一个，其余用限定名）
const ORDERINGS = [
    ("A1Ā1B1B̄1", A1Ā1B1B̄1()),
    ("A1B1B̄1Ā1", GTEMPO.A1B1B̄1Ā1()),
    ("Ā2A1B̄2B1", GTEMPO.Ā2A1B̄2B1()),
]

const δτ = 0.1
const SMALL_N = (2, 4, 8)          # 第一件事：与 GTEMPO 逐 bond 对比
const LARGE_N = (16, 32, 64)       # 第二件事：快速 OSEE 大时间步数扫描

# ---------------------------------------------------------------- 通用工具 --

"""
    chain_modes(lattice, corr) -> (modes, n)

`modes[pos]` 是链上第 `pos` 个格点对应的 Gaussian 模式指标（`a` 为 `p`、
`ā` 为 `n+p`；空向量表示该格点为恒等张量，无自由度）。`n` 为 correlation
function 的 contour 节点数。
"""
function chain_modes(lattice, corr)
    nodes = contour_ordering(corr)
    node_index = Dict(node => k for (k, node) in enumerate(nodes))
    n = length(nodes)
    positions = collect(values(indexmappings(lattice)))
    (length(positions) == length(lattice) && sort!(positions) == collect(1:length(lattice))) ||
        error("ordering $(typeof(OrderingStyle(lattice))) 在单带情形下不是合法排列")
    key_at = Dict(pos => key for (key, pos) in indexmappings(lattice))
    modes = Vector{Vector{Int}}(undef, length(lattice))
    for pos in 1:length(lattice)
        i, conj, br, _ = key_at[pos]
        p = br === :τ ? i - 1 : i          # τ 分支：lattice 比 corr 多一个边界节点
        k = get(node_index, (br, p), 0)
        modes[pos] = k == 0 ? Int[] : (conj ? [n + k] : [k])
    end
    return modes, n
end

"""每个 bond 的快速 OSEE：bond b 的割 = 链前缀格点 1..b 的模式集合（b = 1..L-1）。"""
function fast_profile(modes, s)
    S = Float64[]
    acc = Int[]
    for b in 1:length(modes)-1
        append!(acc, modes[b])
        push!(S, osee(s, acc))
    end
    return S
end

"""GTEMPO 实际构造的 GMPS-IF 的逐 bond 熵（NoTruncation canonicalize 后的精确谱）。"""
function gmps_profile(lattice, corr; D=4096, ϵ=1e-14)
    trunc = truncdimcutoff(D=D, ϵ=ϵ)
    mps = hybriddynamics(lattice, corr, trunc=trunc, band=1)
    canonicalize!(mps, alg=Orthogonalize(trunc=GTEMPO.TK.NoTruncation(), normalize=false))
    M = length(mps)
    S = [b <= M - 1 ? entanglement_entropy(mps; bond=b) : 0.0 for b in 1:M]
    return S[1:M-1], M
end

function make_case(N)
    bath = fermionicbath(semicircular(t=1), β=N * δτ, μ=0)
    lat = GrassmannLattice(N=N, δτ=δτ, contour=:imag, bands=1)
    return lat, correlationfunction(bath, lat)
end

"""按链序指纹合并完全相同的 ordering（同一条链 → 同一个 OSEE 分布）。"""
function group_orderings(N)
    groups = NamedTuple{(:names, :sig, :modes),Tuple{Vector{String},Vector{Tuple},Vector{Vector{Int}}}}[]
    for (oname, ordering) in ORDERINGS
        lat = GrassmannLattice(N=N, δτ=δτ, contour=:imag, ordering=ordering, bands=1)
        corr = correlationfunction(fermionicbath(semicircular(t=1), β=N * δτ, μ=0), lat)
        modes, _ = chain_modes(lat, corr)
        sig = Tuple[Tuple(m) for m in modes]
        k = findfirst(g -> g.sig == sig, groups)
        k === nothing ? push!(groups, (names=[oname], sig=sig, modes=modes)) :
                        push!(groups[k].names, oname)
    end
    return groups
end

# ------------------------------------------------------- 轻量 SVG 绘图（无依赖） --
const _COLORS = ("#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#17becf")
const _DASHES = ("", "7,4", "2,3", "12,4,2,4", "5,3")

_xesc(s) = replace(replace(replace(s, '&' => "&amp;"), '<' => "&lt;"), '>' => "&gt;")

function _nicestep(raw)
    raw <= 0 && return 1.0
    e = 10.0^floor(log10(raw))
    for m in (1.0, 2.0, 2.5, 5.0, 10.0)
        m * e >= raw && return m * e
    end
    return 10.0 * e
end

function _axis(hi, n)
    hi <= 0 && return (1.0, [0.0, 1.0])
    step = _nicestep(hi / n)
    top = step * ceil(hi / step)
    return (top, collect(0.0:step:top))
end

"""
    svg_plot(path, title, subtitle, series, names)

`series` = `(label, xs, ys, style)` 的向量（style 选颜色/线型）；`names` 是
图例下方 “O1 = …” 的对照。写一个独立 SVG 图（不依赖任何绘图包）。
"""
function svg_plot(path, title, subtitle, series, names;
                  xlabel="site (bond)", ylabel="OSEE", width=1000, height=560)
    ml, mr, mt, legend_h = 80, 30, 96, 132
    pw = width - ml - mr
    ph = height - mt - legend_h
    ymax = maximum(maximum(ys) for (_, _, ys, _) in series)
    ytop, yt = _axis(max(ymax, 1e-12) * 1.05, 5)
    xmax = maximum(maximum(xs) for (_, xs, _, _) in series)
    xtop, xt = _axis(Float64(xmax), 8)
    x0, y0 = ml, mt
    px(x) = x0 + pw * x / xtop
    py(y) = y0 + ph - ph * y / ytop

    io = IOBuffer()
    print(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"$width\" height=\"$height\" ",
              "viewBox=\"0 0 $width $height\" font-family=\"monospace\">")
    print(io, "<rect width=\"$width\" height=\"$height\" fill=\"white\"/>")
    print(io, "<text x=\"$(width/2)\" y=\"28\" text-anchor=\"middle\" font-size=\"18\">",
              _xesc(title), "</text>")
    print(io, "<text x=\"$(width/2)\" y=\"$(mt-16)\" text-anchor=\"middle\" font-size=\"13\">",
              _xesc(subtitle), "</text>")
    for y in yt
        print(io, "<line x1=\"$x0\" y1=\"$(py(y))\" x2=\"$(x0+pw)\" y2=\"$(py(y))\" ",
                  "stroke=\"#e6e6e6\" stroke-width=\"1\"/>")
    end
    for (_, xs, ys, st) in series
        pts = join(("$(px(x)),$(py(y))" for (x, y) in zip(xs, ys)), " ")
        dash = _DASHES[mod1(st, length(_DASHES))]
        da = isempty(dash) ? "" : " stroke-dasharray=\"$dash\""
        print(io, "<polyline points=\"$pts\" fill=\"none\" ",
                  "stroke=\"$(_COLORS[mod1(st, length(_COLORS))])\" stroke-width=\"1.8\"$da/>")
    end
    print(io, "<rect x=\"$x0\" y=\"$y0\" width=\"$pw\" height=\"$ph\" fill=\"none\" ",
              "stroke=\"black\" stroke-width=\"1\"/>")
    for y in yt
        print(io, "<text x=\"$(x0-8)\" y=\"$(py(y)+4)\" text-anchor=\"end\" font-size=\"11\">",
                  @sprintf("%.3g", y), "</text>")
    end
    for x in xt
        print(io, "<text x=\"$(px(x))\" y=\"$(y0+ph+18)\" text-anchor=\"middle\" font-size=\"11\">",
                  @sprintf("%.0f", x), "</text>")
    end
    print(io, "<text x=\"$(x0+pw/2)\" y=\"$(y0+ph+42)\" text-anchor=\"middle\" font-size=\"12\">",
              _xesc(xlabel), "</text>")
    print(io, "<text x=\"22\" y=\"$(y0+ph/2)\" text-anchor=\"middle\" font-size=\"12\" ",
              "transform=\"rotate(-90 22 $(y0+ph/2))\">", _xesc(ylabel), "</text>")
    # 图例
    ly = height - legend_h + 22
    lx = ml
    for (lab, _, _, st) in series
        dash = _DASHES[mod1(st, length(_DASHES))]
        da = isempty(dash) ? "" : " stroke-dasharray=\"$dash\""
        print(io, "<line x1=\"$lx\" y1=\"$ly\" x2=\"$(lx+34)\" y2=\"$ly\" ",
                  "stroke=\"$(_COLORS[mod1(st, length(_COLORS))])\" stroke-width=\"2.2\"$da/>")
        print(io, "<text x=\"$(lx+40)\" y=\"$(ly+4)\" font-size=\"12\">", _xesc(lab), "</text>")
        lx += 92
    end
    for (k, nm) in enumerate(names)
        print(io, "<text x=\"76\" y=\"$(ly+26+15*(k-1))\" font-size=\"11\" fill=\"#333\">",
                  _xesc("O$k = $nm"), "</text>")
    end
    print(io, "</svg>")
    write(path, String(take!(io)))
    return path
end

# ---------------------------------------------------------------- 主流程 --

function part1_compare()
    println("\n", "="^96)
    println("1) GTEMPO GMPS bond 熵  vs  QuAPI+GaussianStates 快速 OSEE（虚时，小 N）")
    println("="^96)
    csv = joinpath(OUTDIR, "imag_compare.csv")
    ok = true
    open(csv, "w") do io
        println(io, "N,ordering,bond,key,S_gmps,S_fast,abs_diff")
        for N in SMALL_N
            for (oname, ordering) in ORDERINGS
                lat = GrassmannLattice(N=N, δτ=δτ, contour=:imag, ordering=ordering, bands=1)
                corr = correlationfunction(fermionicbath(semicircular(t=1), β=N * δτ, μ=0), lat)
                key_at = Dict(pos => key for (key, pos) in indexmappings(lat))
                modes, _ = chain_modes(lat, corr)
                s = FermionicGaussianState(fermionic_bdg_matrix(corr))
                S_fast = fast_profile(modes, s)
                S_gmps, M = gmps_profile(lat, corr)
                d = maximum(abs.(S_fast .- S_gmps))
                ok &= d < 1e-6
                @printf("N=%-2d %-12s M=%-3d  max|S_fast-S_gmps| = %.3e   %s\n",
                        N, oname, M, d, d < 1e-6 ? "OK" : "MISMATCH")
                for b in 1:M-1
                    println(io, "$N,\"$oname\",$b,\"$(key_at[b])\",$(S_gmps[b]),$(S_fast[b]),$(abs(S_fast[b]-S_gmps[b]))")
                end
            end
        end
    end
    println(ok ? ">>> 全部一致（< 1e-6）" : ">>> 存在不一致，见上表")
    println("逐 bond 数据 -> ", csv)
    return ok
end

function part2_profiles()
    println("\n", "="^96)
    println("2) 单带虚时 IF：各 ordering 的 OSEE-vs-site（快速 OSEE，大时间步数）")
    println("="^96)
    csv = joinpath(OUTDIR, "imag_profiles.csv")
    open(csv, "w") do io
        println(io, "N,ordering,bond,key,S")
        for N in LARGE_N
            lat, corr = make_case(N)
            key_at = Dict(pos => key for (key, pos) in indexmappings(lat))
            s = FermionicGaussianState(fermionic_bdg_matrix(corr))
            groups = group_orderings(N)
            println("-"^96)
            @printf("N = %d  （链长 %d，%d 个 bond；%d 个声明 ordering → %d 条不同链）\n",
                    N, length(lat), length(lat) - 1, length(ORDERINGS), length(groups))
            series = Tuple{String,Vector{Int},Vector{Float64},Int}[]
            for (k, g) in enumerate(groups)
                S = fast_profile(g.modes, s)
                println("  O$k = ", join(g.names, " ≡ "))
                @printf("       max S = %.6f @ bond %d,  mean S = %.6f\n",
                        maximum(S), argmax(S), sum(S) / length(S))
                for b in 1:length(S)
                    println(io, "$N,\"$(join(g.names, "|"))\",$b,\"$(key_at[b])\",$(S[b])")
                end
                push!(series, ("O$k", collect(1:length(S)), S, k))
            end
            svg = joinpath(OUTDIR, "imag_osee_profile_N$(N).svg")
            svg_plot(svg, "虚时单带 IF 的 OSEE 分布（N=$N, δτ=$δτ, β=$(N*δτ)）",
                     "共 $(length(ORDERINGS)) 个声明 ordering，实际 $(length(groups)) 条不同链",
                     series, [join(g.names, " ≡ ") for g in groups])
            println("  图 -> ", svg)
        end
    end
    println("数据 -> ", csv)
end

if abspath(PROGRAM_FILE) == @__FILE__
    part1_compare()
    part2_profiles()
end
