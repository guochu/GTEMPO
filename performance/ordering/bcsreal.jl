# ============================================================================
# performance/ordering/bcsreal.jl
#
# **BCS bath（单轨道 = 2 band，↑↓ 耦合）** 实时 IF 的 OSEE 随链上 site 的分布
#
# 两件事（与 real.jl 相同的流程）：
#  1. 小时间步数下，GTEMPO `hybriddynamics(..., orbital=1)` 构造的 GrassmannMPS
#     的 bond 熵  vs  GaussianStates 协方差快速 OSEE。
#  2. 用快速 OSEE 在更大时间步数下扫出每种 Grassmann ordering 的
#     OSEE-vs-site 分布，写 CSV + SVG。
#
# IF 指数与配对矩阵 W 的组装见 bcsimag.jl 头部说明（一般复反对称 W，
# 模式序 [a↑; a↓; ā↑; ā↓]）。
#
# **实时分支的节点约定**：GTEMPO 的 BCS PartialIF 对 :+/− 分支只用
# correlation function 的前 Nt 个节点（lattice 的最后一个 +/− 节点是恒等
# 张量），故本脚本的节点表 / W 也取前 Nt 个节点（与 GTEMPO 构造出的 IF 一致）。
# 单带（非 BCS）的 real.jl 用满 kt = Nt+1 个节点，这里不同。
#
# 运行（GTEMPO 项目环境；GaussianStates 由脚本挂到 LOAD_PATH）：
#   JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. \
#       performance/ordering/bcsreal.jl
# ============================================================================

using GTEMPO
using LinearAlgebra
using Printf

const GS_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", "GaussianStates"))
GS_ROOT in LOAD_PATH || push!(LOAD_PATH, GS_ROOT)
using GaussianStates
using GaussianStates: FermionicGaussianState, osee

const OUTDIR = joinpath(@__DIR__, "results")
isdir(OUTDIR) || mkpath(OUTDIR)

const ORDERINGS = [
    ("A1Ā1a1ā1B1B̄1b1b̄1", A1Ā1a1ā1B1B̄1b1b̄1()),
    ("A1Ā1B1B̄1a1ā1b1b̄1", A1Ā1B1B̄1a1ā1b1b̄1()),
    ("time-ascending", a1ā1A1Ā1b1b̄1B1B̄1a2ā2A2Ā2b2b̄2B2B̄2()),
    ("A1B1ā1b̄1Ā1B̄1a1b1", GTEMPO.A1B1ā1b̄1Ā1B̄1a1b1()),
    ("Ā2A1ā1a2B̄2B1b̄1b̄2", GTEMPO.Ā2A1ā1a2B̄2B1b̄1b̄2()),
]

const δt = 0.1
const βbath = 4.0
const GAP = 0.2 + 0.3im
const SMALL_N = (1, 2, 3)          # 第一件事：与 GTEMPO 逐 bond 对比
const LARGE_N = (2, 4, 8)          # 第二件事（Stein 解 n⁶，N 上限 ~8）

# ---------------------------------------------------------------- BCS 工具 --

"""与 GTEMPO BCS PartialIF 一致的 contour 节点表（实时只用前 Nt 个节点）。"""
function bcs_nodes(lattice)
    nodes = Tuple{Symbol,Int}[]
    if :τ in GTEMPO.branches(lattice)
        for j in 1:lattice.Nτ
            push!(nodes, (:τ, j))
        end
    end
    if :+ in GTEMPO.branches(lattice)
        for j in 1:lattice.Nt
            push!(nodes, (:+, j))
            push!(nodes, (:-, j))
        end
    end
    return nodes
end

function bcs_block(x, nodes)
    M = zeros(ComplexF64, length(nodes), length(nodes))
    for (p, (b1, i1)) in enumerate(nodes), (q, (b2, i2)) in enumerate(nodes)
        M[p, q] = GTEMPO.branch(x; b1=b1, b2=b2)[i1, i2]
    end
    return M
end

"""4n×4n 复反对称配对矩阵，模式序 `[a↑; a↓; ā↑; ā↓]`。"""
function bcs_pairing(corr, nodes)
    n = length(nodes)
    cc = bcs_block(corr.ud, nodes)
    cn = bcs_block(corr.uu, nodes)
    nc = bcs_block(corr.dd, nodes)
    nn = bcs_block(corr.du, nodes)
    W = zeros(ComplexF64, 4n, 4n)
    ra1 = 1:n; ra2 = n+1:2n; rā1 = 2n+1:3n; rā2 = 3n+1:4n
    W[rā1, ra1] .= cn
    W[rā1, rā2] .= cc
    W[ra2, rā2] .= nc
    W[ra2, ra1] .= nn
    return W - transpose(W)
end

"""链上每个格点的 Gaussian 模式指标（band1=↑，band2=↓；空 = 恒等张量）。"""
function bcs_chain_modes(lattice, nodes)
    n = length(nodes)
    node_index = Dict(nd => k for (k, nd) in enumerate(nodes))
    positions = collect(values(indexmappings(lattice)))
    (length(positions) == length(lattice) && sort!(positions) == collect(1:length(lattice))) ||
        error("ordering $(typeof(OrderingStyle(lattice))) 在双带情形下不是合法排列")
    key_at = Dict(pos => key for (key, pos) in indexmappings(lattice))
    modes = Vector{Vector{Int}}(undef, length(lattice))
    for pos in 1:length(lattice)
        i, conj, br, band = key_at[pos]
        p = br === :τ ? i - 1 : i
        k = get(node_index, (br, p), 0)
        modes[pos] = if k == 0
            Int[]
        elseif band == 1
            conj ? [2n + k] : [k]
        else
            conj ? [3n + k] : [n + k]
        end
    end
    return modes
end

function fast_profile(modes, s)
    S = Float64[]
    acc = Int[]
    for b in 1:length(modes)-1
        append!(acc, modes[b])
        push!(S, osee(s, acc))
    end
    return S
end

function gmps_profile(lattice, corr; D=1024, ϵ=1e-14)
    trunc = truncdimcutoff(D=D, ϵ=ϵ)
    mps = hybriddynamics(lattice, corr, orbital=1, trunc=trunc)
    canonicalize!(mps, alg=Orthogonalize(trunc=GTEMPO.TK.NoTruncation(), normalize=false))
    M = length(mps)
    S = [b <= M - 1 ? entanglement_entropy(mps; bond=b) : 0.0 for b in 1:M]
    return S[1:M-1], M
end

make_lattice(N; ordering=A1Ā1a1ā1B1B̄1b1b̄1()) =
    GrassmannLattice(N=N, δt=δt, contour=:real, bands=2, ordering=ordering)
make_corr(N) = correlationfunction(
    bcsbath(fermionicbath(semicircular(t=1), β=βbath), Δ=GAP), make_lattice(N))

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

function group_orderings(N, nodes)
    groups = NamedTuple{(:names, :sig, :modes),Tuple{Vector{String},Vector{Tuple},Vector{Vector{Int}}}}[]
    skipped = String[]
    for (oname, ordering) in ORDERINGS
        lat = make_lattice(N; ordering=ordering)
        modes = try
            bcs_chain_modes(lat, nodes)
        catch
            push!(skipped, oname)
            continue
        end
        sig = Tuple[Tuple(m) for m in modes]
        k = findfirst(g -> g.sig == sig, groups)
        k === nothing ? push!(groups, (names=[oname], sig=sig, modes=modes)) :
                        push!(groups[k].names, oname)
    end
    return groups, skipped
end

function part1_compare()
    println("\n", "="^96)
    println("1) BCS bath：GTEMPO GMPS bond 熵  vs  快速 OSEE（实时，小 N）")
    println("="^96)
    csv = joinpath(OUTDIR, "bcs_real_compare.csv")
    ok = true
    open(csv, "w") do io
        println(io, "N,ordering,bond,key,S_gmps,S_fast,abs_diff")
        for N in SMALL_N
            lat0 = make_lattice(N)
            corr = make_corr(N)
            nodes = bcs_nodes(lat0)
            s = FermionicGaussianState(bcs_pairing(corr, nodes))
            for (oname, ordering) in ORDERINGS
                lat = similar(lat0, ordering=ordering)
                key_at = Dict(pos => key for (key, pos) in indexmappings(lat))
                modes = bcs_chain_modes(lat, nodes)
                S_fast = fast_profile(modes, s)
                S_gmps, M = gmps_profile(lat, corr)
                d = maximum(abs.(S_fast .- S_gmps))
                ok &= d < 1e-6
                @printf("N=%-2d %-24s M=%-3d  max|S_fast-S_gmps| = %.3e   %s\n",
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
    println("2) BCS bath 实时 IF：各 ordering 的 OSEE-vs-site（快速 OSEE）")
    println("="^96)
    csv = joinpath(OUTDIR, "bcs_real_profiles.csv")
    open(csv, "w") do io
        println(io, "N,ordering,bond,key,S")
        for N in LARGE_N
            lat0 = make_lattice(N)
            corr = make_corr(N)
            nodes = bcs_nodes(lat0)
            s = FermionicGaussianState(bcs_pairing(corr, nodes))
            key_at = Dict(pos => key for (key, pos) in indexmappings(lat0))
            groups, skipped = group_orderings(N, nodes)
            println("-"^96)
            @printf("N = %d  （链长 %d，%d 个 bond；%d 个声明 ordering → %d 条不同链）\n",
                    N, length(lat0), length(lat0) - 1, length(ORDERINGS), length(groups))
            isempty(skipped) || println("  跳过（非法排列）: ", join(skipped, ", "))
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
            svg = joinpath(OUTDIR, "bcs_real_osee_profile_N$(N).svg")
            svg_plot(svg, "BCS bath 实时 IF 的 OSEE 分布（N=$N, δt=$δt, t=$(N*δt), Δ=$(GAP)）",
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
