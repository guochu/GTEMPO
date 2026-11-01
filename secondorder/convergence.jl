# ============================================================================
# secondorder/convergence.jl
#
# 验证：GTEMPO 实时演化对 δt 的收敛阶数，以及等时（equal-time）耦合的处理
# 方式对精度阶数的影响。
#
# ── 结构分析（本脚本要数值验证的命题）─────────────────────────────────
# Grassmann-TEMPO 实时方案：
#   1. 无系统-环境 Trotter 误差：IF 门作用于 bath Grassmann 腿，系统传播子
#      作用于 path 腿，二者是链上独立张量因子，逐块收缩精确。
#   2. kernel：非对角元 = 精确双胞积分；等时元 = 两个时序三角积分之和
#      （CorrelationMatrix 的 ηⱼₖ[1]+ηₖⱼ[1]，单胞 O(δt³)）。
#   3. IF 构造对该离散 kernel 精确（test/api/quapi_crosscheck.jl 对 dense IF
#      交叉验证）；系统传播子 e^{-ihδt} 为精确矩阵指数（时不变 h）。
#   4. D=128 与 D=384 结果逐位一致 → 截断不是误差源。
#
# ── 已在诊断中发现的两个独立现象 ──────────────────────────────────────
#   a. gt(t=0) 与解析参考差一个不随 δt 变化的 O(1) 常数（greater 函数在
#      t=0 端点的约定/解析延拓差异）—— 本脚本将其排除出收敛阶拟合。
#   b. 闭合 Keldysh 轮迹的 Z 恒等式（线性耦合 bath 下 Z=1 精确成立）在
#      现有等时项处理下只以 O(δt) 收敛（Z−1: 0.137→0.065，δt 减半）；
#      把等时项置零后 Z≈0.998。⟹ 等时项是 O(δt) 误差的来源。
#
# ── 本脚本的实验 ──────────────────────────────────────────────────────
#   等时项缩放参数 κ：CorrelationMatrix 对角改为 κ·ηⱼₖ[1]、κ·ηₖⱼ[1]。
#   κ=1 为现状；κ=0 为全丢；κ=1/2 为对半。对每个 κ 拟合
#     order(gt) = gt(t>0) 最大相对误差的 log-log 斜率
#     order(Z)  = |Z−1| 的 log-log 斜率
#   若存在 κ* 使两者同时为 O(δt²)，即为实时演化的二阶等时处理方案。
#
# 运行：
#   JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. \
#       secondorder/convergence.jl
# ============================================================================
using GTEMPO
using ImpurityModelBase
using LinearAlgebra
using Printf

const OUTDIR = joinpath(@__DIR__, "results")
isdir(OUTDIR) || mkpath(OUTDIR)

# ------------------------------------------------------------------ 模型 --
const μ_SYS = 1.25 * π
const β_BATH = 1.0
const T_FINAL = 1.0
const D_BOND = 128        # 已验证 D=128 与 D=384 无差别（截断不是误差源）

_spec = spectrum(ϵ -> sqrt(1 - (ϵ/10)^2) / π, lb=-10, ub=10)
const BATH = fermionicbath(_spec, β=β_BATH, μ=0)
const MODEL = ToulouseIM(ϵ_d=μ_SYS)

"等时耦合缩放：CorrelationMatrix 对角改为 κ·(两个三角积分)"
function scale_equaltime(corr::RealCorrelationFunction, κ::Real)
    _sc(B::CorrelationMatrix) = begin
        a = copy(B.ηⱼₖ); b = copy(B.ηₖⱼ)
        a[1] *= κ; b[1] *= κ
        CorrelationMatrix(a, b)
    end
    return RealCorrelationFunction(_sc(corr.G₊₊), _sc(corr.G₊₋),
                                   _sc(corr.G₋₊), _sc(corr.G₋₋))
end

"""
    gt_scheme(δt, corr; D) -> (ts, gt, Z)

当前 GTEMPO 全 pipeline 的 greater Green's function：
hybriddynamics（PartialIF）+ sysdynamics + boundarycondition! + gf。
"""
function gt_scheme(δt::Float64, corr; D::Int=D_BOND)
    N = round(Int, T_FINAL / δt)
    lat = GrassmannLattice(N=N, δt=δt, contour=:real)
    trunc = truncdimcutoff(D=D, ϵ=1.0e-10)
    mpsI = hybriddynamics(lat, corr, trunc=trunc)
    mpsK = sysdynamics(lat, MODEL, trunc=trunc)
    mpsK = boundarycondition!(mpsK, lat)
    Z = integrate(lat, mpsK, mpsI)
    gt = Vector{ComplexF64}(undef, lat.k)
    for k in 1:lat.k
        gt[k] = gf(lat, (GrassmannContourIndex(k, conj=false, branch=:+, band=1),
                         GrassmannContourIndex(1, conj=true, branch=:+, band=1)),
                   mpsK, mpsI; Z=Z)
    end
    return collect(0:δt:T_FINAL), gt, Z
end

gt_reference(δt::Float64) =
    [im * toulouse_Gt(BATH, tj; ϵ_d=μ_SYS) for tj in collect(0:δt:T_FINAL)]

"排除 t=0（端点约定差异，见文件头 b/a 说明）后的最大相对误差"
function err_gt(ts, gt, ref)
    m = ts .> 0
    maximum(abs.(gt[m] .- ref[m])) / maximum(abs.(ref[m]))
end

fit_order(hs, errs) = begin
    x = log.(hs); y = log.(errs)
    n = length(hs)
    (n * dot(x, y) - sum(x) * sum(y)) / (n * dot(x, x) - sum(x)^2)
end

# ------------------------------------------------------------------ 主流程 --
function main()
    dts = (0.1, 0.05, 0.025)
    csv = joinpath(OUTDIR, "convergence.csv")
    io = open(csv, "w")
    println(io, "scheme,kappa,dt,N,t,gt,ref,abs_err")

    schemes = [("κ=1 (current)", 1.0), ("κ=1/2", 0.5), ("κ=0 (dropped)", 0.0)]
    summary = Tuple{String,Float64,Float64,Float64}[]   # (name, κ, order_gt, order_Z)

    for (name, κ) in schemes
        errs_gt = Tuple{Float64,Float64}[]
        errs_z = Tuple{Float64,Float64}[]
        println("="^92)
        println("scheme = $name")
        for δt in dts
            N = round(Int, T_FINAL / δt)
            lat = GrassmannLattice(N=N, δt=δt, contour=:real)
            corr = scale_equaltime(correlationfunction(BATH, lat), κ)
            ts, gt, Z = gt_scheme(δt, corr)
            ref = gt_reference(δt)
            e = err_gt(ts, gt, ref)
            push!(errs_gt, (δt, e))
            push!(errs_z, (δt, abs(Z - 1)))
            for (tj, g, r) in zip(ts, gt, ref)
                println(io, name, ",", κ, ",", δt, ",", N, ",", tj, ",", g, ",", r, ",", abs(g - r))
            end
            @printf("  δt=%.5f  N=%3d  err(gt,t>0) = %.6e   |Z-1| = %.6e\n",
                    δt, N, e, abs(Z - 1))
        end
        o_gt = fit_order(first.(errs_gt), last.(errs_gt))
        o_z = fit_order(first.(errs_z), last.(errs_z))
        @printf("  >>> order(gt, t>0) = %.3f    order(|Z-1|) = %.3f\n", o_gt, o_z)
        push!(summary, (name, κ, o_gt, o_z))
        println(io, name, ",", κ, ",order,,,,order_gt,", o_gt)
        println(io, name, ",", κ, ",order,,,,order_Z,", o_z)
    end
    close(io)

    svg = joinpath(OUTDIR, "convergence.svg")
    _svg(summary, dts)
    println("="^92)
    println("数据 -> ", csv)
    println("图   -> ", svg)
    return summary
end

# ------------------------------------------------------- 轻量 SVG（log-log）--
_xesc(s) = replace(replace(replace(s, '&' => "&amp;"), '<' => "&lt;"), '>' => "&gt;")

function _svg(summary, dts)
    svg = joinpath(OUTDIR, "convergence.svg")
    width, height = 960, 640
    ml, mr, mt, mb = 90, 40, 70, 110
    pw, ph = width - ml - mr, height - mt - mb
    colors = Dict(1.0 => "#d62728", 0.5 => "#2ca02c", 0.0 => "#1f77b4")
    dashes = Dict(1.0 => "", 0.5 => "7,4", 0.0 => "2,3")
    labels = Dict(1.0 => "κ=1 (current full triangles)",
                  0.5 => "κ=1/2 (halved)",
                  0.0 => "κ=0 (dropped)")

    # 从 CSV 重读（简单起见直接重算坐标范围用固定框）
    elo, ehi, hlo, hhi = 1e-4, 1.0, 0.020, 0.12
    lx0, lx1 = log10(hlo), log10(hhi)
    ly0, ly1 = log10(elo), log10(ehi)
    px(x) = ml + pw * (x - lx0) / (lx1 - lx0)
    py(y) = mt + ph * (1 - (y - ly0) / (ly1 - ly0))

    io = IOBuffer()
    print(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"$width\" height=\"$height\" ",
              "viewBox=\"0 0 $width $height\" font-family=\"monospace\">")
    print(io, "<rect width=\"$width\" height=\"$height\" fill=\"white\"/>")
    print(io, "<text x=\"$(width/2)\" y=\"30\" text-anchor=\"middle\" font-size=\"17\">",
              "Real-time convergence (gt, t&gt;0), Toulouse T=$T_FINAL β=$β_BATH</text>")
    print(io, "<rect x=\"$ml\" y=\"$mt\" width=\"$pw\" height=\"$ph\" fill=\"none\" stroke=\"black\"/>")
    for e in (1e-4, 1e-3, 1e-2, 1e-1)
        print(io, "<line x1=\"$ml\" y1=\"$(py(log10(e)))\" x2=\"$(ml+pw)\" y2=\"$(py(log10(e)))\" stroke=\"#e6e6e6\"/>")
        print(io, "<text x=\"$(ml-8)\" y=\"$(py(log10(e))+4)\" text-anchor=\"end\" font-size=\"11\">1e",
                  Int(round(log10(e))), "</text>")
    end
    for h in dts
        print(io, "<text x=\"$(px(log10(h)))\" y=\"$(mt+ph+20)\" text-anchor=\"middle\" font-size=\"11\">", h, "</text>")
    end
    print(io, "<text x=\"$(ml+pw/2)\" y=\"$(mt+ph+44)\" text-anchor=\"middle\" font-size=\"13\">δt</text>")
    print(io, "<text x=\"24\" y=\"$(mt+ph/2)\" text-anchor=\"middle\" font-size=\"13\" ",
              "transform=\"rotate(-90 24 $(mt+ph/2))\">max relative error (t&gt;0)</text>")
    for (name, κ, o_gt, _) in summary
        c = colors[κ]; da = dashes[κ]
        da_attr = isempty(da) ? "" : " stroke-dasharray=\"$da\""
        # 重新计算各 δt 的误差：从 CSV 读
        errs = Tuple{Float64,Float64}[]
        csv = joinpath(OUTDIR, "convergence.csv")
        for line in eachline(csv)
            parts = split(line, ',')
            length(parts) >= 8 || continue
            parts[1] == name && parts[3] != "order" || continue
            δt = tryparse(Float64, parts[3])
            δt === nothing && continue
            err = tryparse(Float64, parts[8])
            err === nothing && continue
            push!(errs, (δt, err))
        end
        sort!(errs)
        isempty(errs) && continue
        pts = join(("$(px(log10(h))),$(py(log10(e)))" for (h, e) in errs), " ")
        print(io, "<polyline points=\"$pts\" fill=\"none\" stroke=\"$c\" stroke-width=\"2\"$da_attr/>")
        for (h, e) in errs
            print(io, "<circle cx=\"$(px(log10(h)))\" cy=\"$(py(log10(e)))\" r=\"4\" fill=\"$c\"/>")
        end
        (h1, e1) = errs[1]
        h2 = h1 / 4
        e2 = e1 * (h2 / h1)^o_gt
        print(io, "<line x1=\"$(px(log10(h1)))\" y1=\"$(py(log10(e1)))\" ",
                  "x2=\"$(px(log10(h2)))\" y2=\"$(py(log10(e2)))\" stroke=\"$c\" stroke-dasharray=\"5,4\" opacity=\"0.6\"/>")
        print(io, "<text x=\"$(px(log10(h2))-6)\" y=\"$(py(log10(e2))-10)\" text-anchor=\"end\" ",
                  "font-size=\"13\" fill=\"$c\">order = $(@sprintf("%.2f", o_gt))</text>")
        idx = κ == 1.0 ? 0 : (κ == 0.5 ? 1 : 2)
        ly = mt + 18 + 22 * idx
        print(io, "<line x1=\"$(ml+pw-360)\" y1=\"$ly\" x2=\"$(ml+pw-326)\" y2=\"$ly\" stroke=\"$c\" stroke-width=\"2\"$da_attr/>")
        print(io, "<text x=\"$(ml+pw-318)\" y=\"$(ly+4)\" font-size=\"12\">", _xesc(labels[κ]),
                  "  (order($(Int(κ*100))%) gt=$( @sprintf("%.2f", o_gt)))", "</text>")
    end
    print(io, "</svg>")
    write(svg, String(take!(io)))
    return svg
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
