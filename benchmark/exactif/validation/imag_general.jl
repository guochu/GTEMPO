# =============================================================================
# imag_general.jl —— 虚轴严格 IF 通用管线（参数化 N/order/model）
#
# 用法（分进程避内存击杀，行数据落盘 CSV）：
#   julia ... imag_general.jl <N> h <part 1|2> [order]   # h 行网格（18 εd 分两半）
#   julia ... imag_general.jl <N> g <part 1|2> [order]   # g_raw 行网格（8 εd 分两半）
#   julia ... imag_general.jl <N> closure [order]        # 载入行 → I* → 端到端 Z/n
#
# 配方（N=2 已机器精度闭环）：
#   K 链 = vacuum + 门 exp(GTerm(slice pair, c=e^{−τε_g}−1)) @ 切片 g+1
#         + 杂质分辨率 exp(GTerm(bar_{p+1}, eta_{p+2 环闭合}, −1))
#         + BC exp(GTerm(bar_1, eta_k, −1))
#   I* = 联合 Z+n 行（单项式基，h/g 均 I 侧 Möbius；g 经 _mult_A 直测 raw）
#        min-norm 解 → 序贯消元 exp(GTerm) 因子实现
#
# 运行目录：benchmark/exactif/validation/
# =============================================================================
using Printf
using LinearAlgebra
using Combinatorics: combinations
using GTEMPO
using QuAPI: AbstractCorrelationFunction
using Z2Tensors: NoTruncation
using Z2Tensors: Z2Irrep, blockdim, codomain, domain, fusiontrees, TruncateDimCutoff

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(GT_ROOT, "benchmark", "exactif", "validation", "formb_dense.jl"))
include(joinpath(GT_ROOT, "benchmark", "exactif", "src", "NativeExact.jl"))
using .NativeExact

# ------------------------------------------------------------------ 模型 --
# env 可切换：IMAG_WB / IMAG_AB / IMAG_MODELTAG（文件名标签）
const BETA = 1.0
const WB = parse(Float64, get(ENV, "IMAG_WB", "1.0"))
const AB = parse(Float64, get(ENV, "IMAG_AB", "0.5"))
const MODELTAG = get(ENV, "IMAG_MODELTAG", "a")
HB = ComplexF64[WB;;]
VB = ComplexF64[sqrt(AB) + 0.0im]
hcoupled() = NativeExact._coupled_h(VB, HB)

annih_op_f(dim::Int, mode::Int) = begin
    op = zeros(ComplexF64, dim, dim)
    for n in 0:(dim - 1)
        (n >> mode) & 1 == 1 || continue
        sgn = isodd(count_ones(n & ((1 << mode) - 1))) ? -1 : 1
        op[(n ⊻ (1 << mode)) + 1, n + 1] = sgn
    end
    op
end
function fock_anchor(cells, gate_durs, eps_d)
    m = length(VB); nt = 1 + m; dim = 1 << nt
    ops = [annih_op_f(dim, k) for k in 0:(nt-1)]
    Hq = zeros(ComplexF64, dim, dim)
    h = hcoupled()
    for i in 1:nt, j in 1:nt
        iszero(h[i, j]) || (Hq .+= h[i, j] .* (ops[i]' * ops[j]))
    end
    Him = eps_d .* (ops[1]' * ops[1])
    Z = Matrix{ComplexF64}(I, dim, dim); Zn = Matrix{ComplexF64}(I, dim, dim)
    for p in length(cells):-1:1
        W = exp(-cells[p] .* Hq) * exp(-gate_durs[p] .* Him)
        Z = W * Z; Zn = W * Zn
    end
    Zn = ops[1]' * ops[1] * Zn
    return tr(Z), tr(Zn) / tr(Z)
end

# ------------------------------------------------------------ 格点几何 --
function lat_geom(NL::Int)
    LAT = GTEMPO.GrassmannLattice(N=NL, δτ=BETA / NL, contour=:imag, bands=1)
    L = length(LAT)
    EVENS = Vector{Vector{Int}}([Int[]])
    for k in 2:2:L, comb in combinations(0:(L-1), k)
        push!(EVENS, collect(comb))
    end
    CFGIDX = [1 + sum(1 << p for p in A; init=0) for A in EVENS]
    return LAT, L, EVENS, CFGIDX
end

pairs_of(A) = Tuple{Int,Int}[(A[2i-1], A[2i]) for i in 1:(length(A) ÷ 2)]
function mk_all_subsets(L::Int)
    all_subsets = function(x)
        isempty(x) ? [Tuple{Int,Int}[]] :
            [x[[m == 1 for m in msk]] for msk in Iterators.product(fill(0:1, length(x))...)]
    end
    return all_subsets
end

function mk_injected(LAT)
    injected = function(pairs::Vector{Tuple{Int,Int}})
        mps = GTEMPO.vacuumstate(ComplexF64, LAT)
        for (p1, p2) in sort(pairs)
            mps = GTEMPO.apply!(exp(GTEMPO.GTerm(p1 + 1, p2 + 1; coeff=ComplexF64(1.0))), mps)
        end
        return mps
    end
    return injected
end

"""pipeline-v2 K 链（门 c=e^{−τε}−1 @切片 g=区间起点，含缝位合并门；分辨率 s=−1；BC）。"""
function mk_kchain(LAT, nl::Int, order::Int, ring_links0::Vector{Float64})
    kchain = function(eps_d)
        # dense 环序交错：G₁ R₁ G₂ R₂ … G_nl R_nl，BC 最后
        mps = GTEMPO.vacuumstate(ComplexF64, LAT)
        for p in 1:nl
            iszero(eps_d) || begin
                sg = p + 1
                gb = LAT[GrassmannContourIndex(sg, :τ, 1, true)]
                ge = LAT[GrassmannContourIndex(sg, :τ, 1, false)]
                mps = GTEMPO.apply!(exp(GTEMPO.GTerm(gb, ge;
                                     coeff=ComplexF64(exp(-ring_links0[p] * eps_d) - 1.0))), mps)
            end
            s1 = p + 1
            s2 = p == nl ? 2 : p + 2
            rb = LAT[GrassmannContourIndex(s1, :τ, 1, true)]
            re = LAT[GrassmannContourIndex(s2, :τ, 1, false)]
            mps = GTEMPO.apply!(exp(GTEMPO.GTerm(rb, re; coeff=ComplexF64(-1.0))), mps)
        end
        p1 = LAT[GrassmannContourIndex(1, :τ, 1, true)]
        p2 = LAT[GrassmannContourIndex(LAT.k, :τ, 1, false)]
        cbc = parse(Float64, get(ENV, "IMAG_CBC", "-1.0"))
        mps = GTEMPO.apply!(exp(GTEMPO.GTerm(p1, p2; coeff=ComplexF64(cbc))), mps)
        return mps, ring_links0
    end
    return kchain
end

"""Mobius(I 侧)：pair-乘积注入值 → 单项式系数。"""
function mk_mobiusI(EVENS, pairs_of, all_subsets, idxof)
    mobiusI = function(raw::Vector{ComplexF64})
        out = zeros(ComplexF64, length(EVENS))
        for (ib, B) in enumerate(EVENS)
            val = 0.0 + 0.0im
            for Bsub in all_subsets(pairs_of(B))
                sb = sort(collect(Iterators.flatten(Bsub)))
                sb = isempty(sb) ? Int[] : collect(sb)
                sgn = (-1.0)^((length(B) - length(sb)) ÷ 2)
                val += sgn * raw[idxof(sb)]
            end
            out[ib] = val
        end
        return out
    end
    return mobiusI
end

# ================================================================== 主程序 --
function main()
    NL = parse(Int, ARGS[1])
    mode = ARGS[2]
    order = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 2
    part = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 0
    tag = "n$(NL)o$(order)m$(MODELTAG)"
    outdir = joinpath(GT_ROOT, "benchmark", "exactif", "validation")

    # 先解调度：环区间数 nl（order 4 = 3N ≠ N！lattice/门循环全用 nl）
    # env 覆盖：IMAG_CELLS="c1,c2,..." + IMAG_RLINKS="g1,g2,..."（自定义调度隔离实验）
    if haskey(ENV, "IMAG_CELLS") && haskey(ENV, "IMAG_RLINKS")
        cells_env = [parse(Float64, x) for x in split(ENV["IMAG_CELLS"], ",")]
        rlinks_env = [parse(Float64, x) for x in split(ENV["IMAG_RLINKS"], ",")]
        sched0 = (cells = cells_env, links = rlinks_env)
        ring_links0 = rlinks_env
        nl = length(rlinks_env)
    else
        sched0 = NativeExact.native_schedule(order, NL, BETA / NL)
        ring_links0 = vcat(sched0.links[1] + sched0.links[end], sched0.links[2:end-1])
        nl = length(ring_links0)
    end
    LAT, L, EVENS, CFGIDX = lat_geom(nl)
    pairs_of_l = pairs_of
    all_subsets = mk_all_subsets(L)
    idxof(A) = findfirst(==(A), EVENS)
    injected = mk_injected(LAT)
    kchain = mk_kchain(LAT, nl, order, ring_links0)
    mobiusI = mk_mobiusI(EVENS, pairs_of, all_subsets, idxof)

    hgrid = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0, 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875]
    ggrid = [0.0, 0.125, 0.25, 0.5, 0.75, 1.0, 1.5, 1.875]
    cells = fill(BETA / NL, NL)

    if mode == "h"
        eds = hgrid[part == 1 ? (1:8) : (9:16)]
        file = joinpath(outdir, "imag_$(tag)_h$part.csv")
        open(file, "w") do io
            println(io, "ed,ib,re,im")
            for ed in eds
                mpsK, _ = kchain(ed)
                gB = [GTEMPO.integrate(LAT, mpsK, injected(pairs_of(B))) for B in EVENS]
                h = mobiusI(gB)
                for (ib, v) in enumerate(h)
                    println(io, "$ed,$ib,$(real(v)),$(imag(v))")
                end
                GC.gc(); flush(stdout)
            end
        end
        println("$file 写入完成")
        return
    end

    if mode == "g"
        eds = ggrid[part == 1 ? (1:4) : (5:8)]
        file = joinpath(outdir, "imag_$(tag)_g$part.csv")
        b1p = LAT[GrassmannContourIndex(1, conj=true, branch=:τ, band=1)]
        open(file, "w") do io
            println(io, "ed,ib,z_re,z_im,ge_re,ge_im,g1_re,g1_im")
            for ed in eds
                mpsK, _ = kchain(ed)
                t_end = GTEMPO.GTerm(LAT[GrassmannContourIndex(LAT.Nτ, conj=false, branch=:τ, band=1)], b1p; coeff=1)
                t_1 = GTEMPO.GTerm(LAT[GrassmannContourIndex(1, conj=false, branch=:τ, band=1)], b1p; coeff=1)
                mpsK_end = GTEMPO._mult_A(t_end, mpsK)
                mpsK_1 = GTEMPO._mult_A(t_1, mpsK)
                rawz = [GTEMPO.integrate(LAT, mpsK, injected(pairs_of(B))) for B in EVENS]
                rawe = [GTEMPO.integrate(LAT, mpsK_end, injected(pairs_of(B))) for B in EVENS]
                raw1 = [GTEMPO.integrate(LAT, mpsK_1, injected(pairs_of(B))) for B in EVENS]
                mz = mobiusI(rawz); me = mobiusI(rawe); m1 = mobiusI(raw1)
                for ib in 1:length(EVENS)
                    println(io, "$ed,$ib,$(real(mz[ib])),$(imag(mz[ib])),$(real(me[ib])),$(imag(me[ib])),$(real(m1[ib])),$(imag(m1[ib]))")
                end
                GC.gc(); flush(stdout)
            end
        end
        println("$file 写入完成")
        return
    end

    if mode == "closure"
        cells = sched0.cells
        length(cells) == length(ring_links0) || error("调度表长度不配")
        # 载入行
        hrows = Dict{Float64, Vector{ComplexF64}}()
        for part in (1, 2)
            file = joinpath(outdir, "imag_$(tag)_h$part.csv")
            for (i, ln) in enumerate(eachline(file))
                i == 1 && continue
                f = split(strip(ln), ",")
                ed = parse(Float64, f[1]); ib = parse(Int, f[2])
                get!(hrows, ed, zeros(ComplexF64, length(EVENS)))[ib] =
                    ComplexF64(parse(Float64, f[3]), parse(Float64, f[4]))
            end
        end
        grows = Dict{Float64, Vector{ComplexF64}}()
        for part in (1, 2)
            file = joinpath(outdir, "imag_$(tag)_g$part.csv")
            for (i, ln) in enumerate(eachline(file))
                i == 1 && continue
                f = split(strip(ln), ",")
                ed = parse(Float64, f[1]); ib = parse(Int, f[2])
                g = get!(grows, ed, zeros(ComplexF64, 3 * length(EVENS)))
                g[ib] = ComplexF64(parse(Float64, f[3]), parse(Float64, f[4]))
                g[length(EVENS) + ib] = ComplexF64(parse(Float64, f[5]), parse(Float64, f[6]))
                g[2length(EVENS) + ib] = ComplexF64(parse(Float64, f[7]), parse(Float64, f[8]))
            end
        end
        rows = Vector{Vector{ComplexF64}}(); rhs = ComplexF64[]
        for ed in sort(collect(keys(hrows)))
            push!(rows, hrows[ed]); push!(rhs, fock_anchor(cells, ring_links0, ed)[1])
        end
        for ed in sort(collect(keys(grows)))
            hz = hrows[ed]
            vend = grows[ed][length(EVENS) .+ (1:length(EVENS))]
            nanc = real(fock_anchor(cells, ring_links0, ed)[2])
            push!(rows, hz .- vend .- nanc .* hz); push!(rhs, 0.0 + 0.0im)
        end
        A = hcat(rows...)
        flush(stdout)
        Ist = A * ((A' * A + 1e-13 * I) \ rhs)
        res = norm(A' * Ist - rhs)
        @printf("NL=%d order=%d 联合残差 = %.3e\n", NL, order, res)
        flush(stdout)
        ν_eff = Ist[1]
        println("ν_eff = ", ν_eff); flush(stdout)
        ivec = zeros(ComplexF64, 1 << L)
        for (i, Aset) in enumerate(EVENS)
            ivec[1 + sum(1 << p for p in Aset; init=0)] = Ist[i] / ν_eff
        end
        # 序贯消元实现
        factors = Vector{Tuple{Vector{Int},ComplexF64}}()
        tcur = zeros(ComplexF64, 1 << L); tcur[1] = 1.0
        for Aset in sort(EVENS; by=length)
            isempty(Aset) && continue
            want = ivec[1 + sum(1 << p for p in Aset; init=0)]
            c = want - tcur[1 + sum(1 << p for p in Aset; init=0)]
            abs(c) < 1e-9 && continue
            push!(factors, (Aset, c))
            monos = Tuple{Vector{Int}, ComplexF64}[(Int[], 1.0 + 0.0im), (Aset, c)]
            tcur = mul_monomials!(zeros(ComplexF64, 1 << L), tcur, L, Monomials(monos))
        end
        err = maximum(abs.(tcur - ivec))
        @printf("因子数 = %d，重构误差 = %.3e\n", length(factors), err); flush(stdout)
        err < 1e-8 || error("实现失败")
        # 端到端
        println("端到端（真实 GMPS 链 + integrate + gf）："); flush(stdout)
        maxerr_z = 0.0; maxerr_n = 0.0
        t0 = time()
        # 因子 → 小 MPS，再 mult! 截断乘法链（键维受控）
        mpsI = GTEMPO.vacuumstate(ComplexF64, LAT)
        maxtr = 0.0
        for (fi, (Aset, c)) in enumerate(factors)
            fm = GTEMPO.vacuumstate(ComplexF64, LAT)
            fm = GTEMPO.apply!(exp(GTEMPO.GTerm([p + 1 for p in Aset]...; coeff=ComplexF64(c))), fm)
            mpsI, err = GTEMPO.mult!(mpsI, fm; trunc=TruncateDimCutoff(2048, 1e-13, 1))
            maxtr = max(maxtr, err)
            fi % 10 == 0 && (@printf("  factor %d/%d (%.1fs, err=%.1e)\n", fi, length(factors), time() - t0, maxtr); flush(stdout))
        end
        @printf("mpsI 构建完成（%.1fs，最大截断误差 %.2e）\n", time() - t0, maxtr); flush(stdout)
        for ed in [0.0, 0.35, 0.7, 1.05, 1.4]
            mpsK, _ = kchain(ed)
            z = GTEMPO.integrate(LAT, mpsK, mpsI) * ν_eff
            flush(stdout)
            b1c = LAT[GrassmannContourIndex(1, conj=true, branch=:τ, band=1)]
            tend = GTEMPO.GTerm(LAT[GrassmannContourIndex(LAT.Nτ, conj=false, branch=:τ, band=1)], b1c; coeff=1)
            graw = GTEMPO.integrate(LAT, GTEMPO._mult_A(tend, mpsK), mpsI)
            npipe = 1 - real(graw / (z / ν_eff))
            za, na = fock_anchor(cells, ring_links0, ed)
            ez = abs(real(z) - real(za)); en = abs(npipe - real(na))
            maxerr_z = max(maxerr_z, ez); maxerr_n = max(maxerr_n, en)
            @printf("εd=%.2f : Z = %12.8f（err %9.2e）  n = %9.6f（锚 %9.6f, err %9.2e）\n",
                    ed, real(z), ez, npipe, real(na), en)
        end
        @printf("NL=%d order=%d 最大误差：Z %.2e / n %.2e  %s\n", NL, order, maxerr_z, maxerr_n,
                (maxerr_z < 1e-10 && maxerr_n < 1e-10) ? "✅ PASS" : "✗ FAIL")
    end
end
main()
