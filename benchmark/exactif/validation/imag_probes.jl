# =============================================================================
# imag_probes.jl —— Plan D：I 侧探针族扩充（随机符号注入 + mult! 乘积）
#
# 原理：复合对 I-MPS 张量线性（jwcalib 确认）。注入族 {inj(B)} 的复合张成
#   秩 ~25，Z 不可达（差 0.2）。扩族：
#   (a) 随机符号模式注入 inj±(B, σ)——span 外的新张量方向
#   (b) mult! 乘积探针 mult(inj(B1), inj(B2))
# 行 = composite(K(εd), probe_j)，落盘 CSV。拟合在 imag_dfit.jl 做。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_probes.jl <s|p> <part>
#   s = 符号注入族，p = mult! 乘积族；part = εd 分块（1/2）
# =============================================================================
using Printf
using LinearAlgebra
using Random
using Combinatorics: combinations
using GTEMPO
using Z2Tensors: NoTruncation, TruncateDimCutoff

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(GT_ROOT, "benchmark", "exactif", "src", "NativeExact.jl"))
using .NativeExact

const BETA = 1.0
const WB, AB = 1.0, 0.5
HB = ComplexF64[WB;;]
VB = ComplexF64[sqrt(AB) + 0.0im]

function lat_geom(nl::Int)
    LAT = GTEMPO.GrassmannLattice(N=nl, δτ=BETA / nl, contour=:imag, bands=1)
    L = length(LAT)
    EVENS = Vector{Vector{Int}}([Int[]])
    for k in 2:2:L, comb in combinations(0:(L-1), k)
        push!(EVENS, collect(comb))
    end
    return LAT, L, EVENS
end

function mk_kchain(LAT, ring_links0::Vector{Float64}, eps_d::Float64)
    nl = length(ring_links0)
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
    mps = GTEMPO.apply!(exp(GTEMPO.GTerm(p1, p2; coeff=ComplexF64(-1.0))), mps)
    return mps
end

function main()
    family = ARGS[1]; part = parse(Int, ARGS[2])
    TAG = get(ENV, "IMAG_TAG", "n2o4ma")
    NLm = parse(Int, get(ENV, "IMAG_NL", "2"))
    channel = length(ARGS) >= 3 ? ARGS[3] : "z"      # z = 普通复合；g = 末位插入通道
    outdir = joinpath(GT_ROOT, "benchmark", "exactif", "validation")
    NL, order = NLm, 4
    sched = NativeExact.native_schedule(order, NL, BETA / NL)
    ring_links = vcat(sched.links[1] + sched.links[end], sched.links[2:end-1])
    nl = length(ring_links)
    LAT, L, EVENS = lat_geom(nl)
    NE = length(EVENS)
    rng = MersenneTwister(20261109)

    # 锚点 εd 集 = h 行的 εd 全集（41 点；holdout 4 点也在内）
    eds = Float64[]
    for file in sort(filter(f -> occursin("imag_$(TAG)_h", f) && endswith(f, ".csv"), readdir(outdir, join=true)))
        for (i, ln) in enumerate(eachline(file))
            i == 1 && continue
            push!(eds, parse(Float64, split(strip(ln), ",")[1]))
        end
    end
    eds = sort(unique(eds))
    half = part == 1 ? (1:div(length(eds), 2)) : (div(length(eds), 2)+1:length(eds))
    eds_part = eds[half]
    println("εd 点（part $part）：$(length(eds_part))"); flush(stdout)

    # 探针族定义
    probes = Vector{Any}()
    if family == "c"
        # (c) 随机复系数注入（幅度+相位随机，张量采样更丰富）
        for B in EVENS
            np2 = length(B) ÷ 2
            σ = [ComplexF64(randn(rng), randn(rng)) for _ in 1:np2]
            push!(probes, (B, σ, :sign))
        end
    elseif family == "i"
        # (i) 普通注入族（coeff +1，与 h 行同基）
        for B in EVENS
            push!(probes, (B, ones(Float64, length(B) ÷ 2), :sign))
        end
    elseif family == "s"
        # (a) 8192 随机符号模式注入：σ ∈ {±1}^{npairs(B)} 每B一个随机模式
        for (ib, B) in enumerate(EVENS)
            np2 = length(B) ÷ 2
            σ = [rand(rng) < 0.5 ? -1.0 : 1.0 for _ in 1:np2]
            push!(probes, (B, σ, :sign))
        end
    elseif family == "p"
        # (b) 8192 mult! 乘积：随机 (B1, B2) 对
        for _ in 1:NE
            B1 = EVENS[rand(rng, 1:NE)]; B2 = EVENS[rand(rng, 1:NE)]
            push!(probes, (B1, B2, :mult))
        end
    else
        error("family ∈ {s, p}")
    end
    println("探针数：$(length(probes))"); flush(stdout)

    inj = function(A::Vector{Int}, σ)
        mps = GTEMPO.vacuumstate(ComplexF64, LAT)
        for (i, c) in zip(1:2:length(A), σ)
            mps = GTEMPO.apply!(exp(GTEMPO.GTerm(A[i] + 1, A[i+1] + 1; coeff=ComplexF64(c))), mps)
        end
        return mps
    end

    parttag = part
    file = joinpath(outdir, "imag_$(TAG)_$(family)$(parttag)$(channel == "g" ? "g" : "").csv")
    open(file, "w") do io
        println(io, "ed,jp,re,im")
        for (ipt, ed) in enumerate(eds_part)
            mpsK = mk_kchain(LAT, ring_links, ed)
            if channel == "g"
                b1c = LAT[GrassmannContourIndex(1, conj=true, branch=:τ, band=1)]
                tend = GTEMPO.GTerm(LAT[GrassmannContourIndex(LAT.Nτ, conj=false, branch=:τ, band=1)], b1c; coeff=1)
                mpsK = GTEMPO._mult_A(tend, mpsK)
            end
            t0 = time()
            for (jp, pr) in enumerate(probes)
                if pr[3] === :sign
                    mpsI = inj(pr[1], pr[2])
                else
                    m1 = inj(pr[1], ones(Float64, length(pr[1]) ÷ 2))
                    m2 = inj(pr[2], ones(Float64, length(pr[2]) ÷ 2))
                    mpsI = GTEMPO.mult!(m1, m2; trunc=TruncateDimCutoff(256, 1e-13, 1))[1]
                end
                v = GTEMPO.integrate(LAT, mpsK, mpsI)
                println(io, ed, ",", jp, ",", real(v), ",", imag(v))
            end
            println("  εd=$ed done ($(length(probes)) probes, $(round(time()-t0, digits=1))s)"); flush(stdout)
        end
    end
    println("$file 写入完成")
end
main()
