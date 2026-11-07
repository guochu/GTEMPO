# =============================================================================
# imag_h_grid.jl —— h 向量网格：h(b) = Z_real(真K链(εd,order), 单项式 b)
#
# 恒等式 Z(K, I) = Σ_b I_b·h(b)（τ2 双线性的 I 侧边际；I 侧 Möbius 读全）。
# 不需要 τ2。输出 imag_h.csv 供 imag_calib_solve.jl 拟合。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_h_grid.jl
# =============================================================================
using Printf
using LinearAlgebra
using Combinatorics: combinations
using GTEMPO
using QuAPI: AbstractCorrelationFunction
using Z2Tensors: NoTruncation

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(GT_ROOT, "benchmark", "exactif", "src", "NativeExact.jl"))
using .NativeExact

const BETA = 1.0
const N = 2
const NL = 2
const L = 2 * NL + 2

mutable struct ImagDelta <: AbstractCorrelationFunction
    Δ::Matrix{ComplexF64}
end
GTEMPO.TK.scalartype(::Type{<:ImagDelta}) = ComplexF64
Base.size(x::ImagDelta, args...) = size(x.Δ, args...)
Base.getindex(x::ImagDelta, i::GrassmannContourIndex, j::GrassmannContourIndex) = x.Δ[i.j, j.j]

function injected_mps(lat, pairs::Vector{Tuple{Int,Int}})
    mps = GTEMPO.vacuumstate(ComplexF64, lat)
    for (p1, p2) in sort(pairs)
        mps = GTEMPO.apply!(exp(GTEMPO.GTerm(p1 + 1, p2 + 1; coeff=ComplexF64(1.0))), mps)
    end
    return mps
end

evensets = Vector{Vector{Int}}([Int[]])
for k in (2, 4, 6), comb in combinations(0:(L-1), k)
    push!(evensets, collect(comb))
end
const NS = length(evensets)
idxof(A) = findfirst(==(A), evensets)
pairs_of(A) = Tuple{Int,Int}[(A[2i-1], A[2i]) for i in 1:(length(A) ÷ 2)]
all_subsets(x) = isempty(x) ? [Tuple{Int,Int}[]] :
    [x[[m == 1 for m in msk]] for msk in Iterators.product(fill(0:1, length(x))...)]

"""h(b) = Z_real(K_real, 单项式 b)，I 侧 Möbius 读全。"""
function h_readout(eps_d; order=2)
    sched = NativeExact.native_schedule(order, N, BETA / N)
    ring_links = vcat(sched.links[1] + sched.links[end], sched.links[2:end-1])
    nl = length(ring_links)
    lat = GTEMPO.GrassmannLattice(N=nl, δτ=BETA / nl, contour=:imag, bands=1)
    mpsK = GTEMPO.vacuumstate(ComplexF64, lat)
    model = ToulouseIM(ϵ_d=eps_d)
    for g in 1:nl
        fm = GTEMPO.fock_propagator(model, :τ, ring_links[g])
        sparse = GTEMPO._propagator_sparsegmps(lat, fm, g, :τ, nothing, nothing)
        mpsK = GTEMPO.mult!(mpsK, sparse, trunc=NoTruncation())[1]
    end
    mpsK = GTEMPO.boundarycondition!(mpsK, lat; band=1, trunc=NoTruncation())
    gB = [GTEMPO.integrate(lat, mpsK, injected_mps(lat, pairs_of(B))) for B in evensets]
    h = zeros(ComplexF64, NS)
    for (ib, B) in enumerate(evensets)
        val = 0.0 + 0.0im
        for Bsub in all_subsets(pairs_of(B))
            sb = sort(collect(Iterators.flatten(Bsub)))
            sb = isempty(sb) ? Int[] : collect(sb)
            sgn = (-1.0)^((length(B) - length(sb)) ÷ 2)
            val += sgn * gB[idxof(sb)]
        end
        h[ib] = val
    end
    return h
end

function main()
    grid = collect(0.0:0.125:1.875)         # 16 点：解 16×16 线性系统求 I*
    holdout = [0.525, 1.225]
    outdir = joinpath(GT_ROOT, "benchmark", "exactif", "validation")
    open(joinpath(outdir, "imag_h.csv"), "w") do io
        println(io, "order,epsd,ib,re,im")
        for order in (1, 2)
            for ed in vcat(grid, holdout)
                hv = h_readout(ed; order=order)
                for (ib, v) in enumerate(hv)
                    println(io, "$order,$ed,$ib,$(real(v)),$(imag(v))")
                end
                @printf("order=%d εd=%.3f done (|h|=%.4f)\n", order, ed, norm(hv))
                GC.gc(); flush(stdout)
            end
        end
    end
    println("imag_h.csv 已写入")
end
main()
