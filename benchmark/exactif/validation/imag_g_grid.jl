# =============================================================================
# imag_g_grid.jl —— g-泛函网格：g(b) = gf 读出（I 侧 = 注入单项式 P_b）
#
# n-约束行（线性化）：Σ_b I_b·g(b) − n_dense·Σ_b I_b·h(b) = 0。
# 输出 imag_g.csv（order,epsd,ib,variant,re,im）。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_g_grid.jl
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

const L = 6
evensets = Vector{Vector{Int}}([Int[]])
for k in (2, 4, 6), comb in combinations(0:(L-1), k)
    push!(evensets, collect(comb))
end
idxof(A) = findfirst(==(A), evensets)
pairs_of(A) = Tuple{Int,Int}[(A[2i-1], A[2i]) for i in 1:(length(A) ÷ 2)]

"""g 变体读出：mpsI = 注入 P_b。返回 (z, g_end, g_1)。"""
function g_readout(eps_d; order=2)
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
    out = Vector{Tuple{ComplexF64,ComplexF64,ComplexF64}}()
    for B in evensets
        mpsI = injected_mps(lat, pairs_of(B))
        z = GTEMPO.integrate(lat, mpsK, mpsI)
        b1 = GrassmannContourIndex(1, conj=true, branch=:τ, band=1)
        g_end = GTEMPO.gf(lat, (GrassmannContourIndex(lat.Nτ, conj=false, branch=:τ, band=1), b1), mpsK, mpsI, Z=z)
        g_1 = GTEMPO.gf(lat, (GrassmannContourIndex(1, conj=false, branch=:τ, band=1), b1), mpsK, mpsI, Z=z)
        push!(out, (z, g_end, g_1))
    end
    return out
end

function main()
    grid = [0.0, 0.125, 0.25, 0.5, 0.75, 1.0, 1.5, 1.875]
    outdir = joinpath(GT_ROOT, "benchmark", "exactif", "validation")
    open(joinpath(outdir, "imag_g.csv"), "w") do io
        println(io, "order,epsd,ib,z_re,z_im,ge_re,ge_im,g1_re,g1_im")
        for order in (2,)
            for ed in grid
                rows = g_readout(ed; order=order)
                for (ib, (z, ge, g1)) in enumerate(rows)
                    println(io, "$order,$ed,$ib,$(real(z)),$(imag(z)),$(real(ge)),$(imag(ge)),$(real(g1)),$(imag(g1))")
                end
                @printf("order=%d εd=%.3f done\n", order, ed)
                GC.gc(); flush(stdout)
            end
        end
    end
    println("imag_g.csv 已写入")
end
main()
