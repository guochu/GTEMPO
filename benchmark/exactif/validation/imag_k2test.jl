# =============================================================================
# imag_k2test.jl —— pipeline-v2 K 链（自建 dense-约定门+分辨率+BC）联合一致性判决
#
# K 链全部 GTerm 注入（不用库门张量）：门 c_g × 分辨率 s × BC(+1)。
# combo（ARGS）：
#   1  c_g = −e^{−τε}   s = −1     （imag_dense kernel 镜像候选）
#   2  c_g = +e^{−τε}   s = −1
#   3  c_g = e^{−τε}−1  s = −1
#   4  c_g = 1−e^{−τε}  s = +1
#   5  c_g = −e^{−τε}   s = +1
#   6  c_g = e^{−τε}−1  s = +1
# 判决：Z-行（rank 3）+ n-行（1−g_end 变体）的 26 行联合系统残差。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_k2test.jl <combo>
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
const WB, AB = 1.0, 0.5
HB = ComplexF64[WB;;]
VB = ComplexF64[sqrt(AB) + 0.0im]

function fock_anchor(cells, gate_durs, eps_d)
    m = length(VB); nt = 1 + m; dim = 1 << nt
    ops = [annih_op_f(dim, k) for k in 0:(nt-1)]
    Hq = zeros(ComplexF64, dim, dim)
    h = NativeExact._coupled_h(VB, HB)
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
annih_op_f(dim::Int, mode::Int) = begin
    op = zeros(ComplexF64, dim, dim)
    for n in 0:(dim - 1)
        (n >> mode) & 1 == 1 || continue
        sgn = isodd(count_ones(n & ((1 << mode) - 1))) ? -1 : 1
        op[(n ⊻ (1 << mode)) + 1, n + 1] = sgn
    end
    op
end

const NL = 2
const L = 2 * NL + 2
const LAT = GTEMPO.GrassmannLattice(N=NL, δτ=BETA / NL, contour=:imag, bands=1)
const EVENS = Vector{Vector{Int}}([Int[]])
for k in (2, 4, 6), comb in combinations(0:(L-1), k)
    push!(EVENS, collect(comb))
end
const IDXOF = Dict(A => i for (i, A) in enumerate(EVENS))
pairs_of(A) = Tuple{Int,Int}[(A[2i-1], A[2i]) for i in 1:(length(A) ÷ 2)]
all_subsets(x) = isempty(x) ? [Tuple{Int,Int}[]] :
    [x[[m == 1 for m in msk]] for msk in Iterators.product(fill(0:1, length(x))...)]
idxof(A) = findfirst(==(A), EVENS)

function injected_mps(pairs; coeffs=nothing)
    mps = GTEMPO.vacuumstate(ComplexF64, LAT)
    for (i, (p1, p2)) in enumerate(sort(pairs))
        c = coeffs === nothing ? 1.0 : coeffs[i]
        mps = GTEMPO.apply!(exp(GTEMPO.GTerm(p1 + 1, p2 + 1; coeff=ComplexF64(c))), mps)
    end
    return mps
end

"""pipeline-v2 K 链（自建）。"""
function kchain_v2(eps_d; c_gate::Function, s_res::Float64)
    sched = NativeExact.native_schedule(2, N, BETA / N)
    ring_links = vcat(sched.links[1] + sched.links[end], sched.links[2:end-1])
    mps = GTEMPO.vacuumstate(ComplexF64, LAT)
    # 门（切片 g+1）
    for g in 1:NL
        iszero(eps_d) && continue
        s = g + 1
        pb = LAT[GrassmannContourIndex(s, :τ, 1, true)]
        pe = LAT[GrassmannContourIndex(s, :τ, 1, false)]
        c = c_gate(ring_links[g] * eps_d)
        mps = GTEMPO.apply!(exp(GTEMPO.GTerm(pb, pe; coeff=ComplexF64(c))), mps)
    end
    # 杂质分辨率（bar_{p+1}, eta_{p+2}，环闭合）
    for p in 1:N
        s1 = p + 1
        s2 = p == N ? 2 : p + 2
        pb = LAT[GrassmannContourIndex(s1, :τ, 1, true)]
        pe = LAT[GrassmannContourIndex(s2, :τ, 1, false)]
        mps = GTEMPO.apply!(exp(GTEMPO.GTerm(pb, pe; coeff=ComplexF64(s_res))), mps)
    end
    # BC
    p1 = LAT[GrassmannContourIndex(1, :τ, 1, true)]
    p2 = LAT[GrassmannContourIndex(LAT.k, :τ, 1, false)]
    mps = GTEMPO.apply!(exp(GTEMPO.GTerm(p1, p2; coeff=ComplexF64(-1.0))), mps)
    return mps
end

function main()
    combo = parse(Int, ARGS[1])
    gate_c, s_res = [(τ -> -exp(-τ), -1.0),
                     (τ -> +exp(-τ), -1.0),
                     (τ -> exp(-τ) - 1.0, -1.0),
                     (τ -> 1.0 - exp(-τ), +1.0),
                     (τ -> -exp(-τ), +1.0),
                     (τ -> exp(-τ) - 1.0, +1.0)][combo]
    println("combo $combo: c_gate(τε) 形式 = $(gate_c(0.245))（τε=0.245）, s_res = $s_res")

    hgrid = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0, 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875]
    ggrid = [0.0, 0.125, 0.25, 0.5, 0.75, 1.0, 1.5, 1.875]

    # h 网格
    hrows = Vector{Vector{ComplexF64}}()
    for ed in hgrid
        mpsK = kchain_v2(ed; c_gate=gate_c, s_res=s_res)
        gB = [GTEMPO.integrate(LAT, mpsK, injected_mps(pairs_of(B))) for B in EVENS]
        h = zeros(ComplexF64, length(EVENS))
        for (ib, B) in enumerate(EVENS)
            val = 0.0 + 0.0im
            for Bsub in all_subsets(pairs_of(B))
                sb = sort(collect(Iterators.flatten(Bsub)))
                sb = isempty(sb) ? Int[] : collect(sb)
                sgn = (-1.0)^((length(B) - length(sb)) ÷ 2)
                val += sgn * gB[idxof(sb)]
            end
            h[ib] = val
        end
        push!(hrows, h)
        GC.gc(); flush(stdout)
    end
    # g 网格（z, g_end, g_1）
    grows = Vector{Vector{ComplexF64}}()
    for ed in ggrid
        mpsK = kchain_v2(ed; c_gate=gate_c, s_res=s_res)
        row = zeros(ComplexF64, 3 * length(EVENS))
        # g_raw 经 _mult_A 直接测（GTerm 插在 K 链上；无归一化无守卫）
        b1 = LAT[GrassmannContourIndex(1, conj=true, branch=:τ, band=1)]
        t_end = GTEMPO.GTerm(LAT[GrassmannContourIndex(LAT.Nτ, conj=false, branch=:τ, band=1)], b1; coeff=1)
        t_1 = GTEMPO.GTerm(LAT[GrassmannContourIndex(1, conj=false, branch=:τ, band=1)], b1; coeff=1)
        mpsK_end = GTEMPO._mult_A(t_end, mpsK)
        mpsK_1 = GTEMPO._mult_A(t_1, mpsK)
        for (ib, B) in enumerate(EVENS)
            mpsI = injected_mps(pairs_of(B))
            row[ib] = hrows[findfirst(==(ed), hgrid)][ib]   # 对齐 h
            row[length(EVENS) + ib] = GTEMPO.integrate(LAT, mpsK_end, mpsI)
            row[2length(EVENS) + ib] = GTEMPO.integrate(LAT, mpsK_1, mpsI)
        end
        # I 侧 Möbius：pair-乘积注入值 → 纯单项式系数（与 h-行同基）
        gmono = zeros(ComplexF64, length(EVENS))
        for (ib, B) in enumerate(EVENS)
            val = 0.0 + 0.0im
            for ch in (:z, :ge, :g1)
                col = ch === :z ? 0 : ch === :ge ? length(EVENS) : 2length(EVENS)
                v = 0.0 + 0.0im
                for Bsub in all_subsets(pairs_of(B))
                    sb = sort(collect(Iterators.flatten(Bsub)))
                    sb = isempty(sb) ? Int[] : collect(sb)
                    sgn = (-1.0)^((length(B) - length(sb)) ÷ 2)
                    v += sgn * row[col + idxof(sb)]
                end
                # 写回时只对当前通道操作——改为分别处理
            end
        end
        # （逐通道 Möbius）
        for (col, off) in ((0, 0), (length(EVENS), length(EVENS)), (2length(EVENS), 2length(EVENS)))
            for (ib, B) in enumerate(EVENS)
                v = 0.0 + 0.0im
                for Bsub in all_subsets(pairs_of(B))
                    sb = sort(collect(Iterators.flatten(Bsub)))
                    sb = isempty(sb) ? Int[] : collect(sb)
                    sgn = (-1.0)^((length(B) - length(sb)) ÷ 2)
                    v += sgn * row[col + idxof(sb)]
                end
                row[off + ib] = v
            end
        end
        push!(grows, row)
        GC.gc(); flush(stdout)
    end

    # 联合一致性（n-行用 1−g_end 变体）
    rows = Vector{Vector{ComplexF64}}(); rhs = ComplexF64[]
    for (r, ed) in enumerate(hgrid)
        push!(rows, hrows[r]); push!(rhs, fock_anchor(fill(BETA / N, N), fill(BETA / N, N), ed)[1])
    end
    for (r, ed) in enumerate(ggrid)
        hz = hrows[findfirst(==(ed), hgrid)]
        ge = grows[r][length(EVENS) .+ (1:length(EVENS))]
        nanc = real(fock_anchor(fill(BETA / N, N), fill(BETA / N, N), ed)[2])
        push!(rows, hz .- ge .- nanc .* hz); push!(rhs, 0.0 + 0.0im)
    end
    A = hcat(rows...)
    sol = A * ((A' * A + 1e-12 * I) \ rhs)
    res = norm(A' * sol - rhs)
    @printf("combo %d 联合一致性残差 = %.3e%s\n", combo, res,
            res < 1e-8 ? "  ✅ 一致！" : "  ✗")
end
main()
