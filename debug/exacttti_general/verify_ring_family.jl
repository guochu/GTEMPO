# 探针 4：广义虚时环族候选公式的因子级验证
# 候选（任意 (ba,ca,bb,cb) 味对，M = 块矩阵，n = kτ，B = lattice.bands）：
#   lower[m] = M[m+1, 1]，fit f（y[d] = Σ w q^{d-1}）
#   ch1 = ESS(α₁=w./q, q;   conja=ca, conjb=cb, cell=B, asite=ba, bsite=bb)
#   ch2 = ESS(α₂=w.*q^{n-1}, 1./q; conja=cb, conjb=ca, cell=B, asite=bb, bsite=ba)
#   contact = M[1,1] − Σα₁ + Σα₂, ch3 = ESS([contact],[0]; ca, cb, asite=ba, bsite=bb)
# 目标 = fullfactor(单块核, lattice) 物化 bitwise。
using GTEMPO
using LinearAlgebra: norm
import GrassmannMPS as GM
using Random

"候选环族因子组（单味对）"
function ring_factors(M, lattice, ba, ca, bb, cb, alg)
    n = lattice.kτ
    B = lattice.bands
    lower = ComplexF64[M[j+1, 1] for j in 1:n-1]
    f = GTEMPO.fit_exponentials(alg, lower)
    isempty(f.weights) && return GM.GrassmannFactor[]
    α1 = f.weights ./ f.nodes
    α2 = f.weights .* f.nodes .^ (n - 1)
    sums = [GM.ExpDecaySum(α1, f.nodes; conja=ca, conjb=cb, cell=B, asite=ba, bsite=bb),
            GM.ExpDecaySum(α2, 1 ./ f.nodes; conja=cb, conjb=ca, cell=B, asite=bb, bsite=ba)]
    contact = M[1, 1] - sum(α1) + sum(α2)
    iszero(contact) || push!(sums, GM.ExpDecaySum([contact], [0.0];
        conja=ca, conjb=cb, cell=B, asite=ba, bsite=bb))
    return [GM.ExpDecaySumSum(sums)]
end

"物化因子组 → dense 向量"
function mat(factors, lat)
    spv = GTEMPO.site_spaces(lat)
    acc = GM.togmps(ComplexF64, first(factors), first(spv), length(spv))
    for f in factors[2:end]
        rhs = GM.togmps(ComplexF64, f, first(spv), length(spv))
        acc = GM.mult(acc, rhs)
    end
    return GM.Dense.todense(acc).data
end

"单块核：只含一个键的 HybridizationKernel（虚时）"
function single_block_kernel(M, ba, ca, bb, cb)
    blocks = Dict{Tuple{ContourBlockIndex,ContourBlockIndex},Union{Matrix{ComplexF64},CorrelationMatrix{ComplexF64}}}(
        (ContourBlockIndex(:τ, ba, ca), ContourBlockIndex(:τ, bb, cb)) => M)
    bands = sort(unique([ba, bb]))
    return HybridizationKernel(blocks, 0.0, bands)
end

"环 Toeplitz 热核（反周期闭合）：lower 给定，upper[m] = −lower[n−m]，对角 c"
function ring_kernel(lower, c)
    n = length(lower) + 1
    M = zeros(ComplexF64, n, n)
    for d in 1:n-1
        for k in 1:n-d
            M[k+d, k] = lower[d]
            M[k, k+d] = -lower[n-d]
        end
    end
    M[diagind(M)] .= c
    return M
end

alg = GTEMPO.Defaults.alg_expan()
using LinearAlgebra: diagind

gd(x, y) = maximum(abs.(x .- y))

# ---- ① 跨带 (T,F)：ba=1 → bb=2，环热核 ----
let n = 5
    lat = GTEMPO.GrassmannLattice(N=n, δτ=0.25, bands=2, contour=:imag, symmetry=Z2Symmetry())
    lower = ComplexF64[0.3 * 0.7^(m-1) * cis(0.4 * m) for m in 1:n-1]
    M = ring_kernel(lower, 0.05 - 0.02im)
    k = single_block_kernel(M, 1, true, 2, false)
    cand = ring_factors(M, lat, 1, true, 2, false, alg)
    tgt = GM.togmps(ComplexF64, GTEMPO.fullfactor(k, lat), first(GTEMPO.site_spaces(lat)), length(lat) ÷ 2)
    println("① 跨带 (T,F) 1→2: dev = ", gd(mat(cand, lat), GM.Dense.todense(tgt).data))
    # 反方向 ba=2 → bb=1
    k2 = single_block_kernel(M, 2, true, 1, false)
    cand2 = ring_factors(M, lat, 2, true, 1, false, alg)
    tgt2 = GM.togmps(ComplexF64, GTEMPO.fullfactor(k2, lat), first(GTEMPO.site_spaces(lat)), length(lat) ÷ 2)
    println("① 跨带 (T,F) 2→1: dev = ", gd(mat(cand2, lat), GM.Dense.todense(tgt2).data))
end

# ---- ② 同味 (T,T)：ud 反常（BCS 形）----
let n = 5
    lat = GTEMPO.GrassmannLattice(N=n, δτ=0.25, bands=2, contour=:imag, symmetry=Z2Symmetry())
    lower = ComplexF64[0.25 * 0.6^(m-1) for m in 1:n-1]
    M = ring_kernel(lower, 0.04im)
    k = single_block_kernel(M, 1, true, 2, true)
    cand = ring_factors(M, lat, 1, true, 2, true, alg)
    tgt = GM.togmps(ComplexF64, GTEMPO.fullfactor(k, lat), first(GTEMPO.site_spaces(lat)), length(lat) ÷ 2)
    println("② 同味 (T,T) 1→2: dev = ", gd(mat(cand, lat), GM.Dense.todense(tgt).data))
    # (F,F) du
    k2 = single_block_kernel(M, 2, false, 1, false)
    cand2 = ring_factors(M, lat, 2, false, 1, false, alg)
    tgt2 = GM.togmps(ComplexF64, GTEMPO.fullfactor(k2, lat), first(GTEMPO.site_spaces(lat)), length(lat) ÷ 2)
    println("② 同味 (F,F) 2→1: dev = ", gd(mat(cand2, lat), GM.Dense.todense(tgt2).data))
end

# ---- ③ band-diagonal (T,F) 回归锚：ba==bb==1 应逐位复现现行公式 ----
let n = 5
    lat = GTEMPO.GrassmannLattice(N=n, δτ=0.25, bands=2, contour=:imag)
    lower = ComplexF64[0.3 * 0.7^(m-1) for m in 1:n-1]
    M = ring_kernel(lower, 0.05)
    k = single_block_kernel(M, 1, true, 1, false)
    cand = ring_factors(M, lat, 1, true, 1, false, alg)
    tgt = GM.togmps(ComplexF64, GTEMPO.fullfactor(k, lat), first(GTEMPO.site_spaces(lat)), length(lat) ÷ 2)
    println("③ 带对角 (T,F): dev = ", gd(mat(cand, lat), GM.Dense.todense(tgt).data))
end

# ---- ④ 端到端：BCS 虚时核（3 味对齐全）ExactTTI 候选 ≡ PartialIF ----
let n = 4
    lat = GTEMPO.GrassmannLattice(N=n, δτ=0.25, bands=2, contour=:imag, symmetry=Z2Symmetry())
    bath = GTEMPO.bcsbath(GTEMPO.fermionicbath(GTEMPO.DiracDelta(ω=1.0, α=0.5), β=lat.β), Δ=0.4)
    kernel = correlationfunction(bath, lat)
    println("④ BCS 核键数 = ", length(kernel.blocks), "  anomalous = ", kernel.anomalous)
    allfactors = GM.GrassmannFactor[]
    for (keyt, blk) in kernel.blocks
        (keyt[1].branch === :τ && keyt[2].branch === :τ) || continue
        M = Matrix{ComplexF64}(blk)
        append!(allfactors, ring_factors(M, lat, keyt[1].band, keyt[1].conj, keyt[2].band, keyt[2].conj, alg))
    end
    spv = GTEMPO.site_spaces(lat)
    acc, _ = GM.expmult!(GTEMPO.vacuumstate(lat), first(allfactors), GM.PolarExp(rtol=0.0))
    for f in allfactors[2:end]
        acc, _ = GM.expmult!(acc, f, GM.PolarExp(rtol=0.0))
    end
    full, _ = hybriddynamics(lat, kernel, PartialIF())
    d = maximum(abs.(GM.Dense.todense(acc).data - GM.Dense.todense(full).data))
    println("④ BCS 端到端（候选因子 PolarExp 装配 vs PartialIF）: dev = ", d)
end
