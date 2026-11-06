# ============================================================================
# validation/ed.jl — 单模 Anderson 模型的精确对角化参考
#
# 费米 bit 基（mode m ⟺ bit m−1，JW 符号），与 GTEMPO test/util.jl 的
# singlemode_ed 同一构造：
#   H = μ Σ_imp n_i + U n₁n₂ + Σ_m [ω_m b†_m b_m + √α_m (a†_m b_m + b†_m a_m)]
# 只比较数值（⟨n⟩、G^>、G^<、Gτ），与 K 的 Fock 约定无耦合。
# ============================================================================

module SingleModeED

using LinearAlgebra: eigen, Diagonal, tr

_fermion_sign(state::Int, m::Int) =
    iseven(count_ones(state & ((1 << (m - 1)) - 1))) ? 1 : -1

function _apply_seq(state::Int, ops::Vector{Tuple{Int,Bool}})
    sign = 1; s = state
    for (m, create) in ops
        bit = 1 << (m - 1)
        if create
            (s & bit) != 0 && return 0, 0
            sign *= _fermion_sign(s, m); s |= bit
        else
            (s & bit) == 0 && return 0, 0
            sign *= _fermion_sign(s, m); s &= ~bit
        end
    end
    return sign, s
end

"coeff × c†_i c_j"
function tunneling!(H, i, j, coeff = 1.0)
    for s in 0:size(H, 1)-1
        sign, t = _apply_seq(s, [(j, false), (i, true)])
        sign != 0 && (H[t+1, s+1] += coeff * sign)
    end
    return H
end
density!(H, i, coeff = 1.0) = tunneling!(H, i, i, coeff)

"""
    model(; μ, U, bathspecs, nbands) -> (H, d, H0)

`bathspecs = [(ω, α), …]`（每带一个单模浴）。返回完整 H、去耦 H0、
band 1 杂质湮灭算符 d。
"""
function model(; μ::Real, U::Real, bathspecs::Vector{<:Tuple{Real,Real}},
               nbands::Int = length(bathspecs))
    nmodes = nbands + length(bathspecs)
    dim = 1 << nmodes
    H = zeros(Float64, dim, dim)
    for i in 1:nbands
        density!(H, i, μ)
    end
    if nbands == 2 && U != 0
        for s in 0:dim-1
            ((s & 1) == 1) && ((s >> 1) & 1 == 1) && (H[s+1, s+1] += U)
        end
    end
    H0 = copy(H)
    for (m, (ω, α)) in enumerate(bathspecs)
        b = nbands + m
        imp = ((m - 1) % nbands) + 1
        density!(H, b, ω); density!(H0, b, ω)
        v = sqrt(α)
        tunneling!(H, b, imp, v)   # b† a
        tunneling!(H, imp, b, v)   # a† b
    end
    d = zeros(Float64, dim, dim)
    for s in 0:dim-1
        sign, t = _apply_seq(s, [(1, false)])
        sign != 0 && (d[t+1, s+1] = sign)
    end
    return H, d, H0
end

"热密度矩阵 exp(−βH₀)/tr（β=Inf 给出基态投影子）"
function thermal(H, β)
    vals, vecs = eigen(H)
    if isinf(β)
        return vecs[:, 1] * vecs[:, 1]'
    end
    w = exp.(.-β .* (vals .- minimum(vals)))
    return vecs * Diagonal(w ./ sum(w)) * vecs'
end

"""
    observables(H, d, ts, ρ₀) -> (occ, gt, lt)

初态 ρ₀（含浴）下全空间 Heisenberg 关联：
⟨n(t)⟩、G^>(t,0) = −i⟨d(t)d†(0)⟩、G^<(t,0) = +i⟨d†(0)d(t)⟩。
"""
function observables(H, d, ts, ρ₀)
    vals, vecs = eigen(H)
    U0 = vecs' * d * vecs
    ρe = vecs' * ρ₀ * vecs
    n = U0' * U0
    occ = ComplexF64[]; gt = ComplexF64[]; lt = ComplexF64[]
    for t in ts
        Dmid = Diagonal(exp.(im .* t .* vals))
        dnt = Dmid * U0 * Dmid'                       # d(t)
        push!(occ, tr(Dmid * n * Dmid' * ρe))
        push!(gt, -1.0im * tr(dnt * U0' * ρe))
        push!(lt, 1.0im * tr(U0' * dnt * ρe))
    end
    return occ, gt, lt
end

"虚时热平均：⟨n⟩ 与 Matsubara G(τ,0) = ⟨d(τ)d†(0)⟩（费米时间序）"
function imag_observables(H, d, τs, β)
    vals, vecs = eigen(H)
    U0 = vecs' * d * vecs
    Z = sum(exp.(.-β .* vals))
    n = U0' * U0
    occ = tr(Diagonal(exp.(.-β .* vals)) * n) / Z
    gτ = ComplexF64[]
    for τ in τs
        # Σ_{ab} e^{(τ−β)D_a} e^{−τD_b} |U0_{ab}|² / Z
        s = sum(exp((τ - β) * vals[a]) * exp(-τ * vals[b]) * abs2(U0[a, b])
                for a in eachindex(vals), b in eachindex(vals))
        push!(gτ, s / Z)
    end
    return occ, gτ
end

end # module
