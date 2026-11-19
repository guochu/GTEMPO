# =============================================================================
# R7 排查第四步：反常块幅度 λ-二分（order=2 greater 的地板 vs λ 结构）
# Δ=0.6 核，ud/du 反常块值缩放 λ ∈ {0, 0.25, 0.5, 1}；ED 参照 = λ=1 物理。
# e(λ) = gt_λ − gt_ED(λ=1)：λ=0 应复现臂 B（无地板收敛）；地板随 λ 的阶
# （线性 = O(uv) 机制 / 二次 = O(uv²) 机制）。
# =============================================================================
using GTEMPO, LinearAlgebra, Printf

include(joinpath(@__DIR__, "..", "test", "util.jl"))
include(joinpath(@__DIR__, "..", "test", "normalbath", "util.jl"))

const U = 1.0; const ϵ_ed = 0.8; const ω = 1.0; const α = 0.5
const β = 1.0; const T = 0.4; const Δ_bcs = 0.6
const NMAX = 64

p1 = spin_half_matrices()
n̂ = p1["n"]
I2 = one(n̂)
n_imp = kron(n̂, I2) + kron(I2, n̂)
Himpbare = -ϵ_ed * n_imp + U * kron(n̂, n̂)
I4 = Matrix{ComplexF64}(I, 4, 4)
n_bath = kron(I2, n̂) + kron(n̂, I2)

ts = (1:NMAX) .* (T / NMAX)
model = AndersonIM(U=U, ϵ_d=-ϵ_ed)
ρ0 = thermal_impurity_state(model, β)
H, Ab, Bb, _ = bcs_ed(U, ϵ_ed; ω₀=ω, α=α, Δ=Δ_bcs)
H0dec = kron(Himpbare, I4) + ω * kron(I4, n_bath)
gt_ed64, lt_ed64 = greater_lesser_ed(H, Ab, Bb, H0dec, ts, β)
bath = bcsbath(fermionicbath(DiracDelta(ω=ω, α=α), β=β), Δ=Δ_bcs)
cn = coupling_norm(bath) * 2

"反常块（两腿同味键）缩放 λ"
function scaled_kernel(bath, lat, λ)
    k = correlationfunction(bath, lat)
    blocks = Dict{Tuple{ContourBlockIndex,ContourBlockIndex},eltype(values(k.blocks))}()
    for (key, v) in k.blocks
        if key[1].conj == key[2].conj          # 同味 = 反常块
            blocks[key] = _scale_block(v, λ)
        else
            blocks[key] = v
        end
    end
    return HybridizationKernel(blocks, k.coupling_norm, k.bands)
end
_scale_block(v::CorrelationMatrix, λ) = CorrelationMatrix(λ .* v.ηⱼₖ, λ .* v.ηₖⱼ)
_scale_block(v::AbstractMatrix, λ) = λ .* v

for Nt in (16, 32)
    δt = T / Nt
    lat = GrassmannLattice(N=Nt, δt=δt, bands=2, contour=:real, symmetry=Z2Symmetry())
    stride = NMAX ÷ Nt
    idx = (1:Nt) .* stride
    println("\n======== Nt=$Nt, order=2 ========")
    println("  λ      gt_err(end)    lt_err(end)    gt_err(t=T/2)   |gt_λ(end) − gt_1(end)|")
    gt_prev_end = nothing
    for λ in (0.0, 0.25, 0.5, 1.0)
        kernel = scaled_kernel(bath, lat, λ)
        Is, _ = hybriddynamics(lat, kernel, PartialIF())
        K = sysdynamics(lat, model; order=2, coupling_norm=cn, initial_state=ρ0)
        Z = integrate(K, Is)
        gt = [greater(lat, model, kernel, Is; i=i, j=0, Z=Z, order=2,
                      coupling_norm=cn, initial_state=ρ0) for i in 1:Nt]
        lt = [lesser(lat, model, kernel, Is; i=i, j=0, Z=Z, order=2,
                     coupling_norm=cn, initial_state=ρ0) for i in 1:Nt]
        eend_g = abs(gt[end] - gt_ed64[idx[end]])
        eend_l = abs(lt[end] - lt_ed64[idx[end]])
        emid = abs(gt[Nt÷2] - gt_ed64[idx[Nt÷2]])
        dprev = gt_prev_end === nothing ? NaN : abs(gt[end] - gt_prev_end)
        @printf("  %.2f   %.3e       %.3e       %.3e        %.3e\n", λ, eend_g, eend_l, emid, dprev)
        gt_prev_end = gt[end]
    end
end
