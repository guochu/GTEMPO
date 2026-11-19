# =============================================================================
# R7 排查第五步：U=0 解析 BdG 参照（独立于 bcs_ed 多体 ED）
# 四模（d↑, d↓, b↑, b↓）二次型：H = Σ h c†c + (−Δ b†↑b†↓ + Δ b↓b↑) + h.c.
# BdG 8×8，G⁺(t) = U(t)·(I−f0)，gt = −i·G⁺[d↑,d↑]，lt = +i·(f0·U(t))[d↑,d↑]。
# 对比：① bcs_ed（U=0）② TN order=2 大 Nt 的收敛极限。
# 裁决：bcs_ed 与解析不一致 ⟹ ED 错；TN 极限 ≠ 解析 ⟹ TN 核错。
# =============================================================================
using GTEMPO, LinearAlgebra, Printf

include(joinpath(@__DIR__, "..", "test", "util.jl"))
include(joinpath(@__DIR__, "..", "test", "normalbath", "util.jl"))

const ω = 1.0; const α = 0.5; const β = 1.0; const T = 0.4
const ϵ_d = -0.8          # 模型 μ（= bcs_ed 的 -ϵ_ed）
const Δ_bcs = 0.6

# ---- 解析 BdG（4 模）----
# 模序：(d↑, d↓, b↑, b↓)；h = 对角 + √α 杂化；Π_{b↑,b↓} = −Δ（配对）
function bdg_gt(Δ; ts, β, ϵ_d, ω, α)
    n = 4
    h = zeros(ComplexF64, n, n)
    h[1, 1] = h[2, 2] = ϵ_d           # 单粒子能量 = μ = ϵ_d（模型 H = μ·n）
    h[3, 3] = h[4, 4] = ω
    h[1, 3] = h[3, 1] = sqrt(α)       # d↑–b↑ 杂化
    h[2, 4] = h[4, 2] = sqrt(α)       # d↓–b↓
    Π = zeros(ComplexF64, n, n)
    Π[3, 4] = -Δ                      # b↑,b↓ 配对（H_pair = Π₃₄ b†↑b†↓ + c.c.，k<l 约定）
    Π[4, 3] = +Δ
    Hbdg = [h Π; Π -conj(h)]          # 8×8（h 实对称 ⟹ −h* = −h）
    f0 = Diagonal(1 ./ (exp.(β .* (diag(h) .- ϵ_d)) .+ 1))   # μ = ϵ_d 全模（浴 μ=0）
    # 浴模 μ = 0：f_b = 1/(e^{βω}+1)
    f0[3, 3] = f0[4, 4] = 1 / (exp(β * ω) + 1)
    I4 = Matrix{ComplexF64}(I, n, n)
    gt = [(-im * (exp(-im * Matrix(Hbdg) * t)[1:n, 1:n] * (I4 - Matrix(f0))))[1, 1] for t in ts]
    lt = [( im * (Matrix(f0) * exp(-im * Matrix(Hbdg) * t)[1:n, 1:n]))[1, 1] for t in ts]
    return gt, lt
end

ts = (1:16) .* (T / 16)
gt_an, lt_an = bdg_gt(Δ_bcs; ts=ts, β=β, ϵ_d=ϵ_d, ω=ω, α=α)

# ---- bcs_ed（U=0 多体 ED，16 维）----
H, Ab, Bb, _ = bcs_ed(0.0, -ϵ_d; ω₀=ω, α=α, Δ=Δ_bcs)     # U=0, bcs_ed(ϵ_d=−μ)
n̂ = spin_half_matrices()["n"]
I2 = one(n̂)
n_imp = kron(n̂, I2) + kron(I2, n̂)
Himpbare = ϵ_d * n_imp
I4 = Matrix{ComplexF64}(I, 4, 4)
n_bath = kron(I2, n̂) + kron(n̂, I2)
H0dec = kron(Himpbare, I4) + ω * kron(I4, n_bath)
gt_ed, lt_ed = greater_lesser_ed(H, Ab, Bb, H0dec, ts, β)

println("U=0 Δ=", Δ_bcs, " 解析 BdG vs bcs_ed（16 维 ED）")
println("   t      gt(解析)          gt(ED)          |dev|          lt|dev|")
for (i, t) in enumerate(ts)
    @printf("  %.3f  %+.9f%+9.2eim  %+.9f%+9.2eim  %.2e  %.2e\n", t,
            real(gt_an[i]), imag(gt_an[i]), real(gt_ed[i]), imag(gt_ed[i]),
            abs(gt_an[i] - gt_ed[i]), abs(lt_an[i] - lt_ed[i]))
end
println("max dev: gt = ", maximum(abs.(gt_an .- gt_ed)), "   lt = ", maximum(abs.(lt_an .- lt_ed)))

# ---- TN 收敛极限（U=0，order=2）----
println("\nTN order=2 收敛（vs 解析 BdG）：")
model = AndersonIM(U=0.0, ϵ_d=ϵ_d)
ρ0 = thermal_impurity_state(model, β)
bath = bcsbath(fermionicbath(DiracDelta(ω=ω, α=α), β=β), Δ=Δ_bcs)
cn = coupling_norm(bath) * 2
for Nt in (16, 32)
    δt = T / Nt
    lat = GrassmannLattice(N=Nt, δt=δt, bands=2, contour=:real, symmetry=Z2Symmetry())
    kernel = correlationfunction(bath, lat)
    Is, _ = hybriddynamics(lat, kernel, PartialIF())
    K = sysdynamics(lat, model; order=2, coupling_norm=cn, initial_state=ρ0)
    Z = integrate(K, Is)
    stride = 16 ÷ Nt
    idx = (1:Nt) .* stride
    gt = [greater(lat, model, kernel, Is; i=i, j=0, Z=Z, order=2,
                  coupling_norm=cn, initial_state=ρ0) for i in 1:Nt]
    dev = maximum(abs.(gt .- gt_an[idx]))
    @printf("  Nt=%d: max|gt_TN − gt_解析| = %.3e（端点 %.3e）\n", Nt, dev,
            abs(gt[end] - gt_an[idx[end]]))
end
