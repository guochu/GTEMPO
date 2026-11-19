# =============================================================================
# R7 排查第二步：臂 A（BCS Δ=0.6）误差的逐时刻分解 + order=1 对照
# 问题：地板误差集中在哪些时刻？端点 vs 均匀分布 ⟹ 定位接触门/尾处理/离散化
# =============================================================================
using GTEMPO, LinearAlgebra, Printf

include(joinpath(@__DIR__, "..", "test", "util.jl"))
include(joinpath(@__DIR__, "..", "test", "normalbath", "util.jl"))

const U = 1.0; const ϵ_ed = 0.8; const ω = 1.0; const α = 0.5
const β = 1.0; const T = 0.4; const Δ_bcs = 0.6
const NMAX = 64

p1 = spin_half_matrices()
n̂, σ₊, σ₋, JW = p1["n"], p1["+"], p1["-"], -p1["z"]
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

for order in (1, 2)
    println("\n================ order=$order ================")
    for Nt in (8, 16, 32)
        δt = T / Nt
        lat = GrassmannLattice(N=Nt, δt=δt, bands=2, contour=:real, symmetry=Z2Symmetry())
        kernel = correlationfunction(bath, lat)
        Is, _ = hybriddynamics(lat, kernel, PartialIF())
        K = sysdynamics(lat, model; order=order, coupling_norm=cn, initial_state=ρ0)
        Z = integrate(K, Is)
        stride = NMAX ÷ Nt
        idx = (1:Nt) .* stride
        gt = [greater(lat, model, kernel, Is; i=i, j=0, Z=Z, order=order,
                      coupling_norm=cn, initial_state=ρ0) for i in 1:Nt]
        lt = [lesser(lat, model, kernel, Is; i=i, j=0, Z=Z, order=order,
                     coupling_norm=cn, initial_state=ρ0) for i in 1:Nt]
        println("-- Nt=$Nt（δt=$(round(δt, digits=5))）--")
        println("   i   t        gt_err        lt_err       gt_val")
        for i in 1:Nt
            @printf("  %3d  %.3f   %.3e    %.3e   %.6f\n", i, i*δt,
                    abs(gt[i] - gt_ed64[idx[i]]), abs(lt[i] - lt_ed64[idx[i]]), real(gt[i]))
        end
    end
end
