# =============================================================================
# R7 排查第三步：Z 恒等式探针（U=0 线性耦合 ⟹ Keldysh 轨迹 Z=1 精确成立）
# BCS bath Δ≠0 vs Δ=0，order=1/2，Z−1 的收敛阶。
# 判定：若 Δ≠0 的 Z−1 撞 δt⁰ 地板而 Δ=0 干净 O(δt²) ⟹ 地板在反常 IF 收缩，
# 且有解析锚（Z=1）可迭代定位。
# =============================================================================
using GTEMPO, LinearAlgebra, Printf

include(joinpath(@__DIR__, "..", "test", "util.jl"))
include(joinpath(@__DIR__, "..", "test", "normalbath", "util.jl"))

const ω = 1.0; const α = 0.5; const β = 1.0; const T = 0.4
const ϵ_d = 0.7; const U = 1.0

println("U=1 AndersonIM（Keldysh 轨迹恒等式 Z=1 对线性耦合任意 H_S 精确成立）— BCS bath Z−1 收敛")
println("=" ^ 70)
for Δ in (0.6, 0.0)
    bath = bcsbath(fermionicbath(DiracDelta(ω=ω, α=α), β=β), Δ=Δ)
    cn = coupling_norm(bath) * 2
    model = AndersonIM(U=U, ϵ_d=ϵ_d)
    ρ0 = thermal_impurity_state(model, β)
    for order in (1, 2)
        print("Δ=$Δ order=$order: ")
        for Nt in (4, 8, 16, 32)
            δt = T / Nt
            lat = GrassmannLattice(N=Nt, δt=δt, bands=2, contour=:real, symmetry=Z2Symmetry())
            kernel = correlationfunction(bath, lat)
            Is, _ = hybriddynamics(lat, kernel, PartialIF())
            K = sysdynamics(lat, model; order=order, coupling_norm=cn, initial_state=ρ0)
            Z = integrate(K, Is)
            @printf("Nt=%d: |Z−1|=%.3e   ", Nt, abs(Z - 1))
        end
        println()
    end
end
# normal bath 对照（单带，Z=1 同样精确）
println("\nnormal bath 单带对照（Z=1）：")
bath = discretefermionicbath([ω], [α]; β=β)
model = ToulouseIM(ϵ_d=ϵ_d)
ρ0 = thermal_impurity_state(model, β)
for order in (1, 2)
    print("order=$order: ")
    for Nt in (4, 8, 16, 32)
        δt = T / Nt
        lat = GrassmannLattice(N=Nt, δt=δt, contour=:real)
        kernel = correlationfunction(bath, lat)
        Is, _ = hybriddynamics(lat, kernel, PartialIF())
        K = sysdynamics(lat, model; order=order, coupling_norm=α, initial_state=ρ0)
        Z = integrate(K, Is)
        @printf("Nt=%d: |Z−1|=%.3e   ", Nt, abs(Z - 1))
    end
    println()
end
