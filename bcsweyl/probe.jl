# bcsweyl/probe.jl — 快速探针：imag 轮廓 BCS IF × order 1/2 K 是否可跑、量级合理
# 复刻 test/bcs/fewmodesbath 的 ED 约定（bcs_ed 自带算符、bathZ 归一化分母）
using GTEMPO
using ImpurityModelBase
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "..", "test", "util.jl"))   # bcs_ed / spin_half_matrices

const U, ϵ_ed, ω, α, β = 1.0, 0.5, 1.0, 0.5, 0.8
const Δ = 0.3

# ED 参照（与 test/bcs 同款：bcs_ed 哈密顿量 + 浴自身配分归一化）
H, Ab, Bb, H0 = bcs_ed(U, ϵ_ed; ω₀=ω, α=α, Δ=Δ)
p1 = spin_half_matrices()
n̂, σ₊, σ₋, JW = p1["n"], p1["+"], p1["-"], -p1["z"]
I2 = one(n̂)
Hbath = ω * (kron(n̂, I2) + kron(I2, n̂)) -
        (Δ * kron(JW * σ₊, σ₊) - conj(Δ) * kron(JW * σ₋, σ₋))
bathZ = tr(exp(-β * Hbath))
Zed = real(tr(exp(-β * H))) / real(bathZ)
println("ED: bathZ = $(real(bathZ)),  Z_exact = tr(e^{-βH})/bathZ = $Zed")

for order in (1, 2)
    Nτ = 8
    δτ = β / Nτ
    lat = GrassmannLattice(N=Nτ, δτ=δτ, bands=2, contour=:imag, symmetry=Z2Symmetry())
    model = AndersonIM(U=U, ϵ_d=-ϵ_ed)
    bath = bcsbath(fermionicbath(DiracDelta(ω=ω, α=α), β=β), Δ=Δ)
    Is = bcs_hybriddynamics(lat, bath)
    K = order == 2 ?
        sysdynamics_weyl2order(lat, model; coupling_norm=α) :
        sysdynamics(lat, model)
    Z = integrate(K, Is)
    @printf("order=%d  N=%2d  Z(K,IF) = %.8f   (Z−Zed)/Zed = %+.3e\n",
            order, Nτ, real(Z), (real(Z) - Zed) / Zed)
end
