# 验证：order=1 虚时 ⟨n⟩ 的非单调 = num/Z 误差对消残余
using GTEMPO
import GrassmannMPS as GM
using Printf
using LinearAlgebra: kron, I

include(joinpath(@__DIR__, "..", "test", "util.jl"))
include(joinpath(@__DIR__, "..", "test", "normalbath", "util.jl"))

const U = 1.0; const ϵ_ed = 0.8; const ω = 1.0; const α = 0.5
const β = 1.0; const Δ_bcs = 0.6
const CN = 2α
const bath = bcsbath(fermionicbath(DiracDelta(ω=ω, α=α), β=β), Δ=Δ_bcs)
const model = AndersonIM(U=U, ϵ_d=-ϵ_ed)
p1 = spin_half_matrices()
n̂, σ₊, σ₋, JW = p1["n"], p1["+"], p1["-"], -p1["z"]
I2 = one(n̂)
n_imp = kron(n̂, I2) + kron(I2, n̂)
H, Ab, Bb, H0 = bcs_ed(U, ϵ_ed; ω₀=ω, α=α, Δ=Δ_bcs)
I4 = Matrix{ComplexF64}(I, 4, 4)
Hbath = ω * (kron(n̂, I2) + kron(I2, n̂)) -
        (Δ_bcs * kron(JW * σ₊, σ₊) - conj(Δ_bcs) * kron(JW * σ₋, σ₋))
const Z_exact = real(tr(exp(-β * H)) / tr(exp(-β * Hbath)))
const n_eq_ed = real(tr((x = exp(-β * H); x /= tr(x); x) * kron(n_imp, I4)))
println("Z_exact = $Z_exact   n_eq = $n_eq_ed   num_exact = ", n_eq_ed * Z_exact)

for N in (4, 8, 16, 32)
    h = β / N
    lat = GrassmannLattice(N=N, δτ=h, bands=2, contour=:imag, symmetry=Z2Symmetry())
    kernel = correlationfunction(bath, lat)
    Is, _ = hybriddynamics(lat, kernel, PartialIF())
    Is = Is isa Tuple ? Is : (Is,)
    ρ0 = thermal_impurity_state(model, β)
    K = sysdynamics(lat, model; initial_state=ρ0)
    Z = integrate(K, Is...)
    n = real(occupation(lat, model, nothing, Is...; j=N ÷ 2, branch=:τ, band=nothing,
                        Z=Z, order=1, coupling_norm=CN, initial_state=ρ0))
    num = n * Z
    @printf("N=%2d: reldev(Z)   = %+.4e\n", N, (Z - Z_exact) / Z_exact)
    @printf("      reldev(num) = %+.4e   ⟨n⟩ 残差 = %+.4e\n",
            (num - n_eq_ed * Z_exact) / (n_eq_ed * Z_exact), n - n_eq_ed)
end
