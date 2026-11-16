include("/home/guochu/Documents/Missile/GTEMPO/test/util.jl")
import GTEMPO; using GTEMPO, ImpurityModelBase
using LinearAlgebra
U = 0.0; μ = 0.5; ω = 1.0; α = 0.5
β = 0.8; Nτ = 6; δτ = β / Nτ
norb = 2; nimp = 4
# ED = 对角 1:1（带 b ↔ 浴模 b）
nm = 2nimp
H = zeros(Float64, 1 << nm, 1 << nm)
for b in 1:nimp; ed_density!(H, b, μ); end
for m in 1:nimp; ed_density!(H, nimp + m, ω); end
for b in 1:nimp
    v = sqrt(α)
    ed_tunneling!(H, nimp + b, b, v); ed_tunneling!(H, b, nimp + b, v)
end
Zed = real(tr(exp(-β * H))) / (1 + exp(-β * ω))^nimp

lat = GrassmannLattice(N=Nτ, δτ=δτ, bands=4, contour=:imag, symmetry=U1Symmetry())
# Σ 核 = 4×(α/4)g = α·g：每带恰见自身通道强度（对角混杂）
bath = discretefermionicbath(fill(ω, nimp), fill(α / nimp, nimp); β=β)
model = ImpurityHamiltonian([tunneling(b, b, coeff=μ) for b in 1:4]; bands=4)
Is = (hybriddynamics(lat, bath, PartialIF(); band=1),)
K = sysdynamics(lat, model)
Z = integrate(K, Is...)
occ = occupation(lat, model, nothing, Is...; j=1, branch=:τ, band=1, Z=Z)
println(":u1 bands=4 合并 IF（2 链）: Z/bathZ = ", real(Z)/((1+exp(-β*ω))^nimp),
        "  ED = ", Zed, "  occ(band1) = ", real(occ),
        "  ED n = ", real(tr(fermion_annihilators(8)[1]'fermion_annihilators(8)[1]*exp(-β*H)))/real(tr(exp(-β*H))))
