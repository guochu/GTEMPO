include("/home/guochu/Documents/Missile/GTEMPO/test/util.jl")
import GTEMPO; using GTEMPO, ImpurityModelBase
using LinearAlgebra
U = parse(Float64, get(ENV, "NU", "2.0")); J = parse(Float64, get(ENV, "NJ", "0.4")); μ = 0.5; ω = 1.0; α = 0.5
β = 0.8; Nτ = parse(Int, get(ENV, "NT", "12")); δτ = β / Nτ
norb = 2; nimp = 4

# ED = kanamori_ed 对角 1:1（带 m ↔ 浴模 m）
if U == 0
    # 二次型对角 ED：带 b ↔ 浴模 b（√α）
    nm = 2nimp
    H = zeros(Float64, 1 << nm, 1 << nm)
    for b in 1:nimp; ed_density!(H, b, μ); end
    for m in 1:nimp; ed_density!(H, nimp + m, ω); end
    for b in 1:nimp
        v = sqrt(α)
        ed_tunneling!(H, nimp + b, b, v); ed_tunneling!(H, b, nimp + b, v)
    end
else
    H, _, _, _ = kanamori_ed(U=U, J=J, norb=norb, μ=μ, ω=ω, α=α)
end
cs = fermion_annihilators(8)
Zed = real(tr(exp(-β * H)))
bathZ = (1 + exp(-β * ω))^nimp
n1 = real(tr(cs[1]'cs[1]*exp(-β*H)))/Zed

lat = GrassmannLattice(N=Nτ, δτ=δτ, bands=2norb, contour=:imag, symmetry=:u1u1)
# Σ 核 = Σ_m (α/4) g = α·g：每带恰见自身通道强度（对角混杂 ⟹ 乘积 IF 精确）
bath = discretefermionicbath(fill(ω, nimp), fill(α / nimp, nimp); β=β)
model = KanamoriIM(U=U, J=J, norb=norb, μ=μ)
Is = ntuple(b -> hybriddynamics(lat, bath, PartialIF(); band=b), 2norb)
K = sysdynamics(lat, model)
Z = integrate(K, Is...)
occ = occupation(lat, model, nothing, Is...; j=1, branch=:τ, band=1, Z=Z)
println("norb=2 :u1u1 对角耦合: Z/bathZ = ", real(Z)/bathZ, "  ED/bathZ = ", Zed/bathZ,
        "  occ = ", real(occ), "  ED = ", n1)
