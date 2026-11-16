using Printf
include("/home/guochu/Documents/Missile/GTEMPO/test/util.jl")
import GTEMPO; using GTEMPO, ImpurityModelBase
using LinearAlgebra
μ = 0.5; ω = 1.0
β = 0.8; Nτ = 6; δτ = β / Nτ

function probe(nbands_bath, α, tag)
    nm = 1 + nbands_bath
    H = zeros(Float64, 1 << nm, 1 << nm)
    ed_density!(H, 1, μ)
    for m in 1:nbands_bath; ed_density!(H, 1 + m, ω); end
    for m in 1:nbands_bath
        v = sqrt(α)
        ed_tunneling!(H, 1 + m, 1, v); ed_tunneling!(H, 1, 1 + m, v)
    end
    cs = fermion_annihilators(nm)
    n_ed = real(tr(cs[1]'cs[1]*exp(-β*H)))/real(tr(exp(-β*H)))
    Zed = real(tr(exp(-β*H))) / (1 + exp(-β*ω))^nbands_bath
    lat = GrassmannLattice(N=Nτ, δτ=δτ, bands=1, contour=:imag, symmetry=Z2())
    bath = discretefermionicbath(fill(ω, nbands_bath), fill(α, nbands_bath); β=β)
    Is = (hybriddynamics(lat, bath, PartialIF(); band=1),)
    K = sysdynamics(lat, ToulouseIM(ϵ_d=μ))
    Z = integrate(K, Is...)
    occ = occupation(lat, ToulouseIM(ϵ_d=μ), nothing, Is...; j=1, branch=:τ, Z=Z)
    @printf("%s: Z = %.5f (ED %.5f)  occ = %.5f (ED %.5f)\n", tag, real(Z), Zed, real(occ), n_ed)
end
probe(1, 0.5, "C1 单模 α=0.5")
probe(1, 2.0, "C2 单模 α=2.0")
probe(4, 0.5, "D  4模 α=0.5")
