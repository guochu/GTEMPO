# 新栈（GrassmannMPS 后端）alignment sweep 数据采集
# 用法：julia --project=. validation/alignment/run_new.jl out_new.json
using GTEMPO
using ImpurityModelBase
using LinearAlgebra: Diagonal, tr
using Serialization

μ = 0.4; β = 0.8; ω = 1.0; α = 0.5
out = Dict{String,Any}()

# S1: 虚时单带（δ-bath）
lat = GrassmannLattice(N=4, δτ=0.2, bands=1, contour=:imag)
bath = discretefermionicbath([ω], [α]; β=β)
Is = (hybriddynamics(lat, bath, PartialIF(); band=1),)
K = sysdynamics(lat, ToulouseIM(ϵ_d=μ))
Z = integrate(K, Is...)
out["S1_Z"] = ComplexF64(Z)
out["S1_occ"] = [real(occupation(lat, ToulouseIM(ϵ_d=μ), nothing, Is...; j=i, branch=:τ, Z=Z)) for i in 1:lat.k]

# S2: 实时单带 + 热初态
lat = GrassmannLattice(N=4, δt=0.05, bands=1, contour=:real)
bath = discretefermionicbath([ω], [α]; β=β)
Is = (hybriddynamics(lat, bath, PartialIF(); band=1),)
model = ToulouseIM(ϵ_d=μ)
himρ = fock_thermalstate(model, β)
K = sysdynamics(lat, model; initial_state=Matrix{ComplexF64}(himρ.data))
Z = integrate(K, Is...)
kernel = hybridization(lat, bath)
out["S2_Z"] = ComplexF64(Z)
out["S2_greater"] = ComplexF64[greater(lat, model, kernel, Is...; i=k, j=0, Z=Z, initial_state=Matrix{ComplexF64}(himρ.data)) for k in 1:lat.k]

# S3: 双带 U=1（旧栈 M=2 带互换 bug 场景——预期分歧，注释在案）
lat = GrassmannLattice(N=4, δτ=0.2, bands=2, contour=:imag)
bath = discretefermionicbath([ω], [α]; β=β)
Is = ntuple(b -> hybriddynamics(lat, bath, PartialIF(); band=b), 2)
K = sysdynamics(lat, AndersonIM(U=1.0, ϵ_d=μ))
Z3 = integrate(K, Is...)
out["S3_Z"] = ComplexF64(Z3)
out["S3_occ"] = [real(occupation(lat, AndersonIM(U=1.0, ϵ_d=μ), nothing, Is...; j=i, branch=:τ, band=1, Z=Z3)) for i in 1:lat.k]

open(ARGS[1], "w") do io
    serialize(io, out)
end
println("written ", ARGS[1])
