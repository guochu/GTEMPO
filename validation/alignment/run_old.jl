# 旧栈（backup/GTEMPO）alignment sweep 数据采集
# 用法：julia --project=../backup/GTEMPO validation/alignment/run_old.jl out_old.json
using GTEMPO
using ImpurityModelBase
using LinearAlgebra: Diagonal, tr
using Serialization

μ = 0.4; β = 0.8; ω = 1.0; α = 0.5
out = Dict{String,Any}()

# S1: 虚时单带（δ-bath）
lat = GrassmannLattice(N=4, δτ=0.2, contour=:imag)
bath = fermionicbath(DiracDelta(ω=ω, α=α), β=β)
corr = correlationfunction(bath, lat)
trunc = truncdimcutoff(D=64, ϵ=1e-10)
mpsI = hybriddynamics(lat, corr, band=1, trunc=trunc)
mpsK = sysdynamics(lat, ToulouseIM(ϵ_d=μ), trunc=trunc)
mpsK = boundarycondition!(mpsK, lat)
Z = integrate(lat, mpsK, mpsI)
out["S1_Z"] = ComplexF64(Z)
out["S1_occ"] = [real(occupation(lat, i, mpsK, mpsI, Z=Z)) for i in 1:lat.k]

# S2: 实时单带 + 热初态
lat = GrassmannLattice(N=4, δt=0.05, contour=:real)
corr = correlationfunction(bath, lat)
mpsI = hybriddynamics(lat, corr, band=1, trunc=trunc)
mpsK = sysdynamics(lat, ToulouseIM(ϵ_d=μ), trunc=trunc)
mpsK = boundarycondition!(mpsK, lat)
mpsK = systhermalstate!(mpsK, lat, ToulouseIM(ϵ_d=μ), trunc=trunc, β=β)
Z = integrate(lat, mpsK, mpsI)
gt = [-im * gf(lat, (GrassmannContourIndex(k, conj=false, branch=:+, band=1),
                      GrassmannContourIndex(1, conj=true, branch=:+, band=1)), mpsK, mpsI, Z=Z) for k in 1:lat.k]
out["S2_Z"] = ComplexF64(Z)
out["S2_greater"] = ComplexF64.(gt)

# S3: 双带 U=1（旧栈 M=2 带互换 bug 场景——预期分歧，注释在案）
lat = GrassmannLattice(N=4, δτ=0.2, bands=2, contour=:imag)
bath = fermionicbath(DiracDelta(ω=ω, α=α), β=β)
corr = correlationfunction(bath, lat)
mpsI1 = hybriddynamics(lat, corr, band=1, trunc=trunc)
mpsI2 = hybriddynamics(lat, corr, band=2, trunc=trunc)
mpsK = sysdynamics(lat, AndersonIM(U=1.0, ϵ_d=μ), trunc=trunc)
for band in 1:2
    global mpsK = boundarycondition!(mpsK, lat, band=band, trunc=trunc)
end
Z3 = integrate(lat, mpsK, mpsI1, mpsI2)
out["S3_Z"] = ComplexF64(Z3)
out["S3_occ"] = [real(occupation(lat, i, mpsK, mpsI1, mpsI2, Z=Z3, band=1)) for i in 1:lat.k]

open(ARGS[1], "w") do io
    serialize(io, out)
end
println("written ", ARGS[1])
