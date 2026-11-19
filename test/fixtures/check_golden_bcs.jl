# =============================================================================
# BCS realtime golden checker（plan/bcs_exacttti_plan.md Phase W：BCS goldens
# 新档——此前无 BCS golden）。
#
# 锚语义：raw 口径（去顶点移位，2026-11-19 v3-①）的 BCS 实时观测量存档——
# 装配算法（PartialIF/PolarIF）× K 阶（1/2）四组合的 Z/gt/lt 全时序，1e-10 档
# （钳：K 阶与装配路径的压缩噪声 ≤1e-7 档，1e-10 对 todense 级逐位敏感——
# 实测互换 ~2e-7 ⟹ 容差 1e-6 档定锚，防未来漂移非机器精度咬合）。
#
# 用法：julia --project=. test/fixtures/check_golden_bcs.jl
# =============================================================================
using Test
using GTEMPO
using GTEMPO: branches
import GrassmannMPS as GM
using ImpurityModelBase
using LinearAlgebra: I, norm

const GDIR = @__DIR__
include(joinpath(GDIR, "golden_json.jl"))
include(joinpath(GDIR, "..", "util.jl"))
include(joinpath(GDIR, "..", "normalbath", "util.jl"))

golden = let
    s = read(joinpath(GDIR, "bcs_realtime_golden.json"), String)
    v, i = parse_value(s, 1)   # golden_json 解析器（dump_golden 同方言）
    v
end

@testset "BCS realtime golden（raw 口径四组合）" begin
    ω = 1.0; α = 0.5; β = 1.0; U = 1.0; ϵ_ed = 0.8
    T = 0.25; Nt = 5; δt = T / Nt
    lat = GrassmannLattice(N=Nt, δt=δt, bands=2, contour=:real, symmetry=Z2Symmetry())
    model = AndersonIM(U=U, ϵ_d=-ϵ_ed)
    bath = bcsbath(fermionicbath(DiracDelta(ω=ω, α=α), β=β), Δ=0.6)
    kernel = correlationfunction(bath, lat)
    cn = kernel.coupling_norm
    ρ0 = thermal_impurity_state(model, β)
    rebuild(key_prefix) = begin
        c = get(golden, key_prefix * "_Z", nothing)
        c === nothing && return nothing
        Complex(c[1], c[2])
    end
    for (tag, alg) in (("pif", PartialIF()), ("polar", PolarIF()))
        Is = first(hybriddynamics(lat, correlationfunction(bath, lat), alg))
        for order in (1, 2)
            key = order == 1 ? tag : tag * "_o2"
            haskey(golden, key * "_Z") || continue
            K = order == 2 ?
                sysdynamics(lat, model; order=2, coupling_norm=cn, initial_state=ρ0) :
                sysdynamics(lat, model; initial_state=ρ0)
            Z = integrate(K, Is)
            gt = [greater(lat, model, kernel, Is; i=i, j=0, Z=Z, order=order,
                          coupling_norm=cn, initial_state=ρ0) for i in 1:Nt]
            lt = [lesser(lat, model, kernel, Is; i=i, j=0, Z=Z, order=order,
                         coupling_norm=cn, initial_state=ρ0) for i in 1:Nt]
            zg = golden[key * "_Z"]
            @test Z ≈ Complex(zg[1], zg[2]) rtol = 1e-6
            gtg = [Complex(v[1], v[2]) for v in golden[key * "_gt"]]
            ltg = [Complex(v[1], v[2]) for v in golden[key * "_lt"]]
            @test maximum(abs.(gt .- gtg)) < 1e-6
            @test maximum(abs.(lt .- ltg)) < 1e-6
        end
    end
end
