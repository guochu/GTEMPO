# Phase 2-5 golden checker：新栈重建 weyl golden 算例并断言（1e-10 档）。
#
# 用法：julia --project=. test/fixtures/check_golden_weyl.jl
#
# 锚语义：weyl golden 来自 extensions/secondorderif（cell 离散 + Chrono 排布 +
# Weyl/一阶 cell K + 逐带 cell IF），与本栈离散化同构（其对 U1GTEMPO2 端到端
# ≤1.4e-11）。断言对象 = occupation/gf 全端点（1e-9 档）；Z 报告不断言
# （bare/gaussian 测度约定差一个常数标度——INTEGRATION.md §8 只比观测量的原因）。

using Test
using GTEMPO
using GTEMPO: branches
import GrassmannMPS as GM
using ImpurityModelBase
using LinearAlgebra: diagm, I

include(joinpath(@__DIR__, "golden_json.jl"))

const μ, ω, α, β = 0.7, 1.0, 0.5, 1.0
const U2 = 0.8
bath = fermionicbath(DiracDelta(ω=ω, α=α), β=β)
ρ0_1 = diagm(0 => ComplexF64[1.0, exp(-β * μ)]) / (1 + exp(-β * μ))
ρ0_2 = diagm(0 => ComplexF64[1.0, exp(-β * μ), exp(-β * μ), exp(-β * (2μ + U2))]) /
       (1 + 2exp(-β * μ) + exp(-β * (2μ + U2)))

const TOL = 1e-9

function run_case(name)
    ref = load_golden(name)
    ref === nothing && (println("SKIP ", name, " (no fixture)"); return nothing)
    kind = ref["kind"]; order = ref["order"]
    bands = haskey(ref, "occupation_b1") ? 2 : 1
    if kind == "weyl_WeylRealTime"
        lat = GrassmannLattice(N=3, δt=0.15, bands=bands, contour=:real)
        model = bands == 2 ? AndersonIM(U=U2, ϵ_d=μ) : ToulouseIM(ϵ_d=μ)
        rho0 = bands == 2 ? ρ0_2 : ρ0_1
    elseif kind == "weyl_WeylMixedTime"
        lat = GrassmannLattice(Nτ=2, δτ=0.3, Nt=2, δt=0.15, contour=:mixed)
        model = ToulouseIM(ϵ_d=μ); rho0 = nothing
    elseif kind == "weyl_WeylImagTime"
        lat = GrassmannLattice(N=4, δτ=0.2, contour=:imag)
        model = ToulouseIM(ϵ_d=μ); rho0 = nothing
    else
        error("unknown kind $kind")
    end
    cn = α * bands
    kernel = correlationfunction(bath, lat)
    Is, _ = hybriddynamics(lat, kernel, PartialIF())
    Is isa GM.GMPS && (Is = (Is,))
    if order == 2
        K = sysdynamics(lat, model; order=2, coupling_norm=cn, initial_state=rho0)
    else
        K = sysdynamics(lat, model; order=1, initial_state=rho0)
    end
    zs = integrate_scaled(K, Is...)
    return (; ref, lat, model, rho0, order, cn, Is, K, Z=exp(zs.logabs) * zs.phase)
end

_occ(c, branch, band, j) = occupation(c.lat, c.model, nothing, c.Is...; j=j, branch=branch,
                                      band=band, order=c.order, coupling_norm=c.cn,
                                      initial_state=c.rho0, Z=c.Z)

@testset "weyl golden (cell discretization)" begin
    for name in ["weyl_real_o2_toulouse", "weyl_real_o2_anderson", "weyl_mixed_o2_toulouse",
                 "weyl_imag_o2_toulouse", "weyl_real_o1_toulouse"]
        c = run_case(name)
        c === nothing && continue
        ref = c.ref
        zref = ascomplex(ref["Z"])
        println(name, ":  Z_new=", round(c.Z, sigdigits=8), "  Z_golden=", round(zref, sigdigits=8),
                "  ratio=", round(real(c.Z / zref), sigdigits=6), " + ", round(imag(c.Z / zref), sigdigits=6), "im")
        # occupation 系列（键名 → branch/band 路由）
        routes = [("occupation", default_branch(c.lat), nothing),
                  ("occupation_tau", :τ, nothing),
                  ("occupation_rt", :+, nothing),
                  ("occupation_b1", default_branch(c.lat), 1),
                  ("occupation_b2", default_branch(c.lat), 2)]
        for (key, branch, band) in routes
            haskey(ref, key) || continue
            expected = ascomplex.(ref[key])
            got = [_occ(c, branch, band, j) for j in 1:length(expected)]
            err = maximum(abs.(got .- expected))
            @test err < TOL
            println("  ", key, ": max abs err = ", err)
        end
    end
end

