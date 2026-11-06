# rt_fastpath_check.jl — RevKeldyshOrdering 反向直铺 rt 路径校验
include(joinpath(ARGS[1], "benchmark", "secondorderif", "src", "SecondOrderWeyl.jl"))
using .SecondOrderWeyl
using .SecondOrderWeyl: occupation, greater, lesser
include(joinpath(ARGS[1], "benchmark", "secondorderif", "validation", "ed.jl"))
using .SingleModeED
using GTEMPO
using Printf
using LinearAlgebra: diagm

const μ, ω, α, β, U2 = 0.7, 1.0, 0.5, 1.0, 0.8

println("① 单带实时 ED 回归（Toulouse，order 2）：")
bath1 = fermionicbath(DiracDelta(ω = ω, α = α), β = β)
model1 = ToulouseIM(ϵ_d = μ)
H1, dd1, H01 = SingleModeED.model(μ = μ, U = 0.0, bathspecs = [(ω, α)])
ρ0_imp = diagm(0 => ComplexF64[1.0, exp(-β*μ)]) / (1 + exp(-β*μ))
ρfull1 = SingleModeED.thermal(H01, β)
occ_ed, gt_ed, lt_ed = SingleModeED.observables(H1, dd1, [0.4], ρfull1)
for N in (2, 4, 8)
    c = WeylRealTime(Nt = N, dt = 0.4 / N, bands = 1)
    sol = solve(Problem(c, model1; bath = bath1, order = 2, rho0 = ρ0_imp))
    occ = real(occupation(sol, N))
    g = greater(sol, N, 0)
    l = lesser(sol, N, 0)
    @printf("  N=%d  ⟨n⟩ err=%.3e  gt err=%.3e  lt err=%.3e\n",
            N, abs(occ - occ_ed[1]), abs(g - gt_ed[1]), abs(l - lt_ed[1]))
end

println("② 双带 Anderson 实时（order 2）vs ED + 计时：")
bath2 = fermionicbath(DiracDelta(ω = ω, α = α), β = β)
model2 = AndersonIM(ϵ_d = μ, U = U2)
H2, dd2, H02 = SingleModeED.model(μ = μ, U = U2, bathspecs = [(ω, α), (ω, α)])
ρ0_imp2 = diagm(0 => ComplexF64[1.0, exp(-β*μ), exp(-β*μ), exp(-β*(2μ + U2))]) /
          (1 + 2exp(-β*μ) + exp(-β*(2μ + U2)))
ρfull2 = SingleModeED.thermal(H02, β)
occ2_ed, _, _ = SingleModeED.observables(H2, dd2, [0.4], ρfull2)
ref = 2 * real(occ2_ed[1])
for N in (2, 3, 4, 6, 8)
    c = WeylRealTime(Nt = N, dt = 0.4 / N, bands = 2)
    t = @elapsed sol = solve(Problem(c, model2; bath = bath2, order = 2, rho0 = ρ0_imp2))
    occ = real(occupation(sol, N))
    @printf("  N=%d  %.2fs  ⟨n₁+n₂(0.4)⟩=%.8f  ED=%.8f  err=%.3e\n", N, t, occ, ref, abs(occ - ref))
    flush(stdout)
end
