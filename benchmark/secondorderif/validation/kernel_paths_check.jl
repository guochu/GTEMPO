# ============================================================================
# kernel_paths_check.jl — 实时核 QuAPI 统一路径 + 一般浴端到端校验（完整版模块）
#
# 实时核自 2026-11-07 统一走 QuAPI `Δt` 谱路径（QuAPI 对 Dirac-δ /
# DiscreteSpectrum 离散谱取精确极点装配，数值精确）。旧 Dirac-δ 闭式浴迹
# 特例 `realtime_kernel_trace` 已删除——删除前双路径互验 ≤4.1e-17
#（HISTORY 第八轮 #4 与第十二轮补遗 3）。
# ① Dirac-δ：QuAPI 路径端到端（Problem 默认分派）vs ED；
# ② 双 δ（DiscreteSpectrum，真两模浴）端到端 vs ED——一般浴路径的总验证：
#    τ ⟨n⟩ 收敛 + 实时 ⟨n(t)⟩。
# 运行：julia --project=. benchmark/secondorderif/validation/kernel_paths_check.jl <repo>
# ============================================================================
include(joinpath(ARGS[1], "benchmark", "secondorderif", "src", "SecondOrderWeyl.jl"))
using .SecondOrderWeyl
using .SecondOrderWeyl: occupation, greater, lesser
include(joinpath(ARGS[1], "benchmark", "secondorderif", "validation", "ed.jl"))
using .SingleModeED
using GTEMPO
using ImpurityModelBase: DiscreteSpectrum, spectrumvalues
using Printf
using LinearAlgebra: diagm

const μ, ω, α, β = 0.7, 1.0, 0.5, 1.0

println("① Dirac-δ QuAPI 路径端到端 vs ED（实时 Toulouse）：")
bath1 = GTEMPO.fermionicbath(GTEMPO.DiracDelta(ω = ω, α = α), β = β)
model1 = ToulouseIM(ϵ_d = μ)
H1, dd1, H01 = SingleModeED.model(μ = μ, U = 0.0, bathspecs = [(ω, α)])
ρ0_imp = diagm(0 => ComplexF64[1.0, exp(-β*μ)]) / (1 + exp(-β*μ))
ρfull1 = SingleModeED.thermal(H01, β)
occ_ed, _, _ = SingleModeED.observables(H1, dd1, [0.4], ρfull1)
for N in (2, 4)
    c = WeylRealTime(Nt = N, dt = 0.4 / N, bands = 1)
    prob = Problem(c, model1; bath = bath1, order = 2, rho0 = ρ0_imp)
    sol = solve(prob)
    occ = real(occupation(sol, N))
    @printf("   N=%d  ⟨n(0.4)⟩=%.8f  ED=%.8f  err=%.3e\n", N, occ, real(occ_ed[1]), abs(occ - occ_ed[1]))
end

println("② 双 δ 浴（DiscreteSpectrum）端到端 vs ED——一般浴路径总验证：")
ωs = [0.8, 1.4]; αs = [0.3, 0.2]
spec = DiscreteSpectrum(ωs, αs)
bath2 = GTEMPO.fermionicbath(spec, β = β)
model2 = ToulouseIM(ϵ_d = μ)
H2, dd2, H02 = SingleModeED.model(μ = μ, U = 0.0, bathspecs = [(ωs[1], αs[1]), (ωs[2], αs[2])],
                                  nbands = 1)
occ2_τ_ed, _ = SingleModeED.imag_observables(H2, dd2, [0.0], β)
@printf("   coupling_norm = %.4f (ED Σα = %.4f)   ED ⟨n⟩_τ = %.7f\n",
        α_weight(bath2), sum(spectrumvalues(spec)), real(occ2_τ_ed))
for N in (8, 16)
    cτ = WeylImagTime(Nτ = N, dτ = β / N, bands = 1)
    sol = solve(Problem(cτ, model2; bath = bath2, order = 2))
    n = real(occupation(sol, N ÷ 2))
    @printf("   τ N=%2d  ⟨n⟩=%.7f  err=%.3e\n", N, n, abs(n - real(occ2_τ_ed)))
end
ρfull2 = SingleModeED.thermal(H02, β)
occ2_rt_ed, _, _ = SingleModeED.observables(H2, dd2, [0.4], ρfull2)
for N in (2, 4)
    c = WeylRealTime(Nt = N, dt = 0.4 / N, bands = 1)
    sol = solve(Problem(c, model2; bath = bath2, order = 2, rho0 = ρ0_imp))
    occ = real(occupation(sol, N))
    @printf("   实时 N=%d  ⟨n(0.4)⟩=%.7f  ED=%.7f  err=%.3e\n", N, occ, real(occ2_rt_ed[1]),
            abs(occ - occ2_rt_ed[1]))
end
