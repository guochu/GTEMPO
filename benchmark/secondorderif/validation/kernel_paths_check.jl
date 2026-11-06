# ============================================================================
# kernel_paths_check.jl — 实时核双路径 + 一般浴端到端校验（完整版模块）
#
# ① Dirac-δ：精确路径（闭式浴迹）vs QuAPI 谱路径 逐块对照（应 ≤1e-10）；
# ② Dirac-δ：QuAPI 谱路径端到端（Problem 直接喂谱路径核）vs ED；
# ③ 双 δ（DiscreteSpectrum，真两模浴）端到端 vs ED——一般浴路径的总验证：
#    τ ⟨n⟩ 收敛 + 实时 ⟨n(t)⟩，核走 QuAPI 谱路径（δ-精确求积）。
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

kquapi_kernel(bath, c) =
    HybridizationKernel(c, realtime_kernel_quapi(bath, c), Float64(α_weight(bath)))

const μ, ω, α, β = 0.7, 1.0, 0.5, 1.0

println("① Dirac-δ 双路径逐块对照：")
bath1 = GTEMPO.fermionicbath(GTEMPO.DiracDelta(ω = ω, α = α), β = β)
c4 = WeylRealTime(Nt = 4, dt = 0.1, bands = 1)
kexact = hybridization(bath1, c4)
# 直接调谱路径（绕过分派）
kquapi = HybridizationKernel(c4, realtime_kernel_quapi(bath1, c4), Float64(α_weight(bath1)))
maxd = 0.0
for (b1, b2) in ((:+, :+), (:+, :-), (:-, :+), (:-, :-))
    d = maximum(abs.(kexact.blocks[(b1, b2)] .- kquapi.blocks[(b1, b2)]))
    global maxd = max(maxd, d)
    @printf("   块 (%s,%s): max|Δ| = %.2e\n", b1, b2, d)
end
@printf("   >>> 最大偏差 = %.2e %s\n\n", maxd, maxd < 1e-10 ? "OK" : "✗")

println("② Dirac-δ 谱路径端到端 vs ED（实时 Toulouse）：")
model1 = ToulouseIM(ϵ_d = μ)
H1, dd1, H01 = SingleModeED.model(μ = μ, U = 0.0, bathspecs = [(ω, α)])
ρ0_imp = diagm(0 => ComplexF64[1.0, exp(-β*μ)]) / (1 + exp(-β*μ))
ρfull1 = SingleModeED.thermal(H01, β)
occ_ed, _, _ = SingleModeED.observables(H1, dd1, [0.4], ρfull1)
for N in (2, 4)
    c = WeylRealTime(Nt = N, dt = 0.4 / N, bands = 1)
    prob = Problem(c, model1; kernel = kquapi_kernel(bath1, c), order = 2, rho0 = ρ0_imp)
    sol = solve(prob)
    occ = real(occupation(sol, N))
    @printf("   N=%d  ⟨n(0.4)⟩=%.8f  ED=%.8f  err=%.3e\n", N, occ, real(occ_ed[1]), abs(occ - occ_ed[1]))
end

println("③ 双 δ 浴（DiscreteSpectrum）端到端 vs ED——一般浴路径总验证：")
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
