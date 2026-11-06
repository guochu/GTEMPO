# ============================================================================
# u0_factor_mixed.jl — 混合轮廓 U=0 因子化端到端检验（第十轮，2026-11-06）
#
# 用法：julia --project=<repo> benchmark/secondorderif/validation/u0_factor_mixed.jl <repo>
#
# 本征锚（无跨框架比较）：U=0 双带 = 两个独立单带 ⟹
#   ① ⟨n₁⟩ = ⟨n₂⟩（同浴对称，机器精度）；
#   ② ⟨n₁+n₂⟩ = 2×单带 mixed（逐位）；
#   ③ Z₂/Z₁² = 1（K/IF 因子化；IF 标量在 Z 里 vacuum 归一相消）；
#   ④ order 2 斜率 ≈ 2（vs ED 热值，时间不变）。
# 配套：k_factor_mixed.jl（K₂ vs K₁∧K₁ 逐支撑楔积判决）、
#       u0_factor_rt.jl（实时同款）、if_band_invariance.jl（IF₁ 带嵌入探针）。
# ============================================================================
include(joinpath(ARGS[1], "benchmark", "secondorderif", "src", "SecondOrderWeyl.jl"))
using .SecondOrderWeyl
using .SecondOrderWeyl: occupation
include(joinpath(ARGS[1], "benchmark", "secondorderif", "validation", "ed.jl"))
using .SingleModeED

using GTEMPO
using Printf
using LinearAlgebra: dot

μ, ω, α, β, T = 0.7, 1.0, 0.5, 1.0, 0.4
bath = fermionicbath(DiracDelta(ω = ω, α = α), β = β)
model2 = AndersonIM(ϵ_d = μ, U = 0.0)
H2, dd2, _ = SingleModeED.model(μ = μ, U = 0.0, bathspecs = [(ω, α), (ω, α)])
ρc = SingleModeED.thermal(H2, β)
occ_ref, _, _ = SingleModeED.observables(H2, dd2, [T], ρc)
ref_total = 2 * real(occ_ref[1])
@printf("ED 热值（U=0，时间不变）：⟨n₁+n₂⟩ = %.8f（逐带 %.8f）\n\n", ref_total, real(occ_ref[1]))

fit_slope(hs, errs) = begin
    x = log.(hs); y = log.(errs); n = length(hs)
    (n * dot(x, y) - sum(x) * sum(y)) / (n * dot(x, x) - sum(x)^2)
end

function run_order(order)
    @printf("--- order = %d ---\n", order)
    hs = Float64[]; errs = Float64[]
    worst = (0.0, "")
    for N in (2, 4, 8)
        h = T / N
        c = WeylMixedTime(Nτ = N, dτ = β / N, Nt = N, dt = h, bands = 2)
        sol = solve(Problem(c, model2; bath = bath, order = order))
        n1 = real(occupation(sol, N; band = 1)); n2 = real(occupation(sol, N; band = 2))
        ntot = real(occupation(sol, N))
        c1m = WeylMixedTime(Nτ = N, dτ = β / N, Nt = N, dt = h, bands = 1)
        sol1 = solve(Problem(c1m, ToulouseIM(ϵ_d = μ); bath = bath, order = order))
        nsingle = real(occupation(sol1, N))
        zr = real(sol.Z) / real(sol1.Z)^2
        @printf("  N=%2d h=%.4f ⟨n₁⟩=%.7f ⟨n₂⟩=%.7f |n₁−n₂|=%.1e total=%.7f err(ED)=%.3e err(2×单带)=%.1e Z₂/Z₁²=%.8f\n",
                N, h, n1, n2, abs(n1 - n2), ntot, abs(ntot - ref_total),
                abs(ntot - 2 * nsingle), zr)
        push!(hs, h); push!(errs, abs(ntot - ref_total))
        for (tag, v) in (("对称", abs(n1 - n2)), ("vs2×单带", abs(ntot - 2 * nsingle)),
                         ("Z₂/Z₁²", abs(zr - 1)))
            v > first(worst) && (worst = (v, "N$N $tag"))
        end
    end
    @printf("  >>> order=%d err(ED) 斜率 = %.3f；最大内在偏差 %.1e（%s）\n\n",
            order, fit_slope(hs, errs), first(worst), last(worst))
end
run_order(1)
run_order(2)
println("混合 U=0 因子化检验完成")
