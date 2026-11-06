# ============================================================================
# twopoint_m2_frontier.jl — M≥2 跨格奇插入两点的 α→0 层判决工具（第十轮）
#
# 判决史（2026-11-06）：
#   ① 初判（修复前）：mixed α→0 精确而 α≠0 错 ~1e-1、rt「α→0 也错 ~3e-1」
#     ——后查明 rt 的「错」是 ED 参考初态用错（rt 物理初态 = rho0 杂质真空 ⊗
#     热浴；mixed 有 τ 热化段 ⟹ 耦合热态）。α→0 值只测 K 的 mask-0 常数项
#     （slot 变量补偿不进入），不构成对补偿机制的检验。
#   ② 真根因（修复后确认）：K 装配的 ysv↔site 带映射互换（rt/mixed，M≥2
#     带交换；τ 路径走 _code_from_ys 无此病）+ 补偿缺双向跨带交叉项
#     ⊕_{b≠b'} c_b·p_{b'}（楔积 inv 的逐元素实现）。修复后：
#     mixed 全网格 vs ED 闭合（mixed_twopoint_check.jl ③）、rt α→0 机器精度。
# 用法：julia --project=<repo> benchmark/secondorderif/validation/twopoint_m2_frontier.jl <repo>
# ============================================================================
include(joinpath(ARGS[1], "benchmark", "secondorderif", "src", "SecondOrderWeyl.jl"))
using .SecondOrderWeyl
using .SecondOrderWeyl: raw_correlation
include(joinpath(ARGS[1], "benchmark", "secondorderif", "validation", "ed.jl"))
using .SingleModeED
using GTEMPO
using Printf
using LinearAlgebra: tr, kron

const μ, ω, β, T = 0.7, 1.0, 1.0, 0.4
const bath0 = fermionicbath(DiracDelta(ω = ω, α = 1e-12), β = β)

# ED（decoupled，α≈0 与 mine 的浴一致）
H, d, _ = SingleModeED.model(μ = μ, U = 0.8, bathspecs = [(ω, 0.0), (ω, 0.0)])
ρc_mixed = SingleModeED.thermal(H, β)      # mixed：τ 段热化 ⟹ 耦合热态
hb0 = ComplexF64[0 0; 0 ω]
ρbath0 = exp(-β * hb0) / tr(exp(-β * hb0))
ρimp0 = Matrix{ComplexF64}([1.0 0.0; 0.0 0.0])
ρc_rt = real.(kron(ρbath0, ρbath0, ρimp0, ρimp0))   # rt：rho0 真空 ⊗ 热浴

for (label, mk, ρED) in (("mixed", (N, dt) -> WeylMixedTime(Nτ = N, dτ = β / N, Nt = N, dt = dt, bands = 2), ρc_mixed),
                         ("rt", (N, dt) -> WeylRealTime(Nt = N, dt = dt, bands = 2), ρc_rt))
    println("== $label M=2 α→0 ==")
    N = 4; dt = T / N
    sol = solve(Problem(mk(N, dt), AndersonIM(ϵ_d = μ, U = 0.8); bath = bath0, order = 2))
    for i in (1, 2, 4)
        nΔ, gt, lt = SingleModeED.observables(H, d, [i * dt], ρED)
        gΔ = gt[1] / (-1.0im)
        gmine = raw_correlation(sol, i, 0; kind = :greater)
        @printf("  greater(%d,0): %+.6f%+.6fi vs ED %+.6f%+.6fi |%.2e\n",
                i, real(gmine), imag(gmine), real(gΔ), imag(gΔ), abs(gmine - gΔ))
    end
    println()
end
println("判决完成")
