# ============================================================================
# u0_factor_rt.jl — 实时轮廓 U=0 因子化端到端检验（第十轮，2026-11-06）
#
# 用法：julia --project=<repo> benchmark/secondorderif/validation/u0_factor_rt.jl <repo>
#
# 同 u0_factor_mixed.jl 的本征锚，实时版（ChronoKeldysh K 路径）：
#   ① ⟨n₁⟩ = ⟨n₂⟩；② ⟨n₁+n₂⟩ = 2×单带 rt；③ Z₂/Z₁² = 1；
#   ④ 绝对值 vs ED（初态 = 杂质真空 ⊗ 热浴 β=1 —— 真空浴会被 Pauli 阻断
#      给出恒零占据，不能做参考）。
# ============================================================================
include(joinpath(ARGS[1], "benchmark", "secondorderif", "src", "SecondOrderWeyl.jl"))
using .SecondOrderWeyl
using .SecondOrderWeyl: occupation
include(joinpath(ARGS[1], "benchmark", "secondorderif", "validation", "ed.jl"))
using .SingleModeED

using GTEMPO
using Printf
using LinearAlgebra: dot, diagm, kron, tr, Diagonal

μ, ω, α, β, T = 0.7, 1.0, 0.5, 1.0, 0.4
bath = fermionicbath(DiracDelta(ω = ω, α = α), β = β)
model2 = AndersonIM(ϵ_d = μ, U = 0.0)

# ED：ρ₀ = 杂质真空 ⊗ 热浴（bit 序 = band1, band2, bath1, bath2，band1 = 最高位）
H2, dd2, _ = SingleModeED.model(μ = μ, U = 0.0, bathspecs = [(ω, α), (ω, α)])
hb = ComplexF64[0 0; 0 ω]
ρbath = exp(-β * hb) / tr(exp(-β * hb))
ρimp = Matrix{ComplexF64}([1.0 0.0; 0.0 0.0])
# ED：ρ₀ = 杂质真空 ⊗ 热浴。ed.jl 约定 mode m ⟺ bit m−1：杂质 = 低位
# （mode 1,2）、浴 = 高位（mode 3,4）；kron 第一因子在最高位 ⟹ 顺序
# [bath2, bath1, band2, band1]。
ρ0_ed = real.(kron(ρbath, ρbath, ρimp, ρimp))
occ_ref, _, _ = SingleModeED.observables(H2, dd2, [T], ρ0_ed)
ref_total = 2 * real(occ_ref[1])
@printf("ED（杂质真空⊗热浴 β=1，t=%.2f）：⟨n₁+n₂⟩ = %.8f（逐带 %.8f）\n\n", T, ref_total, real(occ_ref[1]))

fit_slope(hs, errs) = begin
    x = log.(hs); y = log.(errs); n = length(hs)
    (n * dot(x, y) - sum(x) * sum(y)) / (n * dot(x, x) - sum(x)^2)
end

for order in (1, 2)
    @printf("--- order = %d ---\n", order)
    hs = Float64[]; errs = Float64[]
    for N in (2, 4)
        h = T / N
        c = WeylRealTime(Nt = N, dt = h, bands = 2)
        sol = solve(Problem(c, model2; bath = bath, order = order))
        n1 = real(occupation(sol, N; band = 1)); n2 = real(occupation(sol, N; band = 2))
        ntot = real(occupation(sol, N))
        c1m = WeylRealTime(Nt = N, dt = h, bands = 1)
        sol1 = solve(Problem(c1m, ToulouseIM(ϵ_d = μ); bath = bath, order = order))
        nsingle = real(occupation(sol1, N))
        zr = real(sol.Z) / real(sol1.Z)^2
        @printf("  Nt=%d h=%.3f ⟨n₁⟩=%.7f ⟨n₂⟩=%.7f |n₁−n₂|=%.1e total=%.7f err(ED)=%.3e err(2×单带)=%.1e Z₂/Z₁²=%.8f\n",
                N, h, n1, n2, abs(n1 - n2), ntot, abs(ntot - ref_total),
                abs(ntot - 2 * nsingle), zr)
        push!(hs, h); push!(errs, abs(ntot - ref_total))
    end
    @printf("  >>> order=%d err(ED) 斜率 = %.3f\n\n", order, fit_slope(hs, errs))
end
println("实时 U=0 因子化检验完成")
