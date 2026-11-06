# ============================================================================
# continuous_bath_mine.jl — 连续谱浴 Anderson 模型（SecondOrderWeyl 侧）
#
# 与 continuous_bath_u1.jl 完全同一模型/参数（逐带对比 U1 的逐自旋）。
# 输出：results/continuous_bath_mine.csv（case,order,N,h,n_up）
# 运行：julia --project=<GTEMPO仓库> continuous_bath_mine.jl <仓库根> <out.csv>
# ============================================================================
const REPO = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "..", "..")
include(joinpath(REPO, "benchmark", "secondorderif", "src", "SecondOrderWeyl.jl"))
using .SecondOrderWeyl
using .SecondOrderWeyl: occupation
include(joinpath(REPO, "benchmark", "secondorderif", "validation", "ed.jl"))
using .SingleModeED   # （本脚本不用 ED，仅为保持与其它脚本一致的加载顺序）

using GTEMPO
using Printf
using LinearAlgebra: diagm

const μ, β, U2 = 0.7, 1.0, 0.8
const W, ε₀, A₀ = 2.0, 0.2, 0.25          # ∫A = A₀·W = 0.5
A(ε) = abs(ε - ε₀) >= W ? 0.0 : A₀ * cos(π * (ε - ε₀) / (2 * W))^2
spec = bounded(A, ε₀ - W, ε₀ + W)          # 连续有界谱
bath = fermionicbath(spec, β = β)
@printf("coupling_norm（quadgk）= %.8f（解析 A₀·W = %.8f）\n",
        coupling_norm(hybridization(bath, WeylImagTime(Nτ = 4, dτ = β / 4, bands = 1))), A₀ * W)

const model = AndersonIM(ϵ_d = μ, U = U2)
const ρ0_imp2 = diagm(0 => ComplexF64[1.0, exp(-β*μ), exp(-β*μ), exp(-2β*μ)]) /
                (1 + 2exp(-β*μ) + exp(-2β*μ))

out = open(ARGS[2], "w")
println(out, "case,order,N,h,n_up")

println("== SecondOrderWeyl 连续浴 Anderson：虚时 ⟨n↑⟩ ==")
for (order, Ns) in ((1, (4, 8, 16)), (2, (4, 8, 16)))
    for N in Ns
        c = WeylImagTime(Nτ = N, dτ = β / N, bands = 2)
        sol = solve(Problem(c, model; bath = bath, order = order))
        n = real(occupation(sol, N ÷ 2; band = 1))
        println(out, "imag,$order,$N,$(β/N),$n")
        @printf("  order=%d N=%2d  ⟨n₁⟩=%.8f\n", order, N, n); flush(stdout)
        sol = nothing; GC.gc()
    end
end

println("== SecondOrderWeyl 连续浴 Anderson：实时 ⟨n↑(0.4)⟩ ==")
for order in (1, 2)
    for N in (2, 4, 8)
        c = WeylRealTime(Nt = N, dt = 0.4 / N, bands = 2)
        sol = solve(Problem(c, model; bath = bath, order = order, rho0 = ρ0_imp2))
        n = real(occupation(sol, N; band = 1))
        println(out, "real,$order,$N,$(0.4/N),$n")
        @printf("  order=%d N=%2d  ⟨n₁(0.4)⟩=%.8f\n", order, N, n); flush(stdout)
        sol = nothing; GC.gc()
    end
end
close(out)
println("done → ", ARGS[2])
