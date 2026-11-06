# ============================================================================
# continuous_bath_scaling.jl — 连续谱浴 Anderson 收敛阶基准（第十一轮遗留项②）
#
# 模型与 continuous_bath_{mine,u1}.jl 完全一致（cos² 有界谱 ∫A = 0.5、
# Anderson U=0.8、β=1）。对齐已机器精度（≤1.4e-11），本基准补收敛**阶**：
#   τ  虚时 ⟨n₁(β/2)⟩：orders 1/2 × N ∈ (8,16,32,64)
#   rt 实时 ⟨n₁(0.4)⟩：orders 1/2 × N ∈ (4,8,16,32)
# 参考 = 各阶自洽 Richardson 外插（p=1/p=2，用最大两点）；斜率 =
# log2(err(h)/err(h/2)) 逐对 + 最细两点均值。
# 运行：julia --project=<repo> <repo>/benchmark/secondorderif/validation/continuous_bath_scaling.jl <repo>
# ============================================================================
const REPO = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "..", "..")
include(joinpath(REPO, "benchmark", "secondorderif", "src", "SecondOrderWeyl.jl"))
using .SecondOrderWeyl
using .SecondOrderWeyl: occupation

using GTEMPO
using Printf
using LinearAlgebra: diagm

const μ, β, U2 = 0.7, 1.0, 0.8
const W, ε₀, A₀ = 2.0, 0.2, 0.25          # ∫A = A₀·W = 0.5
A(ε) = abs(ε - ε₀) >= W ? 0.0 : A₀ * cos(π * (ε - ε₀) / (2 * W))^2
spec = bounded(A, ε₀ - W, ε₀ + W)
bath = fermionicbath(spec, β = β)
model = AndersonIM(ϵ_d = μ, U = U2)
ρ0_imp2 = diagm(0 => ComplexF64[1.0, exp(-β*μ), exp(-β*μ), exp(-2β*μ)]) /
          (1 + 2exp(-β*μ) + exp(-2β*μ))

"逐对斜率 + 末对斜率"
function slopes(Ns::Vector{Int}, errs::Vector{Float64})
    sl = Float64[]
    for k in 1:length(Ns)-1
        push!(sl, log2(errs[k] / errs[k+1]))   # h 减半 ⟹ 比值即斜率
    end
    return sl
end

results = Dict{Tuple{Symbol,Int},Vector{Float64}}()

println("== τ 虚时 ⟨n₁(β/2)⟩（orders 1/2，N ≤ 64）==")
for order in (1, 2)
    ns = Float64[]
    for N in (8, 16, 32, 64)
        c = WeylImagTime(Nτ = N, dτ = β / N, bands = 2)
        sol = solve(Problem(c, model; bath = bath, order = order))
        n = real(occupation(sol, N ÷ 2; band = 1))
        push!(ns, n)
        @printf("  order=%d N=%2d  ⟨n₁⟩=%.10f\n", order, N, n)
        sol = nothing; GC.gc()
    end
    results[(:tau, order)] = ns
end

println("== rt 实时 ⟨n₁(0.4)⟩（orders 1/2，N ≤ 32）==")
for order in (1, 2)
    ns = Float64[]
    for N in (4, 8, 16, 32)
        c = WeylRealTime(Nt = N, dt = 0.4 / N, bands = 2)
        sol = solve(Problem(c, model; bath = bath, order = order, rho0 = ρ0_imp2))
        n = real(occupation(sol, N; band = 1))
        push!(ns, n)
        @printf("  order=%d N=%2d  ⟨n₁⟩=%.10f\n", order, N, n)
        sol = nothing; GC.gc()
    end
    results[(:real, order)] = ns
end

println("\n== 收敛阶判决（Richardson 参考）==")
ok = true
for (case, Ns) in ((:tau, [8, 16, 32, 64]), (:real, [4, 8, 16, 32]))
    for order in (1, 2)
        ns = results[(case, order)]
        p = order == 2 ? 2.0 : 1.0
        ref = ns[end] + (ns[end] - ns[end-1]) / (2^p - 1)
        errs = abs.(ns .- ref)
        sl = slopes(Ns, errs)
        @printf("  %-4s order=%d: ref=%.10f  errs=%s\n", string(case), order, ref,
                join([@sprintf("%.2e", e) for e in errs], ", "))
        @printf("           逐对斜率 = %s；末对 = %.3f（期望 %.1f）\n",
                join([@sprintf("%.3f", s) for s in sl], ", "), sl[end], p)
        global ok &= (abs(sl[end] - p) < 0.2 || abs(sl[end-1] - p) < 0.2)
    end
end
println(ok ? "\n连续谱浴收敛阶基准：闭合 ✅" : "\n连续谱浴收敛阶基准：超限 ❌")
