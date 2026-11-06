# ============================================================================
# continuous_bath_u1.jl — 连续谱浴 Anderson 模型（U1GTEMPO2 侧）
#
# 谱：A(ε) = A₀·cos²(π(ε−ε₀)/(2W))，ε ∈ [ε₀−W, ε₀+W]，∫A = A₀·W
# 模型：两自旋 SIAM（energy=(μ,μ), U），两自旋同浴；观测量 = 逐自旋 ⟨n↑⟩
# 输出：results/continuous_bath_u1.csv（case,order,N,h,value）
# 运行：julia --project=~/Documents/Missile/U1GTEMPO2 continuous_bath_u1.jl <out.csv>
# ============================================================================
using U1GTEMPO2
using Printf

const μ, β, U2 = 0.7, 1.0, 0.8
const W, ε₀, A₀ = 2.0, 0.2, 0.25          # ∫A = A₀·W = 0.5
A(ε) = abs(ε - ε₀) >= W ? 0.0 : A₀ * cos(π * (ε - ε₀) / (2 * W))^2
bath = FermionicBath(A; support = (ε₀ - W, ε₀ + W), beta = β, mu = 0.0)

# 因子化热初态（实时）：基序 n = 2·n_up + n_down
Z0 = 1 + 2exp(-β * μ) + exp(-2β * μ)
p_th = (1.0 / Z0, exp(-β * μ) / Z0, exp(-β * μ) / Z0, exp(-2β * μ) / Z0)

out = open(ARGS[1], "w")
println(out, "case,order,N,h,n_up")

println("== U1GTEMPO2 连续浴 Anderson：虚时 ⟨n↑⟩ ==")
for (order, Ns) in ((1, (4, 8, 16)), (2, (4, 8, 16)))
    for N in Ns
        c = ImagTime(Nτ = N, dτ = β / N)
        prob = SIAMProblem(c, SIAM(energy = (μ, μ), U = U2); bath = bath, order = order)
        sol = solve(prob)
        n = real(occupation(sol, N ÷ 2; spin = :up))
        println(out, "imag,$order,$N,$(β/N),$n")
        @printf("  order=%d N=%2d  ⟨n↑⟩=%.8f\n", order, N, n); flush(stdout)
        sol = nothing; prob = nothing; GC.gc()
    end
end

println("== U1GTEMPO2 连续浴 Anderson：实时 ⟨n↑(0.4)⟩ ==")
for order in (1, 2)
    for N in (2, 4, 8)
        c = RealTime(Nt = N, dt = 0.4 / N)
        prob = SIAMProblem(c, SIAM(energy = (μ, μ), U = U2); bath = bath, order = order,
                           initial_probabilities = p_th)
        sol = solve(prob)
        n = real(occupation(sol, N; spin = :up))
        println(out, "real,$order,$N,$(0.4/N),$n")
        @printf("  order=%d N=%2d  ⟨n↑(0.4)⟩=%.8f\n", order, N, n); flush(stdout)
        sol = nothing; prob = nothing; GC.gc()
    end
end
close(out)
println("done → ", ARGS[1])
