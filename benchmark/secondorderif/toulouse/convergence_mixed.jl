# ============================================================================
# toulouse/convergence_mixed.jl — 混合时间轴（Kadanoff–Baym L 形）收敛基准
#
# Toulouse 模型 + 单模浴：混合轮廓 ⟨n(t=0.4)⟩，τ 段与 rt 段同步加密
# （Nτ = Nt = N）。ED 参考 = 耦合 H 热态演化（混合轮廓的单算符期望
# = 热平衡值，时间不变——ED 数值验证到 5e-18）。
# 输出：results/toulouse_mixed.csv
# 运行：julia --project=<repo> benchmark/secondorderif/toulouse/convergence_mixed.jl <repo>
# ============================================================================
const REPO = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "..", "..")
include(joinpath(REPO, "benchmark", "secondorderif", "src", "SecondOrderWeyl.jl"))
using .SecondOrderWeyl
using .SecondOrderWeyl: occupation
include(joinpath(REPO, "benchmark", "secondorderif", "validation", "ed.jl"))
using .SingleModeED

using GTEMPO
using Printf
using LinearAlgebra: diagm, dot

const OUT = joinpath(@__DIR__, "results")
isdir(OUT) || mkpath(OUT)

const μ, ω, α, β, T = 0.7, 1.0, 0.5, 1.0, 0.4
const bath = fermionicbath(DiracDelta(ω = ω, α = α), β = β)
const model = ToulouseIM(ϵ_d = μ)
const H1, dd1, H01 = SingleModeED.model(μ = μ, U = 0.0, bathspecs = [(ω, α)])
const ρc = SingleModeED.thermal(H1, β)          # 耦合 H 热态（混合轮廓参考）
const occ_ref, _, _ = SingleModeED.observables(H1, dd1, [T], ρc)
const ref = real(occ_ref[1])

fit_slope(hs, errs) = begin
    x = log.(hs); y = log.(errs); n = length(hs)
    (n * dot(x, y) - sum(x) * sum(y)) / (n * dot(x, x) - sum(x)^2)
end

io = open(joinpath(OUT, "toulouse_mixed.csv"), "w")
println(io, "case,order,N,h,value,reference,abs_err")
@printf("混合轮廓 ⟨n(t=0.4)⟩（ED 参考 = %.7f，时间不变）\n", ref)
data = Dict{Int,Vector{Tuple{Float64,Float64}}}()
for order in (1, 2)
    errs = Float64[]; hs = Float64[]
    for N in (4, 8, 16)
        h = T / N
        c = WeylMixedTime(Nτ = N, dτ = β / N, Nt = N, dt = h, bands = 1)
        sol = solve(Problem(c, model; bath = bath, order = order))
        n = real(occupation(sol, N))
        e = abs(n - ref)
        push!(errs, e); push!(hs, h)
        println(io, "mixed-toulouse,$order,$N,$h,$(n),$ref,$e")
        @printf("  order=%d N=%2d h=%.4f ⟨n⟩=%.7f err=%.3e\n", order, N, h, n, e)
    end
    data[order] = collect(zip(hs, errs))
    @printf("  >>> order=%d 斜率 = %.3f\n\n", order, fit_slope(hs, errs))
end
close(io)

open(joinpath(OUT, "toulouse_mixed_slopes.csv"), "w") do f
    println(f, "case,order,slope")
    for (order, v) in data
        println(f, "mixed,$order,$(fit_slope(first.(v), last.(v)))")
    end
end
println("数据 → $OUT")
