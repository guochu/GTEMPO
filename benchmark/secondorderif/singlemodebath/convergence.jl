# ============================================================================
# singlemodebath/convergence.jl — 两带 Anderson 模型（U≠0，每带一个单模浴）
#                                 二阶 Weyl-GTEMPO 收敛基准
#
# ① 虚时平衡 ⟨n₁+n₂⟩（β=1）：order 1/2，N ∈ {4,8,16}，ED 参考（2×逐带）
# ② 实时 ⟨n₁+n₂(t=0.4)⟩（因子化热初态）：order 1/2，N ∈ {8,16,32}，ED 参考
#    （ChronoKeldysh 时序直铺后 rt K 免 permute!；U≠0 实时离散误差常数较大，
#     纯二阶渐近区 N≥8：err/h² → 1.9e-2）
# 输出：results/anderson_tau.csv、results/anderson_realtime.csv
# 运行：julia --project=<repo> benchmark/secondorderif/singlemodebath/convergence.jl <repo>
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

const μ, ω, α, β, U2 = 0.7, 1.0, 0.5, 1.0, 0.8
const bath = fermionicbath(DiracDelta(ω = ω, α = α), β = β)
const model = AndersonIM(ϵ_d = μ, U = U2)
const H2, dd2, H02 = SingleModeED.model(μ = μ, U = U2, bathspecs = [(ω, α), (ω, α)])
const occ2_ed, _ = SingleModeED.imag_observables(H2, dd2, [0.0], β)
const ed_total_τ = 2 * occ2_ed
const ρ0_imp2 = diagm(0 => ComplexF64[1.0, exp(-β*μ), exp(-β*μ), exp(-β*(2μ + U2))]) /
                (1 + 2exp(-β*μ) + exp(-β*(2μ + U2)))

fit_slope(hs, errs) = begin
    x = log.(hs); y = log.(errs); n = length(hs)
    (n * dot(x, y) - sum(x) * sum(y)) / (n * dot(x, x) - sum(x)^2)
end

function run_case(io, case::String, orders, Ns; realtime::Bool)
    data = Dict{Int,Vector{Tuple{Float64,Float64}}}()
    for order in orders
        errs = Float64[]; hs = Float64[]
        for N in Ns
            h = realtime ? 0.4 / N : β / N
            c = realtime ? WeylRealTime(Nt = N, dt = h, bands = 2) :
                           WeylImagTime(Nτ = N, dτ = h, bands = 2)
            sol = realtime ? solve(Problem(c, model; bath = bath, order = order, rho0 = ρ0_imp2)) :
                             solve(Problem(c, model; bath = bath, order = order))
            j = realtime ? N : N ÷ 2
            occ = real(occupation(sol, j))
            if realtime
                ρfull = SingleModeED.thermal(H02, β)
                occ_ed, _, _ = SingleModeED.observables(H2, dd2, [0.4], ρfull)
                ref = 2 * real(occ_ed[1])
            else
                ref = 2 * real(occ2_ed)
            end
            e = abs(occ - ref)
            push!(errs, e); push!(hs, h)
            println(io, "$case,$order,$N,$h,$(occ),$ref,$e")
            @printf("  %s order=%d N=%2d h=%.4f ⟨n₁+n₂⟩=%.7f ref=%.7f err=%.3e\n",
                    case, order, N, h, occ, ref, e)
        end
        data[order] = collect(zip(hs, errs))
        @printf("  >>> %s order=%d 斜率 = %.3f\n\n", case, order, fit_slope(hs, errs))
    end
    return data
end

io = open(joinpath(OUT, "anderson_tau.csv"), "w")
println(io, "case,order,N,h,value,reference,abs_err")
@printf("① 虚时 ⟨n₁+n₂⟩（β=1，U=%.1f，ED = %.7f）\n", U2, 2real(occ2_ed))
τ_data = run_case(io, "imag-anderson2", (1, 2), (4, 8, 16); realtime = false)
close(io)

io = open(joinpath(OUT, "anderson_realtime.csv"), "w")
println(io, "case,order,N,h,value,reference,abs_err")
println("② 实时 ⟨n₁+n₂(0.4)⟩（因子化热初态）")
rt_data = run_case(io, "real-anderson2", (1, 2), (8, 16, 32); realtime = true)
close(io)

open(joinpath(OUT, "anderson_slopes.csv"), "w") do f
    println(f, "case,order,slope")
    for (tag, d) in (("tau", τ_data), ("realtime", rt_data)), (order, v) in d
        println(f, "$tag,$order,$(fit_slope(first.(v), last.(v)))")
    end
end
println("数据 → $OUT")
