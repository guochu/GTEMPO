# ============================================================================
# toulouse/convergence.jl — Toulouse 模型（无相互作用，U=0，单带 + 单模浴）
#                          二阶 Weyl-GTEMPO 收敛基准
#
# ① 虚时平衡 ⟨n⟩（β=1）：order 1/2，N ∈ {5,10,20,40}，ED 参考
# ② 实时 ⟨n(t=0.4)⟩（因子化热初态）：order 1/2，N ∈ {2,4,8}，ED 参考
# 输出：results/toulouse_tau.csv、results/toulouse_realtime.csv
# 运行：julia --project=<repo> benchmark/secondorderif/toulouse/convergence.jl <repo>
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

const μ, ω, α, β = 0.7, 1.0, 0.5, 1.0
const bath = fermionicbath(DiracDelta(ω = ω, α = α), β = β)
const model = ToulouseIM(ϵ_d = μ)
const H1, dd1, H01 = SingleModeED.model(μ = μ, U = 0.0, bathspecs = [(ω, α)])
const occ1_ed, _ = SingleModeED.imag_observables(H1, dd1, [0.0], β)
const ρ0_imp = diagm(0 => ComplexF64[1.0, exp(-β*μ)]) / (1 + exp(-β*μ))

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
            c = realtime ? WeylRealTime(Nt = N, dt = h, bands = 1) :
                           WeylImagTime(Nτ = N, dτ = h, bands = 1)
            sol = realtime ? solve(Problem(c, model; bath = bath, order = order, rho0 = ρ0_imp)) :
                             solve(Problem(c, model; bath = bath, order = order))
            j = realtime ? N : N ÷ 2
            occ = real(occupation(sol, j))
            if realtime
                ρfull = SingleModeED.thermal(H01, β)
                occ_ed, _, _ = SingleModeED.observables(H1, dd1, [0.4], ρfull)
                ref = real(occ_ed[1])
            else
                ref = real(occ1_ed)
            end
            e = abs(occ - ref)
            push!(errs, e); push!(hs, h)
            println(io, "$case,$order,$N,$h,$(occ),$ref,$e")
            @printf("  %s order=%d N=%2d h=%.4f ⟨n⟩=%.7f ref=%.7f err=%.3e\n",
                    case, order, N, h, occ, ref, e)
        end
        data[order] = collect(zip(hs, errs))
        @printf("  >>> %s order=%d 斜率 = %.3f\n\n", case, order, fit_slope(hs, errs))
    end
    return data
end

io = open(joinpath(OUT, "toulouse_tau.csv"), "w")
println(io, "case,order,N,h,value,reference,abs_err")
@printf("① 虚时 ⟨n⟩（β=1，ED = %.7f）\n", real(occ1_ed))
τ_data = run_case(io, "imag-toulouse", (1, 2), (5, 10, 20, 40); realtime = false)
close(io)

io = open(joinpath(OUT, "toulouse_realtime.csv"), "w")
println(io, "case,order,N,h,value,reference,abs_err")
println("② 实时 ⟨n(0.4)⟩（因子化热初态）")
rt_data = run_case(io, "real-toulouse", (1, 2), (2, 4, 8); realtime = true)
close(io)

open(joinpath(OUT, "toulouse_slopes.csv"), "w") do f
    println(f, "case,order,slope")
    for (tag, d) in (("tau", τ_data), ("realtime", rt_data)), (order, v) in d
        h = real(occ1_ed)
        println(f, "$tag,$order,$(fit_slope(first.(v), last.(v)))")
    end
end
println("数据 → $OUT")
