# Phase 0 golden fixtures dump（secondorderif / cell 离散化锚点）。
#
# 为什么需要第二组 golden：旧 GTEMPO 主栈的实时/混合轮廓是 **vertex 离散化**
# （传播子耦合相邻顶点、IF 用 Δt 节点值），新栈（cell automaton + Chrono）是
# **cell 离散化**（每格全宽源、单元积分核）——两者物理近似差 O(δ)，不可 1e-10
# 对拍。extensions/secondorderif（SecondOrderWeyl）已经是 cell 离散 + Chrono
# 排布 + Weyl/一阶 cell K + 逐带 cell IF，且与 U1GTEMPO2 端到端对齐 ≤1.4e-11
#（INTEGRATION.md §1/§8）——它是新栈的离散化同构锚。
#
# 本脚本依赖旧栈（SecondOrderWeyl include 旧 GTEMPO），必须在换轨接线前运行：
#   julia --project=. test/fixtures/dump_golden_weyl.jl
#
# 断言对象（check_golden Phase 2+ 用新栈重建同算例，1e-10 档）：
#   Z、occupation 全端点（逐支逐带）、greater/lesser/matsubara 全端点。

using GTEMPO
using GTEMPO: branches
using ImpurityModelBase
using LinearAlgebra: diagm, I
using Random

const SOW = include(joinpath(@__DIR__, "..", "..", "extensions", "secondorderif", "src", "SecondOrderWeyl.jl"))
using .SOW
using .SOW: WeylRealTime, WeylImagTime, WeylMixedTime, Problem, solve,
    occupation, greater, lesser, matsubara

Random.seed!(12354)

const OUTDIR = joinpath(@__DIR__, "port_golden")
mkpath(OUTDIR)

jsonnum(x::Float64) = isnan(x) ? "\"NaN\"" : isinf(x) ? (x > 0 ? "\"Inf\"" : "\"-Inf\"") : string(x)
jstr(s) = "\"" * s * "\""
jcplx(z) = "[" * jsonnum(real(z)) * "," * jsonnum(imag(z)) * "]"
jcplxarr(v) = "[" * join(jcplx.(collect(v)), ",") * "]"
jpair(k, v) = Pair{String,String}(k, v)
jobj(items::Pair{String,<:Any}...) = "{" * join((p.first * ":" * string(p.second) for p in items), ",") * "}"
write_json(path, pairs::Pair{String,<:Any}...) = open(path, "w") do io
    write(io, jobj(pairs...))
end

"单一算例 dump"
function dump_weyl_case(name; contour, model, bath, order::Int, rho0, series)
    prob = Problem(contour, model; bath=bath, order=order, rho0=rho0)
    sol = solve(prob)
    fields = Pair{String,<:Any}[]
    push!(fields, jpair("kind", jstr("weyl_" * String(nameof(typeof(contour))))))
    push!(fields, jpair("order", string(order)))
    push!(fields, jpair("Z", jcplx(sol.Z)))
    for (label, values) in series(sol)
        push!(fields, jpair(label, jcplxarr(values)))
    end
    path = joinpath(OUTDIR, "port_golden_$name.json")
    open(path, "w") do io
        write(io, jobj(fields...))
    end
    println("WROTE ", path, "  (", filesize(path), " bytes)")
end

μ, ω, α, β = 0.7, 1.0, 0.5, 1.0
bath = fermionicbath(DiracDelta(ω=ω, α=α), β=β)

# ρ₀：杂质热态对角（validation rt_fastpath 同式）
ρ0_1 = diagm(0 => ComplexF64[1.0, exp(-β * μ)]) / (1 + exp(-β * μ))
U2 = 0.8
ρ0_2 = diagm(0 => ComplexF64[1.0, exp(-β * μ), exp(-β * μ), exp(-β * (2μ + U2))]) /
       (1 + 2exp(-β * μ) + exp(-β * (2μ + U2)))

# ---- 实时 order 2：Toulouse 单带 ----
dump_weyl_case("weyl_real_o2_toulouse"; contour=WeylRealTime(Nt=3, dt=0.15, bands=1),
               model=ToulouseIM(ϵ_d=μ), bath=bath, order=2, rho0=ρ0_1,
               series = sol -> (
                   ("occupation", [occupation(sol, j) for j in 1:3]),
                   ("greater", [greater(sol, j, 0) for j in 1:3]),
                   ("lesser", [lesser(sol, j, 0) for j in 1:3]),
               ))

# ---- 实时 order 2：Anderson 双带（逐带占据） ----
dump_weyl_case("weyl_real_o2_anderson"; contour=WeylRealTime(Nt=3, dt=0.15, bands=2),
               model=AndersonIM(U=U2, ϵ_d=μ), bath=bath, order=2, rho0=ρ0_2,
               series = sol -> (
                   ("occupation", [occupation(sol, j) for j in 1:3]),
                   ("occupation_b1", [occupation(sol, j; band=1) for j in 1:3]),
                   ("occupation_b2", [occupation(sol, j; band=2) for j in 1:3]),
               ))

# ---- 混合 order 2：Toulouse（τ 占据 + rt 占据） ----
dump_weyl_case("weyl_mixed_o2_toulouse"; contour=WeylMixedTime(Nτ=2, dτ=0.3, Nt=2, dt=0.15, bands=1),
               model=ToulouseIM(ϵ_d=μ), bath=bath, order=2, rho0=ρ0_1,
               series = sol -> (
                   ("occupation_tau", [occupation(sol, j; branch=:τ) for j in 1:2]),
                   ("occupation_rt", [occupation(sol, j; branch=:+) for j in 1:2]),
                   ("greater", [greater(sol, j, 0) for j in 1:2]),
                   ("lesser", [lesser(sol, j, 0) for j in 1:2]),
               ))

# ---- 虚时 order 2：Toulouse（matsubara） ----
dump_weyl_case("weyl_imag_o2_toulouse"; contour=WeylImagTime(Nτ=4, dτ=0.2, bands=1),
               model=ToulouseIM(ϵ_d=μ), bath=bath, order=2, rho0=nothing,
               series = sol -> (
                   ("occupation", [occupation(sol, j) for j in 1:4]),
                   ("matsubara", [matsubara(sol, j, 0) for j in 1:4]),
               ))

# ---- 实时 order 1：Toulouse（anti-Wick + Lie 锚） ----
dump_weyl_case("weyl_real_o1_toulouse"; contour=WeylRealTime(Nt=3, dt=0.15, bands=1),
               model=ToulouseIM(ϵ_d=μ), bath=bath, order=1, rho0=ρ0_1,
               series = sol -> (
                   ("occupation", [occupation(sol, j) for j in 1:3]),
                   ("greater", [greater(sol, j, 0) for j in 1:3]),
               ))

println("weyl golden dump complete: ", length(readdir(OUTDIR)), " files total in ", OUTDIR)
