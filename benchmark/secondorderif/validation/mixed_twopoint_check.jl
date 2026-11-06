# ============================================================================
# mixed_twopoint_check.jl — 混合轮廓两点函数 ED 验证（第十轮补口）
#
# 用法：julia --project=<repo> benchmark/secondorderif/validation/mixed_twopoint_check.jl <repo>
#
# ① M=1 mixed（Toulouse）：greater/lesser (i,j) 全网格 vs ED —— 全绿（闭合）。
#    含 j=0（ddag 在 (+,0)=τ→+ 转折点 / − 支闭合点）与 j≥1（− 支插入）——
#    2026-11-06 修复 − 支切位映射 off-by-one（Nτ+2Nt−boundary+1 → Nτ+2Nt−boundary）
#    后全网格通过。
# ② M=2 mixed（Anderson U≠0）：等时两点（i==j，同切位 f_k⊕2=0 无偏移）vs ED。
# ③ M=2 跨格奇插入（i≠j）：确认为已知限制——raw_correlation 应显式抛错
#    （判决链见 HISTORY 第十轮 39-41 / twopoint_m2_frontier.jl）。
# ED 参考（平衡态平移不变）：C(t_i,t_j) = C(t_i−t_j, 0)：
#   greater核 = ⟨d(t_i)ddag(t_j)⟩ = gt(Δ)/(−i)；lesser核 = lt(Δ)/(+i)；
#   等时（i==j）：greater核 = 1−n_th、lesser核 = n_th（0⁺ 约定）。
# ============================================================================
include(joinpath(ARGS[1], "benchmark", "secondorderif", "src", "SecondOrderWeyl.jl"))
using .SecondOrderWeyl
using .SecondOrderWeyl: raw_correlation
include(joinpath(ARGS[1], "benchmark", "secondorderif", "validation", "ed.jl"))
using .SingleModeED

using GTEMPO
using Printf

const μ, ω, α, β, T = 0.7, 1.0, 0.5, 1.0, 0.4
const bath = fermionicbath(DiracDelta(ω = ω, α = α), β = β)

function ed_ref(nbands, U, Δ::Real)
    H, d, _ = SingleModeED.model(μ = μ, U = U,
                                 bathspecs = [(ω, α) for _ in 1:nbands])
    ρc = SingleModeED.thermal(H, β)
    n_th, gt, lt = SingleModeED.observables(H, d, [Δ], ρc)
    return real(n_th[1]), gt[1] / (-1.0im), lt[1] / (1.0im)
end

println("== ① M=1 mixed（Toulouse，N=4 全网格）==")
function check_m1()
    N = 4; dt = T / N
    n_th, _, _ = ed_ref(1, 0.0, 0.0)
    c = WeylMixedTime(Nτ = N, dτ = β / N, Nt = N, dt = dt, bands = 1)
    sol = solve(Problem(c, ToulouseIM(ϵ_d = μ); bath = bath, order = 2))
    worst = 0.0; at = ""
    # API 承诺区间 i ≥ j（greater/lesser 的 push 序 0⁺ 语义只在 d 晚于 ddag
    # 的一侧定义；i<j 返回轮廓时序积类量，不在本文档语义内）
    for i in 1:N, j in 0:i
        Δ = (i - j) * dt
        nΔ, gΔ, lΔ = j == i ? (n_th, 1.0 - n_th, n_th) : ed_ref(1, 0.0, Δ)
        for (k, mine, ref) in ((:g, raw_correlation(sol, i, j; kind = :greater), gΔ),
                               (:l, raw_correlation(sol, i, j; kind = :lesser), lΔ))
            e = abs(mine - ref)
            e > worst && (worst = e; at = "$k($i,$j)")
        end
    end
    @printf("  (i≥j) 全网格 × {g,l} 最大偏差 = %.2e（%s）%s\n\n", worst, at,
            worst < 2e-3 ? "—— 闭合 ✅" : "—— 超限 ❌")
end
check_m1()

println("== ② M=2 mixed（Anderson U≠0，等时 i==j）==")
function check_m2_eq()
    for N in (4, 8)
        dt = T / N
        n_th, _, _ = ed_ref(2, 0.8, 0.0)
        c = WeylMixedTime(Nτ = N, dτ = β / N, Nt = N, dt = dt, bands = 2)
        sol = solve(Problem(c, AndersonIM(ϵ_d = μ, U = 0.8); bath = bath, order = 2))
        worst = 0.0
        for i in 1:N
            gmine = raw_correlation(sol, i, i; kind = :greater)
            lmine = raw_correlation(sol, i, i; kind = :lesser)
            worst = max(worst, abs(gmine - (1.0 - n_th)), abs(lmine - n_th))
        end
        @printf("  N=%d 等时 {greater核=1−n, lesser核=n} 最大偏差 = %.2e%s\n",
                N, worst, worst < 2e-3 ? " ✅" : " ❌")
    end
end
check_m2_eq()

println()
println("== ③ M=2 跨格奇插入（i≠j，i≥j 全网格）—— ysv+双向交叉修复后的端到端判决 ==")
function check_m2_full()
    N = 4; dt = T / N
    c = WeylMixedTime(Nτ = N, dτ = β / N, Nt = N, dt = dt, bands = 2)
    sol = solve(Problem(c, AndersonIM(ϵ_d = μ, U = 0.8); bath = bath, order = 2))
    worst = 0.0; at = ""
    n_th, _, _ = ed_ref(2, 0.8, 0.0)
    for i in 1:N, j in 0:i
        Δ = (i - j) * dt
        if j == i
            gΔ = 1.0 - n_th   # 等时 0⁺ 约定：greater核 = 1−n、lesser核 = n
            lΔ = n_th
        else
            _, gΔ, lΔ = ed_ref(2, 0.8, Δ)
        end
        gmine = raw_correlation(sol, i, j; kind = :greater)
        lmine = raw_correlation(sol, i, j; kind = :lesser)
        eg = abs(gmine - gΔ); el = abs(lmine - lΔ)
        eg > worst && (worst = eg; at = "g($i,$j)")
        el > worst && (worst = el; at = "l($i,$j)")
    end
    @printf("  N=4 (i≥j) 全网格 × {g,l} 最大偏差 = %.2e（%s）%s\n\n", worst, at,
            worst < 2e-3 ? "—— 闭合 ✅" : "—— 超限 ❌")
end
check_m2_full()
println("\n完成")
