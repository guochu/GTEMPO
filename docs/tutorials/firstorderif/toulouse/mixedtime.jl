# =============================================================================
# firstorderif/toulouse/mixedtime.jl —— 一阶影响泛函 GTEMPO 混合 Kadanoff 轴误差标度
#
# 回答的问题：现有 GTEMPO 一阶 IF 管线（QuAPI 核 + Kadanoff 装配 + 精确系统
# 传播子 + boundarycondition）在混合轮廓上，gt/lt/gτ 的误差是否满足一阶
# scaling（relerr ∝ O(δ)¹，δ = δt ∝ δτ）。
#
# 方法：
# - 固定实时段长 T 与 β，Nt = Nτ = N ∈ {3..24} 同步加密（δt = T/N，δτ = β/N）；
# - τ 腿自身演化出相互作用热态 exp(-βH)，参考解 = 同一热态的 ED
#  （greater_lesser_ed_mixed + gτ_ed，见 common.jl）；
# - 误差 = 全时序相对误差（gt、lt、gτ 各一条序列）。
#
# 运行：julia --project=<GTEMPO 根> docs/tutorials/firstorderif/toulouse/mixedtime.jl
# 输出：mixedtime.png（标度图）、mixedtime.csv（数据）、终端误差表
# =============================================================================
include(joinpath(@__DIR__, "common.jl"))

const MODEL_T = 0.2    # 实时段总时长（τ 腿长度 = β）

"""
    mixed_gf(N; trunc) -> [("greater gt", e), ("lesser lt", e), ("gτ", e)]

标准 GTEMPO 一阶 IF 管线的混合 Kadanoff 轮廓：粗格点 + hybriddynamics +
精确系统传播子 + boundarycondition；误差对同一热态的 ED 参考。
"""
function mixed_gf(N::Int; trunc::TruncationScheme)
    δt, δτ = MODEL_T / N, β / N
    ts = collect(0:δt:MODEL_T)
    τs = collect(0:δτ:β)
    lat = GrassmannLattice(Nt=N, δt=δt, Nτ=N, δτ=δτ, contour=:mixed)
    corr = correlationfunction(toulouse_bath(), lat)
    mpsI = hybriddynamics(lat, corr, trunc=trunc)
    mpsK = sysdynamics(lat, ToulouseIM(ϵ_d=ϵ_d), trunc=trunc)
    mpsK = boundarycondition!(mpsK, lat)
    gt, lt, gτ = gtltgτ_series(lat, mpsK, mpsI)
    H, a, adag, H0 = singlemode_ed(μ=ϵ_d, U=0, bathspecs=[(ωb, αb)])
    gt_ed, lt_ed = greater_lesser_ed_mixed(H, a, adag, ts, β)
    gτ_ref = gτ_ed(H, a, adag, τs, β)
    return [("greater gt", relerr(gt, gt_ed)), ("lesser lt", relerr(lt, lt_ed)),
            ("gτ", relerr(gτ, gτ_ref))]
end

function main(; Ns=(3, 6, 12, 24), D=200, ϵ=1.0e-12)
    println("=" ^ 92)
    @printf("GTEMPO 一阶 IF 混合 Kadanoff 误差标度 —— Toulouse 模型 (U=0, ϵ_d=%.2f)，\n", ϵ_d)
    @printf("单模 bath DiracDelta(ω=%.1f, α=%.1f)，β=%.1f，T=%.1f，Nt=Nτ 同步加密，D_bond=%d，ϵ=%.0e\n",
            ωb, αb, β, MODEL_T, D, ϵ)
    println("=" ^ 92)
    return scaling_main(mixed_gf; Ns=Ns, dtof=N -> MODEL_T / N, nlabel="N",
                        csvfile=joinpath(@__DIR__, "mixedtime.csv"),
                        pngfile=joinpath(@__DIR__, "mixedtime.png"),
                        title="GTEMPO first-order IF mixed (Kadanoff) error scaling\n" *
                              "Toulouse model (U=0, ϵ_d=$ϵ_d, β=$β, T=$MODEL_T), single-mode bath",
                        xlabel="δt = T / N  (δτ = β / N)",
                        ylabel="relative error of ĝ>/ĝ</g(τ) series", D=D, ϵ=ϵ)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
