# =============================================================================
# firstorderif/singlemodebath/mixedtime.jl —— 一阶影响泛函 GTEMPO 混合 Kadanoff 轴误差标度
#
# 回答的问题：现有 GTEMPO 一阶 IF 管线（QuAPI 核 + Kadanoff 装配 + 精确系统
# 传播子 + boundarycondition）在相互作用（U≠0、双带）情形下，混合轮廓上
# gt/lt/gτ 的误差是否满足一阶 scaling（relerr ∝ O(δ)¹，δ = δt ∝ δτ）。
#
# 方法：
# - 固定实时段长 T 与 β，Nt = Nτ = N ∈ {3..24} 同步加密（δt = T/N，δτ = β/N）；
# - τ 腿自身演化出相互作用热态 exp(-βH)，参考解 = 同一热态的全模型 ED
#  （greater_lesser_ed_mixed + gτ_ed，见 common.jl）；
# - 误差 = 全时序相对误差（gt、lt、gτ 各一条序列）。
#
# 运行：julia --project=<GTEMPO 根> docs/tutorials/firstorderif/singlemodebath/mixedtime.jl
# 输出：mixedtime.png（标度图）、mixedtime.csv（数据）、终端误差表
# =============================================================================
include(joinpath(@__DIR__, "common.jl"))

const MODEL_T = 0.2    # 实时段总时长（τ 腿长度 = β）

"""
    mixed_gf(N; trunc) -> [("greater gt", e), ("lesser lt", e), ("gτ", e)]

标准 GTEMPO 一阶 IF 管线的混合 Kadanoff 轮廓：粗格点 + 双带 fillband 复制
IF + 精确系统传播子 + boundarycondition；误差对全模型 ED 参考。
"""
function mixed_gf(N::Int; trunc::TruncationScheme)
    δt, δτ = MODEL_T / N, β / N
    ts = collect(0:δt:MODEL_T)
    τs = collect(0:δτ:β)
    lat = GrassmannLattice(Nt=N, δt=δt, Nτ=N, δτ=δτ, contour=:mixed, bands=2)
    bath = fermionicbath(DiracDelta(ω=ωb, α=αb), β=β)
    mpsK, Is = fermionic_setup(lat, bath, AndersonIM(U=U, ϵ_d=ϵ_d), trunc)
    gt, lt, gτ = gtltgτ_series(lat, mpsK, Is...)
    H, a, adag, H0 = singlemode_ed(μ=ϵ_d, U=U, bathspecs=[(ωb, αb), (ωb, αb)])
    gt_ed, lt_ed = greater_lesser_ed_mixed(H, a, adag, ts, β)
    gτ_ref = gτ_ed(H, a, adag, τs, β)
    return [("greater gt", relerr(gt, gt_ed)), ("lesser lt", relerr(lt, lt_ed)),
            ("gτ", relerr(gτ, gτ_ref))]
end

function main(; Ns=(3, 6, 12, 24), D=200, ϵ=1.0e-12)
    println("=" ^ 92)
    @printf("GTEMPO 一阶 IF 混合 Kadanoff 误差标度 —— Anderson 杂质 (U=%.1f, ϵ_d=%.2f)，单模 bath×2，\n", U, ϵ_d)
    @printf("DiracDelta(ω=%.1f, α=%.1f)，β=%.1f，T=%.1f，Nt=Nτ 同步加密，D_bond=%d，ϵ=%.0e\n",
            ωb, αb, β, MODEL_T, D, ϵ)
    println("=" ^ 92)
    return scaling_main(mixed_gf; Ns=Ns, dtof=N -> MODEL_T / N, nlabel="N",
                        csvfile=joinpath(@__DIR__, "mixedtime.csv"),
                        pngfile=joinpath(@__DIR__, "mixedtime.png"),
                        title="GTEMPO first-order IF mixed (Kadanoff) error scaling\n" *
                              "Anderson IM (U=$U, ϵ_d=$ϵ_d, β=$β, T=$MODEL_T), single-mode bath ×2",
                        xlabel="δt = T / N  (δτ = β / N)",
                        ylabel="relative error of ĝ>/ĝ</g(τ) series", D=D, ϵ=ϵ)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
