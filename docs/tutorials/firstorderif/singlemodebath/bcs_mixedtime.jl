# =============================================================================
# firstorderif/singlemodebath/bcs_mixedtime.jl —— 一阶影响泛函 GTEMPO 混合
# Kadanoff 轴误差标度（BCS bath 版）
#
# 模型：BCS bath 下的单模 Anderson 杂质模型（U≠0）。BCS 影响泛函在 Nambu
# 表示下天然是双带，反常通道把两带耦在一起，故单一 IF（orbital=1）携带全部
# bath 影响；谱密度取 DiracDelta(ω, α) 的单模 BCS bath（bcsbath，实能隙 Δ），
# 截断 D=200/ϵ=1e-12 不绑定，测到的是离散化误差。
#
# 回答的问题：现有 GTEMPO 一阶 BCS IF 管线（QuAPI 核 + Kadanoff 装配 + 精确
# 系统传播子 + boundarycondition）在相互作用（U≠0）情形下，混合轮廓上
# gt/lt/gτ 的误差是否满足一阶 scaling（relerr ∝ O(δ)¹，δ = δt ∝ δτ）。
# τ 腿自身演化出相互作用热态 exp(-βH)。
#
# 方法：与 mixedtime.jl 相同，仅 bath/模型/参考解换成 BCS 版——
# - 参考解 = bcs_ed 的 ED（16 维）+ greater_lesser_ed_mixed（全耦合热态）+
#   gτ_ed；注意 bcs_ed 的 H = −ϵ_d(n̂₁+n̂₂)+Un₁n₂+…，模型侧用 AndersonIM
#  (ϵ_d = −ϵ_d) 对齐；
# - 参数沿用 BCS 测试验证过的组合（ϵ_d=0.8, ω=1, α=0.5, β=1, Δ=0.6）。
#
# 运行：julia --project=<GTEMPO 根> docs/tutorials/firstorderif/singlemodebath/bcs_mixedtime.jl
# 输出：bcs_mixedtime.png（标度图）、bcs_mixedtime.csv（数据）、终端误差表
#
# ---------------------------------------------------------------------------
# 实测结论（2026-11-03，ϵ_d=0.8, Δ=0.6, β=1, T=0.2, D=200/ϵ=1e-12，bcs_ed 参考）：
# 一阶 scaling 成立：gt p=1.04（逐点 1.07/1.03/1.02，非常干净）、lt p=1.21、
# gτ p=0.95（细格端 0.78/0.87 的 τ 腿端点假平台，同 bcs_imagtime 的误差曲线
# 过零效应；粗格端 1.25）。Nt=Nτ=3→24 同步加密；N=24 单点约 26 分钟。
# ---------------------------------------------------------------------------
# =============================================================================
include(joinpath(@__DIR__, "common.jl"))

const MODEL_T = 0.2    # 实时段总时长（τ 腿长度 = β）
const BCS_ϵ_d = 0.8    # 杂质能级（BCS 基准，沿测试验证组合；bcs_ed 约定取正号）
const Δ_bcs = 0.6      # BCS 实能隙

"""
    mixed_gf(N; trunc) -> [("greater gt", e), ("lesser lt", e), ("gτ", e)]

标准 GTEMPO 一阶 BCS IF 管线的混合 Kadanoff 轮廓：粗格点 + hybriddynamics
(orbital=1)（单一 IF 携带反常通道）+ 精确系统传播子 + boundarycondition；
误差对全模型 bcs_ed ED 参考。
"""
function mixed_gf(N::Int; trunc::TruncationScheme)
    δt, δτ = MODEL_T / N, β / N
    ts = collect(0:δt:MODEL_T)
    τs = collect(0:δτ:β)
    lat = GrassmannLattice(Nt=N, δt=δt, Nτ=N, δτ=δτ, contour=:Kadanoff, bands=2)
    model = AndersonIM(U=U, ϵ_d=-BCS_ϵ_d)
    mpsK = sysdynamics(lat, model, trunc=trunc)
    for band in 1:2
        mpsK = boundarycondition!(mpsK, lat, band=band, trunc=trunc)
    end
    bath = bcsbath(fermionicbath(DiracDelta(ω=ωb, α=αb), β=β), Δ=Δ_bcs)
    corr = correlationfunction(bath, lat)
    mpsI = hybriddynamics(lat, corr, orbital=1, trunc=trunc)
    gt, lt, gτ = gtltgτ_series(lat, mpsK, mpsI)
    H, a, adag, H0 = bcs_ed(U, BCS_ϵ_d, ω₀=ωb, α=αb, Δ=Δ_bcs)
    gt_ed, lt_ed = greater_lesser_ed_mixed(H, a, adag, ts, β)
    gτ_ref = gτ_ed(H, a, adag, τs, β)
    return [("greater gt", relerr(gt, gt_ed)), ("lesser lt", relerr(lt, lt_ed)),
            ("gτ", relerr(gτ, gτ_ref))]
end

function main(; Ns=(3, 6, 12, 24), D=200, ϵ=1.0e-12)
    println("=" ^ 92)
    @printf("GTEMPO 一阶 IF 混合 Kadanoff 误差标度（BCS bath）—— Anderson 杂质 (U=%.1f, ϵ_d=%.2f)，单模 bath×2，\n", U, BCS_ϵ_d)
    @printf("DiracDelta(ω=%.1f, α=%.1f) + bcsbath(Δ=%.1f)，β=%.1f，T=%.1f，Nt=Nτ 同步加密，D_bond=%d，ϵ=%.0e\n",
            ωb, αb, Δ_bcs, β, MODEL_T, D, ϵ)
    println("=" ^ 92)
    return scaling_main(mixed_gf; Ns=Ns, dtof=N -> MODEL_T / N, nlabel="N",
                        csvfile=joinpath(@__DIR__, "bcs_mixedtime.csv"),
                        pngfile=joinpath(@__DIR__, "bcs_mixedtime.png"),
                        title="GTEMPO first-order IF mixed (Kadanoff) error scaling\n" *
                              "Anderson IM (U=$U, ϵ_d=$BCS_ϵ_d, β=$β, T=$MODEL_T), BCS bath (Δ=$Δ_bcs)",
                        xlabel="δt = T / N  (δτ = β / N)",
                        ylabel="relative error of ĝ>/ĝ</g(τ) series", D=D, ϵ=ϵ)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
