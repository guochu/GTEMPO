# =============================================================================
# firstorderif/toulouse/bcs_mixedtime.jl —— 一阶影响泛函 GTEMPO 混合 Kadanoff 轴
# 误差标度（BCS bath 版）
#
# 模型：BCS bath 下的 Toulouse 杂质（U=0）。BCS 影响泛函在 Nambu 表示下天然
# 是双带，反常通道把两带耦在一起，故单一 IF（orbital=1）携带全部 bath 影响；
# 谱密度取 DiracDelta(ω, α) 的单模 BCS bath（bcsbath，实能隙 Δ），截断
# D=200/ϵ=1e-12 不绑定，测到的是离散化误差。
#
# 回答的问题：现有 GTEMPO 一阶 BCS IF 管线（QuAPI 核 + Kadanoff 装配 + 精确
# 系统传播子 + boundarycondition）在混合轮廓上，gt/lt/gτ 的误差是否满足
# 一阶 scaling（relerr ∝ O(δ)¹，δ = δt ∝ δτ）。τ 腿自身演化出相互作用热态
# exp(-βH)。
#
# 方法：与 mixedtime.jl 相同，仅 bath/模型/参考解换成 BCS 版——
# - 参考解 = bcs_ed 的 ED（16 维）+ greater_lesser_ed_mixed（全耦合热态）+
#   gτ_ed；注意 bcs_ed 的 H = −ϵ_d(n̂₁+n̂₂)+…，模型侧用 AndersonIM(ϵ_d = −ϵ_d)
#   对齐；
# - 参数沿用 BCS 测试验证过的组合（ϵ_d=0.8, ω=1, α=0.5, β=1, Δ=0.6）。
#
# 运行：julia --project=<GTEMPO 根> docs/tutorials/firstorderif/toulouse/bcs_mixedtime.jl
# 输出：bcs_mixedtime.png（标度图）、bcs_mixedtime.csv（数据）、终端误差表
#
# ---------------------------------------------------------------------------
# 实测结论（2026-11-03，ϵ_d=0.8, Δ=0.6, β=1, T=0.2, D=200/ϵ=1e-12，bcs_ed 参考）：
# 一阶 scaling 成立：gt p=1.07（逐点 1.13/1.06/1.03）、lt p=1.29、gτ p=1.59
#（Nt=Nτ=3→24 同步加密）。gτ 的细格端有假平台（N=12→24 阶数 0.28）——τ 腿
# 端点等时 GF 的误差曲线过零效应，与 bcs_imagtime 轴的诊断一致；粗格端
# 2.92/1.57 收敛很快。BCS 混合轴很重：N=24 单点约 21 分钟、14 GB 内存。
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
误差对同一热态的 bcs_ed ED 参考。
"""
function mixed_gf(N::Int; trunc::TruncationScheme)
    δt, δτ = MODEL_T / N, β / N
    ts = collect(0:δt:MODEL_T)
    τs = collect(0:δτ:β)
    lat = GrassmannLattice(Nt=N, δt=δt, Nτ=N, δτ=δτ, contour=:Kadanoff, bands=2)
    model = AndersonIM(U=0.0, ϵ_d=-BCS_ϵ_d)
    mpsK = sysdynamics(lat, model, trunc=trunc)
    for band in 1:2
        mpsK = boundarycondition!(mpsK, lat, band=band, trunc=trunc)
    end
    bath = bcsbath(fermionicbath(DiracDelta(ω=ωb, α=αb), β=β), Δ=Δ_bcs)
    corr = correlationfunction(bath, lat)
    mpsI = hybriddynamics(lat, corr, orbital=1, trunc=trunc)
    gt, lt, gτ = gtltgτ_series(lat, mpsK, mpsI)
    H, a, adag, H0 = bcs_ed(0.0, BCS_ϵ_d, ω₀=ωb, α=αb, Δ=Δ_bcs)
    gt_ed, lt_ed = greater_lesser_ed_mixed(H, a, adag, ts, β)
    gτ_ref = gτ_ed(H, a, adag, τs, β)
    return [("greater gt", relerr(gt, gt_ed)), ("lesser lt", relerr(lt, lt_ed)),
            ("gτ", relerr(gτ, gτ_ref))]
end

function main(; Ns=(3, 6, 12, 24), D=200, ϵ=1.0e-12)
    println("=" ^ 92)
    @printf("GTEMPO 一阶 IF 混合 Kadanoff 误差标度（BCS bath）—— Toulouse 模型 (U=0, ϵ_d=%.2f)，\n", BCS_ϵ_d)
    @printf("单模 DiracDelta(ω=%.1f, α=%.1f) + bcsbath(Δ=%.1f)，β=%.1f，T=%.1f，Nt=Nτ 同步加密，D_bond=%d，ϵ=%.0e\n",
            ωb, αb, Δ_bcs, β, MODEL_T, D, ϵ)
    println("=" ^ 92)
    return scaling_main(mixed_gf; Ns=Ns, dtof=N -> MODEL_T / N, nlabel="N",
                        csvfile=joinpath(@__DIR__, "bcs_mixedtime.csv"),
                        pngfile=joinpath(@__DIR__, "bcs_mixedtime.png"),
                        title="GTEMPO first-order IF mixed (Kadanoff) error scaling\n" *
                              "Toulouse model (U=0, ϵ_d=$BCS_ϵ_d, β=$β, T=$MODEL_T), BCS bath (Δ=$Δ_bcs)",
                        xlabel="δt = T / N  (δτ = β / N)",
                        ylabel="relative error of ĝ>/ĝ</g(τ) series", D=D, ϵ=ϵ)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
