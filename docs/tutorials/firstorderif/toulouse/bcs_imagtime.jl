# =============================================================================
# firstorderif/toulouse/bcs_imagtime.jl —— 一阶影响泛函 GTEMPO 虚时演化误差标度
# （BCS bath 版）
#
# 模型：BCS bath 下的 Toulouse 杂质（U=0）。BCS 影响泛函在 Nambu 表示下天然
# 是双带（带 1/2 = 自旋/粒子-空穴），反常通道把两带耦在一起，故 U=0 也走
# 双带 AndersonIM 管线；谱密度取 DiracDelta(ω, α) 的单模 BCS bath（bcsbath，
# 实能隙 Δ），截断 D=200/ϵ=1e-12 不绑定，测到的是离散化误差。
#
# 回答的问题：现有 GTEMPO 一阶 BCS IF 管线（QuAPI 核 + τ 环装配 + 精确系统
# 传播子 + boundarycondition）的虚时 g(τ)，误差是否满足一阶 scaling
#（relerr ∝ O(δτ)¹）。τ 环自身演化出相互作用热态 exp(-βH)，无需热态初态。
#
# 方法：与 imagtime.jl 相同，仅 bath/模型/参考解换成 BCS 版——
# - 参考解 = bcs_ed 的 ED（16 维）+ gτ_ed；注意 bcs_ed 的 H = −ϵ_d(n̂₁+n̂₂)+…，
#   模型侧用 AndersonIM(ϵ_d = −ϵ_d) 对齐；
# - 参数沿用 BCS 测试验证过的组合（ϵ_d=0.8, ω=1, α=0.5, β=1, Δ=0.6）。
#
# 运行：julia --project=<GTEMPO 根> docs/tutorials/firstorderif/toulouse/bcs_imagtime.jl
# 输出：bcs_imagtime.png（标度图）、bcs_imagtime.csv（数据）、终端误差表
#
# ---------------------------------------------------------------------------
# 实测结论（2026-11-03，ϵ_d=0.8, Δ=0.6, β=1, D=200/ϵ=1e-12，bcs_ed 参考）：
# 一阶 scaling 渐近成立：g(τ) 相对误差 6.6e-4 → 2.1e-5（Nτ=4→128），
# 逐点阶 0.65 → 0.88 → 0.96（Nτ≥32 起）渐近趋近 1；拟合 p=0.87 被 Nτ=8-16
# 的误差曲线过零拉低（误差含异号分量，在 δτ≈0.1 附近相消，使 Nτ=8 意外偏小、
# Nτ=8→16 阶数只剩 0.07）。D=400 对照与 Δ=0 对照均不改变误差（排除截断与
# 反常通道嫌疑）；端点 τ=0/τ=β（等时 GF，离散化敏感处）误差最大，τ 环中部
# 收敛更快（~3 阶）。对照：Anderson（U=1）版本 orders 0.69→1.01 单调干净
# （见 ../singlemodebath/bcs_imagtime.jl）。
# ---------------------------------------------------------------------------
# =============================================================================
include(joinpath(@__DIR__, "common.jl"))

const BCS_ϵ_d = 0.8    # 杂质能级（BCS 基准，沿测试验证组合；bcs_ed 约定取正号）
const Δ_bcs = 0.6      # BCS 实能隙

"""
    imagtime_gf(Nτ; trunc) -> [("gτ", e)]

标准 GTEMPO 一阶 BCS IF 管线的虚时 τ 环：粗格点 + hybriddynamics(orbital=1)
（单一 IF 携带反常通道）+ 精确系统传播子 + boundarycondition；误差对
bcs_ed 的 ED 参考。
"""
function imagtime_gf(Nτ::Int; trunc::TruncationScheme)
    δτ = β / Nτ
    τs = collect(0:δτ:β)
    lat = GrassmannLattice(N=Nτ, δτ=δτ, contour=:imag, bands=2)
    model = AndersonIM(U=0.0, ϵ_d=-BCS_ϵ_d)
    mpsK = sysdynamics(lat, model, trunc=trunc)
    for band in 1:2
        mpsK = boundarycondition!(mpsK, lat, band=band, trunc=trunc)
    end
    bath = bcsbath(fermionicbath(DiracDelta(ω=ωb, α=αb), β=β), Δ=Δ_bcs)
    corr = correlationfunction(bath, lat)
    mpsI = hybriddynamics(lat, corr, orbital=1, trunc=trunc)
    g = gτ_series(lat, mpsK, mpsI)
    H, a, adag, H0 = bcs_ed(0.0, BCS_ϵ_d, ω₀=ωb, α=αb, Δ=Δ_bcs)
    ref = gτ_ed(H, a, adag, τs, β)
    return [("gτ", relerr(g, ref))]
end

function main(; Ns=(4, 8, 16, 32, 64, 128), D=200, ϵ=1.0e-12)
    println("=" ^ 92)
    @printf("GTEMPO 一阶 IF 虚时误差标度（BCS bath）—— Toulouse 模型 (U=0, ϵ_d=%.2f)，\n", BCS_ϵ_d)
    @printf("单模 DiracDelta(ω=%.1f, α=%.1f) + bcsbath(Δ=%.1f)，β=%.1f，D_bond=%d，ϵ=%.0e\n",
            ωb, αb, Δ_bcs, β, D, ϵ)
    println("=" ^ 92)
    return scaling_main(imagtime_gf; Ns=Ns, dtof=N -> β / N, nlabel="Nτ",
                        csvfile=joinpath(@__DIR__, "bcs_imagtime.csv"),
                        pngfile=joinpath(@__DIR__, "bcs_imagtime.png"),
                        title="GTEMPO first-order IF imaginary-time error scaling\n" *
                              "Toulouse model (U=0, ϵ_d=$BCS_ϵ_d, β=$β), BCS bath (Δ=$Δ_bcs)",
                        xlabel="δτ = β / Nτ",
                        ylabel="relative error of g(τ) series", D=D, ϵ=ϵ)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
