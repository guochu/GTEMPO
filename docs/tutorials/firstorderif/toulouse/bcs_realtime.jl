# =============================================================================
# firstorderif/toulouse/bcs_realtime.jl —— 一阶影响泛函 GTEMPO 实时演化误差标度
# （BCS bath 版）
#
# 模型：BCS bath 下的 Toulouse 杂质（U=0）。BCS 影响泛函在 Nambu 表示下
# 天然是双带（带 1/2 = 自旋/粒子-空穴），反常通道把两带耦在一起，故 U=0
# 也走双带 AndersonIM 管线；谱密度取 DiracDelta(ω, α) 的单模 BCS bath
#（bcsbath，实能隙 Δ），截断 D=200/ϵ=1e-12 不绑定，测到的是离散化误差。
#
# 回答的问题：现有 GTEMPO 一阶 BCS IF 管线（QuAPI 核 + 顶点位移装配 + 精确
# 系统传播子 + 可分热态初态）的实时 Keldysh 演化，误差是否满足一阶 scaling
#（relerr ∝ O(δt)¹）。
#
# 方法：与 realtime.jl 相同，仅 bath/模型/参考解换成 BCS 版——
# - 参考解 = 全模型（2 杂质轨道 + 1 BCS 模式，16 维）的 ED（test/util.jl 的
#   bcs_ed / greater_lesser_ed）；注意 bcs_ed 的 H = −ϵ_d(n̂₁+n̂₂)+…，
#   模型侧用 AndersonIM(ϵ_d = −ϵ_d) 对齐；
# - 参数沿用 BCS 测试验证过的组合（ϵ_d=0.8, ω=1, α=0.5, β=1, Δ=0.6）。
#
# 运行：julia --project=<GTEMPO 根> docs/tutorials/firstorderif/toulouse/bcs_realtime.jl
# 输出：bcs_realtime.png（标度图）、bcs_realtime.csv（数据）、终端误差表
#
# ---------------------------------------------------------------------------
# 实测结论（2026-11-03，ϵ_d=0.8, Δ=0.6, β=1, T=0.4, D=200/ϵ=1e-12，bcs_ed 参考）：
# U=0 一阶 scaling 成立：gt/lt 相对误差 1.23e-3 → 7.3e-5（Nt=4→64），
# 逐点阶 1.03→1.01，拟合 p=1.02。BCS IF 比普通 bath 重得多（反常通道使键
# 维数变大：Nt=64 单点约 11 分钟、2.6 GB，独占 CPU 时）。对照：U=1 的同管线
# 曾因 ↓↓ 正常项装配（(a↓,ā↓)·(−Δ)，核转置）出现 δt 无关系统偏差，2026-11-03
# 已修复（QuAPI contourindex.jl 的 BCS getindex 改指数系数语义），
# 修复后见 ../singlemodebath/bcs_realtime.jl（p=1.14/1.04）。
# ---------------------------------------------------------------------------
# =============================================================================
include(joinpath(@__DIR__, "common.jl"))

const MODEL_T = 0.4    # 实时演化总时长
const BCS_ϵ_d = 0.8    # 杂质能级（BCS 基准，沿测试验证组合；bcs_ed 约定取正号）
const Δ_bcs = 0.6      # BCS 实能隙

"""
    realtime_gf(Nt; trunc) -> [("greater gt", e_gt), ("lesser lt", e_lt)]

标准 GTEMPO 一阶 BCS IF 管线的实时 Keldysh 演化：:Keldysh 双带格点 +
hybriddynamics(orbital=1)（单一 IF 携带反常通道）+ 精确系统传播子 +
boundarycondition + 热态初态；误差对全模型 ED 参考（bcs_ed）。
"""
function realtime_gf(Nt::Int; trunc::TruncationScheme)
    δt = MODEL_T / Nt
    ts = collect(0:δt:MODEL_T)
    lat = GrassmannLattice(N=Nt, δt=δt, contour=:Keldysh, bands=2)
    model = AndersonIM(U=0.0, ϵ_d=-BCS_ϵ_d)
    mpsK = sysdynamics(lat, model, trunc=trunc)
    for band in 1:2
        mpsK = boundarycondition!(mpsK, lat, band=band, trunc=trunc)
    end
    mpsK = systhermalstate!(mpsK, lat, model, trunc=trunc, β=β)
    bath = bcsbath(fermionicbath(DiracDelta(ω=ωb, α=αb), β=β), Δ=Δ_bcs)
    corr = correlationfunction(bath, lat)
    mpsI = hybriddynamics(lat, corr, orbital=1, trunc=trunc)
    cache = environments(lat, mpsK, mpsI)
    gt = [-im * cached_greater(lat, i, mpsK, mpsI, cache=cache) for i in 1:lat.k]
    lt = [-im * cached_lesser(lat, i, mpsK, mpsI, cache=cache) for i in 1:lat.k]
    H, a, adag, H0 = bcs_ed(0.0, BCS_ϵ_d, ω₀=ωb, α=αb, Δ=Δ_bcs)
    gt_ed, lt_ed = greater_lesser_ed(H, a, adag, H0, ts, β)
    return [("greater gt", relerr(gt, gt_ed)), ("lesser lt", relerr(lt, lt_ed))]
end

function main(; Ns=(4, 8, 16, 32, 64), D=200, ϵ=1.0e-12)
    println("=" ^ 92)
    @printf("GTEMPO 一阶 IF 实时误差标度（BCS bath）—— Toulouse 模型 (U=0, ϵ_d=%.2f)，\n", BCS_ϵ_d)
    @printf("单模 DiracDelta(ω=%.1f, α=%.1f) + bcsbath(Δ=%.1f)，β=%.1f，T=%.1f，D_bond=%d，ϵ=%.0e\n",
            ωb, αb, Δ_bcs, β, MODEL_T, D, ϵ)
    println("=" ^ 92)
    return scaling_main(realtime_gf; Ns=Ns, dtof=N -> MODEL_T / N, nlabel="Nt",
                        csvfile=joinpath(@__DIR__, "bcs_realtime.csv"),
                        pngfile=joinpath(@__DIR__, "bcs_realtime.png"),
                        title="GTEMPO first-order IF real-time error scaling\n" *
                              "Toulouse model (U=0, ϵ_d=$BCS_ϵ_d, β=$β, T=$MODEL_T), BCS bath (Δ=$Δ_bcs)",
                        xlabel="δt = T / Nt",
                        ylabel="relative error of ĝ>/ĝ< series", D=D, ϵ=ϵ)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
