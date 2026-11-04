# =============================================================================
# firstorderif/singlemodebath/bcs_realtime.jl —— 一阶影响泛函 GTEMPO 实时演化误差标度
# （BCS bath 版）
#
# 模型：BCS bath 下的单模 Anderson 杂质模型（U≠0）。BCS 影响泛函在 Nambu
# 表示下天然是双带（带 1/2 = 自旋/粒子-空穴），反常通道把两带耦在一起，
# 故单一 IF（orbital=1）即携带全部 bath 影响；谱密度取 DiracDelta(ω, α)
# 的单模 BCS bath（bcsbath，实能隙 Δ），截断 D=200/ϵ=1e-12 不绑定，测到的
# 是离散化误差。
#
# 回答的问题：现有 GTEMPO 一阶 BCS IF 管线（QuAPI 核 + 顶点位移装配 + 精确
# 系统传播子 + 可分热态初态）在相互作用（U≠0）情形下，实时 Keldysh 演化
# 的误差是否满足一阶 scaling（relerr ∝ O(δt)¹）。
#
# 方法：与 realtime.jl 相同，仅 bath/模型/参考解换成 BCS 版——
# - 参考解 = 全模型（2 杂质轨道 + 1 BCS 模式，16 维）的 ED（test/util.jl 的
#   bcs_ed / greater_lesser_ed）；注意 bcs_ed 的 H = −ϵ_d(n̂₁+n̂₂)+Un̂₁n̂₂+…，
#   模型侧用 AndersonIM(ϵ_d = −ϵ_d) 对齐；
# - 参数沿用 BCS 测试验证过的组合（ϵ_d=0.8, ω=1, α=0.5, β=1, Δ=0.6）。
#
# 运行：julia --project=<GTEMPO 根> docs/tutorials/firstorderif/singlemodebath/bcs_realtime.jl
# 输出：bcs_realtime.png（标度图）、bcs_realtime.csv（数据）、终端误差表
#
# ---------------------------------------------------------------------------
# 实测结论（2026-11-03，ϵ_d=0.8, Δ=0.6, β=1, T=0.4, D=200/ϵ=1e-12，bcs_ed 参考）：
# 一阶 scaling 成立：gt 5.2e-4 → 2.2e-5（Nt=4→64，逐点阶 1.31→1.03，p=1.14），
# lt 7.9e-4 → 4.4e-5（逐点阶 1.09→1.01，p=1.04）。
#
# 本基准曾揭露并定位一个装配 bug（2026-11-03 已修复）：旧 BCS IF 把 ↓↓ 正常项
# 实现为 (a↓,ā↓)·(−Δ)——对 Grassmann 单项式等价于 ↓↓ 核转置（canonical 应为
# (ā↓,a↓)·(+Δ)，同 QuAPI Dense 的 _delta_terms）。U=0 时 band-2 只贡献标量、
# 在 Z 归约中相消而不可见；U≠0 时给出 δt 无关、随演化时间累积的系统偏差。
# 修复前实测：gt 卡 1.51e-3（p=0.06）、lt 卡 4.90e-4（p=0.28）；D=200/400/800
# 与 Δ=0/0.6 对照均不变（排除截断/反常通道），同参数普通 bath 双 IF 路线干净。
# 修复在 QuAPI contourindex.jl 的 BCS getindex（指数系数语义）。
# 注意 BCS IF 比普通 bath 重得多（反常通道使键维数变大：Nt=64 单点约 11 分钟、
# 2.6 GB，独占 CPU 时）。
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
    model = AndersonIM(U=U, ϵ_d=-BCS_ϵ_d)
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
    H, a, adag, H0 = bcs_ed(U, BCS_ϵ_d, ω₀=ωb, α=αb, Δ=Δ_bcs)
    gt_ed, lt_ed = greater_lesser_ed(H, a, adag, H0, ts, β)
    return [("greater gt", relerr(gt, gt_ed)), ("lesser lt", relerr(lt, lt_ed))]
end

function main(; Ns=(4, 8, 16, 32, 64), D=200, ϵ=1.0e-12)
    println("=" ^ 92)
    @printf("GTEMPO 一阶 IF 实时误差标度（BCS bath）—— Anderson 杂质 (U=%.1f, ϵ_d=%.2f)，\n", U, BCS_ϵ_d)
    @printf("单模 DiracDelta(ω=%.1f, α=%.1f) + bcsbath(Δ=%.1f)，β=%.1f，T=%.1f，D_bond=%d，ϵ=%.0e\n",
            ωb, αb, Δ_bcs, β, MODEL_T, D, ϵ)
    println("=" ^ 92)
    return scaling_main(realtime_gf; Ns=Ns, dtof=N -> MODEL_T / N, nlabel="Nt",
                        csvfile=joinpath(@__DIR__, "bcs_realtime.csv"),
                        pngfile=joinpath(@__DIR__, "bcs_realtime.png"),
                        title="GTEMPO first-order IF real-time error scaling\n" *
                              "Anderson IM (U=$U, ϵ_d=$BCS_ϵ_d, β=$β, T=$MODEL_T), BCS bath (Δ=$Δ_bcs)",
                        xlabel="δt = T / Nt",
                        ylabel="relative error of ĝ>/ĝ< series", D=D, ϵ=ϵ)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
