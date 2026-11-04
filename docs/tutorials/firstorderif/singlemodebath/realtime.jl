# =============================================================================
# firstorderif/singlemodebath/realtime.jl —— 一阶影响泛函 GTEMPO 实时演化误差标度
#
# 回答的问题：现有 GTEMPO 一阶 IF 管线（QuAPI 核 + 顶点位移装配 + 精确系统
# 传播子 + 可分热态初态）在相互作用（U≠0、双带）情形下，实时 Keldysh 演化
# 的误差是否满足一阶 scaling（relerr ∝ O(δt)¹）。
#
# 方法：
# - 固定总时长 T，Nt ∈ {4..64} 加密格点，δt = T/Nt；
# - 参考解 = 全模型的 ED（见 common.jl），初始态 = exp(-βH0) 可分热态，
#   与 GTEMPO 侧 systhermalstate! 同一约定；
# - 误差 = 全时序相对误差 norm(gt − gt_ed)/norm(gt_ed)（lesser 同）。
#
# 运行：julia --project=<GTEMPO 根> docs/tutorials/firstorderif/singlemodebath/realtime.jl
# 输出：realtime.png（标度图）、realtime.csv（数据）、终端误差表
# =============================================================================
include(joinpath(@__DIR__, "common.jl"))

const MODEL_T = 0.5    # 实时演化总时长

"""
    realtime_gf(Nt; trunc) -> [("greater gt", e_gt), ("lesser lt", e_lt)]

标准 GTEMPO 一阶 IF 管线的实时 Keldysh 演化：粗格点 + 双带 fillband 复制
IF + 精确系统传播子 + boundarycondition + 热态初态；误差对全模型 ED 参考。
"""
function realtime_gf(Nt::Int; trunc::TruncationScheme)
    δt = MODEL_T / Nt
    ts = collect(0:δt:MODEL_T)
    lat = GrassmannLattice(N=Nt, δt=δt, contour=:real, bands=2)
    bath = fermionicbath(DiracDelta(ω=ωb, α=αb), β=β)
    mpsK, Is = fermionic_setup(lat, bath, AndersonIM(U=U, ϵ_d=ϵ_d), trunc, β=β)
    gt, lt = gtlt_series(lat, mpsK, Is...)
    H, a, adag, H0 = singlemode_ed(μ=ϵ_d, U=U, bathspecs=[(ωb, αb), (ωb, αb)])
    gt_ed, lt_ed = greater_lesser_ed(H, a, adag, H0, ts, β)
    return [("greater gt", relerr(gt, gt_ed)), ("lesser lt", relerr(lt, lt_ed))]
end

function main(; Ns=(4, 8, 16, 32, 64), D=200, ϵ=1.0e-12)
    println("=" ^ 92)
    @printf("GTEMPO 一阶 IF 实时误差标度 —— Anderson 杂质 (U=%.1f, ϵ_d=%.2f)，单模 bath×2，\n", U, ϵ_d)
    @printf("DiracDelta(ω=%.1f, α=%.1f)，β=%.1f，T=%.1f，D_bond=%d，ϵ=%.0e\n",
            ωb, αb, β, MODEL_T, D, ϵ)
    println("=" ^ 92)
    return scaling_main(realtime_gf; Ns=Ns, dtof=N -> MODEL_T / N, nlabel="Nt",
                        csvfile=joinpath(@__DIR__, "realtime.csv"),
                        pngfile=joinpath(@__DIR__, "realtime.png"),
                        title="GTEMPO first-order IF real-time error scaling\n" *
                              "Anderson IM (U=$U, ϵ_d=$ϵ_d, β=$β, T=$MODEL_T), single-mode bath ×2",
                        xlabel="δt = T / Nt",
                        ylabel="relative error of ĝ>/ĝ< series", D=D, ϵ=ϵ)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
