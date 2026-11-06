# =============================================================================
# exactif/toulouse/realtime.jl —— 严格核 exact_Δt 的实时演化误差标度（Toulouse）
#
# 回答的问题：以严格核 exact_Δt（QuAPI exact_fermionic_Δt，耦合二次型 cell 的
# 精确单粒子传播子）替代 QuAPI 冻结核 Δt 后，GTEMPO 实时 Keldysh 演化对 δt
# 的收敛阶数，以及 Trotter 调度（order 1/2/4）对它的改善。
#
# 各序列：
#   frozen Δt (order 1)   冻结核 + 均匀全步长杂质演化   —— 期望 O(δt)
#   exact_Δt (order 1)    严格核 + 全步长放置            —— 期望 O(δt²)
#   exact_Δt (order 2)    严格核 + Strang 对称放置       —— 期望 O(δt²)
#   exact_Δt (order 4)    严格核 + Suzuki 四阶调度       —— 期望 O(δt⁴)
#
# 模型：Toulouse 杂质（U=0，H = ϵ_d n̂）+ 单模费米 bath（DiracDelta）；
# 参考解 = 单模 bath 的 ED（test/util.jl 的 singlemode_ed / greater_lesser_ed）。
# 误差 = 全时序相对误差 norm(gt − gt_ed)/norm(gt_ed)（lesser 同，取最大）。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/toulouse/realtime.jl
# 输出：realtime.png / realtime.csv / 终端误差表
# =============================================================================
include(joinpath(@__DIR__, "..", "src", "common.jl"))

const MODEL_T = 0.4    # 实时演化总时长
const β = 1.0          # 逆温
const ϵ_d = 0.7        # 杂质能级
const ωb, αb = 1.0, 0.5

toulouse_bath() = fermionicbath(DiracDelta(ω=ωb, α=αb), β=β)
toulouse_model() = ToulouseIM(ϵ_d=ϵ_d)

"""
    realtime_gf(N; trunc) -> [("frozen Δt (order 1)", e), ..., ("exact_Δt (order 4)", e)]

各管线的实时 Keldysh 演化误差（对单模 bath ED 参考）。Toulouse 的 H_imp 为
对角单带，均匀 cell 的 order 1/2 两条链在该模型上重合（IFGTEMPO
test_trotter_orders 的结论），两条曲线应逐点一致。
"""
function realtime_gf(N::Int; trunc::TruncationScheme)
    H, a, adag, H0 = singlemode_ed(μ=ϵ_d, U=0, bathspecs=[(ωb, αb)])
    res = Pair{String,Float64}[]
    # 冻结核参考（QuAPI Δt）：均匀时刻取样
    lat, mpsK, mpsIs = frozen_pipeline(toulouse_bath(), toulouse_model(), N, MODEL_T / N;
                                       trunc=trunc, β=β)
    gt, lt = gtlt_frozen(lat, mpsK, mpsIs)
    ts = collect(0:MODEL_T/N:MODEL_T)
    gt_e, lt_e = greater_lesser_ed(H, a, adag, H0, ts, β)
    push!(res, "frozen Δt (order 1)" => max(relerr(gt, gt_e), relerr(lt, lt_e)))
    # 严格核 + 各阶 Trotter 调度：按各调度的节点时刻取样
    for order in (1, 2, 4)
        pipe = exactif_pipeline(toulouse_bath(), toulouse_model(), N, MODEL_T / N;
                                order=order, trunc=trunc, β=β)
        gt, lt, ts = gtlt_exactif(pipe)
        gt_e, lt_e = greater_lesser_ed(H, a, adag, H0, ts, β)
        push!(res, "exact_Δt (order $order)" => max(relerr(gt, gt_e), relerr(lt, lt_e)))
    end
    return res
end

function main(; Ns=(2, 4, 8, 16), D=200, ϵ=1.0e-12)
    println("=" ^ 100)
    @printf("GTEMPO 严格核 exact_Δt 实时误差标度 —— Toulouse 模型 (U=0, ϵ_d=%.2f, β=%.1f)，\n", ϵ_d, β)
    @printf("单模 bath DiracDelta(ω=%.1f, α=%.1f)，T=%.1f，D_bond=%d，ϵ=%.0e\n",
            ωb, αb, MODEL_T, D, ϵ)
    println("=" ^ 100)
    return scaling_main(realtime_gf; Ns=Ns, T=MODEL_T,
                        csvfile=joinpath(@__DIR__, "realtime.csv"),
                        pngfile=joinpath(@__DIR__, "realtime.png"),
                        title="GTEMPO exact-IF kernel (exact_Δt) real-time error scaling\n" *
                              "Toulouse model (U=0, ϵ_d=$ϵ_d, β=$β, T=$MODEL_T), single-mode bath",
                        D=D, ϵ=ϵ)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
