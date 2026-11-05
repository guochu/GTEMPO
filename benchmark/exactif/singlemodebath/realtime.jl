# =============================================================================
# exactif/singlemodebath/realtime.jl —— 严格核 exact_Δt 的实时演化误差标度
#（单模 bath Anderson 杂质模型，U ≠ 0）
#
# 与 toulouse/realtime.jl 同一标度研究，模型换为两带 Anderson 杂质
# H_imp = μ(n̂₁+n̂₂) + U n̂₁n̂₂ + 单模 bath（每带一个 bath 模）。U ≠ 0 时
# [H_imp(作用项), 耦合] ≠ 0，Trotter 调度的差异真正起效：
#
#   frozen Δt (order 1)   冻结核 + 均匀全步长放置   —— O(δt)
#   exact_Δt (order 1/2)  严格核（均匀 cell，两放置等价）—— 方案误差 O(δt²)
#   exact_Δt (order 4)    严格核 + Suzuki 四阶调度   —— O(δt⁴)
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/singlemodebath/realtime.jl
# =============================================================================
include(joinpath(@__DIR__, "..", "common.jl"))

const MODEL_T = 0.4
const β = 1.0
const μ_d = 0.7       # 杂质化学势（AndersonIM 的 ϵ_d）
const Uint = 1.5      # Hubbard U
const ωb, αb = 1.0, 0.5

anderson_bath() = fermionicbath(DiracDelta(ω=ωb, α=αb), β=β)
anderson_model() = AndersonIM(ϵ_d=μ_d, U=Uint)

function realtime_gf(N::Int; trunc::TruncationScheme)
    H, a, adag, H0 = singlemode_ed(μ=μ_d, U=Uint, bathspecs=[(ωb, αb), (ωb, αb)], nbands=2)
    res = Pair{String,Float64}[]
    # 冻结核参考（QuAPI Δt）：均匀时刻取样
    lat, mpsK, mpsIs = frozen_pipeline(anderson_bath(), anderson_model(), N, MODEL_T / N;
                                       trunc=trunc, β=β, bands=2)
    gt, lt = gtlt_frozen(lat, mpsK, mpsIs)
    ts = collect(0:MODEL_T/N:MODEL_T)
    gt_e, lt_e = greater_lesser_ed(H, a, adag, H0, ts, β)
    push!(res, "frozen Δt (order 1)" => max(relerr(gt, gt_e), relerr(lt, lt_e)))
    # 严格核 + 各阶 Trotter 调度
    for order in (1, 2, 4)
        pipe = exactif_pipeline(anderson_bath(), anderson_model(), N, MODEL_T / N;
                                order=order, trunc=trunc, β=β, bands=2)
        gt, lt, ts = gtlt_exactif(pipe)
        gt_e, lt_e = greater_lesser_ed(H, a, adag, H0, ts, β)
        push!(res, "exact_Δt (order $order)" => max(relerr(gt, gt_e), relerr(lt, lt_e)))
    end
    return res
end

function main(; Ns=(2, 4, 8, 16), D=200, ϵ=1.0e-12)
    println("=" ^ 100)
    @printf("GTEMPO 严格核 exact_Δt 实时误差标度 —— Anderson 模型 (μ=%.2f, U=%.2f, β=%.1f)，\n",
            μ_d, Uint, β)
    @printf("单模 bath DiracDelta(ω=%.1f, α=%.1f)×2，T=%.1f，D_bond=%d，ϵ=%.0e\n",
            ωb, αb, MODEL_T, D, ϵ)
    println("=" ^ 100)
    return scaling_main(realtime_gf; Ns=Ns, T=MODEL_T,
                        csvfile=joinpath(@__DIR__, "realtime.csv"),
                        pngfile=joinpath(@__DIR__, "realtime.png"),
                        title="GTEMPO exact-IF kernel (exact_Δt) real-time error scaling\n" *
                              "Anderson model (μ=$μ_d, U=$Uint, β=$β, T=$MODEL_T), single-mode bath ×2",
                        D=D, ϵ=ϵ)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
