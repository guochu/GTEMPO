# =============================================================================
# exactif/zidentity.jl —— 严格核 exact_Δt 的 Z 恒等式证书
#
# ε_d = 0 时严格核管线的配分函数恒等式（QuAPI exactfermionic.jl 的严格性声明）：
#
#     Z·prefactor == Tr_b[B^N ρ₀⊗ρ_b (B^N)†] = Tr[ρ₀] = 1 （逐位机器精度）
#
# 其中 B = exp(-iδt·(bath+耦合)) 是耦合二次型 cell 的精确演化。该恒等式对
# **任意 Trotter 调度**（order 1/2/4，含 Suzuki 负 cell、零 cell）成立，是
# 区分严格核与冻结核的最锐利证书：
#
#     冻结核 Δt：      |Z − 1| ≈ 5e-3（顶点冻结的 O(δt) 误差，见
#                      secondorder/convergence.jl）
#     严格核 exact_Δt：|Z − 1| ≈ 1e-14（全部 order、全部 N）
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/zidentity.jl
# =============================================================================
using GTEMPO
using Printf
using ImpurityModelBase          # fermionicbath / DiracDelta
using Z2Tensors: truncdimcutoff
include(joinpath(@__DIR__, "..", "src", "NativeExact.jl"))
using .NativeExact

const MODEL_T = 0.4
const β = 1.0
const ωb, αb = 1.0, 0.5

function main(; Ns=(2, 3, 4), D=400, ϵ=1.0e-14)
    bath = fermionicbath(DiracDelta(ω=ωb, α=αb), β=β)
    model0 = ToulouseIM(ϵ_d=0.0)                # ε_d = 0：骨架门 = 恒等传播子
    trunc = truncdimcutoff(D=D, ϵ=ϵ)
    println("=" ^ 92)
    @printf("严格核 Z 恒等式证书 —— Toulouse (ε_d=0, β=%.1f)，单模 bath，T=%.1f\n", β, MODEL_T)
    println("=" ^ 92)
    @printf("%-8s%-8s%-24s%-14s\n", "order", "N", "Z·prefactor", "|Z-1|")
    worst = 0.0
    for order in (1, 2, 4), N in Ns
        pipe = exactif_pipeline(bath, model0, N, MODEL_T / N; order=order, trunc=trunc, β=β)
        Z = integrate(pipe.lat, pipe.mpsK, pipe.mpsIs...) * pipe.kernel.prefactor
        worst = max(worst, abs(Z - 1))
        @printf("%-8d%-8d%-24.14f%-14.3e\n", order, N, real(Z), abs(Z - 1))
        abs(Z - 1) < 1e-10 || error("Z identity violated at order=$order N=$N")
    end
    # 冻结核对照：同一恒等式的 O(δt) 偏差
    println("-" ^ 92)
    for N in Ns
        lat, mpsK, mpsIs = frozen_pipeline(bath, model0, N, MODEL_T / N; trunc=trunc, β=β)
        Z = integrate(lat, mpsK, mpsIs...)
        @printf("%-8s%-8d%-24.14f%-14.3e\n", "frozen", N, real(Z), abs(Z - 1))
    end
    @printf("\n严格核最差 |Z-1| = %.3e（全部 order 通过 1e-10 阈值）\n", worst)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
