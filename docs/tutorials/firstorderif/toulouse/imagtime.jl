# =============================================================================
# firstorderif/toulouse/imagtime.jl —— 一阶影响泛函 GTEMPO 虚时演化误差标度
#
# 回答的问题：现有 GTEMPO 一阶 IF 管线（QuAPI 核 + τ 环装配 + 精确系统传播子
# + boundarycondition）的虚时 g(τ)，误差是否满足一阶 scaling（relerr ∝
# O(δτ)¹）。
#
# 方法：
# - 固定 β，Nτ ∈ {4..64} 加密 τ 环，δτ = β/Nτ；
# - 参考解 = 单模 bath 的 ED（见 common.jl 的 gτ_ed，全耦合 H 的热 trace）；
#   虚时轮廓自身完成到 exp(-βH) 的演化，无需额外热态初态；
# - 误差 = 全时序相对误差 norm(g − g_ed)/norm(g_ed)。
#
# 运行：julia --project=<GTEMPO 根> docs/tutorials/firstorderif/toulouse/imagtime.jl
# 输出：imagtime.png（标度图）、imagtime.csv（数据）、终端误差表
# =============================================================================
include(joinpath(@__DIR__, "common.jl"))

"""
    imagtime_gf(Nτ; trunc) -> [("gτ", e)]

标准 GTEMPO 一阶 IF 管线的虚时 τ 环：粗格点 + hybriddynamics + 精确系统
传播子 + boundarycondition；误差对单模 ED 参考。
"""
function imagtime_gf(Nτ::Int; trunc::TruncationScheme)
    δτ = β / Nτ
    τs = collect(0:δτ:β)
    lat = GrassmannLattice(N=Nτ, δτ=δτ, contour=:imag)
    corr = correlationfunction(toulouse_bath(), lat)
    mpsI = hybriddynamics(lat, corr, trunc=trunc)
    mpsK = sysdynamics(lat, ToulouseIM(ϵ_d=ϵ_d), trunc=trunc)
    mpsK = boundarycondition!(mpsK, lat)
    g = gτ_series(lat, mpsK, mpsI)
    H, a, adag, H0 = singlemode_ed(μ=ϵ_d, U=0, bathspecs=[(ωb, αb)])
    ref = gτ_ed(H, a, adag, τs, β)
    return [("gτ", relerr(g, ref))]
end

function main(; Ns=(4, 8, 16, 32, 64), D=200, ϵ=1.0e-12)
    println("=" ^ 92)
    @printf("GTEMPO 一阶 IF 虚时误差标度 —— Toulouse 模型 (U=0, ϵ_d=%.2f)，\n", ϵ_d)
    @printf("单模 bath DiracDelta(ω=%.1f, α=%.1f)，β=%.1f，D_bond=%d，ϵ=%.0e\n",
            ωb, αb, β, D, ϵ)
    println("=" ^ 92)
    return scaling_main(imagtime_gf; Ns=Ns, dtof=N -> β / N, nlabel="Nτ",
                        csvfile=joinpath(@__DIR__, "imagtime.csv"),
                        pngfile=joinpath(@__DIR__, "imagtime.png"),
                        title="GTEMPO first-order IF imaginary-time error scaling\n" *
                              "Toulouse model (U=0, ϵ_d=$ϵ_d, β=$β), single-mode bath",
                        xlabel="δτ = β / Nτ",
                        ylabel="relative error of g(τ) series", D=D, ϵ=ϵ)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
