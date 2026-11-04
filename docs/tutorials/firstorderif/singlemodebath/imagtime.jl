# =============================================================================
# firstorderif/singlemodebath/imagtime.jl —— 一阶影响泛函 GTEMPO 虚时演化误差标度
#
# 回答的问题：现有 GTEMPO 一阶 IF 管线（QuAPI 核 + τ 环装配 + 精确系统传播子
# + boundarycondition）在相互作用（U≠0、双带）情形下，虚时 g(τ) 的误差是否
# 满足一阶 scaling（relerr ∝ O(δτ)¹）。
#
# 方法：
# - 固定 β，Nτ ∈ {4..64} 加密 τ 环，δτ = β/Nτ；
# - 参考解 = 全模型的 ED（见 common.jl 的 gτ_ed，全耦合 H 的热 trace）；
#   虚时轮廓自身完成到 exp(-βH) 的演化，无需额外热态初态；
# - 误差 = 全时序相对误差 norm(g − g_ed)/norm(g_ed)。
#
# 运行：julia --project=<GTEMPO 根> docs/tutorials/firstorderif/singlemodebath/imagtime.jl
# 输出：imagtime.png（标度图）、imagtime.csv（数据）、终端误差表
# =============================================================================
include(joinpath(@__DIR__, "common.jl"))

"""
    imagtime_gf(Nτ; trunc) -> [("gτ", e)]

标准 GTEMPO 一阶 IF 管线的虚时 τ 环：粗格点 + 双带 fillband 复制 IF + 精确
系统传播子 + boundarycondition；误差对全模型 ED 参考。
"""
function imagtime_gf(Nτ::Int; trunc::TruncationScheme)
    δτ = β / Nτ
    τs = collect(0:δτ:β)
    lat = GrassmannLattice(N=Nτ, δτ=δτ, contour=:imag, bands=2)
    bath = fermionicbath(DiracDelta(ω=ωb, α=αb), β=β)
    mpsK, Is = fermionic_setup(lat, bath, AndersonIM(U=U, ϵ_d=ϵ_d), trunc)
    g = gτ_series(lat, mpsK, Is...)
    H, a, adag, H0 = singlemode_ed(μ=ϵ_d, U=U, bathspecs=[(ωb, αb), (ωb, αb)])
    ref = gτ_ed(H, a, adag, τs, β)
    return [("gτ", relerr(g, ref))]
end

function main(; Ns=(4, 8, 16, 32, 64), D=200, ϵ=1.0e-12)
    println("=" ^ 92)
    @printf("GTEMPO 一阶 IF 虚时误差标度 —— Anderson 杂质 (U=%.1f, ϵ_d=%.2f)，单模 bath×2，\n", U, ϵ_d)
    @printf("DiracDelta(ω=%.1f, α=%.1f)，β=%.1f，D_bond=%d，ϵ=%.0e\n",
            ωb, αb, β, D, ϵ)
    println("=" ^ 92)
    return scaling_main(imagtime_gf; Ns=Ns, dtof=N -> β / N, nlabel="Nτ",
                        csvfile=joinpath(@__DIR__, "imagtime.csv"),
                        pngfile=joinpath(@__DIR__, "imagtime.png"),
                        title="GTEMPO first-order IF imaginary-time error scaling\n" *
                              "Anderson IM (U=$U, ϵ_d=$ϵ_d, β=$β), single-mode bath ×2",
                        xlabel="δτ = β / Nτ",
                        ylabel="relative error of g(τ) series", D=D, ϵ=ϵ)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
