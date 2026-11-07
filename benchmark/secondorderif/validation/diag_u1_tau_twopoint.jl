# ============================================================================
# diag_u1_tau_twopoint.jl — 金标准对照：U1GTEMPO2 的 τ 两点是否带结构尾部
#
# 同一模型三框架对比（双带 Anderson、Dirac-δ 浴 α=0.5、β=1、order 2）：
#   ED（SingleModeED，band-1 参考准确值）
#   U1GTEMPO2（U(1) 电荷管线，matsubara(s,i,0;spin=:up)）
#   SecondOrderWeyl（本基准，Z2 分级管线，matsubara 路由列）
# 判据：G(i,0) 列逐点 |计算 − ED|。若 U1 干净（纯离散化 ~5e-5@N=8）而
# 本基准 ~1e-3，病因在本侧 Z2 构造；若 U1 同带 ~1e-3 尾，则为跨框架架构
# 共性（金标准自身的验收上限）。
# 用法：julia --project=~/Documents/Missile/U1GTEMPO2 <repo>/benchmark/secondorderif/validation/diag_u1_tau_twopoint.jl <GTEMPO repo>
# ============================================================================
using U1GTEMPO2
using Printf

const REPO = length(ARGS) >= 1 ? ARGS[1] : "/home/guochu/Documents/Missile/GTEMPO"
include(joinpath(REPO, "benchmark", "secondorderif", "src", "SecondOrderWeyl.jl"))
using .SecondOrderWeyl
using .SecondOrderWeyl: matsubara as mats_mine
include(joinpath(REPO, "benchmark", "secondorderif", "validation", "ed.jl"))
using .SingleModeED
import GTEMPO

const μ, ω, α, β = 0.7, 1.0, 0.5, 1.0
bath_u1 = DiscreteBath([ω], [α]; beta = β)

for U in (0.0, 0.8)
    println("\n======== U = $U ========")
    H, d1, _ = SingleModeED.model(μ = μ, U = U, bathspecs = [(ω, α), (ω, α)])
    for N in (4, 8)
        # ---- U1 侧 ----
        c_u1 = ImagTime(Nτ = N, dτ = β / N)
        prob_u1 = SIAMProblem(c_u1, SIAM(energy = (μ, μ), U = U); bath = bath_u1, order = 2)
        sol_u1 = U1GTEMPO2.solve(prob_u1)
        g_u1 = [U1GTEMPO2.matsubara(sol_u1, i, 0; spin = :up) for i in 1:N]
        # ---- 本基准侧（matsubara j==0 走平移路由）----
        bath_m = GTEMPO.fermionicbath(GTEMPO.DiracDelta(ω = ω, α = α), β = β)
        mdl_m = U == 0.0 ? GTEMPO.AndersonIM(ϵ_d = μ, U = 0.0) : GTEMPO.AndersonIM(ϵ_d = μ, U = U)
        c_m = SecondOrderWeyl.WeylImagTime(Nτ = N, dτ = β / N, bands = 2)
        sol_m = SecondOrderWeyl.solve(SecondOrderWeyl.Problem(c_m, mdl_m; bath = bath_m, order = 2))
        g_m = [real(mats_mine(sol_m, i, 0)) for i in 1:N]
        # ---- ED 参考与判决 ----
        println("-- N = $N --")
        println("   i     |      ED      |   U1 err    |   mine err")
        for i in 1:N
            _, gτ = SingleModeED.imag_observables(H, d1, [i * β / N], β)
            ref = real(gτ[1])
            # U1 符号约定：取与 ED 一致者
            gu = real(g_u1[i])
            sgn = abs(gu - ref) <= abs(-gu - ref) ? 1 : -1
            eu = abs(sgn * gu - ref)
            em = abs(g_m[i] - ref)
            @printf("  %2d    |  %+.7f  |  %.3e   |  %.3e\n", i, ref, eu, em)
        end
        sol_u1 = nothing; prob_u1 = nothing; sol_m = nothing; GC.gc()
    end
end
println("\n判决完成")
