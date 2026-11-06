# ============================================================================
# diag_j0_repro.jl — matsubara j==0（ddag 折入 fold R）残差复现 + 定位扫描
#
# 第十轮补遗遗留窄口（HISTORY #43）：M≥2 j==0 ~1e-2。本脚本回答三件事：
#   ① 残差的 N 标度（N=4/8/16，U=0.8）；
#   ② U=0 是否同样错（是 ⟹ 补偿/F 路由结构病；否 ⟹ 传播子依赖）；
#   ③ 带不对称性：band-2 ddag 折入是否同样错（带对称掩蔽的反向利用）。
# 用法：julia --project=<repo> <repo>/benchmark/secondorderif/validation/diag_j0_repro.jl <repo>
# ============================================================================
include(joinpath(ARGS[1], "benchmark", "secondorderif", "src", "SecondOrderWeyl.jl"))
using .SecondOrderWeyl
using .SecondOrderWeyl: matsubara
include(joinpath(ARGS[1], "benchmark", "secondorderif", "validation", "ed.jl"))
using .SingleModeED
using GTEMPO
using Printf
using LinearAlgebra: diagm

const μ, ω, α, β = 0.7, 1.0, 0.5, 1.0

"j==0 列扫描：G(τ_i, 0) = ⟨d_b(τ_i) d_b†(0)⟩；b=0 用 matsubara API（band 1），否则直接装配"
function j0scan(N::Int, U::Float64; band::Int = 1)
    bath = fermionicbath(DiracDelta(ω = ω, α = α), β = β)
    mdl = AndersonIM(ϵ_d = μ, U = U)
    H, d1, _ = SingleModeED.model(μ = μ, U = U, bathspecs = [(ω, α), (ω, α)])
    c = WeylImagTime(Nτ = N, dτ = β / N, bands = 2)
    sol = solve(Problem(c, mdl; bath = bath, order = 2))
    dops = SecondOrderWeyl._fock_ops(2)
    ddag = band == 1 ? dops[1]' : dops[2]'
    worst = 0.0; at = ""
    println("== j==0 scan: N=$N U=$U band=$band ==")
    for i in 1:N
        _, gτ = SingleModeED.imag_observables(H, d1, [i * β / N], β)
        # band-2 参考：模型带对称（同 ω/α/μ）⟹ G₂ ≡ G₁，直接用 band-1 ED 值
        ref = real(gτ[1])
        mine = if band == 1
            real(matsubara(sol, i, 0; force = true))
        else
            K = system_factor(c, mdl; order = 2, coupling_norm = α,
                              insertions = [WeylInsertion(:τ, i, dops[2]),
                                            WeylInsertion(:τ, 0, ddag)])
            real(SecondOrderWeyl.expectation(SecondOrderWeyl.weyl_lattice(c), K,
                                             sol.mpsIs...; Z = sol.Z))
        end
        e = abs(mine - ref)
        @printf("  G(%2d,0) err = %.3e   mine %+.7f  ED %+.7f\n", i, e, mine, ref)
        e > worst && (worst = e; at = "($i,0)")
    end
    @printf("  → N=%d U=%.1f band=%d: worst = %.3e at %s\n\n", N, U, band, worst, at)
    return worst
end

for N in (4, 8)
    j0scan(N, 0.8; band = 1)
end
j0scan(8, 0.0; band = 1)
j0scan(8, 0.8; band = 2)
println("diag_j0_repro 完成")
