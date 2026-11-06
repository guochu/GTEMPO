# ============================================================================
# diag_j0_factor.jl — j==0 的 U=0 因子化判决：M=2(U=0) vs M=1 逐点对比
#
# 原理：U=0 ⟹ K₂ = K₁∧K₁（第十轮已证恒等式），同 bath/离散 ⟹ 两轮廓的
# G_j0 应当一致（且都收敛到同一 ED 单带参考）。M=1 j==0 已闭合（3.67e-5@N=16），
# M=2 j==0 差 ~1e-2 ⟹ 差值 = M=2 特有补偿结构（Fb 路由/交叉项）的直接暴露。
# 附：N=1 fold 孤立（链 = 单退化格 Tr[A·R]，无传播）——若 N=1 即错则病在
# fold/close 格本身；若 N=1 对则病随链传播（F 全程路由）。
# 用法：julia --project=<repo> <repo>/benchmark/secondorderif/validation/diag_j0_factor.jl <repo>
# ============================================================================
include(joinpath(ARGS[1], "benchmark", "secondorderif", "src", "SecondOrderWeyl.jl"))
using .SecondOrderWeyl
using .SecondOrderWeyl: matsubara
include(joinpath(ARGS[1], "benchmark", "secondorderif", "validation", "ed.jl"))
using .SingleModeED
using GTEMPO
using Printf

const μ, ω, α, β = 0.7, 1.0, 0.5, 1.0

# 单带 ED 参考（U=0 ⟹ 两带解耦，各带 bath α）
H1, d1, _ = SingleModeED.model(μ = μ, U = 0.0, bathspecs = [(ω, α)])

function gref(τ::Float64)
    _, gτ = SingleModeED.imag_observables(H1, d1, [τ], β)
    return real(gτ[1])
end

"返回 (mine_j0列, N)"
function j0col_m2(N::Int)
    bath = fermionicbath(DiracDelta(ω = ω, α = α), β = β)
    c = WeylImagTime(Nτ = N, dτ = β / N, bands = 2)
    sol = solve(Problem(c, AndersonIM(ϵ_d = μ, U = 0.0); bath = bath, order = 2))
    return [real(matsubara(sol, i, 0; force = true)) for i in 1:N], N
end

function j0col_m1(N::Int)
    bath = fermionicbath(DiracDelta(ω = ω, α = α), β = β)
    c = WeylImagTime(Nτ = N, dτ = β / N, bands = 1)
    sol = solve(Problem(c, ToulouseIM(ϵ_d = μ); bath = bath, order = 2))
    return [real(matsubara(sol, i, 0)) for i in 1:N], N
end

println("== U=0：M=2 vs M=1 j==0 列（同 bath α=$(α), β=$(β)）==")
for N in (1, 2, 4, 8)
    m2, _ = j0col_m2(N)
    m1, _ = N == 1 ? ([m2[1]], 1) : j0col_m1(N)
    println("-- N=$N --")
    for i in 1:N
        τ = i * β / N
        ref = gref(τ)
        @printf("  G(%2d,0): M2 %+.7f (err %.2e) | M1 %+.7f (err %.2e) | M2−M1 %+.2e | ED %+.7f\n",
                i, m2[i], abs(m2[i] - ref), m1[i], abs(m1[i] - ref),
                m2[i] - m1[i], ref)
    end
end
println("diag_j0_factor 完成")
