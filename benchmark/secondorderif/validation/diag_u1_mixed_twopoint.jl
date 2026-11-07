# ============================================================================
# diag_u1_mixed_twopoint.jl — 金标准对照：混合/实时轮廓两点函数（ED / U1 / 本侧）
#
# 同一模型（双带 Anderson U=0.8、Dirac-δ 浴 α=0.5、β=1、T=0.4、order 2）：
# 混合轮廓 WeylMixedTime/MixedTime（Nτ=Nt=N）；实时对照段 WeylRealTime/RealTime。
# ⚠️ 分支约定陷阱：U1 的 lesser 默认 = d@(+,i), d†@(−,j)（与本侧一致）；
#    强制 branch=:+ 会把 lesser 结构性变成 greater（i>j 时 + 支 contour
#    时间序恒给出 greater）——不要覆盖默认分支。
#   greater核 = ⟨d(+,tᵢ) d†(+,tⱼ)⟩（i≥j；j=0 = τ→+ 转折点）
#   lesser核  = ⟨d(+,tᵢ) d†(−,tⱼ)⟩（− 支 ddag）
#   ED 参考（平衡态平移不变）：C(tᵢ,tⱼ) = C(tᵢ−tⱼ)；等时 0⁺：greater=1−n、
#   lesser=n。
# 判决问题：混合两点的失配是本侧特有（修 K）还是跨框架（架构）。
# 用法：julia --project=~/Documents/Missile/U1GTEMPO2 <repo>/benchmark/secondorderif/validation/diag_u1_mixed_twopoint.jl <GTEMPO repo>
# ============================================================================
using U1GTEMPO2
using Printf
using LinearAlgebra: diagm, kron, tr, eigen, Diagonal

const REPO = length(ARGS) >= 1 ? ARGS[1] : "/home/guochu/Documents/Missile/GTEMPO"
include(joinpath(REPO, "benchmark", "secondorderif", "src", "SecondOrderWeyl.jl"))
using .SecondOrderWeyl
using .SecondOrderWeyl: raw_correlation
include(joinpath(REPO, "benchmark", "secondorderif", "validation", "ed.jl"))
using .SingleModeED
import GTEMPO

const μ, ω, α, β, T = 0.7, 1.0, 0.5, 1.0, 0.4
const U = 0.8
bath_u1 = DiscreteBath([ω], [α]; beta = β)

function edref(N, Δ::Real)
    H, d, _ = SingleModeED.model(μ = μ, U = U, bathspecs = [(ω, α), (ω, α)])
    ρc = SingleModeED.thermal(H, β)
    if Δ == 0
        n_th, _, _ = SingleModeED.observables(H, d, [0.0], ρc)
        return real(n_th[1]), 1.0 - real(n_th[1]), real(n_th[1])
    end
    n_th, gt, lt = SingleModeED.observables(H, d, [Δ], ρc)
    return real(n_th[1]), real(gt[1] / (-1.0im)), real(lt[1] / (1.0im))
end

for N in (4, 8)
    println("\n======== MixedTime Nτ=Nt=$N ========")
    c_u1 = MixedTime(Nτ = N, dτ = β / N, Nt = N, dt = T / N)
    prob_u1 = SIAMProblem(c_u1, SIAM(energy = (μ, μ), U = U); bath = bath_u1, order = 2)
    sol_u1 = U1GTEMPO2.solve(prob_u1)

    bath_m = GTEMPO.fermionicbath(GTEMPO.DiracDelta(ω = ω, α = α), β = β)
    c_m = SecondOrderWeyl.WeylMixedTime(Nτ = N, dτ = β / N, Nt = N, dt = T / N, bands = 2)
    sol_m = SecondOrderWeyl.solve(SecondOrderWeyl.Problem(c_m,
        GTEMPO.AndersonIM(ϵ_d = μ, U = U); bath = bath_m, order = 2))

    worst_u1 = 0.0; worst_m = 0.0; worst_d = 0.0
    at_u1 = at_m = at_d = ""
    for i in 1:N, j in 0:i
        Δ = (i - j) * (T / N)
        _, gΔ, lΔ = edref(N, j == i ? 0.0 : Δ)
        gm = real(raw_correlation(sol_m, i, j; kind = :greater))
        lm = real(raw_correlation(sol_m, i, j; kind = :lesser))
        # U1 默认分支：greater = d@(−,i), d†@(+,j)；lesser = d@(+,i), d†@(−,j)
        #（与本侧 raw_correlation 的约定一致；+/+ 强制会把 lesser 结构性变成 greater）
        gu = real(U1GTEMPO2.raw_correlation(sol_u1, i, j; kind = :greater))
        lu = real(U1GTEMPO2.raw_correlation(sol_u1, i, j; kind = :lesser))
        for (tag, v, ref, who) in (("g", gm, gΔ, :m), ("l", lm, lΔ, :m),
                                   ("g", gu, gΔ, :u), ("l", lu, lΔ, :u))
            e = abs(v - ref)
            if who === :m
                e > worst_m && (worst_m = e; at_m = "$tag($i,$j)")
            else
                e > worst_u1 && (worst_u1 = e; at_u1 = "$tag($i,$j)")
            end
        end
        for (tag, vm, vu) in (("g", gm, gu), ("l", lm, lu))
            e = abs(vm - vu)
            e > worst_d && (worst_d = e; at_d = "$tag($i,$j)")
        end
    end
    @printf("  U1 err  = %.3e（%s）\n  mine err = %.3e（%s）\n  |mine−U1| = %.3e（%s）\n",
            worst_u1, at_u1, worst_m, at_m, worst_d, at_d)
    sol_u1 = nothing; prob_u1 = nothing; sol_m = nothing; GC.gc()
end

println("\n======== RealTime 对照（Nt=N，因子化热初态）========")
# ED 参考（关键）：rt 初态 = ρ₀ = e^{−βH₀}/Z 因子化（杂质⊗浴，与 H 不对易），
# 两时关联必须显式演化（平衡态平移不变公式不适用于非对易 ρ₀）：
#   ⟨d(Δ)d†⟩ = Tr[ρ₀ d(Δ)d†]；等时 lesser核 = n(tᵢ)（演化后的占据，非 n(0)！）
H_rt, d_rt, H0_rt = SingleModeED.model(μ = μ, U = U, bathspecs = [(ω, α), (ω, α)])
Z0 = 1 + 2exp(-β * μ) + exp(-β * (2μ + U))
ρ0 = diagm(0 => ComplexF64[1.0 / Z0, exp(-β * μ) / Z0, exp(-β * μ) / Z0, exp(-β * (2μ + U)) / Z0])
p_th = (1.0 / Z0, exp(-β * μ) / Z0, exp(-β * μ) / Z0, exp(-β * (2μ + U)) / Z0)
Hb = diagm(0 => ComplexF64[0.0, ω])
ρb = exp(-β * Hb) / tr(exp(-β * Hb))
ρfull0 = kron(ρb, kron(ρb, ρ0))          # bit 序：bath2,bath1,band2,band1（mode m ⟺ bit m−1）
vals_rt, vecs_rt = eigen(H_rt)
_n_op = d_rt' * d_rt
edn(t::Real) = begin
    D = Diagonal(exp.(im .* t .* vals_rt))
    real(tr(ρfull0 * (vecs_rt * D * vecs_rt') * _n_op * (vecs_rt * D * vecs_rt')'))
end
function edrt_rt(Δ::Real, teq::Real)
    if Δ == 0
        n = edn(teq)
        return 1.0 - n, n
    end
    D = Diagonal(exp.(im .* Δ .* vals_rt))
    dt = (vecs_rt * D * vecs_rt') * d_rt * (vecs_rt * D * vecs_rt')'
    return real(tr(ρfull0 * dt * d_rt')), real(tr(ρfull0 * d_rt' * dt))
end
for N in (4, 8)
    c_u1 = RealTime(Nt = N, dt = T / N)
    prob_u1 = SIAMProblem(c_u1, SIAM(energy = (μ, μ), U = U); bath = bath_u1, order = 2,
                          initial_probabilities = p_th)
    sol_u1 = U1GTEMPO2.solve(prob_u1)
    bath_m = GTEMPO.fermionicbath(GTEMPO.DiracDelta(ω = ω, α = α), β = β)
    c_m = SecondOrderWeyl.WeylRealTime(Nt = N, dt = T / N, bands = 2)
    sol_m = SecondOrderWeyl.solve(SecondOrderWeyl.Problem(c_m,
        GTEMPO.AndersonIM(ϵ_d = μ, U = U); bath = bath_m, order = 2,
        rho0 = ρ0))
    worst_u1 = 0.0; worst_m = 0.0; worst_d = 0.0
    at_u1 = at_m = at_d = ""
    for i in 1:N, j in 0:i
        Δ = (i - j) * (T / N)
        gΔ, lΔ = edrt_rt(j == i ? 0.0 : Δ, i * (T / N))
        gm = real(raw_correlation(sol_m, i, j; kind = :greater))
        lm = real(raw_correlation(sol_m, i, j; kind = :lesser))
        gu = real(U1GTEMPO2.raw_correlation(sol_u1, i, j; kind = :greater))
        lu = real(U1GTEMPO2.raw_correlation(sol_u1, i, j; kind = :lesser))
        for (tag, v, ref, who) in (("g", gm, gΔ, :m), ("l", lm, lΔ, :m),
                                   ("g", gu, gΔ, :u), ("l", lu, lΔ, :u))
            e = abs(v - ref)
            if who === :m
                e > worst_m && (worst_m = e; at_m = "$tag($i,$j)")
            else
                e > worst_u1 && (worst_u1 = e; at_u1 = "$tag($i,$j)")
            end
        end
        for (tag, vm, vu) in (("g", gm, gu), ("l", lm, lu))
            e = abs(vm - vu)
            e > worst_d && (worst_d = e; at_d = "$tag($i,$j)")
        end
    end
    @printf("  Nt=%d: U1 err = %.3e（%s）| mine err = %.3e（%s）| mine−U1 = %.3e（%s）\n",
            N, worst_u1, at_u1, worst_m, at_m, worst_d, at_d)
    sol_u1 = nothing; prob_u1 = nothing; sol_m = nothing; GC.gc()
end
println("\n判决完成")
