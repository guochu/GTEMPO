include("/home/guochu/Documents/Missile/GTEMPO/benchmark/secondorderif/src/SecondOrderWeyl.jl")
using .SecondOrderWeyl
using .SecondOrderWeyl: raw_correlation, matsubara
include("/home/guochu/Documents/Missile/GTEMPO/benchmark/secondorderif/validation/ed.jl")
using .SingleModeED
using GTEMPO
using Printf
using LinearAlgebra: diagm

const μ, ω, α, β, T = 0.7, 1.0, 0.5, 1.0, 0.4

println("== rt M=2（U≠0，因子化热初态 ρ0_imp2，Nt=4 全网格）==")
function rtverdict()
    bath = fermionicbath(DiracDelta(ω = ω, α = α), β = β)
    model = AndersonIM(ϵ_d = μ, U = 0.8)
    H, d, H0 = SingleModeED.model(μ = μ, U = 0.8, bathspecs = [(ω, α), (ω, α)])
    ρ0_imp = diagm(0 => ComplexF64[1.0, exp(-β*μ), exp(-β*μ), exp(-β*(2μ + 0.8))]) /
              (1 + 2exp(-β*μ) + exp(-β*(2μ + 0.8)))
    ρfull = SingleModeED.thermal(H0, β)
    N = 8; dt = T / N
    c = WeylRealTime(Nt = N, dt = dt, bands = 2)
    sol = solve(Problem(c, model; bath = bath, order = 2, rho0 = ρ0_imp))
    worst = 0.0; at = ""
    for i in 1:N, j in 0:i
        Δ = (i - j) * dt
        if j == i
            occ0, _, _ = SingleModeED.observables(H, d, [Δ], ρfull)
            gΔ = 1.0 - real(occ0[1]); lΔ = real(occ0[1])
        else
            _, gt, lt = SingleModeED.observables(H, d, [Δ], ρfull)
            gΔ = gt[1] / (-1.0im); lΔ = lt[1] / (1.0im)
        end
        for (tag, mine, ref) in ((:g, raw_correlation(sol, i, j; kind = :greater), gΔ),
                                 (:l, raw_correlation(sol, i, j; kind = :lesser), lΔ))
            e = abs(mine - ref)
            e > worst && (worst = e; at = "$tag($i,$j)")
        end
    end
    @printf("  (i≥j) 全网格 × {g,l} 最大偏差 = %.2e（%s）%s\n\n", worst, at,
            worst < 2e-3 ? "—— 闭合 ✅" : "—— 超限 ❌")
end
rtverdict()

println("== matsubara（i≥j 全网格，N=8）==")
function matsverdict(nb::Int, mdl, U::Float64)
    bath3 = fermionicbath(DiracDelta(ω = ω, α = α), β = β)
    H, d, _ = SingleModeED.model(μ = μ, U = U, bathspecs = [(ω, α) for _ in 1:nb], nbands = nb)
    N = 16
    c = WeylImagTime(Nτ = N, dτ = β / N, bands = nb)
    sol = solve(Problem(c, mdl; bath = bath3, order = 2))
    worst = 0.0; at = ""
    internal_only = nb >= 2   # M≥2：等时 + 双内部（1≤j≤i≤Nτ−1）已闭合；边界型防守拦截
    for i in 1:N, j in 0:i
        internal_only && j != i && j < 1 && continue   # M≥2：边界型（j==0）防守拦截，跳过
        internal_only && i == N && j != i && continue  # M≥2：close 格（i==Nτ）防守拦截，跳过
        Δτ = (i - j) * (β / N)
        # τ 环 0⁺ 约定（cutop 作用序 ddag·d ⟹ 等时 = ⟨d†d⟩ = n）；τ>0 = ⟨d(τ)d†(0)⟩
        if Δτ == 0
            nτ, _ = SingleModeED.imag_observables(H, d, [β], β)
            ref = real(nτ)
        else
            _, gτ = SingleModeED.imag_observables(H, d, [Δτ], β)
            ref = real(gτ[1])
        end
        mine = real(matsubara(sol, i, j))
        e = abs(mine - ref)
        e > worst && (worst = e; at = "($i,$j)")
    end
    scope = internal_only ? "M=2 等时+双内部" : "M=1 (i≥j) 全网格"
    @printf("  %s最大偏差 = %.2e（%s）%s\n", scope, worst, at,
            worst < 2e-3 ? "—— 闭合 ✅" : "—— 超限 ❌")
end
matsverdict(1, ToulouseIM(ϵ_d = μ), 0.0)
matsverdict(2, AndersonIM(ϵ_d = μ, U = 0.8), 0.8)
