# ============================================================================
# twopoint_band_check.jl — 等时两点 0⁺ 约定 + 逐带 ⟨n_b⟩ 插入验证
#
# 内容：
#   ① 单带实时等时两点核（0⁺）：greater(t,t) = ⟨d(t⁺)d†(t)⟩ = 1−⟨n(t)⟩、
#      lesser(t,t) = ⟨d†(t)d(t⁻)⟩ = ⟨n(t)⟩（含 t=0 边界）；
#   ② 单带不等时回归：gt(t,0)/lt(t,0) vs ED（push 序修正不得破坏）；
#   ③ 两带 τ 逐带 ⟨n₁⟩/⟨n₂⟩ vs ED per-band（对称双带：二者应相等）；
#   ④ 两带实时逐带 ⟨n₁(t)⟩/⟨n₂(t)⟩ vs ED per-band + band-1 等时两点核。
# 运行：julia --project=. benchmark/secondorderif/validation/twopoint_band_check.jl <repo>
# ============================================================================
const REPO = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "..", "..")
include(joinpath(REPO, "benchmark", "secondorderif", "src", "SecondOrderWeyl.jl"))
using .SecondOrderWeyl
using .SecondOrderWeyl: occupation, raw_correlation, greater, lesser
include(joinpath(REPO, "benchmark", "secondorderif", "validation", "ed.jl"))
using .SingleModeED

using GTEMPO
using Printf
using LinearAlgebra: diagm, dot

const μ, ω, α, β, U2 = 0.7, 1.0, 0.5, 1.0, 0.8
const T_OBS = 0.4

function main()

    # ---- ED 参考算符（band-1 与 band-2 的杂质湮灭） ----
    function ed_band_d(nbands::Int, band::Int, bathspecs)
        H, dd, H0 = SingleModeED.model(μ = μ, U = nbands == 2 ? U2 : 0.0, bathspecs = bathspecs,
                                       nbands = nbands)
        dim = size(H, 1)
        d = zeros(Float64, dim, dim)
        for s in 0:dim-1
            sign, t = SingleModeED._apply_seq(s, [(band, false)])
            sign != 0 && (d[t+1, s+1] = sign)
        end
        return H, d, H0
    end

    println("=" ^ 78)
    println("① 单带实时等时两点核（0⁺ 约定）+ ② 不等时回归")
    bath1 = fermionicbath(DiracDelta(ω = ω, α = α), β = β)
    model1 = ToulouseIM(ϵ_d = μ)
    H1, dd1, H01 = SingleModeED.model(μ = μ, U = 0.0, bathspecs = [(ω, α)])
    ρ0_imp = diagm(0 => ComplexF64[1.0, exp(-β*μ)]) / (1 + exp(-β*μ))
    ρfull1 = SingleModeED.thermal(H01, β)
    ts = [0.1, 0.2, 0.3, 0.4]
    occ1_ed, gt1_ed, lt1_ed = SingleModeED.observables(H1, dd1, ts, ρfull1)
    occ0_ed, _, _ = SingleModeED.observables(H1, dd1, [0.0], ρfull1)

    c1 = WeylRealTime(Nt = 4, dt = T_OBS / 4, bands = 1)
    sol1 = solve(Problem(c1, model1; bath = bath1, order = 2, rho0 = ρ0_imp))
    maxerr_eq = 0.0
    for (k, t) in enumerate(ts)
        gk = raw_correlation(sol1, k, k; kind = :greater)      # 核 ⟨d d†⟩
        lk = raw_correlation(sol1, k, k; kind = :lesser)       # 核 ⟨d† d⟩
        e1 = abs(real(gk) - (1 - real(occ1_ed[k])))
        e2 = abs(real(lk) - real(occ1_ed[k]))
        maxerr_eq = max(maxerr_eq, e1, e2)
        @printf("  t=%.1f  greater核 mine %+.6f  ED(1−n) %+.6f |%+.1e   lesser核 mine %+.6f  ED(n) %+.6f |%+.1e\n",
                t, real(gk), 1 - real(occ1_ed[k]), e1, real(lk), real(occ1_ed[k]), e2)
    end
    g00 = raw_correlation(sol1, 0, 0; kind = :greater)
    l00 = raw_correlation(sol1, 0, 0; kind = :lesser)
    e1 = abs(real(g00) - (1 - real(occ0_ed[1]))); e2 = abs(real(l00) - real(occ0_ed[1]))
    maxerr_eq = max(maxerr_eq, e1, e2)
    @printf("  t=0    greater核 mine %+.6f  ED(1−n) %+.6f |%+.1e   lesser核 mine %+.6f  ED(n) %+.6f |%+.1e\n",
            real(g00), 1 - real(occ0_ed[1]), e1, real(l00), real(occ0_ed[1]), e2)

    maxerr_ueq = 0.0
    for (k, t) in enumerate(ts)
        e = max(abs(greater(sol1, k, 0) - gt1_ed[k]), abs(lesser(sol1, k, 0) - lt1_ed[k]))
        maxerr_ueq = max(maxerr_ueq, e)
        @printf("  不等时 t=%.1f  gt(t,0) mine %+.6f%+.6fi  ED %+.6f%+.6fi   lt 同查\n",
                t, real(greater(sol1, k, 0)), imag(greater(sol1, k, 0)), real(gt1_ed[k]), imag(gt1_ed[k]))
        maxerr_ueq = max(maxerr_ueq, abs(lesser(sol1, k, 0) - lt1_ed[k]))
    end
    @printf("  >>> ① 等时最大偏差 = %.2e   ② 不等时回归最大偏差 = %.2e\n\n", maxerr_eq, maxerr_ueq)

    println("=" ^ 78)
    println("③ 两带 τ 逐带 ⟨n_b⟩（order 2）+ ④ 两带实时逐带 ⟨n_b(t)⟩ 与 band-1 等时两点")
    bath2 = fermionicbath(DiracDelta(ω = ω, α = α), β = β)
    model2 = AndersonIM(ϵ_d = μ, U = U2)
    H2, d1_2b, H02 = ed_band_d(2, 1, [(ω, α), (ω, α)])
    _, d2_2b, _ = ed_band_d(2, 2, [(ω, α), (ω, α)])
    occ2_τ_ed, _ = SingleModeED.imag_observables(H2, d1_2b, [0.0], β)
    occ2b_τ_ed, _ = SingleModeED.imag_observables(H2, d2_2b, [0.0], β)
    @printf("  ED per-band τ: ⟨n₁⟩ = %.7f  ⟨n₂⟩ = %.7f\n", real(occ2_τ_ed), real(occ2b_τ_ed))

    errs_b = Float64[]
    for N in (8, 16)
        cτ = WeylImagTime(Nτ = N, dτ = β / N, bands = 2)
        solτ = solve(Problem(cτ, model2; bath = bath2, order = 2))
        n1 = real(occupation(solτ, N ÷ 2; band = 1))
        n2 = real(occupation(solτ, N ÷ 2; band = 2))
        e = max(abs(n1 - real(occ2_τ_ed)), abs(n2 - real(occ2b_τ_ed)))
        push!(errs_b, e)
        @printf("  τ N=%2d: ⟨n₁⟩ = %.7f  ⟨n₂⟩ = %.7f   最大偏差 = %.3e\n", N, n1, n2, e)
    end
    x = log.([β / 8, β / 16]); y = log.(errs_b)
    slope = (2 * dot(x, y) - sum(x) * sum(y)) / (2 * dot(x, x) - sum(x)^2)
    @printf("  >>> ③ τ 逐带收敛斜率 = %.3f\n", slope)

    ts2 = [0.2, 0.4]   # Nt=2 网格（Nt≥4 实时 K 的 permute!（M≥2 预存性能项）显著变慢）
    ρfull2 = SingleModeED.thermal(H02, β)
    occ2_rt_ed, _, _ = SingleModeED.observables(H2, d1_2b, ts2, ρfull2)
    occ2_rt_ed_b2, _, _ = SingleModeED.observables(H2, d2_2b, ts2, ρfull2)
    c2 = WeylRealTime(Nt = 2, dt = T_OBS / 2, bands = 2)
    ρ0_imp2 = diagm(0 => ComplexF64[1.0, exp(-β*μ), exp(-β*μ), exp(-β*(2μ + U2))]) /
              (1 + 2exp(-β*μ) + exp(-β*(2μ + U2)))
    sol2 = solve(Problem(c2, model2; bath = bath2, order = 2, rho0 = ρ0_imp2))
    maxerr = 0.0
    for (k, t) in enumerate(ts2)
        n1t = real(occupation(sol2, k; band = 1))
        n2t = real(occupation(sol2, k; band = 2))
        gt_eq = raw_correlation(sol2, k, k; kind = :greater)
        lt_eq = raw_correlation(sol2, k, k; kind = :lesser)
        e = max(abs(n1t - real(occ2_rt_ed[k])), abs(n2t - real(occ2_rt_ed_b2[k])),
                abs(real(gt_eq) - (1 - real(occ2_rt_ed[k]))), abs(real(lt_eq) - real(occ2_rt_ed[k])))
        maxerr = max(maxerr, e)
        @printf("  实时 t=%.1f: ⟨n₁⟩=%.6f (ED %.6f)  ⟨n₂⟩=%.6f (ED %.6f)  gt核=%.6f (ED %.6f)\n",
                t, n1t, real(occ2_rt_ed[k]), n2t, real(occ2_rt_ed_b2[k]),
                real(gt_eq), 1 - real(occ2_rt_ed[k]))
    end
    @printf("  >>> ④ 实时逐带+等时两点最大偏差 = %.2e\n", maxerr)

end

main()
