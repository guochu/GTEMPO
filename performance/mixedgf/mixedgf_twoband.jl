# ==============================================================================
# Two-band Anderson impurity with U != 0: fully mixed Green's function
# G(τ_i, t_j) = <d(τ_i) d†(t_j)> on the full (τ, t) grid, GTEMPO vs ED.
#
# The ED reference uses the exact 16-dimensional Hamiltonian (two impurity
# bands + one bath mode per band), so the interacting two-band dynamics are
# benchmarked exactly. Mirrors the single-band script mixedgf.jl.
#
# Run from the project root with:
#   julia --project=. performance/mixedgf/mixedgf_twoband.jl
# ============================================================================

using GTEMPO
using Z2Tensors
using ImpurityModelBase
using Random
using LinearAlgebra: tr, exp, norm
using Printf

Random.seed!(1)
include(joinpath(@__DIR__, "..", "..", "test", "util.jl"))
include(joinpath(@__DIR__, "..", "..", "test", "normalbath", "util.jl"))

function main()
    mu = 0.7; U = 1.0; w = 1.0; alpha = 0.5
    beta = 1.0; dtau = 0.1; Ntau = 10
    dt = 0.05; Nt = 6
    trunc = truncdimcutoff(D=100, ϵ=1.0e-10)

    # ED Hamiltonian: 2 impurity bands + 2 bath modes (dim = 16)
    H, a, adag, H0 = singlemode_ed(μ=mu, U=U, bathspecs=[(w, alpha), (w, alpha)])
    rho = exp(-beta * H); rho /= tr(rho)

    spec = DiracDelta(ω=w, α=alpha)
    bath = fermionicbath(spec, β=beta)
    model = AndersonIM(U=U, μ=mu)

    lat = GrassmannLattice(Nt=Nt, δt=dt, Nτ=Ntau, δτ=dtau, contour=:mixed, bands=2)
    mpsK, Is = fermionic_setup(lat, bath, model, trunc)
    cache = environments(lat, mpsK, Is...)
    println("two-band mixed lattice: Nτ=$Ntau, Nt=$Nt, bands=2, U=$U")

    # GTEMPO mixed GF: annihilator on the τ branch (conj=false), creator on a
    # real branch (conj=true); band 1
    function mixed_gtempo(i, br, j)
        cached_gf(lat, (ContourIndex(i, conj=false, branch=:τ, band=1),
                        ContourIndex(j, conj=true, branch=br, band=1)), mpsK, Is...; cache=cache)
    end

    # ED reference: exact circuit on the interacting Hamiltonian
    function mixed_ed(tau_v, t_v)
        op1 = exp(tau_v * H) * a * exp(-tau_v * H)
        op2 = exp(im * t_v * H) * adag * exp(-im * t_v * H)
        return tr(rho * op1 * op2) / tr(rho)
    end

    println("\n=== G(τ_i, t_j): GTEMPO vs ED (interacting, full grid) ===")
    worst = Dict{Symbol, Float64}(:+ => 0.0, :- => 0.0)
    worst_at = Dict{Symbol, Tuple{Int, Int}}(:+ => (0, 0), :- => (0, 0))
    scale = 0.0
    for br in (:+, :-)
        for i in 1:Ntau+1
            tau_v = (i - 1) * dtau
            for j in 1:Nt+1
                t_v = (j - 1) * dt
                g = mixed_gtempo(i, br, j)
                e = mixed_ed(tau_v, t_v)
                scale = max(scale, abs(e))
                if abs(g - e) > worst[br]
                    worst[br] = abs(g - e)
                    worst_at[br] = (i, j)
                end
            end
        end
        @printf("  branch (%s,·): worst |GTEMPO - ED| = %.3e  at (i,j) = %d,%d\n",
                br, worst[br], worst_at[br][1], worst_at[br][2])
    end
    @printf("  max |ED reference| on the grid = %.3f\n", scale)

    println("\n--- sample columns (real part), GTEMPO vs ED ---")
    for i in (1, div(Ntau, 2) + 1, Ntau + 1)
        tau_v = (i - 1) * dtau
        println("τ_i = $(round(tau_v, digits=3)) :")
        for br in (:+, :-)
            print("   branch ", br, "  GTEMPO: ")
            for j in 1:Nt+1
                @printf("%+.4f  ", real(mixed_gtempo(i, br, j)))
            end
            print("\n            ED:     ")
            for j in 1:Nt+1
                @printf("%+.4f  ", real(mixed_ed(tau_v, (j - 1) * dt)))
            end
            println()
        end
    end

    println("\nInterpretation:")
    println("  worst |Δ| at the level of the truncation error (D=100) means the")
    println("  two-band interacting mixed GF is accurate, i.e. the one-site τ shift")
    println("  of the IF does not degrade the cross-branch Green's function.")
end

main()
