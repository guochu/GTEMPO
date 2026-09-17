# Run from the project root with:
#   julia --project=<env> docs/tutorials/bethedmft/imag/bethe_dmft_gtempo.jl
#
# Bethe-lattice DMFT for the Anderson impurity model on the Matsubara axis,
# following the reference workflow of GTEMPOProjects/src/bethedmft/imag/main.jl:
#
#   * initial Δ(iω): analytic semicircular-bath hybridisation `toulouse_Δiw × t²`
#   * Δ(iω) → Δ(τ) kernel:  `Δiw_to_Δτ` (QuAPI, re-exported by GTEMPO)
#   * impurity G(τ):        GTEMPO imaginary-time workflow (two spin bands)
#   * G(τ) → G(iω):         linear interpolation onto a fine τ grid, then
#                           `Gτ_to_Giw`
#   * Bethe self-consistency:
#         G₀⁻¹(iω) = iω + ϵ_d − t² G(iω),   Δ′(iω) = iω + ϵ_d − G₀(iω)
#     with ϵ_d = U/2 (half filling) and t = 1.
#
# The per-iteration Matsubara hybridisations are exported to JSON so that
# the TRIQS cthyb driver (bethe_dmft_ctqmc.py) is fed with the *identical*
# Δ(iω) sequence, starting from the *identical* initial state.

using GTEMPO
using JSON
using Printf
using LinearAlgebra
using Interpolations

# ------------------------------- settings ---------------------------
const smoke = get(ENV, "DMFT_SMOKE", "0") == "1"
const t_bethe = 1.0                       # Bethe hopping, band = [-2t, 2t]
const n_iter = smoke ? 2 : 10             # DMFT iterations
const mix = 0.5                           # linear mixing of Δ
const dτ = smoke ? 0.2 : 0.05             # imaginary-time step
const D_bond = smoke ? 60 : 120           # max bond dimension
const Nω = smoke ? 64 : 1024              # number of Matsubara frequencies (per side)
const δτ_fine = 1.0e-4                    # τ grid for the Gτ → Giw quadrature

const sets = smoke ? [(beta=5.0, U=1.0)] :
			 [(beta=5.0, U=1.0), (beta=5.0, U=5.0),
			  (beta=10.0, U=1.0), (beta=10.0, U=5.0)]
# DMFT_SET=k (1-based) restricts the run to a single parameter set,
# so that different sets can be run as parallel processes
const set_idx = parse(Int, get(ENV, "DMFT_SET", "0"))
const sets_run = (set_idx == 0) ? sets : [sets[set_idx]]

const outdir = joinpath(@__DIR__, "data")
mkpath(outdir)

# semicircular Bethe density of states on the band [-D, D]
ρ₀(ϵ, D) = sqrt(1 - (ϵ / D)^2) * (D / π)
spectrum_func(D) = spectrum(ω -> ρ₀(ω, D), lb=-D, ub=D)

# --------------------------- DMFT loop ------------------------------
function run_dmft(beta::Real, U::Real)
	ϵ_d = U / 2                        # Bethe chemical potential (half filling)
	t = t_bethe
	D_band = 2 * t
	Nτ = round(Int, beta / dτ)
	δτ = beta / Nτ
	println("="^70)
	@printf "Bethe DMFT (GTEMPO): β = %g, U = %g, Nτ = %d, Nω = %d\n" beta U Nτ Nω

	lattice = GrassmannLattice(N=Nτ, δτ=δτ, bands=2, contour=:imag)
	trunc = truncdimcutoff(D=D_bond, ϵ=1.0e-10)

	# impurity model and the Δ-independent part of the path (built once)
	bath = fermionicbath(spectrum_func(D_band), β=beta, μ=0)
	model = AndersonIM(U=U, ϵ_d=-ϵ_d)
	mpsK = sysdynamics(lattice, model, trunc=trunc)
	for band in 1:lattice.bands
		mpsK = boundarycondition!(mpsK, lattice, band=band, trunc=trunc)
	end

	# initial guess: analytic semicircular-bath hybridisation
	Δiw = toulouse_Δiw(bath, n=Nω) .* t^2

	iterations = []
	for it in 1:n_iter
		t0 = time()
		# Δ(iω) → Δ(τ) kernel → influence functional (band 1), then fill band 2
		corr = Δiw_to_Δτ(Δiw, β=beta, N=Nτ)
		mpsI1 = hybriddynamics(lattice, corr, trunc=trunc, band=1)
		mpsI2 = swapband(mpsI1, lattice, 1, 2, trunc=trunc)

		# impurity G(τ) on τ = 0 : δτ : β
		cache = environments(lattice, mpsK, mpsI1, mpsI2)
		gτ = -cached_gf_fast(lattice, mpsK, mpsI1, mpsI2;
							 c1=false, c2=true, b1=:τ, b2=:τ, cache=cache)

		# G(τ) → G(iω): linear interpolation onto a fine τ grid first
		interp = linear_interpolation(0:δτ:beta, gτ)
		gτ′ = interp.(0:δτ_fine:beta)
		Giw = Gτ_to_Giw(gτ′, β=beta, n=Nω)
		ws = ifrequencies(β=beta, n=Nω)

		# Bethe self-consistency:  Δ′(iω) = iω + ϵ_d − G₀(iω)
		G0iw = [1 / (im * w + ϵ_d - t^2 * G) for (w, G) in zip(ws, Giw)]
		Δiw_new = [im * w + ϵ_d - 1 / G0 for (w, G0) in zip(ws, G0iw)]
		err = norm(Δiw_new - Δiw) / norm(Δiw)
		@printf "  iter %2d: norm|Δ′ − Δ|/norm|Δ| = %.3e   (%.1f s)\n" it err (time() - t0)

		push!(iterations, Dict(
			"iter" => it,
			"ws" => ws,                                            # ωₙ (real)
			"Delta_iw" => [[real(d), imag(d)] for d in Δiw],       # Δ used this iteration
			"G_iw" => [[real(g), imag(g)] for g in Giw],
			"G_tau" => gτ, "tau" => collect(0:δτ:beta)))
		Δiw = (1 - mix) * Δiw + mix * Δiw_new
	end

	out = Dict("solver" => "GTEMPO", "beta" => beta, "U" => U, "t" => t,
			   "mu_imp" => -ϵ_d, "nw" => Nω, "ntau" => Nτ, "dtau" => δτ,
			   "mix" => mix, "n_iter" => n_iter, "iterations" => iterations)
	open(joinpath(outdir, "gtempo_b$(round(Int,beta))_U$(round(Int,U)).json"), "w") do f
		JSON.print(f, out)
	end
	println("written gtempo_b$(round(Int,beta))_U$(round(Int,U)).json")
end

for set in sets_run
	run_dmft(set.beta, set.U)
end
println("all done")
