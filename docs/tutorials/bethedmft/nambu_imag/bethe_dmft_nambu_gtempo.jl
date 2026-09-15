# Run from the project root with:
#   julia --project=<env> docs/tutorials/bethedmft/nambu_imag/bethe_dmft_nambu_gtempo.jl
#
# Bethe-lattice DMFT in the superconducting (Nambu) state on the Matsubara
# axis, following the reference workflow of
# GTEMPOProjects/src/bcsbath/siam/main.jl (main_imag_dmft):
#
#   * the bath enters as a BCSCorrelationFunction with four blocks
#     (uu, dd: normal electron/hole channels; ud, du: anomalous channels),
#     each built from Matsubara data via `Δiw_to_Δτ` (QuAPI);
#   * the influence functional is built with `hybriddynamics_naive`
#     (orbital=1, two flavors per orbital) and enters the evaluation as a
#     *single* IF, together with the Δ-independent path mpsK;
#   * observables: normal G(τ) from `cached_gf_fast` (per band), anomalous
#     G(τ) from point-wise `cached_gf` with cross-band ContourIndex pairs;
#   * Bethe self-consistency in Nambu convention (reference code):
#         Δuu′ = t² G_uu ,   Δdd′ = −t² conj(G_dd) ,
#         Δud′ = t² G_ud ,   Δdu′ = t² G_du .
#
# Half filling: ϵ_d = U/2, i.e. `AndersonIM(U=U, μ=−U/2)` (particle–hole
# symmetric for either sign of U).  With the user-specified U = 1, 5 the
# bath probes the anomalous channel; the same script also supports the
# attractive sign (e.g. U = −1 of the reference code) unchanged.
#
# Per-iteration Δ(iωₙ) (all four blocks) is exported to JSON so that the
# TRIQS cthyb driver (bethe_dmft_nambu_ctqmc.py) is fed with the
# *identical* Nambu hybridisation sequence.

using GTEMPO
using JSON
using Printf
using LinearAlgebra
using Interpolations

# ------------------------------- settings ---------------------------
const smoke = get(ENV, "DMFT_SMOKE", "0") == "1"
const t_bethe = 1.0                       # Bethe hopping; half-bandwidth 2t
const n_iter = smoke ? 3 : 10             # DMFT iterations
const mix = 0.5                           # linear mixing of Δ
const dτ = smoke ? 0.25 : 0.05            # imaginary-time step
const D_bond = smoke ? 40 : 64            # max bond dimension
const Nω = smoke ? 64 : 512               # Matsubara frequencies per side
const δτ_fine = 1.0e-4                    # τ grid for Gτ → Giw
const f₀ = 0.3                            # initial anomalous Weiss amplitude
const ω_c = 1.0                           # initial anomalous cutoff scale

const sets = smoke ? [(beta=1.0, U=1.0)] :
			 [(beta=1.0, U=1.0), (beta=1.0, U=5.0)]
# DMFT_SET=k (1-based) restricts the run to a single parameter set,
# so that different sets can be run as parallel processes.
# DMFT_U=... overrides U in all selected sets (U = 0 gives the exactly
# solvable non-interacting cross-check).
const set_idx = parse(Int, get(ENV, "DMFT_SET", "0"))
const sets_run = (set_idx == 0) ? sets : [sets[set_idx]]
const u_env = get(ENV, "DMFT_U", "")
const sets_final = [(beta=s.beta,
					 U=(u_env == "" ? s.U : parse(Float64, u_env)))
					for s in sets_run]

const outdir = joinpath(@__DIR__, "data")
mkpath(outdir)

# --------------------- initial Weiss field --------------------------
"non-interacting Bethe Green's function on the Matsubara axis (Σ = 0)"
function bethe_G0(ωn::Real, t::Real)
	z = im * ωn
	# branch with √(z²−4t²) ~ z as |z| → ∞ (retarded), i.e. the sign of the
	# square root follows the sign of ωₙ on the imaginary axis
	disc = (ωn >= 0 ? im : -im) * sqrt(ωn^2 + 4 * t^2)
	return (z - disc) / (2 * t^2)
end

function initial_weiss(β::Real, Nω::Int, t::Real)
	ws = ifrequencies(β=β, n=Nω)
	Δuu = [t^2 * bethe_G0(w, t) for w in ws]
	Δdd = [-conj(d) for d in Δuu]                  # hole channel of Δuu
	Δud = [f₀ / sqrt(1 + (w / ω_c)^2) for w in ws] # real, even in ωₙ
	Δdu = copy(Δud)
	return ws, Δuu, Δdd, Δud, Δdu
end

# --------------------------- DMFT loop ------------------------------
function run_dmft(beta::Real, U::Real)
	μ_imp = -U / 2                     # particle–hole symmetric filling
	t = t_bethe
	Nτ = round(Int, beta / dτ)
	δτ = beta / Nτ
	println("="^70)
	@printf "Nambu DMFT (GTEMPO): β = %g, U = %g, μ = %g, Nτ = %d, Nω = %d\n" beta U μ_imp Nτ Nω

	trunc = truncdimcutoff(D=D_bond, ϵ=1.0e-10)
	lattice = GrassmannLattice(N=Nτ, δτ=δτ, contour=:imag, bands=2)
	model = AndersonIM(U=U, μ=μ_imp)

	# Δ-independent part of the path (built once)
	mpsK = sysdynamics(lattice, model, trunc=trunc)
	for band in 1:lattice.bands
		mpsK = boundarycondition!(mpsK, lattice, band=band, trunc=trunc)
	end

	ws, Δuu, Δdd, Δud, Δdu = initial_weiss(beta, Nω, t)
	τs = [(i - 1) * δτ for i in 1:Nτ]          # 0 : δτ : β−δτ

	iterations = []
	for it in 1:n_iter
		t0 = time()
		# Δ(iω) → Δ(τ) kernels → BCS influence functional (single IF)
		corr = BCSCorrelationFunction(uu=Δiw_to_Δτ(Δuu, β=beta, N=Nτ),
									  dd=Δiw_to_Δτ(Δdd, β=beta, N=Nτ),
									  ud=Δiw_to_Δτ(Δud, β=beta, N=Nτ),
									  du=Δiw_to_Δτ(Δdu, β=beta, N=Nτ))
		mpsI = hybriddynamics_naive(lattice, corr, orbital=1, trunc=trunc)

		cache = environments(lattice, mpsK, mpsI)

		# normal components (per band): cached_gf_fast returns Nτ+1 points on
		# 0:δτ:β with the τ=β point fixed to 1−C(0); G(τ) = −C(τ)
		Cuu = cached_gf_fast(lattice, mpsK, mpsI; c1=false, c2=true,
							 b1=:τ, b2=:τ, band=1, cache=cache)
		Cdd = cached_gf_fast(lattice, mpsK, mpsI; c1=false, c2=true,
							 b1=:τ, b2=:τ, band=2, cache=cache)
		guu = [-c for c in Cuu]
		gdd = [-c for c in Cdd]

		# anomalous components ⟨d₁(τ) d₂(0)⟩ / ⟨d₂(τ) d₁(0)⟩, point-wise
		gud = [cached_gf(lattice,
						 (ContourIndex(i, conj=false, branch=:τ, band=1),
						  ContourIndex(1, conj=false, branch=:τ, band=2)),
						 mpsK, mpsI; cache=cache) for i in 1:Nτ]
		gdu = [cached_gf(lattice,
						 (ContourIndex(i, conj=false, branch=:τ, band=2),
						  ContourIndex(1, conj=false, branch=:τ, band=1)),
						 mpsK, mpsI; cache=cache) for i in 1:Nτ]

		# τ = β endpoints: the anomalous point-wise series covers
		# 0:δτ:β−δτ, so append F(β) = −F(0) (normal series already
		# contains the τ = β point)
		push!(gud, -gud[1])
		push!(gdu, -gdu[1])

		# G(τ) → G(iω): interpolate onto a fine τ grid, then Gτ_to_Giw
		function to_iw(gτ)
			interp = linear_interpolation(0:δτ:beta, gτ)
			return Gτ_to_Giw(interp.(0:δτ_fine:beta), β=beta, n=Nω)
		end
		Giwuu, Giwdd = to_iw(guu), to_iw(gdd)
		Giwud, Giwdu = to_iw(gud), to_iw(gdu)

		# Bethe self-consistency (Nambu, reference-code convention)
		Δuu_new = [t^2 * g for g in Giwuu]
		Δdd_new = [-t^2 * conj(g) for g in Giwdd]
		Δud_new = [t^2 * g for g in Giwud]
		Δdu_new = [t^2 * g for g in Giwdu]
		err = norm(vcat(Δuu_new - Δuu, Δdd_new - Δdd,
						Δud_new - Δud, Δdu_new - Δdu)) /
			  norm(vcat(Δuu, Δdd, Δud, Δdu))
		@printf "  iter %2d: norm|Δ′ − Δ|/norm|Δ| = %.3e   (%.1f s)\n" it err (time() - t0)

		push!(iterations, Dict(
			"iter" => it,
			"ws" => ws,
			"Delta_uu" => [[real(d), imag(d)] for d in Δuu],
			"Delta_dd" => [[real(d), imag(d)] for d in Δdd],
			"Delta_ud" => [[real(d), imag(d)] for d in Δud],
			"Delta_du" => [[real(d), imag(d)] for d in Δdu],
			"G_tau_uu" => guu, "G_tau_dd" => gdd,
			"G_tau_ud" => gud, "G_tau_du" => gdu,
			"tau" => vcat(τs, beta),
			"G_iw_uu" => [[real(g), imag(g)] for g in Giwuu],
			"G_iw_dd" => [[real(g), imag(g)] for g in Giwdd],
			"G_iw_ud" => [[real(g), imag(g)] for g in Giwud],
			"G_iw_du" => [[real(g), imag(g)] for g in Giwdu]))

		Δuu = (1 - mix) * Δuu + mix * Δuu_new
		Δdd = (1 - mix) * Δdd + mix * Δdd_new
		Δud = (1 - mix) * Δud + mix * Δud_new
		Δdu = (1 - mix) * Δdu + mix * Δdu_new
	end

	out = Dict("solver" => "GTEMPO", "beta" => beta, "U" => U, "t" => t,
			   "mu_imp" => μ_imp, "nw" => Nω, "ntau" => Nτ, "dtau" => δτ,
			   "mix" => mix, "n_iter" => n_iter, "iterations" => iterations)
	open(joinpath(outdir, "gtempo_nambu_b$(round(Int,beta))_U$(round(Int,U)).json"), "w") do f
		JSON.print(f, out)
	end
	println("written gtempo_nambu_b$(round(Int,beta))_U$(round(Int,U)).json")
end

for set in sets_final
	run_dmft(set.beta, set.U)
end
println("all done")
