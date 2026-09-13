# Run from the project root with:
#   julia --project=. docs/tutorials/holstein/holstein.jl
#
# Holstein impurity model: a fermionic impurity coupled to a local phonon
# and a fermionic bath with a semicircular spectral density,
#
#   H = ϵ_d â†â + ω₀ b̂†b̂ + g â†â (b̂† + b̂)
#       + Σₖ ωₖ ĉₖ†ĉₖ + Σₖ Vₖ (â†ĉₖ + ĉₖ†â).
#
# For U = 0 the model is exactly solvable (Lang–Firsov): the phonon only
# couples to the impurity occupation n̂, so the Green's functions factorise
# into a free-impurity part times a phonon decoherence factor. The analytic
# retarded Green's function G^R(t) is obtained from ImpurityModelBase's
# `holstein_G0w_to_Gw`, the continued-fraction expansion of the Holstein
# solution at finite temperature, driven by the free-impurity (Toulouse)
# propagator G₀(ω) = [ω + iδ - ϵ_d - Δ(ω)]⁻¹. For the semicircular bath
# the hybridisation Δ(ω) is known in closed form, which makes the
# evaluation of the reference fast and robust (feeding the raw spectral
# density instead would require a nested numerical Hilbert transform at
# every frequency and leads to pathologically slow adaptive quadratures).
# We benchmark the GTEMPO electron-phonon workflow — the retarded function
# is obtained as G^R(t) = G^>(t) - G^<(t) — against this analytic solution
# and analyse the convergence with the time step δt.

using GTEMPO
using ImpurityModelBase
using LinearAlgebra: norm

# ------------------------- model parameters -------------------------
μ = 0.0            # impurity level ϵ_d = μ (ToulouseIM uses μ n̂)
β = 5.0            # inverse temperature
ω₀ = 0.8           # phonon frequency
α₀ = 0.5           # phonon coupling strength, g = √α₀
g = sqrt(α₀)       # Holstein coupling g n̂ (b̂† + b̂)
th = 1.0           # half-bandwidth of the bath; band = [-t̂, t̂]
spec_f = semicircular(t=th)  # fermionic bath spectral density

# ----------------------- analytic solution --------------------------
# Retarded hybridisation function of the semicircular bath (ρ(ε) =
# 2√(t̂²-ε²)/(πt̂²), normalised to 1): Δ(z) = 2[z - √(z²-t̂²)]/t̂², with the
# branch fixed by Im Δ(ω) ≤ 0 on the real axis. As a check, the free
# path below reproduces ImpurityModelBase's `toulouse_Gt` to ~9 digits.
delta_sem(ω::Real; th=th) = abs(ω) < th ?
	2 * (ω - im * sqrt(th^2 - ω^2)) / th^2 :
	2 * (ω - sign(ω) * sqrt(ω^2 - th^2)) / th^2

# analytic G^R(t) of the Holstein model: continued fraction over the
# free-impurity propagator, followed by one Fourier-type quadrature,
#   G^R(t) = ∫ dω/2π e^{-iωt} [G(ω) - 1/(ω+iδ)] - i .
function holstein_GRt(t; g, ω₀, ϵ_d, β, th=th, δ=1.0e-6, wmax=20.0,
					  maxiter=6)
	G0(w) = 1 / (w + im*δ - ϵ_d - delta_sem(w; th=th))
	integrand = bounded(w -> (holstein_G0w_to_Gw(G0, w; g=g, ω=ω₀, β=β,
						   maxiter=maxiter) - 1 / (w + im*δ)) * exp(-im*w*t),
						-wmax, wmax)
	return quadgkwrapper(integrand) / (2π) - im
end

# evaluate the reference once on the finest grid; the coarser δt grids
# below are sub-grids of it
tmax = 1.0
δt_fine = 0.05
ts_all = collect(0:δt_fine:tmax)
println("computing the analytic reference (holstein continued fraction) ",
		"on ", length(ts_all), " time points ...")
ex_ph = [holstein_GRt(t; g=g, ω₀=ω₀, ϵ_d=μ, β=β) for t in ts_all]
ex_free = [holstein_GRt(t; g=0.0, ω₀=ω₀, ϵ_d=μ, β=β) for t in ts_all]
println("done")

# ----------------------- GTEMPO workflow ----------------------------
# `tmax` is kept fixed while `δt` is refined; returns the retarded
# Green's function on the grid `0:δt:tmax`.
function run_gtempo(δt; tmax=1.0, with_phonon=true)
	Nt = round(Int, tmax / δt)
	trunc = truncdimcutoff(D=100, ϵ=1.0e-10)
	bands = 1
	lattice = GrassmannLattice(N=Nt, δt=δt, contour=:real, order=1, bands=bands)
	flat = FockLattice(N=Nt, δt=δt, contour=:real, order=1, bands=bands)

	# influence functional of the semicircular fermionic bath,
	# built with the TTIIF algorithm on the single-band lattice
	fbath = fermionicbath(spec_f, β=β)
	lattice1 = similar(lattice, bands=1)
	fcorr = correlationfunction(fbath, lattice1)
	alg = ExactTTIIF(algmult=SVDCompression(trunc))
	mpsI_e = hybriddynamics(lattice1, fcorr, alg)
	Is = [fillband(lattice, mpsI_e, band=b) for b in 1:bands]

	# bare impurity dynamics: Toulouse model (single band, U = 0)
	model = ToulouseIM(μ=μ)
	mpsK = sysdynamics(lattice, model, trunc=trunc)
	mpsK = systhermalstate!(mpsK, lattice, model, trunc=trunc, β=β)
	mpsK = boundarycondition!(mpsK, lattice, band=1, trunc=trunc)

	if with_phonon
		# influence functional of the local phonon
		pbath = bosonicbath(DiracDelta(ω=ω₀, α=α₀), β=β)
		pcorr = correlationfunction(pbath, flat)
		mpsI_p = hybriddynamics(flat, pcorr, trunc=trunc)
		# merge the phonon IF into the system state; the fermionic bath IF
		# (`Is`) stays separate and enters the observable evaluation below
		adt = reweighting!(lattice, mpsK, flat, mpsI_p, trunc=trunc)
	else
		adt = mpsK
	end

	cache = environments(lattice, adt, Is...)
	# retarded Green's function on the real-time grid: with
	# G^>(t) = -im * cached_greater(t) and G^<(t) = -im * cached_lesser(t),
	# G^R(t) = G^>(t) - G^<(t) for t > 0
	cg = [cached_greater(lattice, k, adt, Is..., band=1, cache=cache) for k in 1:Nt+1]
	cl = [cached_lesser(lattice, k, adt, Is..., band=1, cache=cache) for k in 1:Nt+1]
	return [-im * (gv - lv) for (gv, lv) in zip(cg, cl)]
end

# -------------------- δt convergence analysis -----------------------
# rel.err = ‖G_TEMPO - G_exact‖ / ‖G_exact‖, accumulated over 0:δt:tmax
println("\nδt        rel.err (phonon)   rel.err (no phonon)")
gts = Dict{Float64, Vector{ComplexF64}}()
for δt in (0.2, 0.1, 0.05)
	gt = run_gtempo(δt)
	gt0 = run_gtempo(δt; with_phonon=false)
	gts[δt] = gt
	# sub-sampling of the fine-grid analytic reference at i·δt
	sel = 1:round(Int, δt / δt_fine):length(ts_all)
	exact_ph = ex_ph[sel]
	exact_free = ex_free[sel]
	err = norm(gt - exact_ph) / norm(exact_ph)
	err0 = norm(gt0 - exact_free) / norm(exact_free)
	println(rpad(δt, 10), rpad(round(err, sigdigits=3), 19), round(err0, sigdigits=3))
end

# detailed comparison at the finest time step
println("\nt        GTEMPO (phonon)      analytic")
for (k, t) in enumerate(ts_all)
	println(rpad(t, 8), rpad(round(gts[δt_fine][k], digits=6), 21),
			round(ex_ph[k], digits=6))
end
