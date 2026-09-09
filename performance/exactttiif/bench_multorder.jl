# Benchmark: effect of the ExactTTIIF `multorder` on accuracy and cost
#
# Run from the repo root:
#   julia --project=. performance/exactttiif/bench_multorder.jl
#   julia --project=. performance/exactttiif/bench_multorder.jl --quick   # smoke run
#
# Setup: single-band noninteracting impurity (ToulouseIM) coupled to a
# multi-mode exponential bath (6 modes). The Prony expansion of the hybridization
# then contains many decay terms (about 20 per branch), so the order in which
# ExactTTIIF multiplies the exact term-MPOs matters through the intermediate
# SVD truncations. The ED solution of the few-mode model is the exact reference.
#
# Outputs (written next to this file):
#   results.csv  one row per (contour, multorder, D)

using GTEMPO
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "..", "..", "test", "util.jl"))
include(joinpath(@__DIR__, "..", "..", "test", "normalbath", "util.jl"))
const QUICK = "--quick" in ARGS

# ----------------------------- setup --------------------------------------
const ws = collect(0.5:0.5:3.0)                                   # 6 bath modes
const αs = [0.5, 0.42, 0.35, 0.28, 0.22, 0.16]
const β = 2.0
const μ = 0.3

const orders = (:λLM, :λSM, :αLM, :αSM, :no)
const Ds = QUICK ? (20,) : (20, 40, 80)
const Dref = 160                                                  # converged reference

bath = fermionicbath(DiscreteSpectrum(ws, αs), β=β)
model = ToulouseIM(μ=μ)
H, a, adag, H0 = singlemode_ed(μ=μ, U=0, bathspecs=collect(zip(ws, αs)))

# imaginary-time contour
const δτ = 0.125
const Nτ = round(Int, β / δτ)
lat_τ = GrassmannLattice(N=Nτ, δτ=δτ, contour=:imag)
corr_τ = correlationfunction(bath, lat_τ)
gτ_ed_ref = gτ_ed(H, a, adag, 0:δτ:(Nτ*δτ), β)

# real-time contour
const δt = 0.125
const Nt = 8
lat_t = GrassmannLattice(N=Nt, δt=δt, contour=:real)
corr_t = correlationfunction(bath, lat_t)
gt_ed_ref, lt_ed_ref = greater_lesser_ed(H, a, adag, H0, 0:δt:(Nt*δt), β)

algexpan = OverDeterminedProny(n=30, tol=1.0e-8, verbosity=0)
algmult(D::Int) = SVDCompression(truncdimcutoff(D=D, ϵ=1.0e-12))
truncscheme(D::Int) = truncdimcutoff(D=D, ϵ=1.0e-12)

# converged (large-D) references for each contour, default multorder
println("building converged D=$Dref references ...")
mpsI_ref_τ = hybriddynamics(lat_τ, corr_τ, ExactTTIIF(algexpan=algexpan, algmult=algmult(Dref), multorder=:αSM))
mpsK_ref_τ = sysdynamics(lat_τ, model, trunc=truncscheme(Dref))
mpsK_ref_τ = boundarycondition!(mpsK_ref_τ, lat_τ)
gτ_conv = gτ_series(lat_τ, mpsK_ref_τ, mpsI_ref_τ)

mpsI_ref_t = hybriddynamics(lat_t, corr_t, ExactTTIIF(algexpan=algexpan, algmult=algmult(Dref), multorder=:αSM))
mpsK_ref_t = sysdynamics(lat_t, model, trunc=truncscheme(Dref))
mpsK_ref_t = boundarycondition!(mpsK_ref_t, lat_t)
mpsK_ref_t = systhermalstate!(mpsK_ref_t, lat_t, model, trunc=truncscheme(Dref), β=β)
gt_conv, lt_conv = gtlt_series(lat_t, mpsK_ref_t, mpsI_ref_t)

# ----------------------------- helpers ------------------------------------
function run_case(contour, order, D)
	multalg = algmult(D)
	alg = ExactTTIIF(algexpan=algexpan, algmult=multalg, multorder=order)
	# K (and the thermal state) always use a generous truncation, so that all
	# measured errors are due to the IF construction alone
	tscheme = truncscheme(Dref)
	if contour == :imag
		t_build = @elapsed mpsI = hybriddynamics(lat_τ, corr_τ, alg)
		bond = maximum(bond_dimension(mpsI))
		mpsK = sysdynamics(lat_τ, model, trunc=tscheme)
		mpsK = boundarycondition!(mpsK, lat_τ)
		g = gτ_series(lat_τ, mpsK, mpsI)
		return (t_build, bond, NaN, NaN, relerr(g, gτ_ed_ref), relerr(g, gτ_conv))
	else
		t_build = @elapsed mpsI = hybriddynamics(lat_t, corr_t, alg)
		bond = maximum(bond_dimension(mpsI))
		mpsK = sysdynamics(lat_t, model, trunc=tscheme)
		mpsK = boundarycondition!(mpsK, lat_t)
		mpsK = systhermalstate!(mpsK, lat_t, model, trunc=tscheme, β=β)
		gt, lt = gtlt_series(lat_t, mpsK, mpsI)
		return (t_build, bond, relerr(gt, gt_ed_ref), relerr(lt, lt_ed_ref),
				relerr(gt, gt_conv), relerr(lt, lt_conv))
	end
end

# ----------------------------- sweep --------------------------------------
path = joinpath(@__DIR__, "results.csv")
open(path, "w") do io
	println(io, "contour,multorder,D,build_time_s,if_bond,relerr_gt_ed,relerr_lt_ed,relerr_gtau_ed,relerr_gt_conv,relerr_lt_conv,relerr_gtau_conv")
	for contour in (:imag, :real)
		for D in Ds
			for order in orders
				r = run_case(contour, order, D)
				if contour == :imag
					println(io, "$contour,$order,$D,$(r[1]),$(r[2]),,,,$(r[5]),,$(r[6])")
					@printf("%-5s multorder=%-4s D=%-3d  t=%7.1fs  bond=%-4d  relerr(gτ) ED=%.3e conv=%.3e\n",
							contour, order, D, r[1], r[2], r[5], r[6])
				else
					println(io, "$contour,$order,$D,$(r[1]),$(r[2]),$(r[3]),$(r[4]),,$(r[5]),$(r[6]),")
					@printf("%-5s multorder=%-4s D=%-3d  t=%7.1fs  bond=%-4d  relerr(gt,lt) ED=%.3e/%.3e conv=%.3e/%.3e\n",
							contour, order, D, r[1], r[2], r[3], r[4], r[5], r[6])
				end
				flush(io)
			end
		end
	end
end
println("results written to ", path)
