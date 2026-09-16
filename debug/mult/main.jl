# Debug: IDMRG1-based mult vs SVDCompression-based mult
#
# Run from the repo root with:
#   JULIA_PROBE_LIBSTDCXX=0 JULIA_DEPOT_PATH="/tmp/jdepot:/home/guochu/.julia" \
#   julia --startup-file=no --project=. debug/mult/main.jl
#
# Questions:
#   1. Are the outputs of mult(x, y, DMRG1(...)) and mult(x, y, SVDCompression(...))
#      strictly identical — even when the set bond dimension is too small and the
#      results are inaccurate?
#   2. Does the IDMRG1 loss function (kvals) decay monotonically throughout?
#
# Reports (also written to debug/mult/report.txt):
#   - per-trial relative errors against the (almost) untruncated product,
#   - the IDMRG1-vs-SVDCompression distance relative to the exact product,
#   - the full per-sweep loss sequences with monotonicity verdicts.

using GTEMPO
using LinearAlgebra: norm
using Dates
using Printf

# ---------------------------------------------------------------- helpers

"""
Run the IDMRG1 mult loop exactly like `GTEMPO.iterativemult`, but return the
loss sequences (per-sweep losses — the new `iterative_compute!` returns the
*last* residual of each sweep — plus the raw per-site residuals of every
sweep) together with the final GMPS.
"""
function idmrg_run(x::GrassmannMPS, y::GrassmannMPS, alg)
	if alg.initguess == :svd
		z = GTEMPO._svd_guess(x, y, alg.D)
	elseif alg.initguess == :rand
		z = randomgmps(promote_type(scalartype(x), scalartype(y)), length(x), D=alg.D)
	else
		error("unsupported initguess")
	end
	cache = GTEMPO.mult_cache(z, x, y)
	sweep_losses = Float64[]             # per-sweep loss (last residual; what iterative_compute! returns)
	sweep_left = Vector{Float64}[]       # per-site residuals of each left sweep  (sites 1…L-1)
	sweep_right = Vector{Float64}[]      # per-site residuals of each right sweep (sites L…2)
	loss_prev, iter, delta = NaN, 0, 2 * alg.tol
	while (iter < alg.maxiter) && (delta >= alg.tol)
		kl = GTEMPO.leftsweep!(cache, alg)
		kr = GTEMPO.rightsweep!(cache, alg)
		push!(sweep_left, collect(kl))
		push!(sweep_right, collect(kr))
		loss_cur = kr[end]
		delta = (iter == 0) ? 2 * alg.tol :
			abs(loss_cur - loss_prev) / max(loss_cur, loss_prev, eps(Float64))
		loss_prev = loss_cur
		push!(sweep_losses, loss_cur)
		iter += 1
	end
	GTEMPO.finalize!(cache, alg)
	setscaling!(cache.z, scaling(x) * scaling(y))
	GTEMPO._rescaling!(cache.z)
	return sweep_losses, sweep_left, sweep_right, cache.z
end

"strict monotone (non-increasing) check; returns (verdict, first violation index)"
function monotone_decreasing(v::Vector{Float64})
	for i in 2:length(v)
		if v[i] > v[i-1]
			return (false, i)
		end
	end
	return (true, 0)
end

"strict monotone (non-decreasing) check up to relative float jitter; returns (verdict, first violation index)"
function monotone_increasing(v::Vector{Float64}; rtol::Float64=1.0e-9)
	for i in 2:length(v)
		if v[i] < v[i-1] * (1 - rtol)
			return (false, i)
		end
	end
	return (true, 0)
end

function run_case(io, name, L, Dx, Dy, D; ntrial=3, maxiter=60)
	trunc = truncdimcutoff(D=D, ϵ=1.0e-14)
	truncbig = truncdimcutoff(D=4096, ϵ=1.0e-16)
	alg_svd = SVDCompression(trunc)
	alg_id = DMRG1(trunc, maxiter=maxiter, tol=1.0e-14, verbosity=-1)
	tag = replace(lowercase(split(name, " (")[1]), " " => "_")
	csv = open(joinpath(@__DIR__, "kvals_$(tag)_D$(D).csv"), "w")
	println(csv, "trial,sweep,seg,pos,site,residual")

	println(io, "="^78)
	println(io, "case: $name  (L=$L, Dx=Dy=$Dx, set bond dimension D=$D, trials=$ntrial)")
	println(io, "="^78)

	for trial in 1:ntrial
		x = randomgmps(L, D=Dx)
		y = randomgmps(L, D=Dy)
		z_exact = mult(x, y, trunc=truncbig)          # almost untruncated reference
		z_svd = mult(x, y, alg_svd)
		losses, lefts, rights, z_id = idmrg_run(x, y, alg_id)
		# dump all kvals: left segment sites 1..L-1, right segment sites L..2
		for (si, kl) in enumerate(lefts), (p, v) in enumerate(kl)
			println(csv, "$trial,$si,left,$p,$p,$v")
		end
		for (si, kr) in enumerate(rights), (p, v) in enumerate(kr)
			println(csv, "$trial,$si,right,$p,$(length(kr) + 1 - p),$v")
		end

		nrm = norm(z_exact)
		e_svd = distance(z_svd, z_exact) / nrm
		e_id = distance(z_id, z_exact) / nrm
		d_id_svd = distance(z_id, z_svd) / nrm
		identical = (z_id == z_svd)                  # bit-wise equality of the GMPS objects
		close01 = d_id_svd < 1.0e-13

		@printf(io, "trial %d: err(SVD)=%9.2e  err(IDMRG1)=%9.2e  |IDMRG1-SVD|/|exact|=%9.2e  bit-identical=%-5s numerically-identical(<1e-13)=%-5s\n",
			trial, e_svd, e_id, d_id_svd, identical, close01)

		# --- loss monotonicity: per-sweep losses (last residual of each sweep,
		#     what iterative_compute! returns); theoretically monotone
		#     non-decreasing (rise to the fixed-point plateau) ---
		mono, bad = monotone_increasing(losses)
		@printf(io, "  per-sweep losses (%d sweeps): %s\n", length(losses), mono ? "monotonically non-decreasing" : "NOT monotone (first decrease at sweep $bad: $(losses[bad-1]) -> $(losses[bad]))")
		@printf(io, "    sweeps: %s\n", join([@sprintf("%.3e", v) for v in losses], "  "))

		# --- within-sweep kvals = norm(mpsj): theoretically monotone non-decreasing,
		#     both within each sweep and across sweeps ---
		linc = count(monotone_increasing(kl)[1] for kl in lefts)
		rinc = count(monotone_increasing(kr)[1] for kr in rights)
		whole = vcat([vcat(kl, kr) for (kl, kr) in zip(lefts, rights)]...)
		winc, wbad = monotone_increasing(whole)
		@printf(io, "  kvals norm(mpsj): left segment non-decreasing in %d/%d sweeps, right in %d/%d, whole run %s%s\n",
			linc, length(lefts), rinc, length(rights), winc ? "non-decreasing" : "DECREASING at step $wbad",
			winc ? "" : @sprintf(" (%.6e -> %.6e)", whole[wbad-1], whole[wbad]))

		# show the residual pattern of the first sweep
		@printf(io, "    sweep 1 left  (sites 1..%d): %s\n", length(lefts[1]),
			join([@sprintf("%.3e", v) for v in lefts[1]], "  "))
		@printf(io, "    sweep 1 right (sites %d..2): %s\n", length(rights[1]) + 1,
			join([@sprintf("%.3e", v) for v in rights[1]], "  "))
	end
	close(csv)
	println(io, "")
end

# ---------------------------------------------------------------- run

open(joinpath(@__DIR__, "report.txt"), "w") do io
	println(io, "# IDMRG1 vs SVDCompression mult — debug report")
	println(io, "# generated by debug/mult/main.jl, ", Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))
	println(io, "")

	# D large enough for the exact product (result should be essentially exact
	# for both algorithms — are they even numerically identical then?)
	run_case(io, "sufficient bond dimension", 6, 4, 4, 64; ntrial=3)

	# D clearly too small: the truncated products are inaccurate; the question is
	# whether the two algorithms still produce the *same* inaccurate state
	run_case(io, "insufficient bond dimension", 6, 4, 4, 2; ntrial=3)
	run_case(io, "insufficient bond dimension", 8, 8, 8, 3; ntrial=3)
end
println("report written to debug/mult/report.txt")
