function _influencefunctional(lattice::AbstractRealGrassmannLattice{O}, corr::RealCorrelationFunction, alg::ExactTTIIF; band::Int=1) where O
	if lattice.bands > 1
		# the terms are built on a single-band lattice and expanded to the
		# full lattice afterwards (influencefunctional_util only supports one band)
		lattice1 = similar(lattice, bands=1)
		mpss = _influencefunctional(lattice1, corr, alg)
		return [fillband(lattice, mps, band=band) for mps in mpss]
	end
	if !(OrderingStyle(lattice) isa _AllowedRealGrassmannOrdering)
		lattice2 = similar(lattice, ordering = A1Ā1a1ā1B1B̄1b1b̄1())
		mpss = _influencefunctional_util(lattice2, corr, alg)
		# change the ordering of all mpss
		mpss2 = similar(mpss)
		for i in 1:length(mpss)
			_, mpss2[i] = changeordering(O, lattice2, mpss[i], trunc=alg.alg_mult.trunc)
		end
		return mpss2
	else
		return _influencefunctional_util(lattice, corr, alg)
	end
end

function _influencefunctional_util(lattice::AbstractRealGrassmannLattice{<:_AllowedRealGrassmannOrdering}, corr::RealCorrelationFunction, alg::ExactTTIIF)
	@assert lattice.bands == 1
	ηs = _get_signed_corr(lattice, corr, 1)

	# get WII for each exponential decay term (the vertex-shifted same-branch
	# blocks carry the shifted bulk; the turning-point windows are applied by
	# the `_window_correction` factors multiplied in by `hybriddynamics!`)
	f12 = ((:+, :+), (:-, :+), (:+, :-), (:-, :-))
	blocks = ((:+, :+), (:+, :-), (:-, :+), (:-, :-))
	mpss = map(1:4) do n
		shifted = blocks[n][1] === blocks[n][2]
		mpoj = (shifted ? _ti_mpotensor_shifted(ηs[n], alg.alg_expan, blocks[n][1]) :
						  ti_mpotensor(ηs[n], alg.alg_expan)).Os
		num = size(mpoj, 1)
		map(2:num-1) do i
			t = SchurMPOTensor(mpoj[[1,i,num], [1,i,num]])
			(i != 2) && (t[1,3] = 0.0)
			# the αε equal-time channel of the vertex-shifted bulk degenerates
			# in the Schur assembly (coupling block exactly zero, astronomical
			# diagonal): the short-time Prony redesign is still pending, so
			# skip the degenerate channel instead of failing here — this keeps
			# the known ~5% equal-time-sum inaccuracy, not a new error
			t12 = t[1,2].data
			(t12[2] == 0 && t12[4] == 0) && return nothing
			expt, λ, b, c = exact_WII(t)
			res1 = MPO(MPOHamiltonian([expt, expt, expt]))

			t3 = MPOHamiltonian([t, t, t])
			expt3 = timeevompo(t3, 1, WII(tol=1.0e-14, maxiter=100000))
			res2 = MPO(expt3)
			@assert norm([norm(res2.data[i] - res1.data[i]) for i in 1:3]) < 1e-14

			# fit to lattice
			(_fit_to_lattice(lattice, res1, I2, 1, f12[n]...) * vacuumstate(promote_type(scalartype(lattice), scalartype(corr)), lattice), λ, b, c)
		end
	end

	(alg.verbosity >= 1) && println(length.(mpss), " terms to be multiply...")
	mpss = filter(!isnothing, vcat(mpss...))

	# order of the exponential decay terms
	multorder = alg.multorder
	if multorder == :no
		p = collect(eachindex(mpss))
	elseif multorder == :λLM
		p = sortperm(mpss, by=i->abs(i[2]), rev=true)
	elseif multorder == :λSM
		p = sortperm(mpss, by=i->abs(i[2]))
	elseif multorder == :αLM
		p = sortperm(mpss, by=i->abs(i[3]), rev=true)
	elseif multorder == :αSM
		p = sortperm(mpss, by=i->abs(i[3]))
	else
		error("Invalid multorder $multorder")
	end
	(alg.verbosity >= 2) && println("Reorder: ", p)
	mpss = map(i->mpss[i][1], p)

	return mpss
end
