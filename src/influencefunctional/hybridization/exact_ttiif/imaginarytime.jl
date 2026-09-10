function _influencefunctional(lattice::ImagGrassmannLattice1Order{O}, corr::ImagCorrelationFunction, alg::ExactTTIIF; band::Int=1) where O
	if lattice.bands > 1
		# the terms are built on a single-band lattice and expanded to the
		# full lattice afterwards (influencefunctional_util only supports one band)
		lattice1 = similar(lattice, bands=1)
		mpss = _influencefunctional(lattice1, corr, alg)
		return [fillband(lattice, mps, band=band) for mps in mpss]
	end
	if !(LayoutStyle(lattice) isa TimeLocalLayout)
		lattice2 = similar(lattice, ordering = A1Ā1B1B̄1())
		mpss = _influencefunctional_util(lattice2, corr, alg)
		# change the ordering of all mpss
		mpss2 = similar(mpss)
		for i in 1:length(mpss)
			_, mpss2[i] = changeordering(O, lattice2, mpss[i], trunc=alg.algmult.trunc)
		end
		return mpss2
	else
		return _influencefunctional_util(lattice, corr, alg)
	end
end

function _influencefunctional_util(lattice::ImagGrassmannLattice1Order, corr::ImagCorrelationFunction, alg::ExactTTIIF)
	@assert lattice.bands == 1

	# get WII for each exponential decay term
	mpoj = ti_mpotensor(corr.data, alg.algexpan).Os
	num = size(mpoj, 1)
	mpss = map(2:num-1) do i
		t = SchurMPOTensor(mpoj[[1,i,num], [1,i,num]])
		(i != 2) && (t[1,3] = 0.0)
		expt, λ, b, c = exact_WII(t)
		res1 = MPO(MPOHamiltonian([expt, expt, expt]))

		t3 = MPOHamiltonian([t, t, t])
		expt3 = timeevompo(t3, 1, WII(tol=1.0e-14, maxiter=100000))
		res2 = MPO(expt3)
		@assert norm([norm(res2.data[i] - res1.data[i]) for i in 1:3]) < 1e-14

		# fit to lattice
		(_fit_to_lattice(lattice, res1, I2, 1) * vacuumstate(lattice), λ, b, c)
	end

	(alg.verbosity >= 1) && println(length(mpss), " terms to be multiply...")

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
