
function differentialinfluencefunctional(lattice::RealGrassmannLattice{O}, corr::RealCorrelationFunction, alg::ExactTranslationInvariantIF) where O
	if !(OrderingStyle(lattice) isa _AllowedRealGrassmannOrdering)
		lattice2 = similar(lattice, ordering = A1Ā1a1ā1B1B̄1b1b̄1())
		mps = _differentialinfluencefunctional(lattice2, corr, alg)
		_, mps2 = changeordering(O, lattice2, mps, trunc=alg.algmult.trunc)
		return mps2
	else
		return _differentialinfluencefunctional(lattice, corr, alg)
	end
end


function _differentialinfluencefunctional(lattice::RealGrassmannLattice{<:_AllowedRealGrassmannOrdering}, corr::RealCorrelationFunction, alg::ExactTranslationInvariantIF)
	@assert lattice.bands == 1
	ηs = _get_signed_corr(lattice, corr, 1)

	# get WII for each exponential decay term
	f12 = ((:+, :+), (:-, :+), (:+, :-), (:-, :-))
	mpss = map(1:4) do n
		mpoj = ti_mpotensor(ηs[n], alg.algexpan).Os
		num = size(mpoj, 1)
		map(2:num-1) do i
			t = SchurMPOTensor(mpoj[[1,i,num], [1,i,num]])
			(i != 2) && (t[1,3] = 0.0)
			expt, λ, b, c = exact_WII(t)
			res1 = MPO(MPOHamiltonian([expt, expt, expt]))

			t3 = MPOHamiltonian([t, t, t])
			expt3 = timeevompo(t3, 1, WII(tol=1.0e-14, maxiter=100000))
			res2 = MPO(expt3)
			@assert norm([norm(res2.data[i] - res1.data[i]) for i in 1:3]) < 1e-14
			
			# fit to lattice
			(_fit_to_lattice(lattice, res1, I2, 1, f12[n]...) * vacuumstate(lattice), λ, b, c)
		end
	end

	(alg.verbosity >= 1) && println(length.(mpss), " terms to be multiply...")
	mpss = vcat(mpss...)

	# multiply together
	multorder = alg.multorder == :default ? :λSM : alg.multorder
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

	mps = mpss[1]
	for i in 2:length(mpss)
		t = @elapsed mps = mult(mps, mpss[i], alg.algmult)
		(alg.verbosity >= 2) && println("$i of $(length(mpss)) takes $t seconds, result mps of bond dimension: ", bond_dimension(mps))
	end
	return mps
end
