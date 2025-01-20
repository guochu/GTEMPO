function accsysdynamics_fast(lattice::RealGrassmannLattice{O}, model::AbstractImpurityHamiltonian; 
							 scaling::Int=10, trunc::TruncationScheme=DefaultKTruncation, kwargs...) where {O}
	if LayoutStyle(lattice) isa TimeLocalLayout
		lattice2 = similar(lattice, ordering=A1B1ā1b̄1Ā1B̄1a1b1())
	else
		lattice2 = similar(lattice, ordering=A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2())
	end
	x = _accsysdynamics_fast(lattice2, model; scaling=scaling, trunc=trunc, kwargs...)
	return changeordering(O, lattice2, x, trunc=trunc)[2]
end

_accsysdynamics_fast(lattice::RealGrassmannLattice{<:A1B1ā1b̄1Ā1B̄1a1b1}, model::AbstractImpurityHamiltonian; kwargs...) = _accsysdynamics_fast_timelocal(lattice, model; kwargs...)

# accsysdynamics_fast(lattice::RealGrassmannLattice, model::AbstractImpurityHamiltonian; kwargs...) = error("accsysdynamics_fast only support BranchLocalLayout for RealGrassmannLattice")
function _accsysdynamics_fast(lattice::RealGrassmannLattice{<:A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2}, model::AbstractImpurityHamiltonian; 
								scaling::Int=10, trunc::TruncationScheme=DefaultKTruncation, kwargs...)
	lattice1 = similar(lattice, N=1)
	lattice_n = zoomin(lattice1, scaling=scaling)
	lattice_n_1 = similar(lattice_n, N=1)
	gmps_n_1 = sysdynamics(lattice_n_1, model; trunc=trunc, kwargs...)
	# println("here: ", bond_dimensions(gmps_n_1))
	# println(bond_dimensions(gmps_n_1))
	posa_f = get_left(lattice_n_1, 2, branch=:+)
	posb_f = get_right(lattice_n_1, 1, branch=:+)
	posa_b = get_left(lattice_n_1, 1, branch=:-)
	posb_b = get_right(lattice_n_1, 2, branch=:-)
	# println(bond_dimensions(gmps_n_1))
	# println("posa_f=", posa_f, " posb_f=", posb_f, " posa_b=", posa_b, " posb_b=", posb_b)
	function check_gmps1(x, tol)
		# println( space_l(x[posa_f]), " ", space_r(x[posb_f]), " ", space_l(x[posa_b]), " ", space_r(x[posb_b]))
		@assert isoneunit(space_l(x[posa_f]))
		@assert isoneunit(space_r(x[posb_f]))
		@assert isoneunit(space_l(x[posa_b]))
		@assert isoneunit(space_r(x[posb_b]))
		s = 1.
		for i in Iterators.flatten((1:posa_f-1, posb_f+1:posa_b-1, posb_b+1:length(x)))
			@assert isoneunit(space_l(x[i])) && isoneunit(space_r(x[i]))
			xi = only(x[i][(Z2Irrep(0), Z2Irrep(0), Z2Irrep(0))])
			s *= xi
			@assert abs(abs(xi) - 1) < tol
			# @assert abs(TK.scalar(x[i]) - 1) < tol
		end
		return s
	end
	s0 = check_gmps1(gmps_n_1, 1.0e-10)
	_s = _scaling(gmps_n_1)^(length(gmps_n_1) / (posb_f-posa_f+1 + posb_b - posa_b+1))
	for i in Iterators.flatten((posa_f:posb_f, posa_b:posb_b))
		gmps_n_1[i] *= _s
	end
	gmps_n_1[posa_f] *= s0

	gmps_n = vacuumstate(lattice_n)
	for i in 1:lattice_n.N
		posa_n = get_left(lattice_n, i+1, branch=:+)
		posb_n = get_right(lattice_n, i, branch=:+)		
		for (idx, idx_n) in zip(posa_f:posb_f, posa_n:posb_n)
			gmps_n[idx_n] = gmps_n_1[idx] 
		end
		posa_n = get_left(lattice_n, i, branch=:-)
		posb_n = get_right(lattice_n, i+1, branch=:-)		
		for (idx, idx_n) in zip(posa_b:posb_b, posa_n:posb_n)
			gmps_n[idx_n] = gmps_n_1[idx] 
		end		
	end
	gmps1 = zoomout(gmps_n, lattice_n, scaling=scaling)
	@assert length(gmps1) == length(lattice_n_1)
	s0 = check_gmps1(gmps1, 1.0e-10)
	_s = _scaling(gmps1)^(length(gmps1) / (posb_f-posa_f+1 + posb_b - posa_b+1))
	for i in Iterators.flatten((posa_f:posb_f, posa_b:posb_b))
		gmps1[i] *= _s
	end
	gmps1[posa_f] *= s0
	gmps = vacuumstate(lattice)
	for i in 1:lattice.N
		posa_n = get_left(lattice, i+1, branch=:+)
		posb_n = get_right(lattice, i, branch=:+)		
		for (idx, idx_n) in zip(posa_f:posb_f, posa_n:posb_n)
			gmps[idx_n] = gmps1[idx]
		end		
		posa_n = get_left(lattice, i, branch=:-)
		posb_n = get_right(lattice, i+1, branch=:-)		
		for (idx, idx_n) in zip(posa_b:posb_b, posa_n:posb_n)
			gmps[idx_n] = gmps1[idx]
		end		
	end

	# return canonicalize!(gmps, alg=Orthogonalize(trunc=trunc))
	return gmps
end

get_left(lattice::RealGrassmannLattice{<:A1B1ā1b̄1Ā1B̄1a1b1}, j::Int) = index(lattice, j, conj=true, band=1, branch=:+)
get_right(lattice::RealGrassmannLattice{<:A1B1ā1b̄1Ā1B̄1a1b1}, j::Int) = index(lattice, j, conj=true, band=lattice.bands, branch=:-)

get_left(lattice::RealGrassmannLattice{<:A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2}, j::Int; branch::Symbol) = index(lattice, j, conj=true, band=lattice.bands, branch=branch)
get_right(lattice::RealGrassmannLattice{<:A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2}, j::Int; branch::Symbol) = index(lattice, j, conj=false, band=lattice.bands, branch=branch)

