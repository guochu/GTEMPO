function accsysdynamics_fast(lattice::ImagGrassmannLattice{O}, model::AbstractImpurityHamiltonian; scaling::Int=10, trunc::TruncationScheme=DefaultKTruncation) where O
	lattice2 = similar(lattice, ordering=A1B1B̄1Ā1())
	x = _accsysdynamics_fast(lattice2, model, scaling=scaling, trunc=trunc)
	return changeordering(O, lattice2, x, trunc=trunc)[2]
end

_accsysdynamics_fast(lattice::ImagGrassmannLattice{<:A1B1B̄1Ā1}, model::AbstractImpurityHamiltonian; kwargs...) = _accsysdynamics_fast_timelocal(lattice, model; kwargs...)

# accsysdynamics_fast(lattice::ImagGrassmannLattice, model::AbstractImpurityHamiltonian; kwargs...) = error("accsysdynamics_fast only support TimeLocalLayout for ImagGrassmannLattice")
function _accsysdynamics_fast_timelocal(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian; 
										scaling::Int=10, trunc::TruncationScheme=DefaultKTruncation, kwargs...)
	lattice1 = similar(lattice, N=1)
	lattice_n = zoomin(lattice1, scaling=scaling)
	lattice_n_1 = similar(lattice_n, N=1)
	gmps_n_1 = sysdynamics(lattice_n_1, model; trunc=trunc, kwargs...)
	posa = get_left(lattice_n_1, 2)
	posb = get_right(lattice_n_1, 1)
	# println(bond_dimensions(gmps_n_1))
	# println("posa=", posa, " posb=", posb)
	function check_gmps1(x, tol)
		# println( space_l(x[posa]), " ", space_r(x[posb]))
		@assert isoneunit(space_l(x[posa]) )
		@assert isoneunit(space_r(x[posb]))
		s = 1.
		for i in 1:posa-1
			@assert isoneunit(space_l(x[i])) && isoneunit(space_r(x[i]))
			xi = only(x[i][(Z2Irrep(0), Z2Irrep(0), Z2Irrep(0))])
			s *= xi
			@assert abs(abs(xi) - 1) < tol
			# @assert abs(TK.scalar(x[i]) - 1) < tol
		end
		for i in posb+1:length(x)
			@assert isoneunit(space_l(x[i])) && isoneunit(space_r(x[i]))
			xi = only(x[i][(Z2Irrep(0), Z2Irrep(0), Z2Irrep(0))])
			s *= xi
			@assert abs(abs(xi) - 1) < tol
			# @assert abs(TK.scalar(x[i]) - 1) < tol
		end
		return s
	end
	s0 = check_gmps1(gmps_n_1, 1.0e-10)
	_s = _scaling(gmps_n_1)^(length(gmps_n_1) / (posb-posa+1))
	for i in posa:posb
		gmps_n_1[i] *= _s
	end
	gmps_n_1[posa] *= s0

	gmps_n = vacuumstate(lattice_n)
	for i in 1:lattice_n.N
		posa_n = get_left(lattice_n, i+1)
		posb_n = get_right(lattice_n, i)		
		for (idx, idx_n) in zip(posa:posb, posa_n:posb_n)
			gmps_n[idx_n] = gmps_n_1[idx] 
		end
	end
	gmps1 = zoomout(gmps_n, lattice_n, scaling=scaling)
	@assert length(gmps1) == length(lattice_n_1)
	s0 = check_gmps1(gmps1, 1.0e-10)
	_s = _scaling(gmps1)^(length(gmps1) / (posb-posa+1))
	for i in posa:posb
		gmps1[i] *= _s
	end
	gmps1[posa] *= s0
	gmps = vacuumstate(lattice)
	for i in 1:lattice.N
		posa_n = get_left(lattice, i+1)
		posb_n = get_right(lattice, i)		
		for (idx, idx_n) in zip(posa:posb, posa_n:posb_n)
			gmps[idx_n] = gmps1[idx]
		end		
	end

	# return canonicalize!(gmps, alg=Orthogonalize(trunc=trunc))
	return gmps
end

get_left(lattice::ImagGrassmannLattice{<:A1B1B̄1Ā1}, j::Int) = index(lattice, j, conj=true, band=lattice.bands)
get_right(lattice::ImagGrassmannLattice{<:A1B1B̄1Ā1}, j::Int) = index(lattice, j, conj=false, band=lattice.bands)