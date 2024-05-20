"""
	zoomout(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, scaling::Int = 10)

zoom out gmps by scaling
"""
function zoomout(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, scaling::Int)
	@assert scaling >= 1
	@assert length(gmps) == length(lattice)
	isa(lattice, MixedGrassmannLattice) && throw(ArgumentError("zoomout does not support MixedGrassmannLattice currently"))
	(scaling == 1) && return gmps
	@assert lattice.N % scaling == 0
	N2 = div(lattice.N, scaling)
	lattice_n = zoomout(lattice, scaling) 
	@assert lattice_n.N * scaling == lattice.N
	gmps_n = vacuumstate(lattice_n)
	if LayoutStyle(lattice) isa TimeLocalLayout
		# the 0th band
		posa, posb = band_boundary(lattice, 0)
		posa_n, posb_n = band_boundary(lattice_n, 0)
		for (idx, idx_n) in zip(posa:posb, posa_n:posb_n)
			gmps_n[idx_n] = gmps[idx] * _scaling(gmps)
		end
		posa, posb = band_boundary(lattice, lattice.k)
		posa_n, posb_n = band_boundary(lattice_n, lattice_n.k)
		for (idx, idx_n) in zip(posa:posb, posa_n:posb_n)
			gmps_n[idx_n] = gmps[idx] * _scaling(gmps)
		end
		# println("the 0-th band, posa=$(posa), posb=$(posb), posa_n=$(posa_n), posb_n=$(posb_n)")
		# the last band
		for j in lattice_n.N:-1:1
			tmp2 = _contract_band(gmps, lattice, j*scaling)
			for kk in j*scaling-1:-1:(j-1)*scaling+2
				tmp2 = tmp2 * _contract_band(gmps, lattice, kk)
			end
			posa, posb = band_boundary(lattice, (j-1)*scaling+1)
			posa_n, posb_n = band_boundary(lattice_n, j)
			for (idx, idx_n) in zip(posa:posb, posa_n:posb_n)
				gmps_n[idx_n] = gmps[idx] * _scaling(gmps)
			end
			# println("the $j-th band, posa=$(posa), posb=$(posb), posa_n=$(posa_n), posb_n=$(posb_n)")
			# println(space_r(gmps_n[posa_n-1])', " ", space_l(gmps_n[posa_n]))
			# @assert space_l(gmps_n[posa_n]) == space_r(gmps_n[posa_n-1])'
			# the rightmost one in the band should absord tmp2
			@tensor tmp3[1,3;4] := tmp2[1,2] * GrassmannTensorMap(gmps_n[posa_n])[2,3,4] 
			gmps_n[posa_n] = get_data(tmp3)
		end
	elseif LayoutStyle(lattice) isa BranchLocalLayout
		# the 0th band
		posa, posb = band_boundary(lattice, 0)
		posa_n, posb_n = band_boundary(lattice_n, 0)
		for (idx, idx_n) in zip(posa:posb, posa_n:posb_n)
			gmps_n[idx_n] = gmps[idx] * _scaling(gmps)
		end
		for forward in (:+, :-)
			posa, posb = band_boundary(lattice, lattice.k, branch=forward)
			posa_n, posb_n = band_boundary(lattice_n, lattice_n.k, branch=forward)
			for (idx, idx_n) in zip(posa:posb, posa_n:posb_n)
				gmps_n[idx_n] = gmps[idx] * _scaling(gmps)
			end
		end
		for j in lattice_n.N:-1:1
			tmp2 = _contract_band(gmps, lattice, j*scaling, branch=:+)
			for kk in j*scaling-1:-1:(j-1)*scaling+2
				tmp2 = tmp2 * _contract_band(gmps, lattice, kk, branch=:+)
			end
			posa, posb = band_boundary(lattice, (j-1)*scaling+1, branch=:+)
			posa_n, posb_n = band_boundary(lattice_n, j, branch=:+)
			for (idx, idx_n) in zip(posa:posb, posa_n:posb_n)
				gmps_n[idx_n] = gmps[idx] * _scaling(gmps)
			end
			@tensor tmp3[1,3;4] := tmp2[1,2] * GrassmannTensorMap(gmps_n[posa_n])[2,3,4] 
			gmps_n[posa_n] = get_data(tmp3)
		end	
		for j in 1:lattice_n.N
			posa, posb = band_boundary(lattice, (j-1)*scaling+1, branch=:-)
			posa_n, posb_n = band_boundary(lattice_n, j, branch=:-)
			for (idx, idx_n) in zip(posa:posb, posa_n:posb_n)
				gmps_n[idx_n] = gmps[idx] * _scaling(gmps)
			end
			tmp2 = _contract_band(gmps, lattice, (j-1)*scaling+2, branch=:-)
			for kk in (j-1)*scaling+3:j*scaling
				tmp2 = tmp2 * _contract_band(gmps, lattice, kk, branch=:-)
			end
			@tensor tmp3[1,2;4] := GrassmannTensorMap(gmps_n[posb_n])[1,2,3] * tmp2[3,4]
			gmps_n[posb_n] = get_data(tmp3)
		end					
	else
		throw(ArgumentError("zoomout not implemented for LayoutStyle $(LayoutStyle(lattice))"))
	end
	return gmps_n
end
zoomout(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice; scaling::Int=10) = zoomout(gmps, lattice, scaling)

function zoomout(lattice::RealGrassmannLattice, scaling::Int)
	# (LayoutStyle(lattice) isa TimeLocalLayout) || throw(ArgumentError("zoomout only implemented for TimeLocalLayout"))
	@assert lattice.N % scaling == 0
	return similar(lattice, N=div(lattice.N, scaling), δt=scaling*lattice.δt)
end
function zoomout(lattice::ImagGrassmannLattice, scaling::Int)
	# (LayoutStyle(lattice) isa TimeLocalLayout) || throw(ArgumentError("zoomout only implemented for TimeLocalLayout"))
	@assert lattice.N % scaling == 0
	return similar(lattice, N=div(lattice.N, scaling), δτ=scaling*lattice.δτ)
end 
zoomout(lattice::AbstractGrassmannLattice; scaling::Int=10) = zoomout(lattice, scaling)
function zoomin(lattice::RealGrassmannLattice, scaling::Int)
	# (LayoutStyle(lattice) isa TimeLocalLayout) || throw(ArgumentError("zoomin only implemented for TimeLocalLayout"))
	similar(lattice, N=lattice.N * scaling, δt=lattice.δt/scaling)
end 
function zoomin(lattice::ImagGrassmannLattice, scaling::Int)
	# (LayoutStyle(lattice) isa TimeLocalLayout) || throw(ArgumentError("zoomin only implemented for TimeLocalLayout"))
	similar(lattice, N=lattice.N * scaling, δτ=lattice.δτ/scaling)
end 
zoomin(lattice::AbstractGrassmannLattice; scaling::Int=10) = zoomin(lattice, scaling)

function accsysdynamics(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian; scaling::Int=10, trunc::TruncationScheme=DefaultKTruncation, kwargs...)
	lattice_scaling = zoomin(lattice, scaling)
	gmps = sysdynamics(lattice_scaling, model; trunc = trunc, kwargs...)
	gmps2 = zoomout(gmps, lattice_scaling, scaling=scaling)
	@assert length(gmps2) == length(lattice)
	return canonicalize!(gmps2, alg=Orthogonalize(trunc=trunc))
end