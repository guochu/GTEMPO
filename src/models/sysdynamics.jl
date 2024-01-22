# the functions here are only applicable for specific orderings

# There is a "root" ordering which can be used to build K accurately and efficiently
# TODO: generalize accsysdynamics to BandLocalLayout

"""
	zoomout(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, scaling::Int = 10)

zoom out gmps by scaling
"""
function zoomout(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, scaling::Int)
	@assert scaling >= 1
	@assert length(gmps) == length(lattice)
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
			gmps_n[posa_n] = @tensor tmp3[1,3;4] := tmp2[1,2] * gmps_n[posa_n][2,3,4] 
		end
	elseif LayoutStyle(lattice) isa BranchLocalLayout
		# the 0th band
		posa, posb = band_boundary(lattice, 0)
		posa_n, posb_n = band_boundary(lattice_n, 0)
		for (idx, idx_n) in zip(posa:posb, posa_n:posb_n)
			gmps_n[idx_n] = gmps[idx] * _scaling(gmps)
		end
		for forward in (true, false)
			posa, posb = band_boundary(lattice, lattice.k, forward=forward)
			posa_n, posb_n = band_boundary(lattice_n, lattice_n.k, forward=forward)
			for (idx, idx_n) in zip(posa:posb, posa_n:posb_n)
				gmps_n[idx_n] = gmps[idx] * _scaling(gmps)
			end
		end
		for j in lattice_n.N:-1:1
			tmp2 = _contract_band(gmps, lattice, j*scaling, forward=true)
			for kk in j*scaling-1:-1:(j-1)*scaling+2
				tmp2 = tmp2 * _contract_band(gmps, lattice, kk, forward=true)
			end
			posa, posb = band_boundary(lattice, (j-1)*scaling+1, forward=true)
			posa_n, posb_n = band_boundary(lattice_n, j, forward=true)
			for (idx, idx_n) in zip(posa:posb, posa_n:posb_n)
				gmps_n[idx_n] = gmps[idx] * _scaling(gmps)
			end
			gmps_n[posa_n] = @tensor tmp3[1,3;4] := tmp2[1,2] * gmps_n[posa_n][2,3,4] 
		end	
		for j in 1:lattice_n.N
			posa, posb = band_boundary(lattice, (j-1)*scaling+1, forward=false)
			posa_n, posb_n = band_boundary(lattice_n, j, forward=false)
			for (idx, idx_n) in zip(posa:posb, posa_n:posb_n)
				gmps_n[idx_n] = gmps[idx] * _scaling(gmps)
			end
			tmp2 = _contract_band(gmps, lattice, (j-1)*scaling+2, forward=false)
			for kk in (j-1)*scaling+3:j*scaling
				tmp2 = tmp2 * _contract_band(gmps, lattice, kk, forward=false)
			end
			gmps_n[posb_n] = @tensor tmp3[1,2;4] := gmps_n[posb_n][1,2,3] * tmp2[3,4]
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

function accsysdynamics(lattice::AbstractGrassmannLattice, model::AbstractImpurityModel; scaling::Int=10, trunc::TruncationScheme=DefaultKTruncation, kwargs...)
	lattice_scaling = zoomin(lattice, scaling)
	gmps = sysdynamics(lattice_scaling, model; trunc = trunc, kwargs...)
	gmps2 = zoomout(gmps, lattice_scaling, scaling=scaling)
	@assert length(gmps2) == length(lattice)
	return canonicalize!(gmps2, alg=Orthogonalize(trunc=trunc))
end

# accsysdynamics_fast(lattice::ImagGrassmannLattice, model::AbstractImpurityModel; kwargs...) = error("accsysdynamics_fast only support TimeLocalLayout for ImagGrassmannLattice")
function _accsysdynamics_fast(lattice::ImagGrassmannLattice{<:A1B1B1A1}, model::AbstractImpurityModel; scaling::Int=10, trunc::TruncationScheme=DefaultKTruncation, kwargs...)
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
function accsysdynamics_fast(lattice::ImagGrassmannLattice{O}, model::AbstractImpurityModel; scaling::Int=10, trunc::TruncationScheme=DefaultKTruncation) where O
	lattice2 = similar(lattice, ordering=A1B1B1A1())
	x = _accsysdynamics_fast(lattice2, model, scaling=scaling, trunc=trunc)
	return changeordering(O, lattice2, x, trunc=trunc)[2]
end

get_left(lattice::ImagGrassmannLattice{<:A1B1B1A1}, j::Int) = index(lattice, j, conj=true, band=lattice.bands)
get_right(lattice::ImagGrassmannLattice{<:A1B1B1A1}, j::Int) = index(lattice, j, conj=false, band=lattice.bands)

_scaling(x::GrassmannMPS) = scaling(x)

# accsysdynamics_fast(lattice::RealGrassmannLattice, model::AbstractImpurityModel; kwargs...) = error("accsysdynamics_fast only support BranchLocalLayout for RealGrassmannLattice")
function _accsysdynamics_fast(lattice::RealGrassmannLattice{<:A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2}, model::AbstractImpurityModel; scaling::Int=10, trunc::TruncationScheme=DefaultKTruncation, kwargs...)
	lattice1 = similar(lattice, N=1)
	lattice_n = zoomin(lattice1, scaling=scaling)
	lattice_n_1 = similar(lattice_n, N=1)
	gmps_n_1 = sysdynamics(lattice_n_1, model; trunc=trunc, kwargs...)
	# println(bond_dimensions(gmps_n_1))
	posa_f = get_left(lattice_n_1, 2, forward=true)
	posb_f = get_right(lattice_n_1, 1, forward=true)
	posa_b = get_left(lattice_n_1, 1, forward=false)
	posb_b = get_right(lattice_n_1, 2, forward=false)
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
		posa_n = get_left(lattice_n, i+1, forward=true)
		posb_n = get_right(lattice_n, i, forward=true)		
		for (idx, idx_n) in zip(posa_f:posb_f, posa_n:posb_n)
			gmps_n[idx_n] = gmps_n_1[idx] 
		end
		posa_n = get_left(lattice_n, i, forward=false)
		posb_n = get_right(lattice_n, i+1, forward=false)		
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
		posa_n = get_left(lattice, i+1, forward=true)
		posb_n = get_right(lattice, i, forward=true)		
		for (idx, idx_n) in zip(posa_f:posb_f, posa_n:posb_n)
			gmps[idx_n] = gmps1[idx]
		end		
		posa_n = get_left(lattice, i, forward=false)
		posb_n = get_right(lattice, i+1, forward=false)		
		for (idx, idx_n) in zip(posa_b:posb_b, posa_n:posb_n)
			gmps[idx_n] = gmps1[idx]
		end		
	end

	# return canonicalize!(gmps, alg=Orthogonalize(trunc=trunc))
	return gmps
end
function accsysdynamics_fast(lattice::RealGrassmannLattice{O}, model::AbstractImpurityModel; scaling::Int=10, trunc::TruncationScheme=DefaultKTruncation, kwargs...) where O
	lattice2 = similar(lattice, ordering=A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2())
	x = _accsysdynamics_fast(lattice2, model; scaling=scaling, trunc=trunc, kwargs...)
	return changeordering(O, lattice2, x, trunc=trunc)[2]
end

get_left(lattice::RealGrassmannLattice{<:A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2}, j::Int; forward::Bool) = index(lattice, j, conj=true, band=lattice.bands, forward=forward)
get_right(lattice::RealGrassmannLattice{<:A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2}, j::Int; forward::Bool) = index(lattice, j, conj=false, band=lattice.bands, forward=forward)


function _contract_band(x::GrassmannMPS, lattice::ImagGrassmannLattice{<:A1A1B1B1}, j::Int)
	local tmp2
	for band in 1:lattice.bands
		pos1 = index(lattice, j, conj=false, band=band)
		pos2 = index(lattice, j, conj=true, band=band)
		@assert pos1 + 1 == pos2
		@tensor tmp4[1,2,4;5] := x[pos1][1,2,3] * x[pos2][3,4,5]
		if band == 1
			tmp2 = g_trace(tmp4, 2)
		else
			@assert @isdefined tmp2
			tmp2 = tmp2 * g_trace(tmp4, 2)
		end
		tmp2 = rmul!(tmp2, scaling(x)^2)
	end
	return tmp2
end

function _contract_band(x::GrassmannMPS, lattice::ImagGrassmannLattice{<:A1B1B1A1}, j::Int)
	local tmp2
	for band in lattice.bands:-1:1
		pos1 = index(lattice, j, conj=false, band=band)
		pos2 = index(lattice, j, conj=true, band=band)
		if band == lattice.bands
			@assert pos1 + 1 == pos2
			@tensor tmp4[1,2,4;5] := x[pos1][1,2,3] * x[pos2][3,4,5]
			tmp2 = g_trace(tmp4, 2)
		else
			@assert @isdefined tmp2
			@tensor tmp4[1,2,5;6] := x[pos1][1,2,3] * tmp2[3,4] * x[pos2][4,5,6]
			tmp2 = g_trace(tmp4, 2)
		end
		tmp2 = rmul!(tmp2, scaling(x)^2)
	end
	return tmp2	
end

function _contract_band(x::GrassmannMPS, lattice::RealGrassmannLattice{<:A1A1a1a1B1B1b1b1}, j::Int)
	local tmp2
	for band in 1:lattice.bands
		if j == 0
			pos1 = index(lattice, 0, conj=false, band=band)
			pos2 = index(lattice, 0, conj=true, band=band)
			@tensor tmp4[1,2,4;5] := x[pos1][1,2,3] * x[pos2][3,4,5]
			@assert pos1 + 1 == pos2
			if band == 1
				tmp2 = g_trace(tmp4, 2)
			else
				@assert @isdefined tmp2
				tmp2 = tmp2 * g_trace(tmp4, 2)
			end
			tmp2 = rmul!(tmp2, scaling(x)^2)
		else
			for f in (true, false)
				pos1 = index(lattice, j, conj=false, forward=f, band=band)
				pos2 = index(lattice, j, conj=true, forward=f, band=band)
				@assert pos1 + 1 == pos2
				@tensor tmp4[1,2,4;5] := x[pos1][1,2,3] * x[pos2][3,4,5]
				if (band == 1) && (f == true)
					tmp2 = g_trace(tmp4, 2)
				else
					@assert @isdefined tmp2
					tmp2 = tmp2 * g_trace(tmp4, 2)
				end	
				tmp2 = rmul!(tmp2, scaling(x)^2)		
			end
		end
	end
	return tmp2
end

function _contract_band(x::GrassmannMPS, lattice::RealGrassmannLattice{<:A1a1B1b1b1B1a1A1}, j::Int)
	local tmp2
	for band in lattice.bands:-1:1
		if j == 0
			pos1 = index(lattice, 0, conj=false, band=band)
			pos2 = index(lattice, 0, conj=true, band=band)
			if band == lattice.bands
				@assert pos1 + 1 == pos2
				@tensor tmp4[1,2,4;5] := x[pos1][1,2,3] * x[pos2][3,4,5]
				tmp2 = g_trace(tmp4, 2)
			else
				@assert @isdefined tmp2
				@tensor tmp4[1,2,5;6] := x[pos1][1,2,3] * tmp2[3,4] * x[pos2][4,5,6]
				tmp2 = g_trace(tmp4, 2)
			end
			tmp2 = rmul!(tmp2, scaling(x)^2)
		else
			for f in (false, true)
				pos1 = index(lattice, j, conj=false, forward=f, band=band)
				pos2 = index(lattice, j, conj=true, forward=f, band=band)
				# println("pos1=$(pos1), pos2=$(pos2)")
				if (band == lattice.bands) && (f == false)
					@assert pos1 + 1 == pos2
					@tensor tmp4[1,2,4;5] := x[pos1][1,2,3] * x[pos2][3,4,5]
					tmp2 = g_trace(tmp4, 2)
				else
					@assert @isdefined tmp2
					@tensor tmp4[1,2,5;6] := x[pos1][1,2,3] * tmp2[3,4] * x[pos2][4,5,6]
					tmp2 = g_trace(tmp4, 2)
				end	
				tmp2 = rmul!(tmp2, scaling(x)^2)		
			end
		end
	end
	return tmp2
end



function _contract_band(x::GrassmannMPS, lattice::RealGrassmannLattice{<:A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2}, j::Int; forward::Bool)
	local tmp2
	for band in 1:lattice.bands
		if j == 0
			pos1 = index(lattice, 0, conj=false, band=band)
			pos2 = index(lattice, 0, conj=true, band=band)
			@tensor tmp4[1,2,4;5] := x[pos1][1,2,3] * x[pos2][3,4,5]
			@assert pos1 + 1 == pos2
			if band == 1
				tmp2 = g_trace(tmp4, 2)
			else
				@assert @isdefined tmp2
				tmp2 = tmp2 * g_trace(tmp4, 2)
			end
			tmp2 = rmul!(tmp2, scaling(x)^2)
		else
			pos1 = index(lattice, j, conj=false, forward=forward, band=band)
			pos2 = index(lattice, j, conj=true, forward=forward, band=band)
			@assert pos1 + 1 == pos2
			@tensor tmp4[1,2,4;5] := x[pos1][1,2,3] * x[pos2][3,4,5]
			if band == 1
				tmp2 = g_trace(tmp4, 2)
			else
				@assert @isdefined tmp2
				tmp2 = tmp2 * g_trace(tmp4, 2)
			end	
			tmp2 = rmul!(tmp2, scaling(x)^2)		
		end
	end
	return tmp2
end

function _contract_band(x::GrassmannMPS, lattice::RealGrassmannLattice{<:A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2}, j::Int; forward::Bool)
	local tmp2
	for band in lattice.bands:-1:1
		if j == 0
			pos1 = index(lattice, 0, conj=false, band=band)
			pos2 = index(lattice, 0, conj=true, band=band)
			if band == lattice.bands
				@assert pos1 + 1 == pos2
				@tensor tmp4[1,2,4;5] := x[pos1][1,2,3] * x[pos2][3,4,5]
				tmp2 = g_trace(tmp4, 2)
			else
				@assert @isdefined tmp2
				@tensor tmp4[1,2,5;6] := x[pos1][1,2,3] * tmp2[3,4] * x[pos2][4,5,6]
				tmp2 = g_trace(tmp4, 2)
			end
			tmp2 = rmul!(tmp2, scaling(x)^2)
		else
			pos1 = index(lattice, j, conj=false, forward=forward, band=band)
			pos2 = index(lattice, j, conj=true, forward=forward, band=band)
			# println("pos1=$(pos1), pos2=$(pos2)")
			if (band == lattice.bands) 
				@assert pos1 + 1 == pos2
				@tensor tmp4[1,2,4;5] := x[pos1][1,2,3] * x[pos2][3,4,5]
				tmp2 = g_trace(tmp4, 2)
			else
				@assert @isdefined tmp2
				@tensor tmp4[1,2,5;6] := x[pos1][1,2,3] * tmp2[3,4] * x[pos2][4,5,6]
				tmp2 = g_trace(tmp4, 2)
			end	
			tmp2 = rmul!(tmp2, scaling(x)^2)		
		end
	end
	return tmp2
end

