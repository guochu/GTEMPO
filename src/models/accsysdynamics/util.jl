_scaling(x::GrassmannMPS) = scaling(x)

function _contract_band(x::GrassmannMPS, lattice::ImagGrassmannLattice{<:A1A1B1B1}, j::Int)
	local tmp2
	for band in 1:lattice.bands
		pos1 = index(lattice, j, conj=false, band=band)
		pos2 = index(lattice, j, conj=true, band=band)
		@assert pos1 + 1 == pos2
		@tensor tmp4[1,2,4;5] := GrassmannTensorMap(x[pos1])[1,2,3] * GrassmannTensorMap(x[pos2])[3,4,5]
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
			@tensor tmp4[1,2,4;5] := GrassmannTensorMap(x[pos1])[1,2,3] * GrassmannTensorMap(x[pos2])[3,4,5]
			tmp2 = g_trace(tmp4, 2)
		else
			@assert @isdefined tmp2
			@tensor tmp4[1,2,5;6] := GrassmannTensorMap(x[pos1])[1,2,3] * tmp2[3,4] * GrassmannTensorMap(x[pos2])[4,5,6]
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
			@tensor tmp4[1,2,4;5] := GrassmannTensorMap(x[pos1])[1,2,3] * GrassmannTensorMap(x[pos2])[3,4,5]
			@assert pos1 + 1 == pos2
			if band == 1
				tmp2 = g_trace(tmp4, 2)
			else
				@assert @isdefined tmp2
				tmp2 = tmp2 * g_trace(tmp4, 2)
			end
			tmp2 = rmul!(tmp2, scaling(x)^2)
		else
			for f in (:+, :-)
				pos1 = index(lattice, j, conj=false, branch=f, band=band)
				pos2 = index(lattice, j, conj=true, branch=f, band=band)
				@assert pos1 + 1 == pos2
				@tensor tmp4[1,2,4;5] := GrassmannTensorMap(x[pos1])[1,2,3] * GrassmannTensorMap(x[pos2])[3,4,5]
				if (band == 1) && (f == :+)
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

function _contract_band(x::GrassmannMPS, lattice::RealGrassmannLattice{<:A1B1ā1b̄1A1B1a1b1}, j::Int)
	local tmp2
	if j == 0
		for band in 1:lattice.bands
			pos1 = index(lattice, 0, conj=false, band=band)
			pos2 = index(lattice, 0, conj=true, band=band)
			@tensor tmp4[1,2,4;5] := GrassmannTensorMap(x[pos1])[1,2,3] * GrassmannTensorMap(x[pos2])[3,4,5]
			@assert pos1 + 1 == pos2
			if band == 1
				tmp2 = g_trace(tmp4, 2)
			else
				@assert @isdefined tmp2
				tmp2 = tmp2 * g_trace(tmp4, 2)
			end
			tmp2 = rmul!(tmp2, scaling(x)^2)
		end
	else
		y = copy(x)
		pos_f, pos_e = band_boundary(lattice, j)
		posa, posb = pos_f, pos_f+2*lattice.bands
		while posb < pos_f + 3*lattice.bands
			@assert posa < posb
			for i in 1:posb-posa-1
				naive_swap!(y, posb-i, trunc=DefaultIntegrationTruncation)
			end
			pos1, pos2 = posa, posa + 1
			@tensor tmp4[1,2,4;5] := GrassmannTensorMap(y[pos1])[1,2,3] * GrassmannTensorMap(y[pos2])[3,4,5]
			if @isdefined tmp2
				tmp2 = tmp2 * g_trace(tmp4, 2)
			else
				tmp2 = g_trace(tmp4, 2)
			end
			tmp2 = rmul!(tmp2, scaling(x)^2)	
			posa += 2
			posb += 1
		end
		posa, posb = pos_f+2*lattice.bands, pos_f+3*lattice.bands
		while posb <= pos_e
			@assert posa < posb
			for i in 1:posb-posa
				naive_swap!(y, posb-i, trunc=DefaultIntegrationTruncation)
			end
			pos1, pos2 = posa, posa + 1
			@tensor tmp4[1,2,4;5] := GrassmannTensorMap(y[pos1])[1,2,3] * GrassmannTensorMap(y[pos2])[3,4,5]
			tmp2 = tmp2 * g_trace(tmp4, 2)
			tmp2 = rmul!(tmp2, scaling(x)^2)	
			posa += 2
			posb += 1
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
				@tensor tmp4[1,2,4;5] := GrassmannTensorMap(x[pos1])[1,2,3] * GrassmannTensorMap(x[pos2])[3,4,5]
				tmp2 = g_trace(tmp4, 2)
			else
				@assert @isdefined tmp2
				@tensor tmp4[1,2,5;6] := GrassmannTensorMap(x[pos1])[1,2,3] * tmp2[3,4] * GrassmannTensorMap(x[pos2])[4,5,6]
				tmp2 = g_trace(tmp4, 2)
			end
			tmp2 = rmul!(tmp2, scaling(x)^2)
		else
			for f in (:-, :+)
				pos1 = index(lattice, j, conj=false, branch=f, band=band)
				pos2 = index(lattice, j, conj=true, branch=f, band=band)
				# println("pos1=$(pos1), pos2=$(pos2)")
				if (band == lattice.bands) && (f == :-)
					@assert pos1 + 1 == pos2
					@tensor tmp4[1,2,4;5] := GrassmannTensorMap(x[pos1])[1,2,3] * GrassmannTensorMap(x[pos2])[3,4,5]
					tmp2 = g_trace(tmp4, 2)
				else
					@assert @isdefined tmp2
					@tensor tmp4[1,2,5;6] := GrassmannTensorMap(x[pos1])[1,2,3] * tmp2[3,4] * GrassmannTensorMap(x[pos2])[4,5,6]
					tmp2 = g_trace(tmp4, 2)
				end	
				tmp2 = rmul!(tmp2, scaling(x)^2)		
			end
		end
	end
	return tmp2
end



function _contract_band(x::GrassmannMPS, lattice::RealGrassmannLattice{<:A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2}, j::Int; branch::Symbol)
	local tmp2
	for band in 1:lattice.bands
		if j == 0
			pos1 = index(lattice, 0, conj=false, band=band)
			pos2 = index(lattice, 0, conj=true, band=band)
			@tensor tmp4[1,2,4;5] := GrassmannTensorMap(x[pos1])[1,2,3] * GrassmannTensorMap(x[pos2])[3,4,5]
			@assert pos1 + 1 == pos2
			if band == 1
				tmp2 = g_trace(tmp4, 2)
			else
				@assert @isdefined tmp2
				tmp2 = tmp2 * g_trace(tmp4, 2)
			end
			tmp2 = rmul!(tmp2, scaling(x)^2)
		else
			pos1 = index(lattice, j, conj=false, branch=branch, band=band)
			pos2 = index(lattice, j, conj=true, branch=branch, band=band)
			@assert pos1 + 1 == pos2
			@tensor tmp4[1,2,4;5] := GrassmannTensorMap(x[pos1])[1,2,3] * GrassmannTensorMap(x[pos2])[3,4,5]
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

function _contract_band(x::GrassmannMPS, lattice::RealGrassmannLattice{<:A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2}, j::Int; branch::Symbol)
	local tmp2
	for band in lattice.bands:-1:1
		if j == 0
			pos1 = index(lattice, 0, conj=false, band=band)
			pos2 = index(lattice, 0, conj=true, band=band)
			if band == lattice.bands
				@assert pos1 + 1 == pos2
				@tensor tmp4[1,2,4;5] := GrassmannTensorMap(x[pos1])[1,2,3] * GrassmannTensorMap(x[pos2])[3,4,5]
				tmp2 = g_trace(tmp4, 2)
			else
				@assert @isdefined tmp2
				@tensor tmp4[1,2,5;6] := GrassmannTensorMap(x[pos1])[1,2,3] * tmp2[3,4] * GrassmannTensorMap(x[pos2])[4,5,6]
				tmp2 = g_trace(tmp4, 2)
			end
			tmp2 = rmul!(tmp2, scaling(x)^2)
		else
			pos1 = index(lattice, j, conj=false, branch=branch, band=band)
			pos2 = index(lattice, j, conj=true, branch=branch, band=band)
			# println("pos1=$(pos1), pos2=$(pos2)")
			if (band == lattice.bands) 
				@assert pos1 + 1 == pos2
				@tensor tmp4[1,2,4;5] := GrassmannTensorMap(x[pos1])[1,2,3] * GrassmannTensorMap(x[pos2])[3,4,5]
				tmp2 = g_trace(tmp4, 2)
			else
				@assert @isdefined tmp2
				@tensor tmp4[1,2,5;6] := GrassmannTensorMap(x[pos1])[1,2,3] * tmp2[3,4] * GrassmannTensorMap(x[pos2])[4,5,6]
				tmp2 = g_trace(tmp4, 2)
			end	
			tmp2 = rmul!(tmp2, scaling(x)^2)		
		end
	end
	return tmp2
end

