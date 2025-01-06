"""
	integrateband(lattice::AbstractGrassmannLattice, x::GrassmannMPS; band::Int=1)
"""
function integrateband(lattice::AbstractGrassmannLattice, x::GrassmannMPS; band::Int=1) 
	data = _integrateband(lattice, x, band=band)
	return GrassmannMPS(data, scaling=scaling(x)^(length(x)/length(data)))
end


function _integrateband(lattice::AbstractGrassmannLattice, x; band::Int=1)
	(ConjugationStyle(lattice) isa AdjacentConjugation) || throw(ArgumentError("integrateband only supports AdjacentConjugation style"))
	(1 <= band <= lattice.bands) || throw(BoundsError(1:lattice.bands, band))
	(length(x) == length(lattice)) || throw(DimensionMismatch())
	(lattice.bands == 1) && return integrate(lattice, x)
	lattice2 = similar(lattice, bands=lattice.bands-1)

	r2 = indexmappings(lattice2)
	r1 = indexmappings(lattice)
	mm = Dict(r1[(j, c, b, ifelse(bj<band, bj, bj+1))]=>pos for ((j, c, b, bj), pos) in r2)
	data = similar(x.data, length(lattice2))
	i = 1
	tmp = nothing
	while isnothing(get(mm, i, nothing))
		tmp2 = _integrate_pairsites(x[i], x[i+1])
		tmp = _absorb_transfer(tmp, tmp2)
		i += 2
	end
	pos2 = 0
	while i <= length(lattice)
		pos2 = get(mm, i, nothing)
		if isnothing(pos2)
			tmp2 = _integrate_pairsites(x[i], x[i+1])
			tmp = _absorb_transfer(tmp, tmp2)
			i += 2
		else
			data[pos2] = _absorb_sitetensor(tmp, x[i])
			i += 1
			tmp = nothing
		end		
	end
	if isnothing(tmp)
		(pos2 == length(lattice2)) || error("something wrong")
	else
		@tensor r[1,2;4] := data[end][1,2,3] * tmp[3,4]
		data[end] = r 
	end
	return data
end


function _integrate_pairsites(a, b)
	@tensor tmp[1,2,4;5] := a[1,2,3] * b[3,4,5]
	return _g_trace(tmp, 2)
end

_absorb_transfer(tmp, tmp2) = isnothing(tmp) ? tmp2 : tmp * tmp2
function _absorb_sitetensor(tmp, x)
	isnothing(tmp) && return x
	@tensor tmp2[1 3; 4] := tmp[1,2] * x[2,3,4]
end