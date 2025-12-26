"""
	SchurMPOTensor{S<:ElementarySpace, M<:MPOTensor, T<:Number} 

Upper triangular matrix, with each entry an MPOTensor
The top left and bottom right must be identity
"""
struct SchurMPOTensor{S<:ElementarySpace, M<:MPOTensor, T<:Number} <: AbstractSparseMPOTensor{S}
	Os::Array{Union{M, T}, 2}
	leftspaces::Vector{S}
	rightspaces::Vector{S}
	pspace::S
end

function SchurMPOTensor{S, M, T}(data::AbstractMatrix) where {S<:ElementarySpace, M <:MPOTensor, T<:Number}
	(size(data, 1) == size(data, 2)) || throw(ArgumentError("SchurMPOTensor requires a square matrix"))
	Os, leftspaces, rightspaces, pspace = compute_mpotensor_data(S, M, T, data)
	for i in 1:size(Os, 1)
		for j in 1:i-1
			(Os[i, j] == zero(T)) || throw(ArgumentError("SchurMPOTensor should be upper triangular"))
		end
	end
	return SchurMPOTensor{S, M, T}(Os, leftspaces, rightspaces, pspace)
end
function SchurMPOTensor(data::AbstractMatrix, ::Type{T}, pspace::S) where {T <: Number, S <: ElementarySpace}
	M = mpotensortype(S, T)
	r = SchurMPOTensor{S, M, T}(data)
	(r.pspace == pspace) || throw(SpaceMismatch("physical space mismatch"))
	return r
end

TK.scalartype(::Type{SchurMPOTensor{S,M,T}}) where {S,M,T} = T
Base.copy(x::SchurMPOTensor) = SchurMPOTensor(copy(x.Os), copy(x.leftspaces), copy(x.rightspaces), x.pspace)

# upper triangular form
# the middle diagonal terms may be identity operator or the JW operator,
"""
	SchurMPOTensor(data::Array{Any, 2}) 
"""
function SchurMPOTensor(data::AbstractMatrix{Any}) 
	S, T = compute_generic_mpotensor_types(data)
	M = mpotensortype(S, T)
	return SchurMPOTensor{S, M, T}(data)
end
SchurMPOTensor(data::AbstractMatrix{SparseMPOTensorElement{M, T}}) where {M<:MPOTensor, T<:Number} = SchurMPOTensor{spacetype(M), M, T}(data)


Base.contains(m::SchurMPOTensor{S, M, T}, i::Int, j::Int) where {S, M, T} = (i <= j) && (m.Os[i, j] != zero(T))
function isscal(x::SchurMPOTensor{S,M,T}, i::Int, j::Int) where {S,M,T}
	sj = x.Os[i, j]
	return (sj isa T) && (abs(sj) > 1.0e-14)
end 

Base.convert(::Type{<:SparseMPOTensor}, t::SchurMPOTensor) = SparseMPOTensor(t.Os, t.leftspaces, t.rightspaces, t.pspace)


function Base.setindex!(m::SchurMPOTensor{S, M, T}, v, i::Int, j::Int) where {S, M, T}
	(i > j) && throw(ArgumentError("not allowed to set the low triangular portion."))
	if isa(v, Number)
		m.Os[i, j] = convert(T, v)
	elseif isa(v, MPOTensor)
		((space(v, 1) == m.leftspaces[i]) && (space(v, 3)' == m.rightspaces[j])) || throw(SpaceMismatch())
		m.Os[i, j] = convert(M, v) 
	elseif isa(v, MPSBondTensor)
		b_iden = isomorphism(T, m.leftspaces[i], m.rightspaces[j])
		@tensor tmp[-1 -2; -3 -4] := b_iden[-1, -3] * v[-2, -4]
		m.Os[i, j] = tmp
	else
		throw(ArgumentError("input should be scalar, MPOTensor or MPSBondTensor type"))
	end
end

function compute_generic_mpotensor_types(data::AbstractMatrix)
	m, n = size(data)
	T = Float64
	local S 
	for sj in data
		if isa(sj, Number)
			T = promote_type(T, typeof(sj))
		elseif isa(sj, AbstractTensorMap)
			T = promote_type(T, scalartype(sj))
			if !@isdefined S
				S = spacetype(sj)
			else
				(S == spacetype(sj)) || throw(SpaceMismatch())
			end
		else
			throw(ArgumentError("eltype must be scalar or TensorMap"))
		end		
	end
	(!@isdefined S) && throw(ArgumentError("there must be at least one TensorMap type"))
	return S, T	
end
# the same row should have the same left space
# the same column should have the same right space
function compute_mpotensor_data(::Type{S}, ::Type{M}, ::Type{T}, data::AbstractMatrix) where {S<:ElementarySpace, M<:MPOTensor{S}, T<:Number}
	@assert !isempty(data)
	m, n = size(data)
	new_data = Array{Union{M, T}, 2}(undef, m, n) 
	rightspaces = Vector{Union{Missing, S}}(missing, n)
	leftspaces = Vector{Union{Missing, S}}(missing, m)

	# check spaces
	local pspace::S
	for i in 1:m
		for j in 1:n
			sj = data[i, j]
			if isa(sj, MPOTensor)
				s_l, s_r = space_l(sj), space_r(sj)'
				if !ismissing(leftspaces[i])
					(leftspaces[i] == s_l) || throw(SpaceMismatch("$i-th left space mismatch: $((leftspaces[i], s_l))"))
				else
					leftspaces[i] = s_l
				end
 				if !ismissing(rightspaces[j])
 					(rightspaces[j] == s_r) || throw(SpaceMismatch("$j-th right space mismatch: $((rightspaces[j], s_r))"))
 				else
 					rightspaces[j] = s_r
 				end		
 				(space(sj, 2) == space(sj, 4)') || throw(SpaceMismatch("physical space mismatch"))
				s_p = space(sj, 2)
 				if !@isdefined pspace
 					pspace = s_p
 				else
 					(pspace == s_p) || throw(SpaceMismatch("physical space mismatch"))
 				end
 			elseif isa(sj, MPSBondTensor)
 				(domain(sj) == codomain(sj)) || throw(SpaceMismatch("physical space mismatch"))			 						
 				s_p = space(sj, 1)
 				if !@isdefined pspace
 					pspace = s_p
 				else
 					(pspace == s_p) || throw(SpaceMismatch("physical space mismatch"))
 				end 				
			end
		end
	end
	(@isdefined pspace) || throw(ArgumentError("pspace is missing."))
	for i in 1:m
		# ismissing(leftspaces[i]) && throw(ArgumentError("leftspaces $i not assigned"))
		if ismissing(leftspaces[i])
			leftspaces[i] = oneunit(pspace)
		end
	end
	for j in 1:n
		# ismissing(rightspaces[j]) && throw(ArgumentError("rightspaces $j not assigned"))
		if ismissing(rightspaces[j])
			rightspaces[j] = oneunit(pspace)			
		end
	end		
	for i in 1:m
		for j in 1:n
			sj = data[i, j]
			if isa(sj, MPSBondTensor)
				_is_id, scal = isid(sj)
				if _is_id
					sj = scal
				else
					virtual = isomorphism(T, leftspaces[i], rightspaces[j])
					sj = @tensor tmp[1,3;2,4] := virtual[1, 2] * sj[3,4]
				end
 			elseif isa(sj, MPOTensor)
 				_is_id, scal = isid(sj)
 				if _is_id
 					sj = scal
 				end
 			end
 			if isa(sj, MPOTensor)
 				sj = convert(M, sj)
 			else
 				isa(sj, Number) || throw(ArgumentError("elt should either be a tensor or a scalar"))
 				sj = convert(T, sj)
 			end
 			new_data[i, j] = sj
		end
	end
	return new_data, convert(Vector{S}, leftspaces), convert(Vector{S}, rightspaces), pspace
end

function compute_mpotensor_data(::Type{M}, ::Type{T}, data::AbstractMatrix, leftspaces::Vector{S}, rightspaces::Vector{S}, pspace::S) where {S<:ElementarySpace, M<:MPOTensor{S}, T<:Number}
	@assert !isempty(data)
	m, n = size(data)
	@assert (m == length(leftspaces) ) && (n == length(rightspaces))
	# T = compute_scalartype(data)
	# M = tensormaptype(S, 2, 2, T)
	new_data = Matrix{Union{M, T}}(undef, m, n) 
	for i in 1:m
		for j in 1:n
			sj = data[i, j]
			if isa(sj, MPSBondTensor)
				(space(sj, 1) == space(sj, 2)' == pspace) || throw(SpaceMismatch("physical space mismatch"))
				# isoneunit(leftspaces[i]) || throw(SpaceMismatch("$i-th left space mismatch"))
				# isoneunit(rightspaces[j]) || throw(SpaceMismatch("$j-th right space mismatch"))
				_is_id, scal = isid(sj)
				if _is_id
					sj = scal
				else
					virtual = isomorphism(T, leftspaces[i], rightspaces[j])
					# sj = _add_legs(sj)
					sj = @tensor tmp[1,3;2,4] := virtual[1, 2] * sj[3,4]
				end
 			elseif isa(sj, MPOTensor)
 				# sj = convert(M, sj)
 				s_l = space(sj, 1)
 				s_r = space(sj, 3)'
 				s_p = space(sj, 2)
 				(s_p == space(sj, 4)' == pspace) || throw(SpaceMismatch("physical space mismatch"))
 				(leftspaces[i] == s_l) || throw(SpaceMismatch("$i-th left space mismatch: $((leftspaces[i], s_l))"))
 				(rightspaces[j] == s_r) || throw(SpaceMismatch("$j-th right space mismatch: $((rightspaces[j], s_r))"))
 				_is_id, scal = isid(sj)
 				if _is_id
 					sj = scal
 				end
 			end
 			if isa(sj, MPOTensor)
 				sj = convert(M, sj)
 			else
 				isa(sj, Number) || throw(ArgumentError("elt should either be a tensor or a scalar"))
 				sj = convert(T, sj)
 			end
 			new_data[i, j] = sj
		end
	end
	return new_data
end
