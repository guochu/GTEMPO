"""
	FockMPS{T<:Number, R<:Real} 

The 1-th position is state |0⟩
The 2-th position is state |1⟩
"""
struct FockMPS{T<:Number, R<:Real} <: Dense1DTN{T}
	data::Vector{Array{T, 3}}
	s::Vector{Union{Missing, Vector{R}}}
	scaling::Ref{Float64}

function FockMPS{T, R}(data::AbstractVector, svectors::Vector, scaling::Ref{R}) where {T<:Number, R<:Number}
	(R == real(T)) || throw(ArgumentError("scalar type for singular vectors must be real"))
	(length(data)+1 == length(svectors)) || throw(DimensionMismatch("length of singular vectors must be length of site tensors+1"))
	_check_mps_space(data)
	new{T, R}(convert(Vector{Array{T, 3}}, data), convert(Vector{Union{Missing, Vector{R}}}, svectors), scaling)
end
end

function FockMPS{T, R}(data::Vector, scaling::Ref{R}) where {T<:Number, R<:Number}
	(R == real(T)) || throw(ArgumentError("scalar type for singular vectors must be real"))
	_check_mps_space(data)
	svectors = Vector{Union{Missing, Vector{R}}}(undef, length(data)+1)
	svectors[1] = ones(space_l(data[1]))
	svectors[end] = ones(space_r(data[end]))
	return FockMPS{T, R}(convert(Vector{Array{T, 3}}, data), svectors, scaling)
end

function FockMPS(data::AbstractVector{<:DenseMPSTensor{T}}, svectors::AbstractVector; scaling::Real=1) where {T <: Number}
	R = real(T)
	return FockMPS{T, R}(data, svectors, Ref(convert(R, scaling)))
end 
function FockMPS(data::AbstractVector{<:DenseMPSTensor{T}}; scaling::Real=1) where {T <: Number}
	R = real(T)
	return FockMPS{T, R}(data, Ref(convert(R, scaling)))
end

# function FockMPS(::Type{T}, L::Int) where {T <: Number}
# 	v = zeros(T, 1, 2, 1)
# 	v[1,1,1] = 1
# 	data = [copy(v) for i in 1:L]
# 	return FockMPS(data, scaling=1)
# end
function FockMPS(::Type{T}, L::Int) where {T <: Number}
	v = ones(T, 1, 2, 1)
	data = [copy(v) for i in 1:L]
	return FockMPS(data, scaling=1)
end
FockMPS(L::Int) = FockMPS(Float64, L)


scaling(x::FockMPS) = x.scaling[]
setscaling!(x::FockMPS, scaling::Real) = (x.scaling[] = scaling)

function TK.normalize!(x::FockMPS)
	setscaling!(x, 1)
	return x
end

Base.copy(psi::FockMPS) = FockMPS(copy(psi.data), copy(psi.s), scaling=scaling(psi))

DMRG.svectors_uninitialized(psi::FockMPS) = any(ismissing, psi.s)
function DMRG.unset_svectors!(psi::FockMPS)
	psi.s[2:end-1] .= missing
	return psi
end

function _check_mps_space(mpstensors::Vector)
	L = length(mpstensors)
	for i in 1:L-1
		(space_r(mpstensors[i]) == space_l(mpstensors[i+1])) || throw(DimensionMismatch())
	end
	(space_l(mpstensors[1]) == 1) || throw(DimensionMismatch("left boundary should be size 1"))
	(space_r(mpstensors[L]) == 1) || throw(DimensionMismatch("right boundary should be size 1"))
	all(x->size(x, 2)==2, mpstensors) || throw(ArgumentError("physical dimension of site tensors must be 2"))
	return true
end


# initializers
function randomfockmps(::Type{T}, L::Int; D::Int) where {T <: Number}
	mpstensors = Vector{Array{T, 3}}(undef, L)
	mpstensors[1] = randn(T, 1,2,D)
	mpstensors[end] = randn(T, D, 2, 1)
	for i in 2:L-1
		mpstensors[i] = randn(T, D, 2, D)
	end
	return FockMPS(mpstensors)
end
randomfockmps(L::Int; kwargs...) = randomfockmps(Float64, L; kwargs...)


function increase_bond!(psi::FockMPS; D::Int)
	if bond_dimension(psi) < D
		for i in 1:length(psi)
			sl = max(D, size(psi[i], 1))
			sr = max(D, size(psi[i], 3))
			m = zeros(scalartype(psi), sl, size(psi[i], 2), sr)
			m[1:size(psi[i], 1), :, 1:size(psi[i], 3)] .= psi[i]
			psi[i] = m
		end
	end
	return psi
end


# check is canonical
DMRG.isleftcanonical(a::FockMPS; kwargs...) = all(x->isleftcanonical(x; kwargs...), a.data)
DMRG.isrightcanonical(a::FockMPS; kwargs...) = all(x->isrightcanonical(x; kwargs...), a.data)

"""
	iscanonical(psi::MPS; kwargs...) = is_right_canonical(psi; kwargs...)
check if the state is right-canonical, the singular vectors are also checked that whether there are the correct Schmidt numbers or not
This form is useful for time evolution for stability issue and also efficient for computing observers of unitary systems
"""
function DMRG.iscanonical(psi::FockMPS; kwargs...)
	isrightcanonical(psi) || return false
	# we also check whether the singular vectors are the correct Schmidt numbers
	svectors_uninitialized(psi) && return false
	hold = l_LL(psi)
	for i in 1:length(psi)-1
		hold = updateleft(hold, psi[i], psi[i])
		tmp = psi.s[i+1]
		isapprox(hold, Diagonal(tmp.^2); kwargs...) || return false
	end
	return true
end
