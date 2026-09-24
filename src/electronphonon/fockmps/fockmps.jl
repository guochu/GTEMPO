"""
	FockMPS{T<:Number, R<:Real}

One-dimensional dense tensor network (`Dense1DTN`) representing a finite matrix product state on the Fock lattice.
The storage payload is a FiniteMPSAlgorithms `CanonicalMPS` (field `.parent`): site tensors, Schmidt values
and the per-site scaling are all carried by the payload, on which the FiniteMPSAlgorithms algorithms operate
in place. `.data` delegates to the payload's site-tensor vector; `.s` / `.scaling` delegate to the payload's
Schmidt values / scaling.

Site tensor convention (payload `CanonicalMPS`):
- Dimension 1: left auxiliary (bond) index, with dimension 1 at the leftmost site
- Dimension 2: physical index; entry 1 corresponds to state |0⟩ and entry 2 to state |1⟩
- Dimension 3: right auxiliary (bond) index, with dimension 1 at the rightmost site

# Examples
```julia
julia> psi = FockMPS(5)               # all-ones vacuum state with 5 sites
julia> psi = randomfockmps(5, D=8)    # random FockMPS with 5 sites and bond dimension 8
```
"""
struct FockMPS{T<:Number, R<:Real} <: Dense1DTN{T}
	parent::CanonicalMPS{T, R}
end

# `.parent` 即内层的 CanonicalMPS payload；`.data` 委托到 payload 的站点张量
# 向量；`.s` / `.scaling` 委托到 payload
function Base.getproperty(psi::FockMPS, s::Symbol)
	s === :parent && return getfield(psi, :parent)
	s === :data && return getfield(psi, :parent).data
	s === :s && return getfield(psi, :parent).s
	s === :scaling && return getfield(psi, :parent).scaling
	throw(ArgumentError("FockMPS has no property $s"))
end
Base.propertynames(::FockMPS) = (:parent, :data, :s, :scaling)

function FockMPS(data::AbstractVector{<:DenseMPSTensor{T}}, svectors::AbstractVector; scaling::Real=1) where {T<:Number}
	return FockMPS(CanonicalMPS(convert(Vector{Array{T, 3}}, data), svectors, scaling))
end
function FockMPS(data::AbstractVector{<:DenseMPSTensor{T}}; scaling::Real=1) where {T<:Number}
	return FockMPS(CanonicalMPS(convert(Vector{Array{T, 3}}, data); scaling))
end

function FockMPS(::Type{T}, L::Int) where {T <: Number}
	return FockMPS(CanonicalMPS(T, L; d=2))
end
FockMPS(L::Int) = FockMPS(Float64, L)


scaling(x::FockMPS) = x.scaling[]
setscaling!(x::FockMPS, scaling::Real) = (x.scaling[] = scaling)

function TK.normalize!(x::FockMPS)
	setscaling!(x, 1)
	return x
end

Base.copy(psi::FockMPS) = FockMPS(copy(psi.parent))
Base.copy!(dst::FockMPS, src::FockMPS) = (copy!(dst.parent, src.parent); dst)

svectors_uninitialized(psi::FockMPS) = FMA.svectors_uninitialized(psi.parent)
function unset_svectors!(psi::FockMPS)
	FMA.unset_svectors!(psi.parent)
	return psi
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
isleftcanonical(a::FockMPS; kwargs...) = all(x->isleftcanonical(x; kwargs...), a.data)
isrightcanonical(a::FockMPS; kwargs...) = all(x->isrightcanonical(x; kwargs...), a.data)

"""
	iscanonical(psi::FockMPS; kwargs...) = is_right_canonical(psi; kwargs...)
check if the state is right-canonical, the singular vectors are also checked that whether there are the correct Schmidt numbers or not
This form is useful for time evolution for stability issue and also efficient for computing observers of unitary systems
"""
function iscanonical(psi::FockMPS; kwargs...)
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
