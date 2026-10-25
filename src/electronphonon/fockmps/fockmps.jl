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
"""
	randomfockmps(L::Int; D) -> FockMPS
	randomfockmps(::Type{T}, L::Int; D) -> FockMPS

Random FockMPS of length `L` (physical dimension 2) and bond dimension `D`,
delegating to FiniteMPSAlgorithms' `randommps`（canonical form，Schmidt 值已
初始化；返回态已归一化）。
"""
randomfockmps(::Type{T}, L::Int; D::Int) where {T <: Number} = FockMPS(FMA.randommps(T, L; d=2, D=D))
randomfockmps(L::Int; kwargs...) = randomfockmps(Float64, L; kwargs...)


"""
	changebond!(psi::FockMPS, D::Int; noise::Real=0) -> psi

Bring the bond profile of `psi` to `min(D, feasible)`: bonds larger than the target are
shrunk by slicing the leading bond indices, smaller bonds are grown with random entries
of magnitude `noise` (zero by default, so that the represented state is unchanged).
Delegates to FiniteMPSAlgorithms' `changebond!` on the inner `CanonicalMPS` payload;
the chain is returned in right-canonical form with the Schmidt values initialized.
"""
changebond!(psi::FockMPS, D::Int; noise::Real=0) = (FMA.changebond!(psi.parent; D, noise); psi)


# check is canonical
isleftcanonical(a::FockMPS; kwargs...) = all(x->isleftcanonical(x; kwargs...), a.data)
isrightcanonical(a::FockMPS; kwargs...) = all(x->isrightcanonical(x; kwargs...), a.data)

"""
	iscanonical(psi::FockMPS; atol=1e-8)

Whether `psi` is in canonical form（全部位点张量右正交、Schmidt 值已初始化、
且逐键 `Diagonal(s.^2) ≈` 左环境），委托给 FiniteMPSAlgorithms 的
`iscanonical`。
"""
iscanonical(psi::FockMPS; kwargs...) = FMA.iscanonical(psi.parent; kwargs...)
