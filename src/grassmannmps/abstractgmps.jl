abstract type AbstractGMPS{A<:MPSTensor} end
abstract type AbstractFiniteGMPS{A<:MPSTensor} <: AbstractGMPS{A} end

TK.scalartype(::Type{<:AbstractGMPS{A}}) where {A<:MPSTensor} = scalartype(A)
TK.spacetype(::Type{<:AbstractGMPS{A}}) where {A<:MPSTensor} = spacetype(A)
TK.spacetype(m::AbstractGMPS) = spacetype(typeof(m))
TK.sectortype(A::Type{<:AbstractGMPS}) = sectortype(spacetype(A))
TK.sectortype(a::AbstractGMPS) = sectortype(typeof(a))
DMRG.mpstensortype(::Type{<:AbstractGMPS{A}}) where {A<:MPSTensor} = A
DMRG.mpstensortype(m::AbstractGMPS) = mpstensortype(typeof(m))


DMRG.space_l(a::AbstractGMPS) = space_l(a[1])
DMRG.space_r(a::AbstractGMPS) = space_r(a[end])

DMRG.bond_dimension(a::AbstractGMPS, bond::Int) = begin
	((bond >= 1) && (bond <= length(a))) || throw(BoundsError(storage(a), bond))
	dim(space(a[bond], 3))
end 
DMRG.bond_dimensions(a::AbstractGMPS) = [bond_dimension(a, i) for i in 1:length(a)]
DMRG.bond_dimension(a::AbstractGMPS) = maximum(bond_dimensions(a))

DMRG.physical_space(a::AbstractGMPS, i::Int) = physical_space(a[i])
DMRG.physical_spaces(a::AbstractGMPS) = [physical_space(a[i]) for i in 1:length(a)]
DMRG.left_virtualspace(a::AbstractGMPS, i::Int) = space_l(a[i])
DMRG.right_virtualspace(a::AbstractGMPS, i::Int) = space_r(a[i])
DMRG.left_virtualspaces(a::AbstractGMPS) = [left_virtualspace(a, i) for i in 1:length(a)]
DMRG.right_virtualspaces(a::AbstractGMPS) = [right_virtualspace(a, i) for i in 1:length(a)]


edgespace(::Type{<:AbstractFiniteGMPS}) = Rep[ℤ₂](0=>1)
edgespace(x::AbstractGMPS) = edgespace(typeof(x))

function DMRG.l_LL(f, vspace::ElementarySpace, x::AbstractGMPS, ys::AbstractGMPS...)
	T = promote_type(scalartype(x), map(scalartype, ys)...)
	return GrassmannTensorMap(TensorMap(f, T, vspace, ⊗(space_l(x), map(y->space_l(y), ys)...)))
end
DMRG.l_LL(x::AbstractGMPS, ys::AbstractGMPS...) = l_LL(ones, edgespace(x), x, ys...)

function DMRG.r_RR(f, vspace::ElementarySpace, x::AbstractGMPS, ys::AbstractGMPS...)
	T = promote_type(scalartype(x), map(scalartype, ys)...)
	return GrassmannTensorMap(TensorMap(f, T, ⊗(space_r(x)', map(y->space_r(y)', reverse(ys))...), vspace))
end
DMRG.r_RR(x::AbstractGMPS, ys::AbstractGMPS...) = r_RR(ones, edgespace(x), x, ys...)
