abstract type AbstractGMPS{A<:MPSTensor} end
abstract type AbstractFiniteGMPS{A<:MPSTensor} <: AbstractGMPS{A} end

TK.scalartype(::Type{<:AbstractGMPS{A}}) where {A<:MPSTensor} = scalartype(A)
TK.spacetype(::Type{<:AbstractGMPS{A}}) where {A<:MPSTensor} = spacetype(A)
TK.spacetype(m::AbstractGMPS) = spacetype(typeof(m))
TK.sectortype(A::Type{<:AbstractGMPS}) = sectortype(spacetype(A))
TK.sectortype(a::AbstractGMPS) = sectortype(typeof(a))
mpstensortype(::Type{<:AbstractGMPS{A}}) where {A<:MPSTensor} = A
mpstensortype(m::AbstractGMPS) = mpstensortype(typeof(m))


space_l(a::AbstractGMPS) = space_l(a[1])
space_r(a::AbstractGMPS) = space_r(a[end])

bonddim(a::AbstractGMPS, bond::Int) = begin
	((bond >= 1) && (bond <= length(a))) || throw(BoundsError(storage(a), bond))
	dim(space(a[bond], 3))
end 
bonddims(a::AbstractGMPS) = [bonddim(a, i) for i in 1:length(a)]
bonddim(a::AbstractGMPS) = maximum(bonddims(a))

physpace(a::AbstractGMPS, i::Int) = physpace(a[i])
physpaces(a::AbstractGMPS) = [physpace(a[i]) for i in 1:length(a)]
left_virtualspace(a::AbstractGMPS, i::Int) = space_l(a[i])
right_virtualspace(a::AbstractGMPS, i::Int) = space_r(a[i])
left_virtualspaces(a::AbstractGMPS) = [left_virtualspace(a, i) for i in 1:length(a)]
right_virtualspaces(a::AbstractGMPS) = [right_virtualspace(a, i) for i in 1:length(a)]


edgespace(::Type{<:AbstractFiniteGMPS}) = Z2Space(0=>1)
edgespace(x::AbstractGMPS) = edgespace(typeof(x))

function l_LL(f, vspace::ElementarySpace, x::AbstractGMPS, ys::AbstractGMPS...)
	T = promote_type(scalartype(x), map(scalartype, ys)...)
	return f(T, vspace, ⊗(space_l(x), map(y->space_l(y), ys)...))
end
l_LL(x::AbstractGMPS, ys::AbstractGMPS...) = l_LL(ones, edgespace(x), x, ys...)

function r_RR(f, vspace::ElementarySpace, x::AbstractGMPS, ys::AbstractGMPS...)
	T = promote_type(scalartype(x), map(scalartype, ys)...)
	return f(T, ⊗(space_r(x)', map(y->space_r(y)', reverse(ys))...), vspace)
end
r_RR(x::AbstractGMPS, ys::AbstractGMPS...) = r_RR(ones, edgespace(x), x, ys...)
