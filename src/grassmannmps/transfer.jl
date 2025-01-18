

struct GrassmannTransferMatrix{M<:MPSTensor, N}
	states::NTuple{N, Vector{GrassmannTensorMap{M}}}
	scaling::Float64

	function GrassmannTransferMatrix{M, N}(states::NTuple{N, Vector{M}}, scaling::Real=1) where {M<:MPSTensor, N}
		(N > 0) || throw(ArgumentError("no element"))
		L = length(states[1])
		all(x->length(x)==L, states) || throw(ArgumentError("elements have different lengths"))
		new{M, N}(map(x->GrassmannTensorMap.(x), states), convert(Float64, scaling))
	end
end
GrassmannTransferMatrix(states::NTuple{N, Vector{M}}, scaling::Real=1) where {M<:MPSTensor, N} = GrassmannTransferMatrix{M, N}(states, scaling)
GrassmannTransferMatrix(states::Vararg{Vector{M}, N}; scaling::Real=1) where {M<:MPSTensor, N} = GrassmannTransferMatrix{M, N}(states, scaling)

Base.length(x::GrassmannTransferMatrix) = length(x.states[1])
TK.scalartype(::Type{GrassmannTransferMatrix{M, N}}) where {M, N} = scalartype(M)
scaling(x::GrassmannTransferMatrix) = x.scaling


function update_pair_left end
function update_pair_right end

function Base.:*(left::GrassmannTensorMap{<:AbstractTensorMap{<:Number, S, 1, N}}, m::GrassmannTransferMatrix{M, N}) where {S, M, N}
	@assert length(m) % 2 == 0
	for i in 1:div(length(m), 2)
		left = lmul!(scaling(m), update_pair_left(left, i, m.states...)) 
	end
	return left
end
function Base.:*(m::GrassmannTransferMatrix{M, N}, right::GrassmannTensorMap{<:AbstractTensorMap{<:Number, S, N, 1}}) where {S, M, N}
	@assert length(m) % 2 == 0
	for i in div(length(m), 2):-1:1
		right = lmul!(scaling(m), update_pair_right(right, i, m.states...)) 
	end
	return right
end

function DMRG.l_LL(f, vspace::ElementarySpace, m::GrassmannTransferMatrix)
	return GrassmannTensorMap(f(scalartype(m), vspace, ⊗(map(y->space_l(y[1].data), m.states)...)))
end

function DMRG.r_RR(f, vspace::ElementarySpace, m::GrassmannTransferMatrix)
	return GrassmannTensorMap(f(scalartype(m), ⊗(map(y->space_r(y[end].data)', reverse(m.states))...), vspace))
end


GrassmannTransferMatrix(states::Vararg{M, N}) where {M <: GrassmannMPS, N} = GrassmannTransferMatrix(map(x->x.data, states), scaling(states...)^2)

function GrassmannTransferMatrix(j::Int, states::Vararg{M, N}) where {M <: GrassmannMPS, N}
	posa, posb = 2*j-1, 2*j
	return GrassmannTransferMatrix(map(x->[x[posa], x[posb]], states), scaling(states...)^2)
end 


scaling(x::GrassmannMPS, y::GrassmannMPS, zs::GrassmannMPS...) = scaling(x) * scaling(y) * prod(map(scaling, zs))