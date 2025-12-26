abstract type AbstractMPO{A<:MPOTensor} end
abstract type AbstractFiniteMPO{A<:MPOTensor} <: AbstractMPO{A} end

TK.scalartype(::Type{<:AbstractMPO{A}}) where {A<:MPOTensor} = scalartype(A)
TK.spacetype(::Type{<:AbstractMPO{A}}) where {A<:MPOTensor} = spacetype(A)
TK.spacetype(m::AbstractMPO) = spacetype(typeof(m))

space_l(state::AbstractMPO) = space_l(state[1])
space_r(state::AbstractMPO) = space_r(state[end])
sector(a::AbstractMPO) = _sector(a)

