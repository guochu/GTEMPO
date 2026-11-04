

# coeff * α^n, α must be in [0, 1]
"""
	struct ExponentialDecayTerm{M1<:SiteOperator, M<:SiteOperator, M2, T <:Number}

Exponential decay of the form coeff * [â ⊗ (αm̂)^⊗n ⊗ b̂]
"""
struct ExponentialDecayTerm{M1<:SiteOperator, M<:SiteOperator, M2, T <:Number} <: AbstractLongRangeTerm
    a::M1
    m::M
    b::M2
    α::T
    coeff::T
end

function ExponentialDecayTerm(a::SiteOperator, b::SiteOperator; middle::MPSBondTensor=id(physpace(a)), α::Number=1., coeff::Number=1.)
    T = promote_type(typeof(α), typeof(coeff))
    check_a_b(a, b)
    return ExponentialDecayTerm(a, middle, b, convert(T, α), convert(T, coeff))
end

TK.scalartype(::Type{ExponentialDecayTerm{M1, M, M2, T}}) where {M1, M, M2, T} = promote_type(scalartype(M1), scalartype(M), scalartype(M2), T)
TK.spacetype(::Type{ExponentialDecayTerm{M1, M, M2, T}}) where {M1, M, M2, T} = spacetype(M1)

Base.adjoint(x::ExponentialDecayTerm) = ExponentialDecayTerm(_op_adjoint(x.a, x.m, x.b)..., conj(x.α), conj(coeff(x)))

# coeff * [â ⊗ b̂] at distance 1, nothing beyond
"""
    NearestNeighbourTerm{M1<:SiteOperator, M2, T <:Number}

Finite-range term of the form coeff * [â ⊗ b̂] that contributes only at
distance 1. Unlike `ExponentialDecayTerm` the closing operator carries no
pole factor, so the coefficient lands exactly on the nearest-neighbour pair
and nothing leaks to larger distances. This is what makes an explicit
distance-1 contact representable: a pole α → 0 exponential channel cannot
express it (the α on the closing operator kills the channel while the
compensating 1/α coefficient degenerates the Schur assembly).
"""
struct NearestNeighbourTerm{M1<:SiteOperator, M2, T <:Number} <: AbstractLongRangeTerm
    a::M1
    b::M2
    coeff::T
end

function NearestNeighbourTerm(a::SiteOperator, b::SiteOperator; coeff::Number=1.)
    check_a_b(a, b)
    return NearestNeighbourTerm(a, b, coeff)
end

Base.adjoint(x::NearestNeighbourTerm) = NearestNeighbourTerm(x.a', x.b', conj(coeff(x)))
TK.scalartype(::Type{NearestNeighbourTerm{M1, M2, T}}) where {M1, M2, T} = promote_type(scalartype(M1), scalartype(M2), T)
TK.spacetype(::Type{NearestNeighbourTerm{M1, M2, T}}) where {M1, M2, T} = spacetype(M1)
_op_adjoint(a::MPSBondTensor, m::MPSBondTensor, b::MPSBondTensor) = (a', m', b')
_op_adjoint(a::MPOTensor, m::MPSBondTensor, b::MPOTensor) = (unsafe_mpotensor_adjoint(a), m', unsafe_mpotensor_adjoint(b))


function _longrange_schurmpo_util(h1, h2s::Vector{<:AbstractLongRangeTerm})
    isempty(h2s) && throw(ArgumentError("empty interactions."))
	pspace = physpace(h2s[1].a)
	N = length(h2s)
	T = Float64
	for item in h2s
		T = promote_type(T, scalartype(item))
	end
	cell = Matrix{Any}(undef, N+2, N+2)
	for i in 1:length(cell)
		cell[i] = zero(T)
	end
	# diagonals
	cell[1, 1] = 1
	cell[end, end] = 1
	cell[1, end] = h1
	for i in 1:N
        h2 = h2s[i]
        if isa(h2.a, MPSBondTensor)
            b_iden = id(T, oneunit(spacetype(h2.a)))
        else
            b_iden = id(T, space_r(h2.a)')
        end
        if h2 isa NearestNeighbourTerm
            # distance-1 only: the channel state is dead (a zero tensor keeps
            # the row/col aux-space layout) and the closing operator carries
            # no pole, so coeff lands exactly on the nearest-neighbour pair
            @tensor nn_iden[-1 -2; -3 -4] := b_iden[-1, -3] * id(physpace(h2.a))[-2, -4]
            cell[i+1, i+1] = zero(T) * nn_iden
            cell[1, i+1] = h2.coeff * h2.a
            cell[i+1, end] = h2.b
        elseif h2 isa ExponentialDecayTerm
            m = h2.m
            @tensor iden[-1 -2; -3 -4] := b_iden[-1, -3] * m[-2, -4]
            cell[i+1, i+1] = h2.α * iden
            cell[1, i+1] = h2.coeff * h2.a
            cell[i+1, end] = h2.α * h2.b
        else
            throw(ArgumentError("unsupported long-range term $(typeof(h2)); expand it into exponential channels first"))
        end
	end
	return SchurMPOTensor(cell)
end

function check_a_b(a::SiteOperator, b::SiteOperator)
    if isa(a, MPOTensor)
        isa(b, MPOTensor) || throw(ArgumentError("a and b must both be MPOTensor or MPSBondTensor"))
        s_l = space_l(a)
        (s_l == space_r(b)' == oneunit(s_l)) || throw(ArgumentError("only strict MPOTensor is allowed"))
    else
        isa(b, MPSBondTensor) || throw(ArgumentError("a and b must both be MPOTensor or MPSBondTensor"))
    end    
end

"""
    SchurMPOTensor(h1::ScalarSiteOp, h2s::Vector{<:AbstractLongRangeTerm})
    SchurMPOTensor(h2s::Vector{<:AbstractLongRangeTerm})

Return an SchurMPOTensor, with outer matrix size (N+2)×(N+2) (N=length(h2s))
Algorithm reference: "Time-evolving a matrix product state with long-ranged interactions"
"""
SchurMPOTensor(h1::MPSBondTensor, h2s::Vector{<:AbstractLongRangeTerm}) = _longrange_schurmpo_util(h1, h2s)
SchurMPOTensor(h2s::Vector{<:AbstractLongRangeTerm}) = _longrange_schurmpo_util(0., h2s)
