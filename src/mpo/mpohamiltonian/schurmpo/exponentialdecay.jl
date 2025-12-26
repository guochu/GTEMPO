

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

function ExponentialDecayTerm(a::SiteOperator, b::SiteOperator; middle::MPSBondTensor=id(physical_space(a)), α::Number=1., coeff::Number=1.) 
    T = promote_type(typeof(α), typeof(coeff))
    check_a_b(a, b)
    return ExponentialDecayTerm(a, middle, b, convert(T, α), convert(T, coeff))
end

TK.scalartype(::Type{ExponentialDecayTerm{M1, M, M2, T}}) where {M1, M, M2, T} = promote_type(scalartype(M1), scalartype(M), scalartype(M2), T)
TK.spacetype(::Type{ExponentialDecayTerm{M1, M, M2, T}}) where {M1, M, M2, T} = spacetype(M1)

Base.adjoint(x::ExponentialDecayTerm) = ExponentialDecayTerm(_op_adjoint(x.a, x.m, x.b)..., conj(x.α), conj(coeff(x)))
_op_adjoint(a::MPSBondTensor, m::MPSBondTensor, b::MPSBondTensor) = (a', m', b')
_op_adjoint(a::MPOTensor, m::MPSBondTensor, b::MPOTensor) = (unsafe_mpotensor_adjoint(a), m', unsafe_mpotensor_adjoint(b))


function _longrange_schurmpo_util(h1, h2s::Vector{<:ExponentialDecayTerm})
    isempty(h2s) && throw(ArgumentError("empty interactions."))
	pspace = physical_space(h2s[1].a)
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
        if isa(h2s[i].a, MPSBondTensor)
            b_iden = id(T, oneunit(spacetype(h2s[i].a)))
        else
            b_iden = id(T, space_r(h2s[i].a)')
        end
        m = h2s[i].m
        @tensor iden[-1 -2; -3 -4] := b_iden[-1, -3] * m[-2, -4]
		cell[i+1, i+1] = h2s[i].α * iden
		cell[1, i+1] = h2s[i].coeff * h2s[i].a
		cell[i+1, end] = h2s[i].α * h2s[i].b
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
    SchurMPOTensor(h1::ScalarSiteOp, h2s::Vector{<:ExponentialDecayTerm})
    SchurMPOTensor(h2s::Vector{<:ExponentialDecayTerm})

Return an SchurMPOTensor, with outer matrix size (N+2)×(N+2) (N=length(h2s))
Algorithm reference: "Time-evolving a matrix product state with long-ranged interactions"
"""
SchurMPOTensor(h1::MPSBondTensor, h2s::Vector{<:ExponentialDecayTerm}) = _longrange_schurmpo_util(h1, h2s)
SchurMPOTensor(h2s::Vector{<:ExponentialDecayTerm}) = _longrange_schurmpo_util(0., h2s)
