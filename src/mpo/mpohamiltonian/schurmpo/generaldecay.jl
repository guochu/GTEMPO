

"""
    struct GenericDecayTerm{M1, M<:MPSBondTensor, M2, F, T <: Number}

Generic decay of the form coeff * [â ⊗ f(n)*m̂^⊗n ⊗ b̂]
"""
struct GenericDecayTerm{M1, M<:MPSBondTensor, M2, F, T <: Number} <: AbstractLongRangeTerm
    a::M1
    m::M
    b::M2
    f::F
    coeff::T
end

GenericDecayTerm(a::SiteOperator, b::SiteOperator, f; middle::MPSBondTensor = id(physical_space(a)), coeff::Number=1.) = GenericDecayTerm(a, middle, b, f, coeff)

function GenericDecayTerm(a::SiteOperator, b::SiteOperator; 
                            middle::MPSBondTensor = id(physical_space(a)), f, coeff::Number=1.) 
    check_a_b(a, b)
    GenericDecayTerm(a, middle, b, f, coeff)
end 
TK.scalartype(::Type{GenericDecayTerm{M1, M, M2, F, T}}) where {M1, M, M2, F<:AbstractVector, T} = promote_type(scalartype(M1), scalartype(M), scalartype(M2), eltype(F), T)
TK.scalartype(x::GenericDecayTerm{M1, M, M2, F, T}) where {M1, M, M2, F<:AbstractVector, T} = scalartype(typeof(x))
TK.scalartype(x::GenericDecayTerm{M1, M, M2, F, T}) where {M1, M, M2, F, T} = promote_type(scalartype(M1), scalartype(M), scalartype(M2), T, typeof(x.f(0.)))
Base.adjoint(x::GenericDecayTerm) = GenericDecayTerm(_op_adjoint(x.a, x.m, x.b)..., _conj(x.f), conj(coeff(x)))
TK.spacetype(::Type{GenericDecayTerm{M1, M, M2, F, T}}) where {M1, M, M2, F, T} = spacetype(M1)

_conj(f) = x->conj(x.f(x))
_conj(f::AbstractVector) = conj(f)

""" 
    PowerlawDecayTerm(a::M, b::M; α::Number=1., kwargs...)

coeff * n^α, α should be negative in principle (diverging otherwise)
"""
PowerlawDecayTerm(a::M, b::M; α::Number=1., kwargs...) where {M<:SiteOperator} = GenericDecayTerm(a, b; f=x->x^α, kwargs...)



# L is the number of sites

"""
    exponential_expansion(x::GenericDecayTerm{M, T, F}; len::Int, alg)

Convert a GenericDecayTerm into a list of ExponentialDecayTerm
"""
function exponential_expansion(x::GenericDecayTerm{M1, M, M2, F, T}; len::Union{Int, Nothing}=nothing, alg::ExponentialExpansionAlgorithm=PronyExpansion()) where {M1, M, M2, F, T}
    if F <: AbstractVector
        xs, lambdas = exponential_expansion(x.f, alg=alg)
        isa(len, Int) && println("key len ignored")
    else
        isa(len, Int) || throw(ArgumentError("key len should be Int when F is not a vector"))
        xs, lambdas = exponential_expansion(x.f, len-1, alg=alg)
    end
    r = ExponentialDecayTerm{M1, M, M2, eltype(lambdas)}[]
    for (c, alpha) in zip(xs, lambdas)
        push!(r, ExponentialDecayTerm(x.a, x.b; middle=x.m, α=alpha, coeff=c * coeff(x)))
    end
    return r
end
