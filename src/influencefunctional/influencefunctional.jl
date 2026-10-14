"""
	abstract type InfluenceFunctionalAlgorithm

Supertype of all influence functional (IF) construction algorithms.
"""
abstract type InfluenceFunctionalAlgorithm end
"""
	struct PartialIF

Build the IF as the product of partial MPOs, each with D=2
see [SciPost Phys. Core 7, 063 (2024)]
"""
@kwdef struct PartialIF <: InfluenceFunctionalAlgorithm
	trunc::TruncateDimCutoff = Defaults.alg_trunc_if()
	verbosity::Int = 0
end
"""
	struct XTRGIF

Build the IF as a translational variant MPO
see [SciPost Phys. Core 7, 063 (2024)]
"""
@kwdef struct XTRGIF{T<:ExponentialExpansionAlgorithm, E<:TimeEvoMPOAlgorithm, M<:DMRGAlgorithm} <: InfluenceFunctionalAlgorithm
	alg_expan::T = Defaults.alg_expan()
	alg_evo::E = Defaults.alg_evo()
	alg_mult::M = Defaults.alg_mult()
	k::Int = 5
	fast::Bool = true
	verbosity::Int = 0
end

function Base.getproperty(x::XTRGIF, s::Symbol)
	if s == :trunc
		return x.alg_mult.trunc
	else
		getfield(x, s)
	end
end


"""
	struct ExactTTIIF <: InfluenceFunctionalAlgorithm

Build the IF as a translationally invariant MPO by exactly exponentiating the
Trotter decomposition of the quadratic IF kernel. `multorder` controls the
ordering of the exponential decay terms (`:αSM` by default).
"""
@kwdef struct ExactTTIIF{T<:ExponentialExpansionAlgorithm, M<:DMRGAlgorithm} <: InfluenceFunctionalAlgorithm
	alg_expan::T = Defaults.alg_expan()
	alg_mult::M = Defaults.alg_mult()
	multorder::Symbol = :αSM
	verbosity::Int = 0
end
# allowed order:
# :λLM, λ large first
# :λSM, λ small first
# :αLM, α large first
# :αSM, α small first, the default
# :no, no order

"""
	struct TDVPIF

Construct the influence functional with a second-order single-site TDVP
imaginary-time flow: the influence functional is viewed as the "equilibrium
state" IF = exp(H) of the influence operator H (the MPO form returned by
`influenceoperator`), and is computed by evolving the identity influence
functional (β = 0) along the flow dz/dτ = H·z from τ = 0 to τ = 1.

Each flow step is one forward-backward TDVP sweep: the center tensor AC is
evolved by +δτ/2 (Krylov exponentiation of the local effective map) and the
bond matrix C by -δτ/2, with a full +δτ step at the turning site. The state
is zero-padded to bond dimension `trunc.D` before the flow, so the bond
dimension grows with the correlations up to `trunc.D`.

# Fields
- `alg_expan::ExponentialExpansionAlgorithm`: exponential (Prony) expansion algorithm for the bath correlation function.
- `trunc::TruncateDimCutoff`: bond dimension of the flow manifold / final influence functional.
- `δ::Float64`: imaginary-time step of the flow (0 < δ ≤ 1, adjusted so that 1/δ is an integer).
- `verbosity::Int`: verbosity level of the output.

On real-time lattices the influence operator driving the flow is the sum of
the 4 branch MPOs returned by `influenceoperators` ((+,+), (+,−), (−,+),
(−,−)), since the site-wise product algebra satisfies e^a∘e^b = e^{a+b}.
"""
@kwdef struct TDVPIF{T<:ExponentialExpansionAlgorithm} <: InfluenceFunctionalAlgorithm
	alg_expan::T = Defaults.alg_expan()
	trunc::TruncateDimCutoff = Defaults.alg_trunc_if()   # bond dimension of the flow manifold / final IF
	δ::Float64 = 0.1                                # imaginary-time step of the flow (0 < δ ≤ 1; 1/δ is rounded to an integer number of steps)
	verbosity::Int = 0
end


include("hybridization/hybridization.jl")
include("retardedinteract/retardedinteract.jl")