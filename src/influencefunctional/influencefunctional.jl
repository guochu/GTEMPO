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
struct PartialIF <: InfluenceFunctionalAlgorithm 
	trunc::TruncationDimCutoff
	verbosity::Int
end
PartialIF(; trunc::TruncationDimCutoff=DefaultITruncation, verbosity::Int=0) = PartialIF(trunc, verbosity)
"""
	struct XTRGIF

Build the IF as a translational variant MPO
see [SciPost Phys. Core 7, 063 (2024)]
"""
struct XTRGIF{T<:ExponentialExpansionAlgorithm, E<:TimeEvoMPOAlgorithm, M<:DMRGAlgorithm} <: InfluenceFunctionalAlgorithm 
	algexpan::T
	algevo::E
	algmult::M
	# trunc::TruncationDimCutoff
	k::Int
	fast::Bool
	verbosity::Int
end
XTRGIF(; algexpan::ExponentialExpansionAlgorithm=OverDeterminedProny(n=15, tol=1.0e-4, verbosity=0), 
						 algevo::TimeEvoMPOAlgorithm=WII(), 
						 algmult::DMRGAlgorithm=DefaultMultAlg,
						 k::Int=5, 
						 fast::Bool=true,
						 verbosity::Int=0) = XTRGIF(algexpan, algevo, algmult, k, fast, verbosity)

function Base.getproperty(x::XTRGIF, s::Symbol)
	if s == :trunc
		return x.algmult.trunc
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
struct ExactTTIIF{T<:ExponentialExpansionAlgorithm, M<:DMRGAlgorithm, M2<:DMRGAlgorithm} <: InfluenceFunctionalAlgorithm
	algexpan::T
	algmult::M
	algmult2::M2 # only used in iGTEMPO._differentialinfluencefunctional2
	multorder::Symbol
	verbosity::Int
end
ExactTTIIF(; algexpan::ExponentialExpansionAlgorithm=OverDeterminedProny(n=15, tol=1.0e-4, verbosity=0), 
						 algmult::DMRGAlgorithm=DefaultMultAlg, algmult2::DMRGAlgorithm=algmult,
						 multorder::Symbol = :αSM,
						 verbosity::Int=0) = ExactTTIIF(algexpan, algmult, algmult2, multorder, verbosity)
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
- `algexpan::ExponentialExpansionAlgorithm`: exponential (Prony) expansion algorithm for the bath correlation function.
- `trunc::TruncationDimCutoff`: bond dimension of the flow manifold / final influence functional.
- `δ::Float64`: imaginary-time step of the flow (0 < δ ≤ 1, adjusted so that 1/δ is an integer).
- `verbosity::Int`: verbosity level of the output.
- `callback::Function`: callback function invoked after the flow.

On real-time lattices the influence operator driving the flow is the sum of
the 4 branch MPOs returned by `influenceoperators` ((+,+), (+,−), (−,+),
(−,−)), since the site-wise product algebra satisfies e^a∘e^b = e^{a+b}.
"""
struct TDVPIF <: InfluenceFunctionalAlgorithm
	algexpan::ExponentialExpansionAlgorithm
	trunc::TruncationDimCutoff      # bond dimension of the flow manifold / final IF
	δ::Float64                      # imaginary-time step of the flow (0 < δ ≤ 1, adjusted so that 1/δ is an integer)
	verbosity::Int
	callback::Function
end
"""
	TDVPIF(; algexpan, trunc, δ, verbosity, callback)

Keyword constructor for `TDVPIF`; all parameters have default values and usually need not be passed explicitly.
"""
function TDVPIF(; algexpan::ExponentialExpansionAlgorithm=OverDeterminedProny(n=15, tol=1.0e-4, verbosity=0),
				trunc::TruncationDimCutoff=DefaultITruncation,
				δ::Real=0.1,
				verbosity::Int=0,
				callback::Function=Returns(nothing))
	(0 < δ <= 1) || throw(ArgumentError("δ must be a real number in (0, 1], got $δ"))
	# adjust δ so that 1/δ is an integer number of steps
	n = max(round(Int, 1 / δ), 1)
	return TDVPIF(algexpan, trunc, 1 / n, verbosity, callback)
end

# temporary solution
changetrunc(x::DMRG1; trunc=x.trunc) = similar(x, trunc=trunc)
changetrunc(x::DMRG2; trunc=x.trunc) = similar(x, trunc=trunc)
changetrunc(x::SVDCompression; trunc=x.trunc) = similar(x, D=trunc.D, tol=trunc.ϵ)


include("hybridization/hybridization.jl")
include("retardedinteract/retardedinteract.jl")