abstract type InfluenceFunctionalAlgorithm end
"""
	struct PartialIF

Build the IF as the product of partial MPOs, each with D=2
"""
struct PartialIF <: InfluenceFunctionalAlgorithm end
"""
	struct TranslationInvariantIF

Build the IF as a translational variant MPO
"""
struct TranslationInvariantIF{T<:ExponentialExpansionAlgorithm, E<:TimeEvoMPOAlgorithm} <: InfluenceFunctionalAlgorithm 
	algexpan::T
	algevo::E
	k::Int
end
TranslationInvariantIF(; algexpan::ExponentialExpansionAlgorithm=PronyExpansion(tol=1.0e-5, verbosity=0), algevo::TimeEvoMPOAlgorithm=WII(), k::Int=5) = TranslationInvariantIF(
						algexpan, algevo, k)

include("util.jl")
include("imaginarytime.jl")
include("realtime.jl")

# two algorithms to build the IF
include("partialif.jl")
include("fullif.jl")

