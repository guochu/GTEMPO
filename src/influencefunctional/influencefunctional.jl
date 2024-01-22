abstract type InfluenceFunctionalAlgorithm end
"""
	struct PartialIF
"""
struct PartialIF <: InfluenceFunctionalAlgorithm end
"""
	struct FullIF
"""
struct FullIF{T<:ExponentialExpansionAlgorithm, E<:TimeEvoMPOAlgorithm} <: InfluenceFunctionalAlgorithm 
	algexpan::T
	algevo::E
end
FullIF(; algexpan::ExponentialExpansionAlgorithm=PronyExpansion(), algevo::TimeEvoMPOAlgorithm=WII()) = FullIF(algexpan, algevo)

include("util.jl")
include("imaginarytime.jl")
include("realtime.jl")

# two algorithms to build the IF
include("partialif.jl")
include("fullif.jl")

