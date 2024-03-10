abstract type InfluenceFunctionalAlgorithm end
"""
	struct PartialIF

Build the IF as the product of partial MPOs, each with D=2
"""
struct PartialIF <: InfluenceFunctionalAlgorithm 
	trunc::TruncationDimCutoff
end
PartialIF(; trunc::TruncationDimCutoff=DefaultITruncation) = PartialIF(trunc)
"""
	struct TranslationInvariantIF

Build the IF as a translational variant MPO
"""
struct TranslationInvariantIF{T<:ExponentialExpansionAlgorithm, E<:TimeEvoMPOAlgorithm, M<:DMRGAlgorithm} <: InfluenceFunctionalAlgorithm 
	algexpan::T
	algevo::E
	algmult::M
	# trunc::TruncationDimCutoff
	k::Int
	verbosity::Int
end
TranslationInvariantIF(; algexpan::ExponentialExpansionAlgorithm=PronyExpansion(n=15, tol=1.0e-4, verbosity=0), 
						 algevo::TimeEvoMPOAlgorithm=WII(), 
						 algmult::DMRGAlgorithm=DefaultMultAlg,
						 k::Int=5, verbosity::Int=0) = TranslationInvariantIF(algexpan, algevo, algmult, k, verbosity)

function Base.getproperty(x::TranslationInvariantIF, s::Symbol)
	if s == :trunc
		return x.algmult.trunc
	else
		getfield(x, s)
	end
end

# temporary solution
changetrunc(x::DMRG1; trunc=x.trunc) = similar(x, D=trunc.D, tolgauge=trunc.ϵ, tol=trunc.ϵ)
changetrunc(x::DMRG2; trunc=x.trunc) = similar(x, trunc=trunc)
changetrunc(x::SVDCompression; trunc=x.trunc) = similar(x, D=trunc.D, tol=trunc.ϵ)

include("util.jl")
include("imaginarytime.jl")
include("realtime.jl")
include("mixedtime.jl")

# two algorithms to build the IF
include("partialif.jl")
include("fullif.jl")

