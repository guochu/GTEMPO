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
	struct TranslationInvariantIF

Build the IF as a translational variant MPO
see [SciPost Phys. Core 7, 063 (2024)]
"""
struct TranslationInvariantIF{T<:ExponentialExpansionAlgorithm, E<:TimeEvoMPOAlgorithm, M<:DMRGAlgorithm} <: InfluenceFunctionalAlgorithm 
	algexpan::T
	algevo::E
	algmult::M
	# trunc::TruncationDimCutoff
	k::Int
	fast::Bool
	verbosity::Int
end
TranslationInvariantIF(; algexpan::ExponentialExpansionAlgorithm=PronyExpansion(n=15, tol=1.0e-4, verbosity=0), 
						 algevo::TimeEvoMPOAlgorithm=WII(), 
						 algmult::DMRGAlgorithm=DefaultMultAlg,
						 k::Int=5, 
						 fast::Bool=true,
						 verbosity::Int=0) = TranslationInvariantIF(algexpan, algevo, algmult, k, fast, verbosity)

function Base.getproperty(x::TranslationInvariantIF, s::Symbol)
	if s == :trunc
		return x.algmult.trunc
	else
		getfield(x, s)
	end
end


struct ExactTranslationInvariantIF{T<:ExponentialExpansionAlgorithm, M<:DMRGAlgorithm, M2<:DMRGAlgorithm} <: InfluenceFunctionalAlgorithm 
	algexpan::T
	algmult::M
	algmult2::M2 # only used in iGTEMPO._differentialinfluencefunctional2
	multorder::Symbol
	verbosity::Int
end
ExactTranslationInvariantIF(; algexpan::ExponentialExpansionAlgorithm=PronyExpansion(n=15, tol=1.0e-4, verbosity=0), 
						 algmult::DMRGAlgorithm=DefaultMultAlg, algmult2::DMRGAlgorithm=algmult,
						 multorder::Symbol = :αSM,
						 verbosity::Int=0) = ExactTranslationInvariantIF(algexpan, algmult, algmult2, multorder, verbosity)
# allowed order: 
# :λLM, λ large first
# :λSM, λ small first
# :αLM, α large first
# :αSM, α small first, the default
# :no, no order

# temporary solution
changetrunc(x::DMRGMult1; trunc=x.trunc) = similar(x, trunc=trunc)
changetrunc(x::DMRGMult2; trunc=x.trunc) = similar(x, trunc=trunc)
changetrunc(x::SVDCompression; trunc=x.trunc) = similar(x, D=trunc.D, tol=trunc.ϵ)


include("hybridization/hybridization.jl")
include("retardedinteract/retardedinteract.jl")