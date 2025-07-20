abstract type IntegrationAlgorithm end
struct ExactIntegrate <: IntegrationAlgorithm end
struct BMPSIntegrate{T<:TruncationScheme} <: IntegrationAlgorithm 
	trunc::T
end
BMPSIntegrate(; trunc::TruncationScheme=DefaultIntegrationTruncation) = BMPSIntegrate(trunc)



include("acintegrate/util.jl")
include("acintegrate/ac_integrate.jl")
include("acintegrate/integrate1.jl")
include("acintegrate/integrate2.jl")
include("acintegrate/integrate3.jl")
include("acintegrate/integrate4.jl")
include("acintegrate/integrate5.jl")

# include("acintegrate/integrate6.jl")
# include("acintegrate/integrate7.jl")

include("acintegrate/ac_cached_integrate.jl")
include("acintegrate/ac_bmps_integrate.jl")

include("conversion.jl")

"""
	integrate(alg::IntegrationAlgorithm, lattice::AbstractGrassmannLattice, x0::GrassmannMPS, x1::GrassmannMPS...; trunc)

integrate all the pairs of GVs of the result of the multiplication of x0, x1...
return the result as a single scalar
alg is the algorithm used for the integration, for which two choices are supported currently
ExactIntegrate: do the integration exactly using the zipup algorithm [see PRB 109, 045140 (2024)]
BMPSIntegrate: do the integration approximately using the BoundaryMPS algorithm, in this case the keyword 
trunc can be specified to control the truncation
"""
function integrate(alg::IntegrationAlgorithm, lattice::AbstractGrassmannLattice, x0::GrassmannMPS, x1::Vararg{GrassmannMPS}; kwargs...)
	if ConjugationStyle(lattice) isa AdjacentConjugation
		return _ac_integrate(alg, lattice, x0, x1...; kwargs...)
	else
		r = toadjacentordering(lattice, x0, x1...; kwargs...)
		return integrate(alg, first(r), Base.tail(r)...; kwargs...)
	end
end
integrate(lattice::AbstractGrassmannLattice, x0::GrassmannMPS, x1::Vararg{GrassmannMPS}; 
			alg::IntegrationAlgorithm=ExactIntegrate(), kwargs...) = integrate(alg, lattice, x0, x1...; kwargs...)

function integrate(lattice::AbstractGrassmannLattice, x0::Vector{<:GrassmannMPS}, x1::Vararg{GrassmannMPS}; kwargs...)
	r = zero(scalartype(lattice))
	for item in x0
		r += integrate(lattice, item, x1...; kwargs...)
	end	
	return r
end

# more complicatd integration
include("parallelintegrate.jl")
include("partialintegrate.jl")
include("partialintegrate2.jl")