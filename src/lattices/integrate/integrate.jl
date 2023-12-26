abstract type IntegrationAlgorithm end
struct ExactIntegrate <: IntegrationAlgorithm end
struct BMPSIntegrate <: IntegrationAlgorithm end



include("acintegrate/util.jl")
include("acintegrate/ac_integrate.jl")
include("acintegrate/integrate1.jl")
include("acintegrate/integrate2.jl")
include("acintegrate/integrate3.jl")
include("acintegrate/integrate4.jl")
include("acintegrate/integrate6.jl")
include("acintegrate/integrate7.jl")

include("acintegrate/ac_cached_integrate.jl")
include("acintegrate/ac_bmps_integrate.jl")

include("conversion.jl")


function _integrate(alg::IntegrationAlgorithm, lattice::AbstractGrassmannLattice, x0::GrassmannMPS, x1::GrassmannMPS...; kwargs...)
	if ConjugationStyle(lattice) isa AdjacentConjugation
		return _ac_integrate(alg, lattice, x0, x1...; kwargs...)
	else
		r = toadjacentordering(lattice, x0, x1...; kwargs...)
		return _integrate(alg, first(r), Base.tail(r)...; kwargs...)
	end
end
integrate(lattice::AbstractGrassmannLattice, x0::GrassmannMPS, x1::GrassmannMPS...; 
			alg::IntegrationAlgorithm=ExactIntegrate(), kwargs...) = _integrate(alg, lattice, x0, x1...; kwargs...)

function integrate(lattice::AbstractGrassmannLattice, x0::Vector{<:GrassmannMPS}, x1::GrassmannMPS...; kwargs...)
	r = zero(scalartype(lattice))
	for item in x0
		r += integrate(lattice, item, x1...; kwargs...)
	end	
	return r
end
