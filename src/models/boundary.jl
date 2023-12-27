boundarycondition(x::GrassmannMPS, lattice::AbstractGrassmannLattice; kwargs...) = boundarycondition!(copy(x), lattice; kwargs...)

function boundarycondition!(x::GrassmannMPS, lattice::ImagGrassmannLattice; band::Int=1, trunc::TruncationScheme=DefaultIntegrationTruncation)
	pos1, pos2 = index(lattice, 0, conj=true, band=band), index(lattice, lattice.k, conj=false, band=band)
	apply!(exp(GTerm(pos1, pos2, coeff=-1)), x)
	canonicalize!(x, alg=Orthogonalize(trunc=trunc))
	pos1, pos2 = index(lattice, 1, conj=true, band=band), index(lattice, 0, conj=false, band=band)
	apply!(exp(GTerm(pos1, pos2, coeff=1)), x)
	canonicalize!(x, alg=Orthogonalize(trunc=trunc))
	return x
end
function boundarycondition2(x0::GrassmannMPS, lattice::ImagGrassmannLattice; band::Int=1, trunc::TruncationScheme=DefaultIntegrationTruncation)
	pos1, pos2 = index(lattice, 0, conj=true, band=band), index(lattice, lattice.k, conj=false, band=band)
	t = exp(GTerm(pos1, pos2, coeff=-1))
	x = t * x0
	canonicalize!(x, alg=Orthogonalize(trunc=trunc))
	pos1, pos2 = index(lattice, 1, conj=true, band=band), index(lattice, 0, conj=false, band=band)
	t = GTerm(pos1, pos2, coeff=1)
	x2 = t * x
	canonicalize!(x2, alg=Orthogonalize(trunc=trunc))
	return x, x2
end
function boundarycondition2(x0::Vector{<:GrassmannMPS}, lattice::AbstractGrassmannLattice; kwargs...)
	r = eltype(x0)[]
	for item in x0
		y1, y2 = boundarycondition2(item, lattice; kwargs...)
		push!(r, y1)
		push!(r, y2)
	end
	return r
end
function boundarycondition!(x::GrassmannMPS, lattice::RealGrassmannLattice; band::Int=1, trunc::TruncationScheme=DefaultIntegrationTruncation)
	pos1, pos2 = index(lattice, 0, conj=true, band=band, forward=true), index(lattice, lattice.k, conj=false, band=band, forward=true)
	apply!(exp(GTerm(pos1, pos2, coeff=-1)), x)
	canonicalize!(x, alg=Orthogonalize(trunc=trunc))
	pos1, pos2 = index(lattice, lattice.k, conj=true, band=band, forward=false), index(lattice, 0, conj=false, band=band, forward=false)
	apply!(exp(GTerm(pos1, pos2, coeff=1)), x)
	canonicalize!(x, alg=Orthogonalize(trunc=trunc))
	return x
end
function boundarycondition2(x0::GrassmannMPS, lattice::RealGrassmannLattice; band::Int=1, trunc::TruncationScheme=DefaultIntegrationTruncation) 
	(LayoutStyle(lattice) isa BranchLocalLayout) || throw(ArgumentError("boundarycondition2 only work with BranchLocalLayout for RealGrassmannLattice"))
	pos1, pos2 = index(lattice, 0, conj=true, band=band, forward=true), index(lattice, lattice.k, conj=false, band=band, forward=true)
	t = exp(GTerm(pos1, pos2, coeff=-1))
	x = t * x0
	canonicalize!(x, alg=Orthogonalize(trunc=trunc))
	pos1, pos2 = index(lattice, lattice.k, conj=true, band=band, forward=false), index(lattice, 0, conj=false, band=band, forward=false)
	t = GTerm(pos1, pos2, coeff=1)
	x2 = t * x
	canonicalize!(x2, alg=Orthogonalize(trunc=trunc))
	return x, x2
end
