boundarycondition(x::GrassmannMPS, lattice::AbstractGrassmannLattice; kwargs...) = boundarycondition!(copy(x), lattice; kwargs...)

function boundarycondition_branching(x0::Vector{<:GrassmannMPS}, lattice::AbstractGrassmannLattice; kwargs...)
	r = eltype(x0)[]
	for item in x0
		y1, y2 = boundarycondition_branching(item, lattice; kwargs...)
		push!(r, y1)
		push!(r, y2)
	end
	return r
end

# imaginary-time
function boundarycondition!(x::GrassmannMPS, lattice::ImagGrassmannLattice; band::Int=1, trunc::TruncationScheme=DefaultIntegrationTruncation)
	pos1, pos2 = index(lattice, 0, conj=true, band=band), index(lattice, lattice.k, conj=false, band=band)
	apply!(exp(GTerm(pos1, pos2, coeff=-1)), x)
	canonicalize!(x, alg=Orthogonalize(trunc=trunc))
	pos1, pos2 = index(lattice, 1, conj=true, band=band), index(lattice, 0, conj=false, band=band)
	apply!(exp(GTerm(pos1, pos2, coeff=1)), x)
	canonicalize!(x, alg=Orthogonalize(trunc=trunc))
	return x
end
function boundarycondition_branching(x0::GrassmannMPS, lattice::ImagGrassmannLattice; band::Int=1, trunc::TruncationScheme=DefaultIntegrationTruncation)
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

# real-time
function boundarycondition!(x::GrassmannMPS, lattice::RealGrassmannLattice; band::Int=1, trunc::TruncationScheme=DefaultIntegrationTruncation)
	pos1, pos2 = index(lattice, 0, conj=true, band=band, branch=:+), index(lattice, lattice.k, conj=false, band=band, branch=:+)
	apply!(exp(GTerm(pos1, pos2, coeff=-1)), x)
	canonicalize!(x, alg=Orthogonalize(trunc=trunc))
	pos1, pos2 = index(lattice, lattice.k, conj=true, band=band, branch=:-), index(lattice, 0, conj=false, band=band, branch=:-)
	apply!(exp(GTerm(pos1, pos2, coeff=1)), x)
	canonicalize!(x, alg=Orthogonalize(trunc=trunc))
	return x
end
function boundarycondition_branching(x0::GrassmannMPS, lattice::RealGrassmannLattice; band::Int=1, trunc::TruncationScheme=DefaultIntegrationTruncation) 
	(LayoutStyle(lattice) isa BranchLocalLayout) || throw(ArgumentError("boundarycondition_branching only work with BranchLocalLayout for RealGrassmannLattice"))
	pos1, pos2 = index(lattice, 0, conj=true, band=band, branch=:+), index(lattice, lattice.k, conj=false, band=band, branch=:+)
	t = exp(GTerm(pos1, pos2, coeff=-1))
	x = t * x0
	canonicalize!(x, alg=Orthogonalize(trunc=trunc))
	pos1, pos2 = index(lattice, lattice.k, conj=true, band=band, branch=:-), index(lattice, 0, conj=false, band=band, branch=:-)
	t = GTerm(pos1, pos2, coeff=1)
	x2 = t * x
	canonicalize!(x2, alg=Orthogonalize(trunc=trunc))
	return x, x2
end
# mixed-time
function boundarycondition!(x::GrassmannMPS, lattice::MixedGrassmannLattice; band::Int=1, trunc::TruncationScheme=DefaultIntegrationTruncation)
	pos1, pos2 = index(lattice, 0, conj=true, band=band), index(lattice, lattice.Nτ, conj=false, band=band, branch=:τ)
	apply!(exp(GTerm(pos1, pos2, coeff=-1)), x)
	canonicalize!(x, alg=Orthogonalize(trunc=trunc))
	pos1, pos2 = index(lattice, 1, conj=true, band=band, branch=:+), index(lattice, 0, conj=false, band=band)
	apply!(exp(GTerm(pos1, pos2, coeff=1)), x)
	canonicalize!(x, alg=Orthogonalize(trunc=trunc))
	return x
end
# function boundarycondition_branching(x0::GrassmannMPS, lattice::MixedGrassmannLattice; band::Int=1, trunc::TruncationScheme=DefaultIntegrationTruncation) 
# 	(LayoutStyle(lattice) isa BranchLocalLayout) || throw(ArgumentError("boundarycondition_branching only work with BranchLocalLayout for MixedGrassmannLattice"))
# 	pos1, pos2 = index(lattice, 0, conj=true, band=band), index(lattice, lattice.Nτ, conj=false, band=band, branch=:τ)
# 	t = exp(GTerm(pos1, pos2, coeff=-1))
# 	x = t * x0
# 	canonicalize!(x, alg=Orthogonalize(trunc=trunc))
# 	pos1, pos2 = index(lattice, 1, conj=true, band=band, branch=:+), index(lattice, 0, conj=false, band=band)
# 	t = GTerm(pos1, pos2, coeff=1)
# 	x2 = t * x
# 	canonicalize!(x2, alg=Orthogonalize(trunc=trunc))
# 	return x, x2
# end