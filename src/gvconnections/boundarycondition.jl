"""
	boundarycondition(x::GrassmannMPS, lattice::AbstractGrassmannLattice; kwargs...)

Apply the boundary connection term onto the given GMPS, which connects the GVs on step 1 and N
One could also use the function boundarycondition_branching instead, which return a vector
of GMPSs instead of a single one. The sum of these GMPS will be equal to the result of the 
current function. But the sum may be performed on the fly for efficiency [see PRB 109, 165113 (2024)]

Note that we have different boundary conditions for different contours
"""
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
	pos1, pos2 = index(lattice, lattice.kt, conj=true, band=band, branch=:-), index(lattice, lattice.kt, conj=false, band=band, branch=:+)
	apply!(exp(GTerm(pos1, pos2, coeff=1)), x)

	pos1, pos2 = index(lattice, 1, conj=true, band=band, branch=:τ), index(lattice, 1, conj=false, band=band, branch=:-)
	apply!(exp(GTerm(pos1, pos2, coeff=1)), x)
	
	pos1, pos2 = index(lattice, 0, conj=true, band=band), index(lattice, lattice.kτ, conj=false, band=band, branch=:τ)
	apply!(exp(GTerm(pos1, pos2, coeff=-1)), x)
	canonicalize!(x, alg=Orthogonalize(trunc=trunc))
	pos1, pos2 = index(lattice, 1, conj=true, band=band, branch=:+), index(lattice, 0, conj=false, band=band)
	apply!(exp(GTerm(pos1, pos2, coeff=1)), x)
	canonicalize!(x, alg=Orthogonalize(trunc=trunc))
	return x
end
# function boundarycondition_branching(x0::GrassmannMPS, lattice::MixedGrassmannLattice; band::Int=1, trunc::TruncationScheme=DefaultIntegrationTruncation)
# 	pos1, pos2 = index(lattice, lattice.kt, conj=true, band=band, branch=:-), index(lattice, lattice.kt, conj=false, band=band, branch=:+)
# 	x = exp(GTerm(pos1, pos2, coeff=1)) * x0
	
# 	pos1, pos2 = index(lattice, 1, conj=true, band=band, branch=:τ), index(lattice, 1, conj=false, band=band, branch=:-)
# 	x = exp(GTerm(pos1, pos2, coeff=1)) * x
# 	canonicalize!(x, alg=Orthogonalize(trunc=trunc))

# 	pos1, pos2 = index(lattice, 1, conj=true, band=band, branch=:+), index(lattice, 0, conj=false, band=band)
# 	x = exp(GTerm(pos1, pos2, coeff=1)) * x
# 	canonicalize!(x, alg=Orthogonalize(trunc=trunc))

# 	pos1, pos2 = index(lattice, 0, conj=true, band=band), index(lattice, lattice.kτ, conj=false, band=band, branch=:τ)
# 	x2 = GTerm(pos1, pos2, coeff=-1) * x
# 	canonicalize!(x2, alg=Orthogonalize(trunc=trunc))
# 	return x, x2
# end