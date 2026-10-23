"""
	boundarycondition(x::GrassmannMPS, lattice::AbstractGrassmannLattice; kwargs...)

Apply the boundary connection term onto the given GMPS, which closes the
Grassmann contour directly between the physical end GVs of the lattice.

Note that we have different boundary conditions for different contours
"""
boundarycondition(x::GrassmannMPS, lattice::AbstractGrassmannLattice; kwargs...) = boundarycondition!(copy(x), lattice; kwargs...)

"""
	boundarycondition!(x::GrassmannMPS, lattice::AbstractGrassmannLattice; band=1, trunc=Defaults.alg_trunc_sys())

In-place version of `boundarycondition`: applies the boundary connection term
onto `x` directly and returns it.

Note that we have different boundary conditions for different contours
"""
# imaginary-time
function boundarycondition!(x::GrassmannMPS, lattice::AbstractImagGrassmannLattice; band::Int=1, trunc::TruncationScheme=Defaults.alg_trunc_sys())
	# trace closure (antiperiodic boundary condition) directly between the
	# last and the first physical time slice
	pos1, pos2 = lattice[GrassmannContourIndex(1, :τ, band, true)], lattice[GrassmannContourIndex(lattice.k, :τ, band, false)]
	apply!(exp(GTerm(pos1, pos2, coeff=-1)), x)
	canonicalize!(x, alg=Orthogonalize(trunc=trunc))
	return x
end

# real-time
function boundarycondition!(x::GrassmannMPS, lattice::AbstractRealGrassmannLattice; band::Int=1, trunc::TruncationScheme=Defaults.alg_trunc_sys())
	# connection at the turning point of the Keldysh contour
	pos1, pos2 = lattice[GrassmannContourIndex(lattice.k, :-, band, true)], lattice[GrassmannContourIndex(lattice.k, :+, band, false)]
	apply!(exp(GTerm(pos1, pos2, coeff=-1)), x)
	canonicalize!(x, alg=Orthogonalize(trunc=trunc))
	return x
end
# mixed-time
function boundarycondition!(x::GrassmannMPS, lattice::AbstractMixedGrassmannLattice; band::Int=1, trunc::TruncationScheme=Defaults.alg_trunc_sys())
	# turning point of the real-time part
	pos1, pos2 = lattice[GrassmannContourIndex(lattice.kt, :-, band, true)], lattice[GrassmannContourIndex(lattice.kt, :+, band, false)]
	apply!(exp(GTerm(pos1, pos2, coeff=1)), x)

	# junction between the imaginary-time end (τ=0) and the real-time start
	pos1, pos2 = lattice[GrassmannContourIndex(1, :τ, band, true)], lattice[GrassmannContourIndex(1, :-, band, false)]
	apply!(exp(GTerm(pos1, pos2, coeff=1)), x)

	# closure of the Kadanoff contour: from τ=β back to t=0 on the forward branch
	pos1, pos2 = lattice[GrassmannContourIndex(1, :+, band, true)], lattice[GrassmannContourIndex(lattice.kτ, :τ, band, false)]
	apply!(exp(GTerm(pos1, pos2, coeff=-1)), x)
	canonicalize!(x, alg=Orthogonalize(trunc=trunc))
	return x
end