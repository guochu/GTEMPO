push!(LOAD_PATH, "../../../src")


using JSON
using Serialization

using GTEMPO

function J(D::Real, ω::Real)
	return (D/pi) * sqrt(1 - (ω/D)^2 ) * 0.1
end

spectrum_func(D) = SpectrumFunction(ω -> J(D, ω), lb = -D, ub = D)

function vacuum(lattice::RealGrassmannLattice, i::Int, A::GrassmannMPS, B::GrassmannMPS...; 
            Z::Number = integrate(lattice, A, B...))
    pos1, pos2 = index(lattice, i, conj=false, branch=:+, band=1), index(lattice, i, conj=true, branch=:+, band=1)
    pos1′, pos2′ = index(lattice, i, conj=false, branch=:+, band=2), index(lattice, i, conj=true, branch=:+, band=2)
    t = GTerm(pos1, pos2, pos1′, pos2′, coeff=1)
    A2 = mult(t, A)
    g = integrate(lattice, A2, B...)/Z
    return g
end

function nup(lattice::RealGrassmannLattice, i::Int, A::GrassmannMPS, B::GrassmannMPS...; 
            Z::Number = integrate(lattice, A, B...))
    pos1, pos2 = index(lattice, i, conj=false, branch=:+, band=1), index(lattice, i, conj=true, branch=:-, band=1)
    pos1′, pos2′ = index(lattice, i, conj=false, branch=:+, band=2), index(lattice, i, conj=true, branch=:+, band=2)
    t = GTerm(pos1, pos2, pos1′, pos2′, coeff=1)
    A2 = mult(t, A)
    g = integrate(lattice, A2, B...)/Z
    return g
end

function ndown(lattice::RealGrassmannLattice, i::Int, A::GrassmannMPS, B::GrassmannMPS...; 
            Z::Number = integrate(lattice, A, B...))
    pos1, pos2 = index(lattice, i, conj=false, branch=:+, band=1), index(lattice, i, conj=true, branch=:+, band=1)
    pos1′, pos2′ = index(lattice, i, conj=false, branch=:+, band=2), index(lattice, i, conj=true, branch=:-, band=2)
    t = GTerm(pos1, pos2, pos1′, pos2′, coeff=1)
    A2 = mult(t, A)
    g = integrate(lattice, A2, B...)/Z
    return g
end

function nn(lattice::RealGrassmannLattice, i::Int, A::GrassmannMPS, B::GrassmannMPS...; 
            Z::Number = integrate(lattice, A, B...))
    pos1, pos2 = index(lattice, i, conj=false, branch=:+, band=1), index(lattice, i, conj=true, branch=:-, band=1)
    pos1′, pos2′ = index(lattice, i, conj=false, branch=:+, band=2), index(lattice, i, conj=true, branch=:-, band=2)
    t = GTerm(pos1, pos2, pos1′, pos2′, coeff=1)
    A2 = mult(t, A)
    g = integrate(lattice, A2, B...)/Z
    return g
end

function cached_vacuum(lattice::RealGrassmannLattice, i::Int, A::GrassmannMPS, B::GrassmannMPS...; 
            cache=environments(lattice, A, B...))
    pos1, pos2 = index(lattice, i, conj=false, branch=:+, band=1), index(lattice, i, conj=true, branch=:+, band=1)
    pos1′, pos2′ = index(lattice, i, conj=false, branch=:+, band=2), index(lattice, i, conj=true, branch=:+, band=2)
    t = GTerm(pos1, pos2, pos1′, pos2′, coeff=1)
    A2 = t * A
    a, b = positions(t)[1], positions(t)[end]
    return MPSImpuritySolvers.cached_integrate_util(lattice, a, b, cache, A2, B...)
end

function cached_nup(lattice::RealGrassmannLattice, i::Int, A::GrassmannMPS, B::GrassmannMPS...; 
            cache=environments(lattice, A, B...))
    pos1, pos2 = index(lattice, i, conj=false, branch=:+, band=1), index(lattice, i, conj=true, branch=:-, band=1)
    pos1′, pos2′ = index(lattice, i, conj=false, branch=:+, band=2), index(lattice, i, conj=true, branch=:+, band=2)
    t = GTerm(pos1, pos2, pos1′, pos2′, coeff=1)
    A2 = t * A
    a, b = positions(t)[1], positions(t)[end]
    return MPSImpuritySolvers.cached_integrate_util(lattice, a, b, cache, A2, B...)
end

function cached_ndown(lattice::RealGrassmannLattice, i::Int, A::GrassmannMPS, B::GrassmannMPS...; 
            cache=environments(lattice, A, B...))
    pos1, pos2 = index(lattice, i, conj=false, branch=:+, band=1), index(lattice, i, conj=true, branch=:+, band=1)
    pos1′, pos2′ = index(lattice, i, conj=false, branch=:+, band=2), index(lattice, i, conj=true, branch=:-, band=2)
    t = GTerm(pos1, pos2, pos1′, pos2′, coeff=1)
    A2 = t * A
    a, b = positions(t)[1], positions(t)[end]
    return MPSImpuritySolvers.cached_integrate_util(lattice, a, b, cache, A2, B...)
end

function cached_nn(lattice::RealGrassmannLattice, i::Int, A::GrassmannMPS, B::GrassmannMPS...; 
            cache=environments(lattice, A, B...))
    pos1, pos2 = index(lattice, i, conj=false, branch=:+, band=1), index(lattice, i, conj=true, branch=:-, band=1)
    pos1′, pos2′ = index(lattice, i, conj=false, branch=:+, band=2), index(lattice, i, conj=true, branch=:-, band=2)
    t = GTerm(pos1, pos2, pos1′, pos2′, coeff=1)
    A2 = t * A
    a, b = positions(t)[1], positions(t)[end]
    return MPSImpuritySolvers.cached_integrate_util(lattice, a, b, cache, A2, B...)
end




# Γ = 0.1

function main(δt, order=7)
	β = 20.
	U = 0.25*pi
	ϵ_d=-U/2
	D = 1

	# δt = 0.05
	t = 50.
	N = round(Int, t / δt)

	bath = fermionicbath(spectrum_func(D), β=β, μ=0)
	exact_model = SISB(bath, μ = ϵ_d, U=U)

	bands = 2
	lattice = GrassmannLattice(N=N, δt=δt, contour=:real, order=1, bands=bands)
	println("number of sites, ", length(lattice))


	corr = correlationfunction(exact_model.bath, lattice)

	trunc = truncdimcutoff(D=1024, ϵ=10.0^(-order), add_back=0)
	truncK = truncdimcutoff(D=1024, ϵ=1.0e-10, add_back=0)

	mpspath = "data/tempo1_beta$(β)_N$(N)_order$(order).mps"
	if ispath(mpspath)
		println("load MPS-IF from path ", mpspath)
		mpsI1, mpsI2 = Serialization.deserialize(mpspath)
	else
		println("computing MPS-IF...")
		@time mpsI1 = hybriddynamics(lattice, corr, trunc=trunc, band=1)
		@time mpsI2 = hybriddynamics(lattice, corr, trunc=trunc, band=2)
		# println("Z is ", integrate(mpsI, lattice))
		println("save MPS-IF to path ", mpspath)
		Serialization.serialize(mpspath, (mpsI1, mpsI2))
	end
	bds = bond_dimensions(mpsI1)

	println("mpsI bond dimension is ", bond_dimension(mpsI1), " ", bond_dimension(mpsI2))


	@time mpsK = sysdynamics!(vacuumstate(lattice), lattice, exact_model, trunc=truncK)
	println("mpsK bond dimension is ", bond_dimension(mpsK))
	mpsK = boundarycondition(mpsK, lattice, band=1)
	mpsK = boundarycondition(mpsK, lattice, band=2)

	cache = environments(lattice, mpsK, mpsI1, mpsI2)
	@time nup_ = [cached_nup(lattice, i, mpsK, mpsI1, mpsI2, cache=cache) for i in 1:lattice.N]
	@time ndown_ = [cached_ndown(lattice, i, mpsK, mpsI1, mpsI2, cache=cache) for i in 1:lattice.N]
	@time nn_ = [cached_nn(lattice, i, mpsK, mpsI1, mpsI2, cache=cache) for i in 1:lattice.N]
	@time vacuum_ = [cached_vacuum(lattice, i, mpsK, mpsI1, mpsI2, cache=cache) for i in 1:lattice.N]


	ts = [i*δt for i in 1:N]

	data_path = "result/tempo1_beta$(β)_N$(N)_order$(order).json"
	results = Dict("ts"=>ts, "nup"=>real(nup_), "ndown"=>real(ndown_), "nn"=>real(nn_), "vacuum"=>real(vacuum_), "bd"=>bds)
	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return nup_, ndown_, nn_, vacuum_
end



