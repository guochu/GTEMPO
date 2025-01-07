# push!(LOAD_PATH, "../../../src")

# using GTEMPO, JSON, Serialization

include("../../../src/includes.jl")
using JSON, Serialization

function J(D::Real, ω::Real)
	return (D/(2*pi)) * sqrt(1 - (ω/D)^2 ) * 0.1 / 2
end
spectrum_func(D) = SpectrumFunction(ω -> J(D, ω), lb = -D, ub = D)

function transport_model(Lsys::Int, J::Real=1)
	h = ImpurityHamiltonian(bands = Lsys)
	# for band in 1:h.bands
	# 	push!(h, tunneling(band, band, coeff=ϵ_d))
	# end
	for band in 1:h.bands-1
		t = tunneling(band, band+1, coeff=J)
		push!(h, t)
		push!(h, t')
	end
	return h
end

function main(V::Real, t::Real, Lsys::Int; δt=0.1, chi=60, chi2=500)
	β = Inf
	D = 2.
	N = round(Int, t / δt)

	ts = [i*δt for i in 0:N]


	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)
	trunc2 = truncdimcutoff(D=chi2, ϵ=1.0e-10, add_back=0)
	lattice1 = GrassmannLattice(N=N, δt=t/N, bands=1, contour=:real)

	println("number of sites ", length(lattice1))

	leftmu = V / 2
	rightmu = -V / 2


	leftbath = fermionicbath(spectrum_func(D), β=β, μ=leftmu)
	rightbath = fermionicbath(spectrum_func(D), β=β, μ=rightmu)

	leftcorr = correlationfunction(leftbath, lattice1)
	rightcorr = correlationfunction(rightbath, lattice1)

	mpspath = "data/tempo_V$(V)_t$(t)_dt$(δt)_chi$(chi).mps"
	if ispath(mpspath)
		println("load MPS-IF from path ", mpspath)
		mpsI_left, mpsI_right = Serialization.deserialize(mpspath)
	else
		println("computing MPS-IF...")
		@time mpsI_left = hybriddynamics(lattice1, leftcorr, trunc=trunc, band=1)
		@time mpsI_right = hybriddynamics(lattice1, rightcorr, trunc=trunc, band=1)
		Serialization.serialize(mpspath, (mpsI_left, mpsI_right))
		println("save MPS-IF to path ", mpspath)
	end

	# println("Z is ", integrate(mpsI, lattice))
	println("bond dimension of mpsI is ", bond_dimension(mpsI_left), ", ", bond_dimension(mpsI_right))

	# impurity dynamics
	Lhalf = div(Lsys, 2)
	lattice2 = similar(lattice1, bands=2)
	J = 1

	h = ImpurityHamiltonian(bands = 2)
	term = TunnelingTerm(1, 2, coeff=J)
	push!(h, term)
	push!(h, term')

	# calculate right currents
	mpsK = baresysdynamics(lattice2, h, trunc=trunc)
	mpsK = bulkconnection!(mpsK, lattice2, band=1)
	mpsK = boundarycondition!(mpsK, lattice2, band=1)

	adt_right = mpsI_left
	for band in 1:Lsys-1
		adt_right2 = fillband(lattice2, adt_right, band=1)
		tmp = mult(mpsK, adt_right2, trunc=trunc2)
		adt_right = integrateband(lattice2, tmp, band=1)
		canonicalize!(adt_right, alg = Orthogonalize(trunc=trunc2))
	end

	adt_right = bulkconnection!(adt_right, lattice1, band=1)
	adt_right = boundarycondition!(adt_right, lattice1, band=1)
	
	println("bond dimension of adt_right is ", bond_dimension(adt_right))
	println("adt_right scale is ", scaling(adt_right))

	cache = environments(lattice1, adt_right, mpsI_right)
	@time rightcurrents = [cached_electriccurrent_fast(lattice1, rightcorr, k+1, adt_right, mpsI_right, cache=cache) for k in 1:N]

	# calculate left currents
	adt_left = mpsI_right
	mpsK = baresysdynamics(lattice2, h, trunc=trunc)
	mpsK = bulkconnection!(mpsK, lattice2, band=2)
	mpsK = boundarycondition!(mpsK, lattice2, band=2)
	for band in Lsys:-1:2
		adt_left2 = fillband(lattice2, adt_left, band=2)
		tmp = mult(mpsK, adt_left2, trunc=trunc2)
		adt_left = integrateband(lattice2, tmp, band=2)
		canonicalize!(adt_left, alg = Orthogonalize(trunc=trunc2))
	end
	adt_left = bulkconnection!(adt_left, lattice1, band=1)
	adt_left = boundarycondition!(adt_left, lattice1, band=1)

	cache = environments(lattice1, mpsI_left, adt_left)
	@time leftcurrents = [cached_electriccurrent_fast(lattice1, leftcorr, k+1, mpsI_left, adt_left, cache=cache) for k in 1:N]

	data_path = "result/tempo_V$(V)_t$(t)_dt$(δt)_L$(Lsys)_chi$(chi)_chi2$(chi2).json"

	results = Dict("ts"=>ts, "leftcurrent"=>leftcurrents, "rightcurrent"=>rightcurrents)

	println("save results to ", data_path)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return ts, leftcurrents, rightcurrents

end
