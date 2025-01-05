# push!(LOAD_PATH, "../../../src")

# using GTEMPO, JSON, Serialization

include("../../../src/includes.jl")
using JSON, Serialization

J(D, ε) = sqrt(D^2-ε^2)/pi
spectrum_func(D=1) = SpectrumFunction(ω -> J(D, ω), lb = -D, ub = D)

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
	D = 1.
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

	adt_right = mpsI_left
	for band in 1:Lsys-1
		mpsK = baresysdynamics(lattice2, h, trunc=trunc)
		mpsK = bulkcondition!(mpsK, lattice2, band=1)
		mpsK = boundarycondition!(mpsK, lattice2, band=1)
		adt_right2 = fillband(lattice2, adt_right, band=1)

		tmp = mult(mpsK, adt_right2, trunc=trunc2)
		adt_right = integrateband(lattice2, tmp, band=1)
		canonicalize!(adt_right, alg = Orthogonalize(trunc=trunc2))
	end

	adt_right = bulkcondition!(adt_right, lattice1, band=1)
	adt_right = boundarycondition!(adt_right, lattice1, band=1)
	

	println("bond dimension of adt_right is ", bond_dimension(adt_right))
	println("adt_right scale is ", scaling(adt_right))


	cache = environments(lattice1, adt_right, mpsI_right)
	@time gt = [-im*cached_greater(lattice1, k, adt_right, mpsI_right, cache=cache) for k in 1:lattice1.kt]
	@time lt = [im*cached_lesser(lattice1, k, adt_right, mpsI_right, cache=cache) for k in 1:lattice1.kt]

	data_path = "result/tempo_V$(V)_t$(t)_dt$(δt)_L$(Lsys)_chi$(chi)_chi2$(chi2).json"

	results = Dict("ts"=>ts, "gt" => gt, "lt"=>lt)

	println("save results to ", data_path)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return ts, gt

end
