push!(LOAD_PATH, "../../src")
using GTEMPO

using JSON, Serialization


J(D, ε) = sqrt(D^2-ε^2)/pi
spectrum_func(D=1) = SpectrumFunction(ω -> J(D, ω), lb = -D, ub = D)


function main(t; β=5, δτ=0.1, δt = 0.1, chi=60, chi2=4*chi)
	# β = 5.
	norb = 2
	U = 2.
	J=0.5
	μ = (3*U - 5*J)/2
	# μ = U / 2

	Nt = round(Int, t/δt)
	Nτ = round(Int, β/δτ)
	# chi = 2048

	ts = [i*δt for i in 0:N]

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)
	trunc2 = truncdimcutoff(D=chi2, ϵ=1.0e-10, add_back=0)
	algmult = DMRGMult1(trunc=trunc2)

	lattice = GrassmannLattice(Nτ=Nτ, δτ=δτ, Nt=Nt, δt=δt, bands=2*norb, contour=:Kadanoff)
	lattice1 = similar(lattice, bands=1)

	println("number of sites ", length(lattice))

	bath = fermionicbath(spectrum_func(), β=β, μ=0)
	exact_model = KanamoriIM(U=U, J=J, μ=-μ, norb=norb)

	mpspath = "data/kadanoff_norb$(norb)_beta$(β)_N$(N)_chi$(chi).mps"
	if ispath(mpspath)
		println("load MPS-IF from path ", mpspath)
		mpsI = Serialization.deserialize(mpspath)
	else
		corr = correlationfunction(bath, lattice)
		println("computing MPS-IF...")
		@time mpsI = hybriddynamics(lattice1, corr, trunc=trunc) 

		println("save MPS-IF to path ", mpspath)
		Serialization.serialize(mpspath, mpsI)
	end

	println("bond dimension of mpsI is ", bond_dimension(mpsI))

	@time mpsK = accsysdynamics_fast(lattice, exact_model, trunc=trunc2)
	# mpsK = boundarycondition!(mpsK, lattice, band=1)
	println("bond dimension of mpsK is ", bond_dimension(mpsK))
	mps_adt = mpsK

	lattice_tmp = lattice

	for band in 1:lattice.bands-1
		println("iteration at $(band)-th band...")
		mps_adt = boundarycondition!(mps_adt, lattice_tmp, band=1)
		mpsI1 = fillband(lattice_tmp, mpsI, band=1)

		mps_adt = integrateband(lattice_tmp, mps_adt, mpsI1, algmult, band=1)
		lattice_tmp = similar(lattice_tmp, bands=lattice_tmp.bands-1)
	end
	mps_adt = boundarycondition!(mps_adt, lattice1, band=1)


	
	cache = environments(lattice1, mps_adt, mpsI)

	@time gt = cached_greater_fast(lattice1, mps_adt, mpsI, cache=cache)
	@time lt = cached_lesser_fast(lattice1, mps_adt, mpsI, cache=cache)

	data_path = "result/kadanoff_norb$(norb)_beta$(β)_t$(t)_U$(U)_J$(J)_mu$(μ)_N$(N)_chi$(chi)_chi2$(chi2).json"

	results = Dict("ts"=>ts, "gt" => gt, "lt"=>lt)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end


end