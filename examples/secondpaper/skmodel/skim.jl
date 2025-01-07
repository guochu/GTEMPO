push!(LOAD_PATH, "../../../src")
using GTEMPO

using JSON, Serialization


J(D, ε) = sqrt(D^2-ε^2)/pi
spectrum_func(D=1) = SpectrumFunction(ω -> J(D, ω), lb = -D, ub = D)


function main(β, J=0.5)
	# β = 5.
	norb = 2
	U = 2.
	# J = 0.5
	μ = (3*U - 5*J)/2
	# μ = U / 2

	δτ = 0.1
	N = round(Int, β/δτ)
	chi = 2048

	τs = [i*δτ for i in 0:N]

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-9, add_back=0)
	lattice = GrassmannLattice(N=N, δτ=β/N, bands=2*norb, contour=:imag, ordering=A1A1B1B1())

	println("number of sites ", length(lattice))

	bath = fermionicbath(spectrum_func(), β=β, μ=0)
	exact_model = SKIM(bath, U=U, J=J, μ=-μ, norb=norb)

	mpspath = "data/tempo1_norb$(norb)_beta$(β)_N$(N)_chi$(chi)_b.mps"
	if ispath(mpspath)
		println("load MPS-IF from path ", mpspath)
		mpsI1, mpsI2 = Serialization.deserialize(mpspath)
	else
		corr = correlationfunction(bath, lattice)
		println("computing MPS-IF...")
		# @time mpsIs = [hybriddynamics(lattice, corr, trunc=trunc, band=band) for band in 1:lattice.bands]
		mpsI1 = vacuumstate(lattice)
		@time for band in 1:norb
			mpsI1 = hybriddynamics(mpsI1, lattice, corr, trunc=trunc, band=band)
		end

		mpsI2 = vacuumstate(lattice)
		@time for band in norb+1:2*norb
			mpsI2 = hybriddynamics(mpsI2, lattice, corr, trunc=trunc, band=band)
		end
		println("save MPS-IF to path ", mpspath)
		Serialization.serialize(mpspath, (mpsI1, mpsI2))
	end

	# println("Z is ", integrate(mpsI, lattice))
	# println("mpsI scale is ", scale(mpsI))
	println("bond dimension of mpsI is ", bond_dimension(mpsI1), " ", bond_dimension(mpsI2))

	truncK = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)
	@time mpsK = accsysdynamics_fast(lattice, exact_model, trunc=truncK, scaling=100)
	mpsKs = [mpsK]
	for band in 1:lattice.bands
		mpsKs = boundarycondition_branching(mpsKs, lattice, band=band)
	end
	println("bond dimension of mpsK is ", bond_dimension(mpsKs[1]), ", number of Ks ", length(mpsKs))

	@time g = parallel_Gτ(lattice, mpsKs, mpsI1, mpsI2)

	data_path = "result/anderson_tempo1_norb$(norb)_beta$(β)_U$(U)_J$(J)_mu$(μ)_N$(N)_b.json"

	results = Dict("ts"=>τs, "gf" => g, "bd"=>bond_dimensions(mpsI1))

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end


end