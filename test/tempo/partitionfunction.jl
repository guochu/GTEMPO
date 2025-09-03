println("------------------------------------")
println("|       Partition Function         |")
println("------------------------------------")



@testset "Partition Function" begin
	rtol = 1.0e-5

	w0 = 0.5
	α = 0.7
	w1 = 1

	trunc = truncdimcutoff(D=60, ϵ=1.0e-10, add_back=0)

	# a toy model
	Himp = w0 .* fermiondensityoperator(2, 1)
	Hbath = w1 .* fermiondensityoperator(2, 2)
	Hhyb = sqrt(α) .* (fermionadagoperator(2, 1) * fermionaoperator(2, 2))
	Hhyb = Hhyb + Hhyb'
	H = Himp + Hhyb + Hbath


	# imaginary time
	N = 10
	δτ = 0.1
	β = N * δτ

	model = AndersonIM(U=0, μ=-w0)
	spec = DiracDelta(ω=w1, α=α)
	bath = fermionicbath(spec, β=β)

	lattice = GrassmannLattice(N=N, δτ=δτ, contour=:imag, bands=1)
	mpsK = sysdynamics(lattice, model, trunc=trunc)
	for band in 1:lattice.bands
		mpsK = boundarycondition!(mpsK, lattice, trunc=trunc)
	end

	corr = correlationfunction(bath, lattice)
	mpsI = hybriddynamics(lattice, corr, trunc=trunc)

	Z1 = tr(exp(-β .* H))
	Z2 = integrate(lattice, mpsK, mpsI)

	println("Z1=", Z2, " Z2=", Z2)

end	