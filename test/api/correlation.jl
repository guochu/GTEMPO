@testset "API: correlationfunction" begin
	β = 1.0

	# imaginary-time lattice
	lat = GrassmannLattice(N=4, δτ=0.25, contour=:imag)
	corr = correlationfunction(fermionicbath(spectrum_func(), β=β, μ=0), lat)
	@test corr isa ImagCorrelationFunction
	@test size(corr.data) == (lat.N, lat.N)

	# a bosonic bath works on the Grassmann lattice as well (retarded channel)
	corr = correlationfunction(bosonicbath(DiracDelta(ω=1, α=0.5), β=β), lat)
	@test corr isa ImagCorrelationFunction

	# real-time lattice
	lat = GrassmannLattice(N=4, δt=0.1, contour=:real)
	corr = correlationfunction(fermionicbath(spectrum_func(), β=β, μ=0), lat)
	@test corr isa RealCorrelationFunction
	@test size(corr.G₊₊) == (lat.k, lat.k)   # branch block (+,+)
	@test size(corr.G₋₊) == (lat.k, lat.k)   # branch block (-,+)

	# mixed (Kadanoff) lattice: index(corr, i, j; b1, b2) over three branches
	lat = GrassmannLattice(Nt=3, δt=0.1, Nτ=4, δτ=0.25, contour=:mixed)
	corr = correlationfunction(fermionicbath(spectrum_func(), β=β, μ=0), lat)
	@test corr isa MixedCorrelationFunction
	for b1 in (:+, :-, :τ), b2 in (:+, :-, :τ)
		@test index(corr, 1, 1, b1=b1, b2=b2) isa Number
	end
end

@testset "API: baths" begin
	β = 1.0
	# fermionic bath with a continuous spectrum
	bath = fermionicbath(spectrum_func(), β=β, μ=0.3)
	@test bath.β == β && bath.μ == 0.3
	# delta-spectrum bath
	bath = fermionicbath(DiracDelta(ω=1, α=0.5), β=β)
	@test bath.β == β
	# bosonic bath
	bath = bosonicbath(Leggett(d=3, ωc=1), β=β)
	@test bath.β == β
	# bcs bath wraps a normal fermionic bath with a pairing gap
	bath2 = bcsbath(fermionicbath(spectrum_func(), β=β, μ=0), Δ=0.3)
	@test bath2.Δ == 0.3 && bath2.β == β
	# discretization of a continuous bath
	disbath = discretebath(fermionicbath(spectrum_func(), β=β, μ=0), δw=0.2)
	@test num_sites(disbath) > 0
end
