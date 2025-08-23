println("------------------------------------")
println("|              Models              |")
println("------------------------------------")



@testset "Grassmann Ordering" begin
	N = 5
	δτ = 0.05
	ϵ_d = 1.25*pi
	dw = 0.2
	β = N * δτ
	τs = collect(0:δτ:β)
	rtol = 1.0e-5
	trunc = truncdimcutoff(D=300, ϵ=1.0e-6, add_back=0)

	for μ in (-5, 0, 5)
		bath = fermionicbath(spectrum_func(), β=β, μ=μ)
		exact_model = AndersonIM(μ=ϵ_d, U=1.3)

		Inrms = Float64[]
		Knrms = Float64[]
		for ordering in imag_grassmann_orderings
			lattice = GrassmannLattice(N=N, δτ=β/N, bands=2, contour=:imag, ordering=ordering)
			mpsI = vacuumstate(lattice)
			corr = correlationfunction(bath, lattice)
			for band in 1:lattice.bands
				mpsI = hybriddynamics!(mpsI, lattice, corr, trunc=trunc, band=band)
			end
			push!(Inrms, norm(mpsI))
			mpsK = sysdynamics(lattice, exact_model, trunc=trunc)
			push!(Knrms, norm(mpsK))
		end
		for i in 2:length(Inrms)
			@test abs((Inrms[i] - Inrms[1]) / Inrms[1]) < rtol
			@test abs((Knrms[i] - Knrms[1]) / Knrms[1]) < rtol
		end

		Inrms = Float64[]
		Knrms = Float64[]
		for ordering in real_grassmann_orderings
			lattice = GrassmannLattice(N=N, δt=β/N, bands=2, contour=:real, ordering=ordering)
			mpsI = vacuumstate(lattice)
			corr = correlationfunction(bath, lattice)
			for band in 1:lattice.bands
				mpsI = hybriddynamics!(mpsI, lattice, corr, trunc=trunc, band=band)
			end
			push!(Inrms, norm(mpsI))
			mpsK = sysdynamics(lattice, exact_model, trunc=trunc)
			push!(Knrms, norm(mpsK))
		end
		for i in 2:length(Inrms)
			@test abs((Inrms[i] - Inrms[1]) / Inrms[1]) < rtol
			@test abs((Knrms[i] - Knrms[1]) / Knrms[1]) < rtol
		end
	end
end

@testset "GF-imaginary time: benchmarking with ED, Analytic, TEMPO" begin
	N = 25
	δτ = 0.01
	ϵ_d = 1.25*pi
	dw = 0.2
	β = N * δτ
	τs = collect(0:δτ:β)
	rtol = 1.0e-2

	trunc = truncdimcutoff(D=300, ϵ=1.0e-6, add_back=0)

	for μ in (-5, 0, 5)
		# println("μ = ", μ)
		spec = spectrum_func()
		bath = fermionicbath(spec, β=β, μ=μ)
		b2 = discretebath(bath, δw=dw)
		model = Toulouse(b2, ϵ_d=ϵ_d)
		g₁ = toulouse_Gτ(model, τs)
		g₂ = [toulouse_Gτ(spec, τ, β = β, ϵ_d = ϵ_d, μ = μ) for τ in τs]
		@test norm(g₁ - g₂) / norm(g₁) < rtol

		exact_model = AndersonIM(μ=ϵ_d, U=0)
		for ordering in imag_grassmann_orderings
			lattice = GrassmannLattice(N=N, δτ=β/N, contour=:imag, ordering=ordering)
			mpsI = vacuumstate(lattice)
			corr = correlationfunction(bath, lattice)
			for band in 1:lattice.bands
				mpsI = hybriddynamics!(mpsI, lattice, corr, trunc=trunc, band=band)
			end
			mpsI = boundarycondition(mpsI, lattice)
			mpsK = sysdynamics(lattice, exact_model, trunc=trunc)
			g₃ = Gτ(lattice, mpsK, mpsI)
			@test norm(g₁ - g₃) / norm(g₁) < rtol

			Z = integrate(lattice, mpsK, mpsI)

			n1 = [occupation(lattice, i, mpsK, mpsI, Z=Z) for i in 1:lattice.k-1]
			n2 = [occupation2(lattice, i, mpsK, mpsI, Z=Z) for i in 1:lattice.k-1]
			@test norm(n1 - n2) / norm(n1) < rtol


		end

		# delta spectrum 
		spec = DiracDelta(ω=1, α=0.5)
		bath = fermionicbath(spec, β=β, μ=μ)
		g₂ = [toulouse_Gτ(spec, τ, β = β, ϵ_d = ϵ_d, μ = μ) for τ in τs]

		exact_model = AndersonIM(μ=ϵ_d, U=0)
		for ordering in imag_grassmann_orderings
			lattice = GrassmannLattice(N=N, δτ=β/N, contour=:imag, ordering=ordering)
			mpsI = vacuumstate(lattice)
			corr = correlationfunction(bath, lattice)
			for band in 1:lattice.bands
				mpsI = hybriddynamics!(mpsI, lattice, corr, trunc=trunc, band=band)
			end
			mpsI = boundarycondition(mpsI, lattice)
			mpsK = sysdynamics(lattice, exact_model, trunc=trunc)
			g₃ = Gτ(lattice, mpsK, mpsI)
			@test norm(g₂ - g₃) / norm(g₂) < rtol
		end
	end
end

@testset "GF-real time: benchmarking with Analytic" begin
	N = 10
	δt = 0.01
	t = N * δt
	ϵ_d = 1.25*pi
	β = 1.
	ts = [i*δt for i in 0:N]
	rtol = 5*1.0e-2

	trunc = truncdimcutoff(D=100, ϵ=1.0e-10, add_back=0)
		

	for spec in (spectrum_func(), DiracDelta())
		bath = fermionicbath(spec, β=β, μ=0.)
		gt = [im*toulouse_Gt(spectrum_func(), tj, ϵ_d = ϵ_d, μ = 0.) for tj in ts]

		exact_model = AndersonIM(μ=ϵ_d, U=0)
		corr = Δt(bath, N=N, t=t)
		for ordering in real_ac_grassmann_orderings
			lattice = GrassmannLattice(N=N, δt=δt, contour=:real, ordering=ordering)
			mpsI = hybriddynamics(lattice, corr, trunc=trunc) 
			mpsK = sysdynamics(lattice, exact_model, trunc=trunc)
			mpsK = boundarycondition(mpsK, lattice)
			mpsK = systhermalstate!(mpsK, lattice, exact_model, β=β)
			cache = environments(lattice, mpsK, mpsI)
			g1 = [cached_Gt(lattice, k, 1, mpsK, mpsI, c1=false, c2=true, b1=:+, b2=:+, cache=cache) for k in 1:lattice.k]
			g2 = [cached_Gt(lattice, 1, k, mpsK, mpsI, c1=true, c2=false, b1=:-, b2=:+, cache=cache) for k in 1:lattice.k]
			@test norm(gt - g1) / norm(gt) < rtol
			@test norm(g2) / length(g2) < rtol

			n1 = [cached_occupation(lattice, i, mpsK, mpsI, cache=cache) for i in 1:lattice.k-1]
			n2 = [occupation2(lattice, i, mpsK, mpsI, branch=:+, Z=Zvalue(cache)) for i in 1:lattice.k-1]
			@test norm(n1 - n2) / norm(n1) < rtol

			n1 = [real(cached_nn2(lattice, i, 2, mpsK, mpsI, b1=:+, b2=:-, cache=cache)) for i in 1:lattice.k-1]
			n2 = [real(nn(lattice, i, 2, mpsK, mpsI, b1=:+, b2=:+, Z=Zvalue(cache))) for i in 1:lattice.k-1]
			@test norm(n1 - n2) / norm(n1) < rtol
			n1 = [real(cached_nn(lattice, i, 2, mpsK, mpsI, b1=:+, b2=:+, cache=cache)) for i in 1:lattice.k-1]
			@test norm(n1 - n2) / norm(n1) < rtol
		end

	end
end

@testset "Bare impurity dynamics: benchmarking mixed-time and imaginary GTEMPO with Analytic" begin
	N = 25
	δτ = 0.05
	ϵ_d = 1.25*pi
	β = N * δτ
	τs = collect(0:δτ:β)
	δt = 0.1
	Nt = 20
	t = δt * Nt
	ts = [i*δt for i in 0:Nt]

	rtol = 1.0e-5

	trunc = truncdimcutoff(D=100, ϵ=1.0e-10, add_back=0)
	bath = fermionicbath(spectrum_func(), β=β, μ=0)
	gt1 = [freefermion_greater(tj, β=β, μ=ϵ_d) for tj in ts]
	gt2 = [freefermion_lesser(tj, β=β, μ=ϵ_d) for tj in ts]
	gτ = [freefermion_Gτ(τ, β=β, μ=ϵ_d) for τ in τs]

	exact_model = AndersonIM(μ=-ϵ_d, U=0)
	for ordering in mixed_ac_grassmann_orderings
		lattice = GrassmannLattice(Nt=Nt, δt=δt, Nτ=N, δτ=β/N, contour=:mixed, ordering=ordering)
		mpsK = sysdynamics(lattice, exact_model, trunc=trunc)
		mpsK = boundarycondition(mpsK, lattice)
		cache = environments(lattice, mpsK)
		g1 = [-im*cached_Gm(lattice, k, 1, mpsK, c1=false, c2=true, b1=:+, b2=:+, cache=cache) for k in 1:lattice.kt]
		g2 = [im*cached_Gm(lattice, 1, k, mpsK, c1=true, c2=false, b1=:-, b2=:+, cache=cache) for k in 1:lattice.kt]
		g3 = [cached_Gm(lattice, k, 1, mpsK, c1=false, c2=true, b1=:τ, b2=:τ, cache=cache) for k in 1:lattice.kτ]
		@test norm(gt1 - g1) / norm(gt1) < rtol
		@test norm(gt2 - g2) / norm(gt2) < rtol
		@test norm(gτ - g3) / norm(gτ) < rtol
	end

	for ordering in imag_ac_grassmann_orderings
		lattice = GrassmannLattice(N=N, δτ=δτ, contour=:imag, ordering=ordering)
		mpsK = sysdynamics(lattice, exact_model, trunc=trunc)
		mpsK = boundarycondition(mpsK, lattice)
		cache = environments(lattice, mpsK)
		g3 = [cached_Gτ(lattice, k, 1, mpsK, c1=false, c2=true, cache=cache) for k in 1:lattice.k]
		@test norm(gτ - g3) / norm(gτ) < rtol		
	end
end

@testset "GF-mixed time: benchmarking ED" begin

	N = 25
	δτ = 0.01
	ϵ_d = 1.25*pi
	dw = 0.1
	β = N * δτ
	τs = collect(0:δτ:β)
	δt = 0.01
	Nt = 20
	t = δt * Nt
	ts = [i*δt for i in 0:Nt]

	rtol = 1.0e-2

	trunc = truncdimcutoff(D=100, ϵ=1.0e-10, add_back=0)
		
	bath = fermionicbath(spectrum_func(), β=β, μ=0)
	b2 = discretebath(bath, δw=dw)
	gτ = toulouse_Gτ(Toulouse(b2, ϵ_d=ϵ_d), τs)
	gt1, gt2 = toulouse_greater_lesser(Toulouse(b2, ϵ_d=ϵ_d), ts)
	gt1, gt2 = im*gt1, -im*gt2

	exact_model = AndersonIM(μ=ϵ_d, U=0)
	for ordering in mixed_ac_grassmann_orderings
		lattice = GrassmannLattice(Nt=Nt, δt=δt, Nτ=N, δτ=β/N, contour=:mixed, ordering=ordering)
		corr = correlationfunction(bath, lattice)
		mpsI = hybriddynamics(lattice, corr, trunc=trunc) 
		mpsK = sysdynamics(lattice, exact_model, trunc=trunc)
		mpsK = boundarycondition(mpsK, lattice)
		cache = environments(lattice, mpsK, mpsI)
		g1 = [cached_greater(lattice, k, mpsK, mpsI, cache=cache) for k in 1:lattice.kt]
		g2 = [cached_lesser(lattice, k, mpsK, mpsI, cache=cache) for k in 1:lattice.kt]
		g3 = [cached_Gτ(lattice, k, mpsK, mpsI, cache=cache) for k in 1:lattice.kτ]
		@test norm(gt1 - g1) / norm(gt1) < rtol
		@test norm(gt2 - g2) / norm(gt2) < rtol
		@test norm(gτ - g3) / norm(gτ) < rtol

		n1 = [real(cached_nn2(lattice, i, 2, mpsK, mpsI, b1=:+, b2=:-, cache=cache)) for i in 1:lattice.kt-1]
		n2 = [real(nn2(lattice, i, 2, mpsK, mpsI, b1=:+, b2=:+, Z=Zvalue(cache))) for i in 1:lattice.kt-1]
		@test norm(n1 - n2) / norm(n1) < rtol
		n1 = [real(cached_nn(lattice, i, 2, mpsK, mpsI, b1=:+, b2=:+, cache=cache)) for i in 1:lattice.kt-1]
		@test norm(n1 - n2) / norm(n1) < rtol
	end

end

@testset "GF-real time, thermal initial state" begin

	function _test_AndersonIM(U, ϵ_d)
		trunc = truncdimcutoff(D=200, ϵ=1.0e-8, add_back=0)
		δτ = 0.05
		β = 1.
		N = round(Int, β / δτ)

		bath = fermionicbath(spectrum_func(), β=β, μ=0)
		exact_model = AndersonIM(U=U, μ=ϵ_d)

		bands = (U == 0.) ? 1 : 2
		lattice = GrassmannLattice(δτ=δτ, N=N, bands=bands, contour=:imag)

		mpsKs = [sysdynamics(lattice, exact_model)]
		for band in 1:lattice.bands
			mpsKs = boundarycondition_branching(mpsKs, lattice, band=band)
		end
		

		for ordering in real_grassmann_orderings

			lattice_r = GrassmannLattice(δt=0.1, N=5, bands=bands, contour=:real)
			exact_model = AndersonIM(U=U, μ=ϵ_d)
			mps = sysdynamics(lattice_r, exact_model)
			mps = systhermalstate!(mps, lattice_r, exact_model, β= β)
			for band in 1:lattice.bands
				mps = boundarycondition(mps, lattice_r, band=band)
			end

			for band in 1:lattice.bands
				n1 = 1-Gτ(lattice, 1, mpsKs, band=band)
				for i in 1:lattice_r.N
					n2 = occupation(lattice_r, i, mps, band=band)
					@test abs((n2-n1)/n1) < 1.0e-2
				end
			end

		end
	end

	function _test_AndersonIM_2(U, ϵ_d)
		trunc = truncdimcutoff(D=200, ϵ=1.0e-8, add_back=0)
		δτ = 0.05
		β = 10
		N = round(Int, β / δτ)

		bath = fermionicbath(spectrum_func(), β=β, μ=0)
		exact_model = AndersonIM(U=U, μ=ϵ_d)

		bands = (U == 0.) ? 1 : 2
		lattice = GrassmannLattice(δτ=δτ, N=N, bands=bands, contour=:imag)

		mpsKs = [accsysdynamics_fast(lattice, exact_model)]
		for band in 1:lattice.bands
			mpsKs = boundarycondition_branching(mpsKs, lattice, band=band)
		end
		
		for ordering in [A1Ā1B1B̄1a1ā1b1b̄1(), A1Ā1a1ā1B1B̄1b1b̄1(), A1Ā1B1B̄1b̄1B̄1ā1Ā1(), A2Ā2A1Ā1a2ā2a1ā1B2B̄2B1B̄1b2b̄2b1b̄1(), A2Ā2B2B̄2A1Ā1B1B̄1a1ā1b1b̄1a2ā2b2b̄2()]

			lattice_r = GrassmannLattice(δt=0.1, N=5, bands=bands, contour=:real)
			exact_model = AndersonIM(U=U, μ=ϵ_d)
			mps = accsysdynamics_fast(lattice_r, exact_model)
			mps = systhermalstate!(mps, lattice_r, exact_model, β= β)
			for band in 1:lattice.bands
				mps = boundarycondition(mps, lattice_r, band=band)
			end

			for band in 1:lattice.bands
				n1 = 1-Gτ(lattice, 1, mpsKs, band=band)
				for i in 1:lattice_r.N
					n2 = occupation(lattice_r, i, mps, band=band)
					@test abs((n2-n1)/n1) < 1.0e-2
				end
			end

		end		
	end

	function _test_sk(ϵ_d)
		U = 2.
		J = 0.45
		ϵ_d = 0.7
		# ϵ_d = -(3*U - 5*J)/2
		norb=2
		δτ = 0.05
		β = 10
		D = 1.
		N = round(Int, β / δτ)

		bath = fermionicbath(spectrum_func(D), β=β, μ=0)
		exact_model = KanamoriIM(U=U, μ=ϵ_d, J=J, norb=norb)

		bands = 2 * norb
		lattice = GrassmannLattice(δτ=δτ, N=N, bands=bands, contour=:imag)

		mpsKs = [accsysdynamics_fast(lattice, exact_model)]
		for band in 1:lattice.bands
			mpsKs = boundarycondition_branching(mpsKs, lattice, band=band)
		end


		lattice_r = GrassmannLattice(δt=0.1, N=2, bands=bands, contour=:real, ordering = A2Ā2B2B̄2A1Ā1B1B̄1a1ā1b1b̄1a2ā2b2b̄2())
		mps = accsysdynamics_fast(lattice_r, exact_model)
		mps = systhermalstate!(mps, lattice_r, exact_model, β= β)
		for band in 1:lattice.bands
			mps = boundarycondition(mps, lattice_r, band=band)
		end

		for band in 1:lattice.bands
			n1 = 1-Gτ(lattice, 1, mpsKs, band=band)
			for i in 1:lattice_r.N
				n2 = occupation(lattice_r, i, mps, band=band)
				@test abs((n2-n1)/n1) < 1.0e-2
			end
		end

	end

	_test_AndersonIM(0, 1)
	_test_AndersonIM(0, -1)
	_test_AndersonIM(1, 0.5)
	_test_AndersonIM(1, -0.5)

	_test_AndersonIM_2(0, 1)
	_test_AndersonIM_2(0, -1)
	_test_AndersonIM_2(1, 0.7)
	_test_AndersonIM_2(1, -0.7)

	_test_sk(0.7)
	_test_sk(-1.875)
end

@testset "GF-real time, occupation and particle current: benchmarking ED, TEMPO" begin
	β = 1.0
	D = 1.0
	ϵ_d = -1.0
	N = 10
	δt = 0.02
	t = N*δt	
	dw = 0.02

	rtol = 5*1.0e-2

	trunc = truncdimcutoff(D=200, ϵ=1.0e-9, add_back=0)

	for μ in (-2, 0, 2)
		# println("μ = ", μ)
		bath = fermionicbath(spectrum_func(D), β=β, μ=μ)
		# ED
		b2 = discretebath(bath, δw=dw)
		model = Toulouse(b2, ϵ_d=ϵ_d)
		ρ₀ = separablecdm(model, 0)
		h = cmatrix(model)
		cache = eigencache(transpose(h))
		observer = particlecurrent_cmatrix(model)
		observer2 = heatcurrent_cmatrix(model)
		currents = ComplexF64[]
		heatcurrents = ComplexF64[]
		ns = Float64[]
		for i in 1:N
			ρ = ImpurityModelBase.timeevo(ρ₀, h, -im*i*δt, cache)
			push!(ns, real(ρ[1, 1]))
			push!(currents, sum(observer .* ρ))
			push!(heatcurrents, sum(observer2 .* ρ))
		end
		heatcurrents = real(heatcurrents)

		# TEMPO 1 order
		exact_model = AndersonIM(μ=ϵ_d, U=0)
		for ordering in real_grassmann_orderings
			lattice = GrassmannLattice(N=N, δt=δt, contour=:real, order=1, ordering=ordering)

			corr = correlationfunction(bath, lattice)
			mpsI = hybriddynamics(lattice, corr, trunc=trunc) 
			mpsI = boundarycondition(mpsI, lattice)
			mpsK = sysdynamics(lattice, exact_model, trunc=trunc)

			ns2 = occupation(lattice, mpsK, mpsI)
			currents2 = electriccurrent(lattice, corr, mpsK, mpsI)
			@test norm(ns-ns2) / norm(ns) < rtol
			@test norm(currents - currents2) / norm(currents) < rtol

			heatcurrents2 = real(heatcurrent_fast(lattice, bath, mpsK, mpsI))
			@test norm(heatcurrents - heatcurrents2) / norm(heatcurrents) < rtol
		end
		# TEMPO 1 order 2
		for ordering in (A1Ā1B1B̄1a1ā1b1b̄1(), A1Ā1a1ā1B1B̄1b1b̄1(), A1Ā1B1B̄1b̄1B̄1ā1Ā1())
			lattice_o = GrassmannLattice(N=N, δt=δt, contour=:real, order=1, bands=1)
			corr = correlationfunction(bath, lattice_o)
			lattice = similar(lattice_o, N=0)
			ns2 = zeros(Float64, N)
			currents2 = zeros(ComplexF64, N)
			mpsI = vacuumstate(lattice)
			mpsK = vacuumstate(lattice)
			for k in 2:N+1
				lattice, mpsI, mpsK = makestep(lattice, mpsI, mpsK)
				@test k == timesteps(mpsI, lattice) == lattice.k
				mpsI = hybriddynamicsstepper(mpsI, lattice, corr, trunc=trunc)
				mpsI′ = boundarycondition(mpsI, lattice)
				mpsK = sysdynamicsstepper!(mpsK, lattice, exact_model, trunc=trunc)
				ns2[k-1] = occupation(lattice, k-1, mpsK, mpsI′)
				currents2[k-1] = electriccurrent_fast(lattice, corr, k, mpsK, mpsI′)
			end
			@test norm(ns-ns2) / norm(ns) < rtol
			@test norm(currents - currents2) / norm(currents) < rtol
		end
		# TEMPO 2 order
		for ordering in (A1Ā1B1B̄1a1ā1b1b̄1(), A1Ā1a1ā1B1B̄1b1b̄1(), A1Ā1B1B̄1b̄1B̄1ā1Ā1())
			lattice_o = GrassmannLattice(N=N, δt=δt, contour=:real, order=2, bands=1)
			corr = correlationfunction(bath, lattice_o)
			lattice = similar(lattice_o, N=0)
			ns2 = zeros(Float64, N)
			currents2 = zeros(ComplexF64, N)
			mpsI = vacuumstate(lattice)
			mpsK = vacuumstate(lattice)
			for k in 2:N+1
				lattice, mpsI, mpsK = makestep(lattice, mpsI, mpsK)
				@test k == timesteps(mpsI, lattice) == lattice.k
				mpsI2 = hybriddynamicsstepper(mpsI, lattice, corr, finalize=true, trunc=trunc)
				mpsI2 = boundarycondition(mpsI2, lattice)
				mpsI = hybriddynamicsstepper(mpsI, lattice, corr, trunc=trunc)
				mpsK = sysdynamicsstepper!(mpsK, lattice, exact_model, trunc=trunc)
				cache = environments(lattice, mpsK, mpsI2)
				ns2[k-1] = cached_occupation(lattice, mpsK, mpsI2, cache=cache)
				currents2[k-1] = cached_electriccurrent_fast(lattice, corr, mpsK, mpsI2, cache=cache)
			end
			@test norm(ns-ns2) / norm(ns) < rtol
			@test norm(currents - currents2) / norm(currents) < rtol
		end
	end
end

@testset "Benchmarking interacting systems with different GrassmannOrdering" begin
	epsilon_d = -1.25  * pi
	leftmu = -5
	rightmu = 3
	leftbeta = 1
	rightbeta = 0.5
	U = 3.5*pi
	# U = 0
	Dl = 10.
	Dr = 10.
	rtol = 5*1.0e-2

	δt = 0.02
	N = 5
	t = N*δt
	ts = [i*δt for i in 0:N] 

	leftbath = fermionicbath(spectrum_func(Dl), β=leftbeta, μ=leftmu)
	rightbath = fermionicbath(spectrum_func(Dr), β=rightbeta, μ=rightmu)
	trunc = truncdimcutoff(D=200, ϵ=1.0e-7, add_back=0)

	exact_model = AndersonIM(μ=epsilon_d, U=U)
	# TEMPO, order ABB̄Ā
	lattice = GrassmannLattice(N=N, δt=δt, contour=:real, bands=2, order=1, ordering=A1Ā1B1B̄1b̄1B̄1ā1Ā1())
	leftcorr = correlationfunction(leftbath, lattice)
	rightcorr = correlationfunction(rightbath, lattice)
	corr = leftcorr + rightcorr
	mpsI = vacuumstate(lattice)
	for band in 1:lattice.bands
		mpsI = hybriddynamics!(mpsI, lattice, corr, trunc=trunc, band=band) 
	end
	mpsI = boundarycondition(boundarycondition(mpsI, lattice, band=1), lattice, band=2)
	mpsK = sysdynamics(lattice, exact_model, trunc=trunc)

	ns = occupation(lattice, mpsK, mpsI)
	currents_left = [electriccurrent(lattice, leftcorr, k+1, mpsK, mpsI) for k in 1:N]
	currents_right = [electriccurrent(lattice, rightcorr, k+1, mpsK, mpsI) for k in 1:N]
	# TEMPO, order AĀBB̄
	lattice = GrassmannLattice(N=N, δt=δt, contour=:real, bands=2, order=1, ordering=A1Ā1a1ā1B1B̄1b1b̄1())
	leftcorr = correlationfunction(leftbath, lattice)
	rightcorr = correlationfunction(rightbath, lattice)
	corr = leftcorr + rightcorr
	mpsI1 = hybriddynamics(lattice, corr, trunc=trunc, band=1) 
	mpsI1 = boundarycondition(mpsI1, lattice, band=1)
	mpsI2 = hybriddynamics(lattice, corr, trunc=trunc, band=2) 
	mpsI2 = boundarycondition(mpsI2, lattice, band=2)
	mpsK = sysdynamics(lattice, exact_model, trunc=trunc)

	cache = environments(lattice, mpsK, mpsI1, mpsI2)
	ns2 = cached_occupation(lattice, mpsK, mpsI1, mpsI2, cache=cache)
	@test norm(ns2-ns) / norm(ns) < rtol
	currents_left2 = [cached_electriccurrent_fast(lattice, leftcorr, k+1, mpsK, mpsI1, mpsI2, cache=cache) for k in 1:N]
	@test norm(currents_left2-currents_left) / norm(currents_left) < rtol
	currents_right2 = [cached_electriccurrent_fast(lattice, rightcorr, k+1, mpsK, mpsI1, mpsI2, cache=cache) for k in 1:N]
	@test norm(currents_right2-currents_right) / norm(currents_right) < rtol

end

