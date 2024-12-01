println("------------------------------------")
println("|        Independent Bosons        |")
println("------------------------------------")

using LinearAlgebra: Diagonal, eigen, Hermitian

function spin_half_matrices()
	s_SP = Array{Float64, 2}([0 0; 1 0])
	s_SM = Array{Float64, 2}([0 1; 0 0])
	s_Z = Array{Float64, 2}([-1 0; 0 1])
	s_x = s_SP+s_SM
	s_y = -im*(s_SP-s_SM)
	n = Array{Float64, 2}([0 0; 0 1])
	return Dict("x"=>s_x, "y"=>s_y, "z"=>s_Z, "+"=>s_SP, "-"=>s_SM, "n"=>n)
end

function Aop(d::Int)
	(d <= 1) && error("d must be larger than 1.")
	a = zeros(Float64, d, d)
	for i = 1:(d - 1)
		a[i, i+1] = sqrt(i)
	end
	return a
end

ADAGop(d::Int) = Array(transpose(Aop(d)))

Nop(d::Int) = ADAGop(d) * Aop(d)


function boson(;d::Int=5)
	_N = Nop(d)
	_N2 = _N * _N
	return Dict("a"=>Aop(d),"adag"=>ADAGop(d), "n"=>_N, "n2"=>_N2)
end

# <e^τH A e^-τH B>
function gf_imag(H, A, B, β::Real, n::Int)
	δτ = β / n
	τs = 0:δτ:β
	λs, U = eigen(Hermitian(H))
	ρ = U * Diagonal(exp.(-β .* λs)) * U'
	tr_ρ = tr(ρ)
	g(τ) = tr(U * Diagonal(exp.(τ .* λs)) * U' * A * U * Diagonal(exp.(-τ .* λs)) * U' * B * ρ) / tr_ρ
	return g.(τs)
end

# <e^iHt A e^-iHt B>
function gf_real(H, A, B, β::Real, t::Real, n::Int, ρ=exp(-β*H))
	δt = t / n
	ts = 0:δt:t
	λs, U = eigen(Hermitian(H))
	# ρ = U * Diagonal(exp.(-β .* λs)) * U'
	tr_ρ = tr(ρ)
	gt(tj) = -im*tr(U * Diagonal(exp.(im*tj .* λs)) * U' * A * U * Diagonal(exp.(-im*tj .* λs)) * U' * B * ρ) / tr_ρ
	ls(tj) = im*tr(B * U * Diagonal(exp.(im*tj .* λs)) * U' * A * U * Diagonal(exp.(-im*tj .* λs)) * U' * ρ) / tr_ρ
	return gt.(ts), ls.(ts)
end



# ϵ_d n̂ + α n̂(b̂ + b̂†) + ω₀b̂†b̂ 
function noninteracting_operators(ϵ_d; ω₀=1, α=0.5, d=100)
	p1 = spin_half_matrices()
	n̂, σ₊, σ₋ = p1["n"], p1["+"], p1["-"]
	Is = one(n̂)
	p2 = boson(d=d)
	b̂, b̂′, n̂b = p2["a"], p2["adag"], p2["n"]
	Ib = one(b̂)
	# total Hamiltonian
	Himp = -ϵ_d * kron(n̂, Ib)
	Hbath = ω₀ * kron(Is, n̂b)
	Hhyb = sqrt(α) * kron(n̂, b̂′ + b̂)
	H = Himp + Hhyb + Hbath
	A, B = kron(σ₋, Ib), kron(σ₊, Ib)

	return H, A, B, Himp + Hbath
end

# ϵ_d(n̂↑ + n̂↓) + U n̂↑n̂↓ + α (n̂↑ + n̂↓)(b̂ + b̂†) + ω₀b̂†b̂ 
function interacting_operators(U, ϵ_d=U/2; ω₀=1, α=0.5, d=100)
	p1 = spin_half_matrices()
	n̂, σ₊, σ₋ = p1["n"], p1["+"], p1["-"]
	Is = one(n̂)
	n_ud = kron(n̂, Is) + kron(Is, n̂)
	nn = kron(n̂,n̂)
	p2 = boson(d=d)
	b̂, b̂′, n̂b = p2["a"], p2["adag"], p2["n"]
	Ib = one(b̂)
	# total Hamiltonian
	Himp = kron(-ϵ_d*n_ud + U * nn, Ib)
	Hbath = ω₀ * kron(kron(Is, Is), n̂b)
	Hhyb = sqrt(α) * kron(n_ud, b̂′ + b̂)
	H =  Himp + Hhyb + Hbath

	A, B = kron(kron(σ₋, Is), Ib), kron(kron(σ₊, Is), Ib)

	return H, A, B, Himp + Hbath
end


function noninteracting_real(ϵ_d)
	β = 0.1
	δt=0.01
	N = 10
	t = N * δt

	H, a, adag, H0 = noninteracting_operators(ϵ_d, ω₀=1, α=0.5, d=100)
	ρ = exp(-β*H0)
	g1, g2 = gf_real(H, a, adag, β, t, N, ρ)

	return g1, g2
end

function interacting_real(U, ϵ_d=U/2)
	β = 0.1
	δt=0.01
	N = 10
	t = N * δt

	H, a, adag, H0 = interacting_operators(U, ϵ_d, ω₀=1, α=0.5, d=100)
	ρ = exp(-β*H0)
	g1, g2 = gf_real(H, a, adag, β, t, N, ρ)

	return g1, g2
end

@testset "Independent bosons: imaginary time" begin
	rtol = 1.0e-2
	rtol2 = 1.0e-5
	δτ=0.01
	N = 10
	β = N * δτ
	chi = 100

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)
	truncK = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)

	lattice = GrassmannLattice(N=N, δτ=δτ, contour=:imag, order=1)

	for ϵ_d in (-0.5, 0, 0.5)
		for spec in (Leggett(d=3, ωc=1), DiracDelta(ω₀=1, α=0.5))

			bath = bosonicbath(spec, β=β)
			corr = correlationfunction(bath, lattice)

			mpsI = retardedinteractdynamics(lattice, corr, trunc=trunc)
			mpsI′ = retardedinteractdynamics_naive(lattice, corr, trunc=trunc)
			@test distance(mpsI, mpsI′) / norm(mpsI) <= rtol2

			fbath = fermionicbath(semicircular(), β=β, μ=0)
			exact_model = SISB(fbath, U=0., μ=-ϵ_d)
			mpsK = sysdynamics(lattice, exact_model, trunc=truncK)
			mpsKs = [mpsK]
			for band in 1:lattice.bands
				mpsKs = boundarycondition_branching(mpsKs, lattice, band=band)
			end
			cache = environments(lattice, mpsKs, mpsI)
			g1 = cached_Gτ(lattice, mpsKs, mpsI, cache=cache)
			g2 = independentbosons_Gτ(spec, β=β, ϵ_d=-ϵ_d, N=N)

			@test norm(g1 - g2) / norm(g1) < rtol	
		end	
	end


	ϵ_d = 0.7
	lattice = GrassmannLattice(N=N, δτ=δτ, contour=:imag, order=1, bands=2)
	for U in (0, 1)
		for spec in (Leggett(d=3, ωc=1), DiracDelta(ω₀=1, α=0.5))

			bath = bosonicbath(spec, β=β)
			corr = correlationfunction(bath, lattice)

			mpsI = retardedinteractdynamics(lattice, corr, trunc=trunc)
			mpsI′ = retardedinteractdynamics_naive(lattice, corr, trunc=trunc)
			@test distance(mpsI, mpsI′) / norm(mpsI) <= rtol2

			fbath = fermionicbath(semicircular(), β=β, μ=0)
			exact_model = SISB(fbath, U=U, μ=-ϵ_d)
			mpsK = sysdynamics(lattice, exact_model, trunc=truncK)
			mpsKs = [mpsK]
			for band in 1:lattice.bands
				mpsKs = boundarycondition_branching(mpsKs, lattice, band=band)
			end
			cache = environments(lattice, mpsKs, mpsI)
			g1 = cached_Gτ(lattice, mpsKs, mpsI, cache=cache)
			g2 = independentbosons_Gτ(spec, β=β, ϵ_d=-ϵ_d, N=N, U=U, bands=2)

			@test norm(g1 - g2) / norm(g1) < rtol	
		end	
	end
end

@testset "Independent bosons: real time" begin
	rtol = 1.0e-2
	rtol2 = 1.0e-3
	β = 0.1
	δt=0.01
	Nt = 10
	t = Nt * δt
	chi = 60

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)
	truncK = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)

	lattice = GrassmannLattice(N=Nt, δt=δt, contour=:real, order=1)
	for ϵ_d in (-0.5, 0, 0.5)

		bath = bosonicbath(DiracDelta(ω₀=1, α=0.5), β=β)
		corr = correlationfunction(bath, lattice)
		mpsI = retardedinteractdynamics(lattice, corr, trunc=trunc)
		mpsI′ = retardedinteractdynamics_naive(lattice, corr, trunc=trunc)
		@test distance(mpsI, mpsI′) / norm(mpsI) <= rtol2

		fbath = fermionicbath(semicircular(), β=β, μ=0)
		exact_model = SISB(fbath, U=0., μ=-ϵ_d)
		mpsK = sysdynamics(lattice, exact_model, trunc=truncK)
		for band in 1:lattice.bands
			mpsK = boundarycondition!(mpsK, lattice, band=band)
		end
		mpsK = systhermalstate!(mpsK, lattice, exact_model, trunc=trunc, δτ=0.001)
		cache = environments(lattice, mpsK, mpsI)

		g1 = [-im*cached_greater(lattice, k, mpsK, mpsI, c1=false, c2=true, b1=:+, b2=:+, band=1, cache=cache) for k in 1:Nt+1]
		g2 = [-im*cached_lesser(lattice, k, mpsK, mpsI, c1=true, c2=false, b1=:-, b2=:+, band=1, cache=cache) for k in 1:Nt+1]

		g1′, g2′ = noninteracting_real(ϵ_d)
		@test norm(g1 - g1′) / norm(g1) < rtol	
		@test norm(g2 - g2′) / norm(g2) < rtol	
	end

	lattice = GrassmannLattice(N=Nt, δt=δt, contour=:real, order=1, bands=2)
	ϵ_d = 0.3
	for U in (0, 1)

		bath = bosonicbath(DiracDelta(ω₀=1, α=0.5), β=β)
		corr = correlationfunction(bath, lattice)
		mpsI = retardedinteractdynamics(lattice, corr, trunc=trunc)
		mpsI′ = retardedinteractdynamics_naive(lattice, corr, trunc=trunc)
		@test distance(mpsI, mpsI′) / norm(mpsI) <= rtol2

		fbath = fermionicbath(semicircular(), β=β, μ=0)
		exact_model = SISB(fbath, U=U, μ=-ϵ_d)
		mpsK = sysdynamics(lattice, exact_model, trunc=truncK)
		for band in 1:lattice.bands
			mpsK = boundarycondition!(mpsK, lattice, band=band)
		end
		mpsK = systhermalstate!(mpsK, lattice, exact_model, trunc=trunc, δτ=0.001)
		cache = environments(lattice, mpsK, mpsI)

		g1 = [-im*cached_greater(lattice, k, mpsK, mpsI, c1=false, c2=true, b1=:+, b2=:+, band=1, cache=cache) for k in 1:Nt+1]
		g2 = [-im*cached_lesser(lattice, k, mpsK, mpsI, c1=true, c2=false, b1=:-, b2=:+, band=1, cache=cache) for k in 1:Nt+1]

		g1′, g2′ = interacting_real(U, ϵ_d)
		@test norm(g1 - g1′) / norm(g1) < rtol	
		@test norm(g2 - g2′) / norm(g2) < rtol	

	end

end

@testset "Independent bosons: mixed time" begin
	rtol = 5.0e-2
	rtol2 = 1.0e-3
	β = 0.05
	δτ = 0.01
	Nτ = round(Int, β/δτ)
	δt=0.01
	Nt = 5
	t = Nt * δt
	chi = 80

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)
	truncK = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)


	# noninteracting case
	lattice = GrassmannLattice(Nt=Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed, order=1)

	ϵ_d = 0.5
	spec = Leggett(d=3, ωc=1)
	bath = bosonicbath(spec, β=β)
	corr = correlationfunction(bath, lattice)

	mpsI = retardedinteractdynamics(lattice, corr, trunc=trunc)
	mpsI′ = retardedinteractdynamics_naive(lattice, corr, trunc=trunc)
	@test distance(mpsI, mpsI′) / norm(mpsI) <= rtol2

	fbath = fermionicbath(semicircular(), β=β, μ=0)
	exact_model = SISB(fbath, U=0., μ=-ϵ_d)
	mpsK = sysdynamics(lattice, exact_model, trunc=truncK)
	for band in 1:lattice.bands
		mpsK = boundarycondition!(mpsK, lattice, band=band)
	end
	cache = environments(lattice, mpsK, mpsI)

	g1 = [-im*cached_Gm(lattice, k, 1, mpsK, mpsI, c1=false, c2=true, b1=:+, b2=:+, band=1, cache=cache) for k in 1:Nt+1]
	g2 = [im*cached_Gm(lattice, 1, k, mpsK, mpsI, c1=true, c2=false, b1=:-, b2=:+, band=1, cache=cache) for k in 1:Nt+1]
	g3 = [cached_Gm(lattice, k, 1, mpsK, mpsI, c1=false, c2=true, b1=:τ, b2=:τ, band=1, cache=cache) for k in 1:Nτ+1]

	g1′ = [independentbosons_greater(spec, tj, β=β, ϵ_d=-ϵ_d) for tj in 0:δt:t]
	g2′ = [independentbosons_lesser(spec, tj, β=β, ϵ_d=-ϵ_d) for tj in 0:δt:t]
	g3′ = independentbosons_Gτ(spec, β=β, ϵ_d=-ϵ_d, N=Nτ)

	@test norm(g1 - g1′) / norm(g1) < rtol	
	@test norm(g2 - g2′) / norm(g2) < rtol	
	@test norm(g3 - g3′) / norm(g2) < rtol	

	# interacting case
	U = -0.8
	lattice = GrassmannLattice(Nt=Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed, order=1, bands=2)

	mpsI = retardedinteractdynamics(lattice, corr, trunc=trunc)
	mpsI′ = retardedinteractdynamics_naive(lattice, corr, trunc=trunc)
	@test distance(mpsI, mpsI′) / norm(mpsI) <= rtol2
	
	fbath = fermionicbath(semicircular(), β=β, μ=0)
	exact_model = SISB(fbath, U=U, μ=-ϵ_d)
	mpsK = sysdynamics(lattice, exact_model, trunc=truncK)
	for band in 1:lattice.bands
		mpsK = boundarycondition!(mpsK, lattice, band=band)
	end
	cache = environments(lattice, mpsK, mpsI)

	g1 = [-im*cached_Gm(lattice, k, 1, mpsK, mpsI, c1=false, c2=true, b1=:+, b2=:+, band=1, cache=cache) for k in 1:Nt+1]
	g2 = [im*cached_Gm(lattice, 1, k, mpsK, mpsI, c1=true, c2=false, b1=:-, b2=:+, band=1, cache=cache) for k in 1:Nt+1]
	g3 = [cached_Gm(lattice, k, 1, mpsK, mpsI, c1=false, c2=true, b1=:τ, b2=:τ, band=1, cache=cache) for k in 1:Nτ+1]

	g1′ = [independentbosons_greater(spec, tj, β=β, ϵ_d=-ϵ_d, U=U, bands=2) for tj in 0:δt:t]
	g2′ = [independentbosons_lesser(spec, tj, β=β, ϵ_d=-ϵ_d, U=U, bands=2) for tj in 0:δt:t]
	g3′ = independentbosons_Gτ(spec, β=β, ϵ_d=-ϵ_d, N=Nτ, U=U, bands=2)

	@test norm(g1 - g1′) / norm(g1) < rtol	
	@test norm(g2 - g2′) / norm(g2) < rtol	
	@test norm(g3 - g3′) / norm(g2) < rtol	

end