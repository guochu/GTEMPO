using LinearAlgebra, JSON
using LinearAlgebra: ishermitian, Hermitian, eigen

struct EigenCache{T<:Number, R<:Real}
	h::Matrix{T}
	U::Matrix{T}
	λs::Vector{R}
end

function EigenCache(h::AbstractMatrix)
	@assert ishermitian(h)
	# h2 = convert(Matrix{ComplexF64}, h)
	h2 = convert(Matrix{eltype(h)}, h)
	λs, U = eigen(Hermitian(h2))
	return EigenCache(h2, U, λs)
end


# initializers
# do transpose here

"""
	freefermions_cache(h::AbstractMatrix)

h is the coefficient matrix, namely, 
h = [h₁₁ c₁†c₁, h₁₂ c₁†c₂, h₁₃ c₁†c₃...; h₂₁ c₂†c₁, h₂₂ c₂†c₂, h₂₃ c₂†c₃; ...]

time evolution for the coefficient matrix h of free fermions is dρ/dt = -i [h^t, ρ],
where ρ is the quadratic observables
"""
eigencache(h::AbstractMatrix) = EigenCache(h)


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
function gf_imag(A, B, β::Real, n::Int, cache::EigenCache)
	δτ = β / n
	τs = 0:δτ:β
	λs, U = cache.λs, cache.U
	ρ = U * Diagonal(exp.(-β .* λs)) * U'
	tr_ρ = tr(ρ)
	g(τ) = tr(U * Diagonal(exp.(τ .* λs)) * U' * A * U * Diagonal(exp.(-τ .* λs)) * U' * B * ρ) / tr_ρ
	return g.(τs)
end

# <e^iHt A e^-iHt B>
function gf_real(A, B, β::Real, t::Real, n::Int, cache::EigenCache, ρ=exp(-β*cache.h))
	δt = t / n
	ts = 0:δt:t
	λs, U = cache.λs, cache.U
	# ρ = U * Diagonal(exp.(-β .* λs)) * U'
	tr_ρ = tr(ρ)
	gt(tj) = -im*tr(U * Diagonal(exp.(im*tj .* λs)) * U' * A * U * Diagonal(exp.(-im*tj .* λs)) * U' * B * ρ) / tr_ρ
	ls(tj) = im*tr(B * U * Diagonal(exp.(im*tj .* λs)) * U' * A * U * Diagonal(exp.(-im*tj .* λs)) * U' * ρ) / tr_ρ
	return gt.(ts), ls.(ts)
end
function gf_real_nn(Nop, β::Real, t::Real, n::Int, cache::EigenCache, ρ=exp(-β*cache.h))
	δt = t / n
	ts = 0:δt:t
	λs, U = cache.λs, cache.U
	# ρ = U * Diagonal(exp.(-β .* λs)) * U'
	tr_ρ = tr(ρ)
	gt(tj) = tr(U * Diagonal(exp.(im*tj .* λs)) * U' * Nop * U * Diagonal(exp.(-im*tj .* λs)) * U' * Nop * ρ) / tr_ρ
	return gt.(ts)
	
end


# ϵ_d â†â + α₀ â†â(b̂ + b̂†) + ω₀b̂†b̂ + α₁(â†ĉ + ĉ†â) + ω₁ĉ†ĉ
function noninteracting_operators(ϵ_d; ω₀=1, α₀=0.5, ω₁=1, α₁=1, d=100)
	p1 = spin_half_matrices()
	n̂, σ₊, σ₋, JW = p1["n"], p1["+"], p1["-"], -p1["z"]
	Is = one(n̂)
	p2 = boson(d=d)
	b̂, b̂′, n̂b = p2["a"], p2["adag"], p2["n"]
	Ib = one(b̂)
	# total Hamiltonian
	Nimp = kron(n̂, kron(Is, Ib))
	Himp = -ϵ_d * Nimp
	Hbath0 = ω₀ * kron(Is, kron(Is, n̂b))
	Hbath1 = ω₁ * kron(Is, kron(n̂, Ib))
	Hhyb0 = sqrt(α₀) * kron(n̂, kron(Is, b̂′ + b̂))
	tmp = kron(σ₊, kron(JW*σ₋, Ib))
	Hhyb1 = sqrt(α₁) * (tmp + tmp')
	H = Himp + Hhyb0 + Hhyb1 + Hbath0 + Hbath1
	A, B = kron(σ₋, kron(Is, Ib)), kron(σ₊, kron(Is, Ib))

	return H, A, B, Nimp, Himp + Hbath0 + Hbath1
end

function noninteracting_imag(ϵ_d; β=1, N=100, ω₀=1, α₀=0.5, ω₁=1, α₁=1)
	δτ=β/N

	H, a, adag, Nimp, H0 = noninteracting_operators(ϵ_d, ω₀=ω₀, α₀=α₀, ω₁=ω₁, α₁=α₁, d=100)

	cache = eigencache(H)
	g1 = gf_imag(a, adag, β, N, cache)
	g3 = gf_imag(Nimp, Nimp, β, N, cache)

	data_path = "result/noninteracting_eq_ED_imag_beta$(β)_mu$(ϵ_d)_dtau$(δτ)_omega0$(ω₀)_alpha0$(α₀)_omega1$(ω₁)_alpha1$(α₁).json"

	τs = collect(0:δτ:β)
	results = Dict("taus"=>τs, "gt" => g1, "nn"=>g3)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return g1, g3

end

function noninteracting_neq(ϵ_d; β=1, t=1, N=100, ω₀=1, α₀=0.5, ω₁=1, α₁=1)
	δt=t/N

	H, a, adag, Nimp, H0 = noninteracting_operators(ϵ_d, ω₀=ω₀, α₀=α₀, ω₁=ω₁, α₁=α₁, d=10)
	ρ = exp(-β*H0)
	cache = eigencache(H)
	g1, g2 = gf_real(a, adag, β, t, N, cache, ρ)
	g3 = gf_real_nn(Nimp, β, t, N, cache, ρ)

	data_path = "result/noninteracting_neq_ED_real_beta$(β)_mu$(ϵ_d)_t$(t)_N$(N)_omega0$(ω₀)_alpha0$(α₀)_omega1$(ω₁)_alpha1$(α₁).json"

	ts = collect(0:δt:t)
	results = Dict("ts"=>ts, "gt" => g1, "lt"=>g2, "nn"=>g3)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return g1, g2, g3

end

function noninteracting_eq(ϵ_d; β=1, t=1, N=100, ω₀=1, α₀=0.5, ω₁=1, α₁=1)
	δt=t/N

	H, a, adag, Nimp, H0 = noninteracting_operators(ϵ_d, ω₀=ω₀, α₀=α₀, ω₁=ω₁, α₁=α₁, d=100)
	cache = eigencache(H)
	ρ = exp(-β*cache.h)
	g1, g2 = gf_real(a, adag, β, t, N, cache, ρ)
	g3 = gf_real_nn(Nimp, β, t, N, cache, ρ)

	data_path = "result/noninteracting_eq_ED_real_beta$(β)_mu$(ϵ_d)_t$(t)_N$(N)_omega0$(ω₀)_alpha0$(α₀)_omega1$(ω₁)_alpha1$(α₁).json"

	ts = collect(0:δt:t)
	results = Dict("ts"=>ts, "gt" => g1, "lt"=>g2, "nn"=>g3)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return g1, g2, g3
end
