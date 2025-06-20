
f(D, ϵ) = sqrt(1-(ϵ/D)^2) / π
spectrum_func(D) = SpectrumFunction(ϵ->f(D, ϵ), lb=-D, ub=D)
spectrum_func() = spectrum_func(10)
spectrum_func2(D=10) = SpectrumFunction(ϵ->f(D, ϵ), lb=0, ub=D)

function _error(a::Number, b::Number, tol::Real)
	ab = a - b
	if abs(a) < tol
		return abs(ab)
	else
		return abs(ab / a)
	end
end


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
	# Hbathbare = ω₀ * n̂b
	# return H, A, B, Hbathbare
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
	# Hbathbare = ω₀ * n̂b
	# return H, A, B, Hbathbare
end


# function noninteracting_real(ϵ_d)
# 	β = 0.1
# 	δt=0.01
# 	N = 10
# 	t = N * δt

# 	H, a, adag, H0 = noninteracting_operators(ϵ_d, ω₀=1, α=0.5, d=100)
# 	ρ = exp(-β*H0)
# 	g1, g2 = gf_real(H, a, adag, β, t, N, ρ)

# 	return g1, g2
# end

# function interacting_real(U, ϵ_d=U/2)
# 	β = 0.1
# 	δt=0.01
# 	N = 10
# 	t = N * δt

# 	H, a, adag, H0 = interacting_operators(U, ϵ_d, ω₀=1, α=0.5, d=100)
# 	ρ = exp(-β*H0)
# 	g1, g2 = gf_real(H, a, adag, β, t, N, ρ)

# 	return g1, g2
# end
