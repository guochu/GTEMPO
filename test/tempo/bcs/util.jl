


"""
H = -ϵ_d(n̂₁+n̂₂) + Un̂₁n̂₂ + √α(â₁†ĉ₁ + ĉ₁†â₁ + â₂†ĉ₂ + ĉ₂†â₂) + ω₀(m̂₁+m̂₂) - (Δĉ₁†ĉ†₂ + Δ'ĉ₂ĉ₁)
where n̂ⱼ = âⱼ†âⱼ, m̂ⱼ = ĉⱼ†ĉⱼ
"""
function bcs_operators(U, ϵ_d; ω₀=1, α=0.5, Δ=0.3)
	p1 = spin_half_matrices()
	n̂, σ₊, σ₋, JW = p1["n"], p1["+"], p1["-"], -p1["z"]
	Is = one(n̂)
	n_ud = kron(n̂, Is) + kron(Is, n̂)
	nn = kron(n̂,n̂)
	I_ud = kron(Is, Is)
	# total Hamiltonian
	Himpbare = -ϵ_d*n_ud + U * nn
	Himp = kron(Himpbare, I_ud)
	Hbathbare = ω₀ * n_ud - (Δ * kron(JW*σ₊, σ₊)- conj(Δ) * kron(JW*σ₋, σ₋))
	Hbath = kron(I_ud, Hbathbare)
	tmp = sqrt(α) * (kron(kron(kron(JW*σ₊, JW), σ₋), Is) + kron(Is, kron(JW*σ₊, kron(JW, σ₋))))
	Hhyb = tmp + tmp'
	H = Himp + Hhyb + Hbath
	A, B = kron(kron(σ₋, Is), I_ud), kron(kron(σ₊, Is), I_ud)

	return H, A, B, Himp + Hbath
end

"""
H = -ϵ_d(n̂₁+n̂₂) + Un̂₁n̂₂ + â₁†â₂ + â₂†â₁ + â₁†â₂† + â₂â₁
	√α(â₁†ĉ₁ + ĉ₁†â₁ + â₂†ĉ₂ + ĉ₂†â₂) + ω₀(m̂₁+m̂₂) - (Δĉ₁†ĉ†₂ + Δ'ĉ₂ĉ₁)
where n̂ⱼ = âⱼ†âⱼ, m̂ⱼ = ĉⱼ†ĉⱼ
"""

function bcs_operators2(U, ϵ_d; ω₀=1, α=0.5, Δ=0.3)
	p1 = spin_half_matrices()
	n̂, σ₊, σ₋, JW = p1["n"], p1["+"], p1["-"], -p1["z"]
	Is = one(n̂)
	n_ud = kron(n̂, Is) + kron(Is, n̂)
	nn = kron(n̂,n̂)
	I_ud = kron(Is, Is)
	# total Hamiltonian
	Himpbare = -ϵ_d*n_ud + U * nn
	tmp = kron(JW*σ₊, σ₋)
	tmp = tmp + tmp'
	Himpbare += tmp
	tmp = kron(JW*σ₊, σ₊)
	tmp = tmp + tmp'
	Himpbare += tmp
	Himp = kron(Himpbare, I_ud)
	Hbathbare = ω₀ * n_ud - (Δ * kron(JW*σ₊, σ₊)- conj(Δ) * kron(JW*σ₋, σ₋))
	Hbath = kron(I_ud, Hbathbare)
	tmp = sqrt(α) * (kron(kron(kron(JW*σ₊, JW), σ₋), Is) + kron(Is, kron(JW*σ₊, kron(JW, σ₋))))
	Hhyb = tmp + tmp'
	H = Himp + Hhyb + Hbath
	A, B = kron(kron(σ₋, Is), I_ud), kron(kron(σ₊, Is), I_ud)

	return H, A, B, Hbathbare
end

function bcs_siam(; μ::Real, U::Real)
	h = ImpurityHamiltonian(bands=2)
	push!(h, interaction(1,2,2,1, coeff=U))
	for band in 1:h.bands
		push!(h, tunneling(band, band, coeff=μ))
	end
	t = tunneling(1, 2, coeff=1)
	push!(h, t)
	push!(h, t')
	t = adagadag(1, 2, coeff=1)
	push!(h, t)
	push!(h, t')
	return h
end

# function bcs_gfs(bath::AbstractDiscreteBath, ϵ_d, Δ)
# 	m = bcs_cmatrix(bath, ϵ_d, Δ)
	
# end

# function bcs_cmatrix(bath::AbstractDiscreteBath, ϵ_d, Δ)
# 	N = num_sites(bath)
# 	m = zeros(2*(N+1), 2*(N+1))
# 	m[1,1] = m[2,2] = -ϵ_d
# 	ws, fs = frequencies(bath), spectrumvalues(bath)
# 	for i in 1:n
# 		j = i + 1
# 		m[2*j-1, 2*j-1] = ws[i]
# 		m[2*j, 2*j] = ws[i]
# 		m[2*j-1, 2*j] = -Δ
# 		m[2*j, 2*j-1] = -Δ
# 		m[1, 2*j-1] = fs[i]
# 		m[2*j-1-1] = fs[i]
# 		m[1, 2*j] = fs[i]
# 		m[2*j-1] = fs[i]
# 	end
# 	return m
# end