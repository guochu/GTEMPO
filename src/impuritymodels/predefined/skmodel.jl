# The SK impurity model

"""
	KanamoriIM(; U, J, norb, μ=-U/2) -> ImpurityHamiltonian

Kanamori (multi-orbital SK) impurity Hamiltonian with `norb` orbitals
(2*norb bands, orbital `a` occupies bands 2a-1 and 2a):

H = μ Σᵢ nᵢ
  + U Σₐ n_{a↑} n_{a↓}
  + (U-2J) Σ_{a≠b} n_{a↑} n_{b↓}
  + (U-3J) Σ_{a>b} (n_{a↑} n_{b↑} + n_{a↓} n_{b↓})
  - J Σ_{a≠b} (c†_{a↑}c†_{a↓}c_{b↓}c_{b↑} + c†_{a↑}c†_{b↓}c_{b↑}c_{a↓})

constructed as a generic `ImpurityHamiltonian`.
"""
function KanamoriIM(; U::Real, J::Real, norb::Int, μ::Real=-U/2)
	U, J, μ = float(U), float(J), float(μ)
	h = ImpurityHamiltonian(bands=2*norb)

	# chemical potential
	for band in 1:h.bands
		push!(h, tunneling(band, band, coeff=μ))
	end
	# intra-orbital interaction
	for a in 1:norb
		i, j = 2*a-1, 2*a
		push!(h, interaction(i, j, j, i, coeff=U))
	end
	# inter-orbital density-density interaction
	for a in 1:norb
		for b in 1:norb
			if a != b
				i, j = 2*a-1, 2*b
				push!(h, interaction(i, j, j, i, coeff=U-2*J))
			end
		end
	end
	for b in 1:norb
		for a in (b+1):norb
			i, j = 2*a-1, 2*b-1
			push!(h, interaction(i, j, j, i, coeff=U-3*J))
			i, j = 2*a, 2*b
			push!(h, interaction(i, j, j, i, coeff=U-3*J))
		end
	end
	# spin-flip and pair-hopping terms
	for a in 1:norb
		for b in 1:norb
			if a != b
				i, j, k, l = 2*a-1, 2*a, 2*b-1, 2*b
				push!(h, interaction(i, j, k, l, coeff=-J))
				i, j, k, l = 2*a-1, 2*b, 2*b-1, 2*a
				push!(h, interaction(i, j, k, l, coeff=-J))
			end
		end
	end
	return h
end
