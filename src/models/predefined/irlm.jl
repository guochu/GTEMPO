"""
	IRLM(; μ, J, U) -> ImpurityHamiltonian

Interacting resonant level model, a spinless fermionic model with
three impurities and two baths:

H = (μ-U) n₂ + J (c†₁c₂ + c†₂c₁ + c†₂c₃ + c†₃c₂) + U (n₁n₂ + n₂n₃)

constructed as a generic `ImpurityHamiltonian` with 3 bands.
"""
function IRLM(; μ::Real, J::Real, U::Real)
	μ, J, U = float(μ), float(J), float(U)
	h = ImpurityHamiltonian(bands=3)
	push!(h, tunneling(2, 2, coeff=μ-U))

	t = tunneling(2, 1, coeff=J)
	push!(h, t)
	push!(h, t')

	t = tunneling(2, 3, coeff=J)
	push!(h, t)
	push!(h, t')

	push!(h, interaction(1, 2, 2, 1, coeff=U))
	push!(h, interaction(3, 2, 2, 3, coeff=U))
	return h
end
