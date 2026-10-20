"""
	bond Schmidt 奇异值：键 `bond` 的 Schmidt 谱，归一化为 Σs² = 1。

	`bond` 位于 site `bond` 与 `bond+1` 之间（`1 <= bond <= length(psi)-1`），
	谱存放在 `psi.s[bond+1]`（Z2 分块对角张量）。若 Schmidt 值尚未初始化，
	就地调用 `canonicalize!`（与 FiniteMPSAlgorithms 的语义一致）。
"""
function _bond_spectrum(psi::GrassmannMPS, bond::Integer)
	(1 <= bond <= length(psi) - 1) || throw(BoundsError(1:length(psi)-1, bond))
	svectors_uninitialized(psi) && canonicalize!(psi)
	s = psi.s[bond + 1]
	ismissing(s) && throw(ArgumentError("Schmidt values at bond $bond are not initialized"))
	sv = Float64[]
	for (f1, f2) in fusiontrees(s)
		blk = s[f1, f2]
		for i in 1:size(blk, 1)
			push!(sv, abs(blk[i, i]))
		end
	end
	nrm = sqrt(sum(abs2, sv))
	(nrm > zero(nrm)) || throw(ArgumentError("Schmidt spectrum at bond $bond is zero"))
	return sv ./ nrm
end

"""
	schmidt_values(psi::GrassmannMPS; bond=div(length(psi), 2)) -> Vector{Float64}

返回键 `bond` 的 Schmidt 奇异值（降序不保证）。
语义与 `FiniteMPSAlgorithms.schmidt_values` 一致，只是谱取自 Z2 分块对角张量。
"""
schmidt_values(psi::GrassmannMPS; bond::Integer=div(length(psi), 2)) =
	_bond_spectrum(psi, bond)

"""
	entanglement_spectrum(psi::GrassmannMPS; bond=div(length(psi), 2)) -> Vector{Float64}

返回键 `bond` 的纠缠谱（Schmidt 奇异值平方，归一化概率向量）。
语义与 `FiniteMPSAlgorithms.entanglement_spectrum` 一致。
"""
entanglement_spectrum(psi::GrassmannMPS; bond::Integer=div(length(psi), 2)) =
	abs2.(schmidt_values(psi; bond=bond))

"""
	entanglement_entropy(psi::GrassmannMPS; bond=div(length(psi), 2), α=1) -> Float64

键 `bond` 的 Rényi 纠缠熵（自然对数）：`α=1` 为 von Neumann 熵，
`α=2` 为二阶 Rényi 熵。语义与 `FiniteMPSAlgorithms.entanglement_entropy` 一致。

对由 `hybriddynamics` 构造的影响泛函 `psi`，该值就是影响泛函在相应时间割上的
operator-space entanglement entropy（OSEE）。
"""
entanglement_entropy(psi::GrassmannMPS; bond::Integer=div(length(psi), 2), α::Real=1) =
	FMA.renyi_entropy(abs2.(schmidt_values(psi; bond=bond)); α)
