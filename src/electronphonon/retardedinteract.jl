"""
	get_pair_pos(lattice::AbstractGrassmannLattice, j::Int, band::Int, b::Symbol)

密度算符 `n̂` 在第 `j` 个节点的两条 Grassmann 腿（链上位置）以及 Keldysh
分支的定向符号：`:+`/`:τ` 为 `(ā_{j+1}, a_j, +1)`，`:-` 为
`(a_{j+1}, ā_j, -1)`。
"""
function get_pair_pos(lattice::AbstractGrassmannLattice, j::Int, band::Int, b::Symbol)
	if b == :-
		return lattice[GrassmannContourIndex(j+1, b, band, false)], lattice[GrassmannContourIndex(j, b, band, true)], -1
	else
		return lattice[GrassmannContourIndex(j+1, b, band, true)], lattice[GrassmannContourIndex(j, b, band, false)], 1
	end
end

"""
	_retarded_nodes(lattice::AbstractGrassmannLattice)

`(branch, time)` 节点表；`:τ` 分支有 `Nτ` 个节点，`:+`/`:-` 有 `Nt` 个
（Grassmann lattice 的 `kτ`/`kt` 多出的边界节点不进入 IF）。
"""
_retarded_nodes(lattice::AbstractGrassmannLattice) =
	[(br, i) for br in branches(lattice) for i in 1:(br === :τ ? lattice.Nτ : lattice.Nt)]

"""
	retardedinteractdynamics_naive!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice,
	                                corr::AbstractCorrelationFunction; trunc)

Build the retarded (density–density) interaction `exp(Σᵢⱼ Δᵢⱼ n̂ᵢ n̂ⱼ)` as a
Grassmann MPS by applying one `exp(Δᵢⱼ n̂ᵢ n̂ⱼ)` factor at a time.

`n̂ᵢ = Σ_band āᵢaᵢ` is the total occupation of contour node `i`; the sum runs
over **all** band pairs, so there is no restriction on the number of bands.
The contour branches are read from `lattice`, so the same implementation
covers imaginary, real and mixed time.
"""
function retardedinteractdynamics_naive!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice,
										 corr::AbstractCorrelationFunction;
										 trunc::TruncationScheme=Defaults.alg_trunc_if())
	alg = Orthogonalize(orth=TK.SVD(), trunc=trunc)
	T = promote_type(scalartype(lattice), scalartype(corr))
	nodes = _retarded_nodes(lattice)
	for (br1, i) in nodes, band1 in 1:lattice.bands
		pos1a, pos1b, c1 = get_pair_pos(lattice, i, band1, br1)
		tmp = vacuumstate(T, lattice)
		for (br2, j) in nodes, band2 in 1:lattice.bands
			pos2a, pos2b, c2 = get_pair_pos(lattice, j, band2, br2)
			c = corr[GrassmannContourIndex(i, br1, band1, false),
					 GrassmannContourIndex(j, br2, band1, true)] * c1 * c2
			t = pos1a == pos2a ? exp(GTerm(pos1a, pos1b, coeff=c)) :
								 exp(GTerm(pos1a, pos1b, pos2a, pos2b, coeff=c))
			apply!(t, tmp)
			canonicalize!(tmp, alg=alg)
		end
		mult!(gmps, tmp, trunc=trunc)
	end
	return gmps
end

"""
	retardedinteractdynamics_naive(lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...)

Build the retarded interaction `e^{ΣᵢⱼΔᵢⱼnᵢnⱼ}` as a GMPS using the PartialIF algorithm
(the TTI-IF algorithm for this case is to be developed)
corr: bosonic hybridization function calculated using QuAPI
"""
retardedinteractdynamics_naive(lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...) = retardedinteractdynamics_naive!(
				vacuumstate(promote_type(scalartype(lattice), scalartype(corr)), lattice), lattice, corr; kwargs...)
