"""
    hybriddynamics(lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; band=1, trunc=Defaults.alg_trunc_if())
    hybriddynamics(lattice::AbstractGrassmannLattice, corr, alg::InfluenceFunctionalAlgorithm; band=1)
    hybriddynamics!(gmps, lattice, corr; band=1, trunc=Defaults.alg_trunc_if())

Build the Grassmann MPS representation `I` of the hybridization influence
functional from the discretized bath correlation function `corr` on `lattice`.

- The plain methods (and `hybriddynamics!` acting in place on `gmps`) use the
  default `PartialIF` algorithm (product of bond-dimension-2 partial MPOs,
  [SciPost Phys. Core 7, 063 (2024)]).
- Passing an explicit algorithm (`PartialIF`, `XTRGIF`, `ExactTTIIF`, `TDVPIF`)
  as the third positional argument selects that construction algorithm; see the
  docstrings of the individual algorithm types for the available options.
"""
function hybriddynamics!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, corr::RealCorrelationFunction; band::Int=1, trunc::TruncationScheme=Defaults.alg_trunc_if())
	(1 <= band <= lattice.bands) || throw(BoundsError(1:lattice.bands, band))
	η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = corr.G₊₊, corr.G₊₋, corr.G₋₊, corr.G₋₋
	@assert size(η⁺⁺) == size(η⁺⁻) == size(η⁻⁺) == size(η⁻⁻)
	k = lattice.k
	for i in 1:k
		tmp1 = partialif_hybrid(lattice, i, view(η⁺⁺, i, 1:k), view(η⁺⁻, i, 1:k), b1=:+, band=band)
		tmp3 = partialif_hybrid(lattice, i, view(η⁻⁺, i, 1:k), view(η⁻⁻, i, 1:k), b1=:-, band=band)

		gmps = mult!(gmps, tmp1, trunc=trunc)[1]
		gmps = mult!(gmps, tmp3, trunc=trunc)[1]
	end
	return gmps
end

function hybriddynamics_naive!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, corr::RealCorrelationFunction; band::Int=1, trunc::TruncationScheme=Defaults.alg_trunc_if())
	(1 <= band <= lattice.bands) || throw(BoundsError(1:lattice.bands, band))
	η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = corr.G₊₊, corr.G₊₋, corr.G₋₊, corr.G₋₋
	@assert size(η⁺⁺) == size(η⁺⁻) == size(η⁻⁺) == size(η⁻⁻)
	k = lattice.k
	for i in 1:k
		tmp1 = partialif_hybrid_naive(lattice, i, view(η⁺⁺, i, 1:k), view(η⁺⁻, i, 1:k), b1=:+, band=band, trunc=trunc)
		tmp3 = partialif_hybrid_naive(lattice, i, view(η⁻⁺, i, 1:k), view(η⁻⁻, i, 1:k), b1=:-, band=band, trunc=trunc)

		gmps = mult!(gmps, tmp1, trunc=trunc)[1]
		gmps = mult!(gmps, tmp3, trunc=trunc)[1]
	end
	return gmps		
end

function partialif_hybrid(lattice::AbstractRealGrassmannLattice, i::Int, cols::AbstractVector; b1::Symbol, b2::Symbol, band::Int=1)
	row = lattice[GrassmannContourIndex(i, b1, band, true)]
	col_pos = [lattice[GrassmannContourIndex(j, b2, band, false)] for j in length(cols):-1:1]
	mpo = partialmpo(row, col_pos, reverse(cols))
	return mpo * vacuumstate(promote_type(scalartype(lattice), eltype(cols)), lattice)
end
function partialif_hybrid(lattice::AbstractRealGrassmannLattice, i::Int, cols_f::AbstractVector, cols_b::AbstractVector; b1::Symbol, band::Int=1)
	@assert length(cols_f) == length(cols_b)
	row = lattice[GrassmannContourIndex(i, b1, band, true)]
	cols = eltype(cols_f)[]
	col_pos = Int[]
	for j in length(cols_f):-1:1
		pos1, pos2 = lattice[GrassmannContourIndex(j, :+, band, false)], lattice[GrassmannContourIndex(j, :-, band, false)]
		push!(col_pos, pos1)
		push!(col_pos, pos2)
		push!(cols, cols_f[j])
		push!(cols, cols_b[j])
	end
	mpo = partialmpo(row, col_pos, cols)
	return mpo * vacuumstate(promote_type(scalartype(lattice), eltype(cols)), lattice)
end
function partialif_hybrid_naive(lattice::AbstractRealGrassmannLattice, i::Int, cols_f::AbstractVector, cols_b::AbstractVector; 
							b1::Symbol, band::Int=1, trunc::TruncationScheme=Defaults.alg_trunc_if())
	@assert length(cols_f) == length(cols_b)
	alg = Orthogonalize(orth=TK.SVD(), trunc=trunc)
	gmps = vacuumstate(promote_type(scalartype(lattice), eltype(cols_f), eltype(cols_b)), lattice)
	for j in length(cols_f):-1:1
		pos1, pos2 = lattice[GrassmannContourIndex(i, b1, band, true)], lattice[GrassmannContourIndex(j, :+, band, false)]
		t = exp(GTerm(pos1, pos2, coeff=cols_f[j]))
		apply!(t, gmps)
		canonicalize!(gmps, alg=alg)
		pos1, pos2 = lattice[GrassmannContourIndex(i, b1, band, true)], lattice[GrassmannContourIndex(j, :-, band, false)]
		t = exp(GTerm(pos1, pos2, coeff=cols_b[j]))
		apply!(t, gmps)
		canonicalize!(gmps, alg=alg)
	end
	return gmps
end

# reverse col and row, only used in stepper
function partialif_hybrid(lattice::AbstractRealGrassmannLattice, rows::AbstractVector, j::Int; b1::Symbol, b2::Symbol, band::Int=1)
	col = lattice[GrassmannContourIndex(j, b2, band, false)]
	row_pos = [lattice[GrassmannContourIndex(i, b1, band, true)] for i in length(rows):-1:1]
	mpo = partialmpo(col, row_pos, -reverse(rows))
	return mpo * vacuumstate(promote_type(scalartype(lattice), eltype(rows)), lattice)
end
function partialif_hybrid(lattice::AbstractRealGrassmannLattice, rows_f::AbstractVector, rows_b::AbstractVector, j::Int; b2::Symbol, band::Int=1)
	col = lattice[GrassmannContourIndex(j, b2, band, false)]
	rows = eltype(rows_f)[]
	row_pos = Int[]
	for i in length(rows_f):-1:1
		pos1, pos2 = lattice[GrassmannContourIndex(i, :+, band, true)], lattice[GrassmannContourIndex(i, :-, band, true)]
		push!(row_pos, pos1)
		push!(row_pos, pos2)
		push!(rows, -rows_f[i])
		push!(rows, -rows_b[i])
	end
	mpo = partialmpo(col, row_pos, rows)
	return mpo * vacuumstate(promote_type(scalartype(lattice), eltype(rows)), lattice)
end