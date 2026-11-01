# wrapper for QuAPI package

"""
    correlationfunction(bath::AbstractBath, lattice::AbstractGrassmannLattice)

Compute the discrete hybridization functions using QuAPI
"""
function correlationfunction(bath::AbstractBath, lattice::ImagGrassmannLattice)
    # @assert lattice.β == bath.β
    isapprox(lattice.β, bath.β; rtol=1.0e-12, atol=1.0e-14) || @warn "lattice.β=$(lattice.β), but bath.β=$(bath.β)"
    Δτ(bath, N=lattice.N, δτ=lattice.δτ)
end 
correlationfunction(bath::AbstractBath, lattice::RealGrassmannLattice) = Δt(bath, N=lattice.N, t=lattice.t) 
function correlationfunction(bath::AbstractBath, lattice::MixedGrassmannLattice)
    isapprox(lattice.β, bath.β; rtol=1.0e-12, atol=1.0e-14) || @warn "lattice.β=$(lattice.β), but bath.β=$(bath.β)"
    Δm(bath, Nτ=lattice.Nτ, t=lattice.t, Nt=lattice.Nt)
end  

# ImagCorrelationFunction{T} 的参数 T 是 CorrelationMatrix 的元素类型（QuAPI 2026-10-30 重构后）
TK.scalartype(::Type{<:ImagCorrelationFunction{T}}) where {T} = T
TK.scalartype(::Type{<:RealCorrelationFunction}) = ComplexF64
TK.scalartype(::Type{<:MixedCorrelationFunction}) = ComplexF64
TK.scalartype(::Type{<:BCSCorrelationFunction{M}}) where M = scalartype(M)
# ---------------------------------------------------------------------------
# 2026-10-30 QuAPI unified real-time convention adaptation: the compact
# RealCorrelationFunction / GenericRealCorrelationFunction now embed the
# Keldysh relative minus into the CROSS-branch blocks
# ((+,−) = −(c, b), (−,+) = −(a, d) — same storage as MixedCorrelationFunction),
# while the same-branch blocks are unchanged. GTEMPO's IF exponents and
# observable formulas consume the PHYSICAL hybridization blocks, so every
# cross-branch read must compensate with one minus.
# ---------------------------------------------------------------------------
"""
	_physical_branch(corr::AbstractCorrelationFunction, b1::Symbol, b2::Symbol)

The physical (Keldysh-sign-free) branch block of a real-time correlation
function: undoes the embedded relative minus of the unified QuAPI storage
for cross-branch blocks. Same-branch blocks and non-real-time contours pass
through unchanged.
"""
function _physical_branch(corr, b1::Symbol, b2::Symbol)
	η = branch(corr; b1=b1, b2=b2)
	return _embedded_real_cross(corr, b1, b2) ? -η : η
end

# true when (b1, b2) addresses a cross block of a sign-embedded real-time
# representation — either directly, or through the BCS wrapper on the real
# contour (its uu/ud/du wrap the embedded representation there)
function _embedded_real_cross(corr, b1::Symbol, b2::Symbol)
	(b1 in (:+, :-) && b2 in (:+, :-) && b1 !== b2) || return false
	corr isa AbstractRealCorrelationFunction && return true
	return corr isa BCSCorrelationFunction && corr.uu isa AbstractRealCorrelationFunction
end

"""
	_physical_coef(corr::AbstractCorrelationFunction, row, col)

Contraction coefficient between the `row` (ā-type) and `col` (a-type) legs
with the unified real-time convention compensated: the QuAPI accessor
serves the sign-embedded cross blocks, the IF exponent needs the physical
ones.
"""
function _physical_coef(corr, row, col)
	v = corr[conj(row), conj(col)]
	return _embedded_real_cross(corr, branch(row), branch(col)) ? -v : v
end
