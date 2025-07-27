# wrapper for QuAPI package

"""
    correlationfunction(bath::AbstractBath, lattice::AbstractGrassmannLattice)

Compute the discrete hybridization functions using QuAPI
"""
function correlationfunction(bath::AbstractBath, lattice::ImagGrassmannLattice1Order)
    # @assert lattice.β == bath.β
    (lattice.β == bath.β) || @warn "lattice.β=$(lattice.β), but bath.β=$(bath.β)"
    Δτ(bath, N=lattice.N, δτ=lattice.δτ)
end 
correlationfunction(bath::AbstractBath, lattice::RealGrassmannLattice1Order) = Δt(bath, N=lattice.N, t=lattice.t) 
correlationfunction(bath::AbstractBath, lattice::RealGrassmannLattice2Order) = Δt(bath, N=2*lattice.N, t=lattice.t)
function correlationfunction(bath::AbstractBath, lattice::MixedGrassmannLattice1Order)
    (lattice.β == bath.β) || @warn "lattice.β=$(lattice.β), but bath.β=$(bath.β)"
    Δm(bath, Nτ=lattice.Nτ, t=lattice.t, Nt=lattice.Nt)
end  

TK.scalartype(::Type{<:ImagCorrelationFunction{<:AbstractMatrix{T}}}) where {T} = T
TK.scalartype(::Type{<:RealCorrelationFunction}) = ComplexF64
TK.scalartype(::Type{<:MixedCorrelationFunction}) = ComplexF64
TK.scalartype(::Type{<:BCSCorrelationFunction{M}}) where M = scalartype(M)