

"""
    correlationfunction(bath::AbstractBath, lattice::AbstractGrassmannLattice)

Compute the discrete correlation functions (QUAPI). More details see
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