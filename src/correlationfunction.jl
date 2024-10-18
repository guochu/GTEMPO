

"""
    correlationfunction(bath::AbstractFermionicBath, lattice::AbstractGrassmannLattice)

Compute the discrete correlation functions (QUAPI)
"""
function correlationfunction(bath::AbstractFermionicBath, lattice::ImagGrassmannLattice1Order)
    # @assert lattice.β == bath.β
    (lattice.β == bath.β) || @warn "lattice.β=$(lattice.β), but bath.β=$(bath.β)"
    Cτ(bath, N=lattice.N, δτ=lattice.δτ)
end 
correlationfunction(bath::AbstractFermionicBath, lattice::RealGrassmannLattice1Order) = Ct(bath, N=lattice.N, t=lattice.t) 
correlationfunction(bath::AbstractFermionicBath, lattice::RealGrassmannLattice2Order) = Ct(bath, N=2*lattice.N, t=lattice.t)
function correlationfunction(bath::AbstractFermionicBath, lattice::MixedGrassmannLattice1Order)
    (lattice.β == bath.β) || @warn "lattice.β=$(lattice.β), but bath.β=$(bath.β)"
    Cm(bath, Nτ=lattice.Nτ, t=lattice.t, Nt=lattice.Nt)
end  