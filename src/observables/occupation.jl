# NOTE: occupation2 (the equivalent-time Green's function route, ⟨aᵢ†aᵢ⟩ = 1 + Gτ(i,i)
# or = G⁻(tᵢ,tᵢ)) has been commented out: benchmarked against the Toulouse ED reference,
# it is slightly more accurate than `occupation` on the Keldysh contour (0.35% vs 0.53%,
# both at the discretization/truncation noise level), but wrong on the imaginary-time
# contour for i ≥ 2 (a spurious jump of ~6%; only i = 1 happens to be correct, so the
# boundary Grassmann trace contribution is not properly accounted for).
#
# # real-time first order
# """
#     occupation2(lattice::RealGrassmannLattice1Order, i::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
#             band::Int=1, alg=ExactIntegrate(), Z = integrate(lattice, A, B..., alg=alg))
#
# Return the occupation at time step i on the Keldysh contour
# """
# function occupation2(lattice::RealGrassmannLattice1Order, i::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; band::Int=1, kwargs...)
#     a = ContourIndex(i, conj=false, branch=:+, band=band)
#     b = ContourIndex(i, conj=true, branch=:-, band=band)
#     return real(gf(lattice, (a, b), A, B...; kwargs...))
# end
# occupation2(lattice::RealGrassmannLattice1Order, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; band::Int=1,
#             alg::IntegrationAlgorithm=ExactIntegrate(),
#             Z::Number = integrate(lattice, A, B..., alg=alg)) = [occupation2(lattice, i, A, B...; alg=alg, Z=Z, band=band) for i in 1:lattice.N]
#
# function occupation2(lattice::ImagGrassmannLattice1Order, i::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; band::Int=1, kwargs...)
#     a = ContourIndex(i, conj=true, branch=:τ, band=band)
#     b = ContourIndex(i, conj=false, branch=:τ, band=band)
#     return 1 + gf(lattice, (a, b), A, B...; kwargs...)
# end
#
#
# # real-time second order
# function occupation2(lattice::RealGrassmannLattice2Order, A::GrassmannMPS, B::Vararg{GrassmannMPS}; band::Int=1, kwargs...)
#     i = lattice.k
#     a = ContourIndex(i, conj=false, branch=:+, band=band)
#     b = ContourIndex(i, conj=true, branch=:-, band=band)
#     return real(gf(lattice, (a, b), A, B...; kwargs...))
# end


"""
    occupation(lattice::AbstractGrassmannLattice, i::Int, A, B...; band=1, branch=:τ, alg=ExactIntegrate(), Z=integrate(lattice, A, B...))

Impurity occupation number ⟨n̂(tᵢ)⟩ at the `i`-th time step, computed by
inserting the density operator n̂ = a†a (as a `GTerm`) at position `i` on
`branch` and integrating the product with the influence functional `B...`
normalized by `Z`. For imaginary-time lattices the default `branch = :τ` is
correct; on real-time lattices pass `branch = :+`.
"""
function occupation(lattice::AbstractGrassmannLattice, i::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; band::Int=1,
                     branch::Symbol=:τ, alg::IntegrationAlgorithm=ExactIntegrate(), Z::Number = integrate(lattice, A, B..., alg=alg))
    A′ = _insert_n(lattice, A, i, branch, band)
    return integrate(lattice, A′, B..., alg=alg) / Z
end

_insert_n(lattice, A::GrassmannMPS, i, b1, band) = insert_n(lattice, A, i, branch=b1, band=band)
_insert_n(lattice, A::Vector, i, b1, band) = [insert_n(lattice, Aj, i, branch=b1, band=band) for Aj in A]
