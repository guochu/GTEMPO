# real-time first order
"""
    occupation(lattice::RealGrassmannLattice1Order, i::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
            band::Int=1, alg=ExactIntegrate(), Z = integrate(lattice, A, B..., alg=alg))

Return the occupation at time step i on the Keldysh contour
"""
function occupation(lattice::RealGrassmannLattice1Order, i::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; kwargs...) 
    return real(Gt(lattice, i, i, A, B...; c1=false, c2=true, b1=:+, b2=:-, kwargs...))
end
occupation(lattice::RealGrassmannLattice1Order, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; band::Int=1, 
            alg::IntegrationAlgorithm=ExactIntegrate(), 
            Z::Number = integrate(lattice, A, B..., alg=alg)) = [occupation(lattice, i, A, B...; alg=alg, Z=Z, band=band) for i in 1:lattice.N]

function occupation(lattice::ImagGrassmannLattice1Order, i::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; kwargs...) 
    return 1 + Gτ(lattice, i, i, A, B...; c1=true, c2=false, kwargs...)
end


# real-time second order
function occupation(lattice::RealGrassmannLattice2Order, A::GrassmannMPS, B::Vararg{GrassmannMPS}; kwargs...) 
    return real(Gt(lattice, lattice.k, lattice.k, A, B...; c1=false, c2=true, b1=:+, b2=:-, kwargs...))
end


function occupation2(lattice::AbstractGrassmannLattice, i::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; band::Int=1, 
                     branch::Symbol=:τ, alg::IntegrationAlgorithm=ExactIntegrate(), Z::Number = integrate(lattice, A, B..., alg=alg)) 
    A′ = _insert_n(lattice, A, i, branch, band)
    return integrate(lattice, A′, B..., alg=alg) / Z
end

_insert_n(lattice, A::GrassmannMPS, i, b1, band) = insert_n(lattice, A, i, branch=b1, band=band)
_insert_n(lattice, A::Vector, i, b1, band) = [insert_n(lattice, Aj, i, branch=b1, band=band) for Aj in A]
