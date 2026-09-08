###--------------imaginary time----------------

"""
    gf(lattice::AbstractGrassmannLattice, a::ContourIndex, b::ContourIndex, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
                alg=ExactIntegrate(), Z = integrate(lattice, A, B..., alg=alg))

Calculate the Green's function ⟨x y⟩, where x, y are GVs 
specified by contour indices a and b

If A is type Vector{GrassmannMPS}, it means A is a sum of GMPSs
The multiplication of A and B... is assumed
alg: the algorithm to perform the integration
Z: the partition function values
"""
function gf(lattice::AbstractGrassmannLattice, a::NTuple{N, ContourIndex}, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
            alg::IntegrationAlgorithm=ExactIntegrate(), Z::Number = integrate(lattice, A, B..., alg=alg)) where {N}
    pos = map(x->lattice[x], a)
    t = GTerm(pos, coeff=1)
    A2 = _mult_A(t, A)
    return integrate(lattice, A2, B..., alg=alg)/Z    
end

"""
    contour_ordered_gf(lattice::AbstractGrassmannLattice, a::ContourIndex, b::ContourIndex, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
                alg=ExactIntegrate(), Z = integrate(lattice, A, B..., alg=alg))

similar to gf, but return the contour ordered Green's function
"""
function contour_ordered_gf(lattice::AbstractGrassmannLattice, a::ContourIndex, b::ContourIndex, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
                            alg::IntegrationAlgorithm=ExactIntegrate(), Z::Real = integrate(lattice, A, B..., alg=alg)) 
    ((!a.conj) && (b.conj)) || throw(ArgumentError("conj(a)=false and conj(b)=true should be satisfied"))
    return (a < b) ? -gf(lattice, (b, a), A, B...; alg=alg, Z=Z) : gf(lattice, (a, b), A, B...; alg=alg, Z=Z) 
end

"""
    greater(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
            band::Int=1, alg=ExactIntegrate(), Z = integrate(lattice, A, B..., alg=alg))

Return the greater Green's function ⟨aᵢ bⱼ⟩ on the Keldysh or the Kadanoff-Baym contour
"""
function greater(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; 
                 band::Union{Int, Tuple{Int, Int}}=1, alg::IntegrationAlgorithm=ExactIntegrate(), Z::Number = integrate(lattice, A, B..., alg=alg))
    @assert i >= j
    band isa Int && (band = (band, band))
    a = ContourIndex(i, conj=false, branch=:+, band=band[1])
    b = ContourIndex(j, conj=true, branch=:+, band=band[2])
    return gf(lattice, (a, b), A, B...; alg=alg, Z=Z)
end
greater(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; kwargs...) = greater(lattice, i, 1, A, B...; kwargs...)

"""
    lesser(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
            band::Int=1, alg=ExactIntegrate(), Z = integrate(lattice, A, B..., alg=alg))

Return the lesser Green's function ⟨aᵢ bⱼ⟩ on the Keldysh or the Kadanoff-Baym contour
"""
function lesser(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; 
                 band::Union{Int, Tuple{Int, Int}}=1, alg::IntegrationAlgorithm=ExactIntegrate(), Z::Number = integrate(lattice, A, B..., alg=alg))
    @assert i <= j
    band isa Int && (band = (band, band))
    a = ContourIndex(i, conj=true, branch=:-, band=band[1])
    b = ContourIndex(j, conj=false, branch=:+, band=band[2])
    return gf(lattice, (a, b), A, B...; alg=alg, Z=Z)
end
lesser(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; kwargs...) = lesser(lattice, 1, i, A, B...; kwargs...)

# other real-time observables

