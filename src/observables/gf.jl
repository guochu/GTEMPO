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

# function gf(lattice::AbstractGrassmannLattice, a::ContourIndex, b::ContourIndex, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
#                 alg::IntegrationAlgorithm=ExactIntegrate(), Z::Number = integrate(lattice, A, B..., alg=alg))
#     pos1, pos2 = lattice[a], lattice[b]
#     t = GTerm(pos1, pos2, coeff=1)
#     A2 = _mult_A(t, A)
#     return integrate(lattice, A2, B..., alg=alg)/Z    
# end

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
    Gτ(lattice::ImagGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
        band::Int=1, c1::Bool=false, c2::Bool=true, alg=ExactIntegrate(), Z = integrate(lattice, A, B..., alg=alg))

Return the Matsubara Green's function ⟨aᵢ bⱼ⟩ on the imaginary-time axis
band: the band to calculate the Green's function
c1: whether to take the conjugate of the first GV a
c2: whether to take the conjugate of the second GV b
The other keywords are the same as gf
"""
function Gτ(lattice::ImagGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
            band::Union{Int, Tuple{Int, Int}}=1, c1::Bool=false, c2::Bool=true, alg::IntegrationAlgorithm=ExactIntegrate(), Z::Real = integrate(lattice, A, B..., alg=alg))
    if isa(band, Int)
        band = (band, band)
    end
    band1, band2 = band
    a, b = ContourIndex(i, conj=c1, branch=:τ, band=band1), ContourIndex(j, conj=c2, branch=:τ, band=band2)
    return gf(lattice, (a, b), A, B...; alg=alg, Z=Z)
end
Gτ(lattice::ImagGrassmannLattice, i::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; kwargs...) = Gτ(lattice, i, 1, A, B...; kwargs...)

###--------------real time 1 order----------------

"""
    Gt(lattice::RealGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
        b1, b2, band::Int=1, c1::Bool=true, c2::Bool=false, alg=ExactIntegrate(), Z = integrate(lattice, A, B..., alg=alg))

Return the real-time Green's function ⟨aᵢ bⱼ⟩ on the Keldysh contour
band: the band to calculate the Green's function
c1: whether to take the conjugate of the first GV a
c2: whether to take the conjugate of the second GV b
b1: the branch of the first GV (can be :+ or :-)
b2: the branch of the second GV (can be :+ or :-)
The other keywords are the same as gf
"""
function Gt(lattice::RealGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
            b1::Symbol, b2::Symbol, c1::Bool=true, c2::Bool=false, band::Union{Int, Tuple{Int, Int}}=1, 
            alg::IntegrationAlgorithm=ExactIntegrate(), 
            Z::Number = integrate(lattice, A, B..., alg=alg))
    if isa(band, Int)
        band = (band, band)
    end
    band1, band2 = band
    a, b = ContourIndex(i, conj=c1, branch=b1, band=band1), ContourIndex(j, conj=c2, branch=b2, band=band2)
    return gf(lattice, (a, b), A, B...; alg=alg, Z=Z)
end

###--------------mixed time 1 order----------------

"""
    Gm(lattice::MixedGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
        b1, b2, band::Int=1, c1::Bool=true, c2::Bool=false, alg=ExactIntegrate(), Z = integrate(lattice, A, B..., alg=alg))

Return the Green's function ⟨aᵢ bⱼ⟩ on the L-shaped Kadanoff-Baym contour
band: the band to calculate the Green's function
c1: whether to take the conjugate of the first GV a
c2: whether to take the conjugate of the second GV b
b1: the branch of the first GV (can be :+, :- or :τ)
b2: the branch of the second GV (can be :+, :- or :τ)
The other keywords are the same as gf
"""
function Gm(lattice::MixedGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; 
            b1::Symbol, b2::Symbol, c1::Bool=true, c2::Bool=false, band::Union{Int, Tuple{Int, Int}}=1, 
            alg::IntegrationAlgorithm=ExactIntegrate(), 
            Z::Number = integrate(lattice, A, B..., alg=alg))
    if isa(band, Int)
        band = (band, band)
    end
    band1, band2 = band    
    a, b = ContourIndex(i, conj=c1, branch=b1, band=band1), ContourIndex(j, conj=c2, branch=b2, band=band2)
    return gf(lattice, (a, b), A, B...; alg=alg, Z=Z)
end

########************************#######

"""
    parallel_Gτ(lattice::ImagGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
        band::Int=1, c1::Bool=false, c2::Bool=true, alg=ExactIntegrate(), Z = integrate(lattice, A, B..., alg=alg))

The same as Gτ, but use multi-threading
"""
function parallel_Gτ(lattice::ImagGrassmannLattice, i::Int, j::Int, A::Vector{<:GrassmannMPS}, B::GrassmannMPS...; 
            band::Union{Int, Tuple{Int, Int}}=1, c1::Bool=false, c2::Bool=true, alg::IntegrationAlgorithm=ExactIntegrate(), Z::Real = parallel_integrate(lattice, A, B..., alg=alg) )
    if isa(band, Int)
        band = (band, band)
    end
    band1, band2 = band    
    pos1, pos2 = index(lattice, i, conj=c1, band=band1), index(lattice, j, conj=c2, band=band2)
    t = GTerm(pos1, pos2, coeff=1)
    A2 = _mult_A(t, A)
    return parallel_integrate(lattice, A2, B..., alg=alg)/Z
end
parallel_Gτ(lattice::ImagGrassmannLattice, i::Int, A::Vector{<:GrassmannMPS}, B::GrassmannMPS...; kwargs...) = parallel_Gτ(lattice, i, 1, A, B...; kwargs...)

function Gτ(lattice::ImagGrassmannLattice1Order, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...;
            band::Int=1, alg::IntegrationAlgorithm=ExactIntegrate(), Z::Real = integrate(lattice, A, B..., alg=alg))
    g = zeros(Float64, lattice.k)
    for i in 1:lattice.k-1
        g[i] = Gτ(lattice, i, A, B...; band=band, alg=alg, Z=Z)
    end
    g[end] = 1 - g[1]
    return g
end

function parallel_Gτ(lattice::ImagGrassmannLattice1Order, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; 
                     band::Int=1, alg::IntegrationAlgorithm=ExactIntegrate(), Z::Real = parallel_integrate(lattice, A, B..., alg=alg))
    g = zeros(Float64, lattice.k)
    function _f(ist::Int, ifn::Int, r, lat, args...; kws...)
        for i in ist:ifn
            r[i] = Gτ(lat, i, args...; kws...)
        end
    end
    parallel_run(lattice.k-1, Threads.nthreads(), _f, g, lattice, A, B...; Z=Z, band=band, alg=alg)
    g[end] = 1 - g[1]
    return g
end

function Gt(lattice::MixedGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; b1::Symbol, b2::Symbol, kwargs...)
    (b1 in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
    (b2 in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
    return Gm(lattice, i, j, A, B...; b1=b1, b2=b2, kwargs...)
end

function Gτ(lattice::MixedGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; c1::Bool=false, c2::Bool=true, kwargs...)
    return Gm(lattice, i, j, A, B...; b1=:τ, b2=:τ, c1=c1, c2=c2, kwargs...)
end
Gτ(lattice::MixedGrassmannLattice, i::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; kwargs...) = Gτ(lattice, i, 1, A, B...; kwargs...)

"""
    greater(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
            band::Int=1, alg=ExactIntegrate(), Z = integrate(lattice, A, B..., alg=alg))

Return the greater Green's function ⟨aᵢ bⱼ⟩ on the Keldysh or the Kadanoff-Baym contour
"""
function greater(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; 
                 band::Union{Int, Tuple{Int, Int}}=1, alg::IntegrationAlgorithm=ExactIntegrate(), Z::Number = integrate(lattice, A, B..., alg=alg))
    @assert i >= j
    return Gt(lattice, i, j, A, B...; b1=:+, b2=:+, c1=false, c2=true, band=band, alg=alg, Z=Z)
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
    return Gt(lattice, i, j, A, B...; b1=:-, b2=:+, c1=true, c2=false, band=band, alg=alg, Z=Z)
end
lesser(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; kwargs...) = lesser(lattice, 1, i, A, B...; kwargs...)


# other real-time observables




