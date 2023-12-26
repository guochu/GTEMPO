# ###--------------imaginary time----------------
# function bmps_gf(lattice::ImagGrassmannLattice, i::Int, A::GrassmannMPS, B::GrassmannMPS...; band::Int=1, kwargs...)
# 	pos1, pos2 = index(lattice, i, conj=false, band=band), index(lattice, 1, conj=true, band=band)
# 	t = GTerm(pos1, pos2, coeff=1)
# 	A2 = t * A
#     As = [A2, B...]
#     left = _bmps_integrate_left(lattice, As, i; kwargs...)
#     right = _bmps_integrate_right(lattice, As, 1; kwargs...)
#     center = div(1+i, 2)
#     r = _integrate_center(lattice, As, left, i, right, 1, center; kwargs...)
#     As = [A, B...]
#     Z = _integrate_center(lattice, As, left, i, right, 1, center; kwargs...)
# 	return r/Z
# end

# function bmps_gf2(lattice::ImagGrassmannLattice, i::Int, A::GrassmannMPS, B::GrassmannMPS...; 
#                     band::Int=1, trunc::TruncationScheme=DefaultIntegrationTruncation, Z::Real=bmps_integrate(lattice, A, B..., trunc=trunc))
#     pos1, pos2 = index(lattice, i, conj=false, band=band), index(lattice, 1, conj=true, band=band)
#     t = GTerm(pos1, pos2, coeff=1)
#     A2 = t * A
#     return bmps_integrate(lattice, A2, B..., trunc=trunc)/Z
# end

# function bmps_gf(lattice::ImagGrassmannLattice, i::Int, A::Vector{<:GrassmannMPS}, B::GrassmannMPS...; band::Int=1, kwargs...)
#     pos1, pos2 = index(lattice, i, conj=false, band=band), index(lattice, 1, conj=true, band=band)
#     t = GTerm(pos1, pos2, coeff=1)
#     rs = zeros(eltype(lattice), length(A))
#     Zs = zeros(eltype(lattice), length(A))
#     center = div(1+i, 2)
#     for j in 1:length(A)
#         Aj = A[j]
#         Aj2 = t * Aj
#         As = [Aj2, B...]
#         left = _bmps_integrate_left(lattice, As, i; kwargs...)
#         right = _bmps_integrate_right(lattice, As, 1; kwargs...)
#         rs[j] = _integrate_center(lattice, As, left, i, right, 1, center; kwargs...)
#         As = [Aj, B...]
#         Zs[j] = _integrate_center(lattice, As, left, i, right, 1, center; kwargs...)
#     end
#     return sum(rs)/sum(Zs)
# end

# function bmps_gf2(lattice::ImagGrassmannLattice, i::Int, A::Vector{<:GrassmannMPS}, B::GrassmannMPS...; 
#                     band::Int=1, trunc::TruncationScheme=DefaultIntegrationTruncation, Z::Real=bmps_integrate(lattice, A, B..., trunc=trunc))
#     pos1, pos2 = index(lattice, i, conj=false, band=band), index(lattice, 1, conj=true, band=band)
#     t = GTerm(pos1, pos2, coeff=1)
#     A2 = [t * item for item in A]
#     return bmps_integrate(lattice, A2, B..., trunc=trunc)/Z
# end

# ###--------------real time 1 order----------------
# function bmps_gf(lattice::RealGrassmannLattice, i::Int, j::Int, A::GrassmannMPS, B::GrassmannMPS...; 
#                 f1::Bool, f2::Bool, c1::Bool=true, c2::Bool=false, band::Int=1, kwargs...)
#     pos1, pos2 = index(lattice, i, conj=c1, forward=f1, band=band), index(lattice, j, conj=c2, forward=f2, band=band)
#     t = GTerm(pos1, pos2, coeff=1)
#     A2 = t * A

#     a, b = max(i, j), min(i, j)
#     As = [A2, B...]
#     left = _bmps_integrate_left(lattice, As, a; kwargs...)
#     right = _bmps_integrate_right(lattice, As, b; kwargs...)
#     center = div(a+b, 2)
#     r = _integrate_center(lattice, As, left, a, right, b, center; kwargs...)
#     As = [A, B...]
#     Z = _integrate_center(lattice, As, left, a, right, b, center; kwargs...)
#     return r/Z
# end

# function bmps_gf2(lattice::RealGrassmannLattice, i::Int, j::Int, A::GrassmannMPS, B::GrassmannMPS...; 
#                 f1::Bool, f2::Bool, c1::Bool=true, c2::Bool=false, band::Int=1, 
#                 trunc::TruncationScheme=DefaultIntegrationTruncation, Z::Number=bmps_integrate(lattice, A, B..., trunc=trunc))
#     pos1, pos2 = index(lattice, i, conj=c1, forward=f1, band=band), index(lattice, j, conj=c2, forward=f2, band=band)
#     t = GTerm(pos1, pos2, coeff=1)
#     A2 = t * A
#     g = bmps_integrate(lattice, A2, B..., trunc=trunc)/Z
#     return g
# end

# function bmps_gf(lattice::RealGrassmannLattice, i::Int, j::Int, A::Vector{<:GrassmannMPS}, B::GrassmannMPS...; 
#                 f1::Bool, f2::Bool, c1::Bool=true, c2::Bool=false, band::Int=1, kwargs...)
#     pos1, pos2 = index(lattice, i, conj=c1, forward=f1, band=band), index(lattice, j, conj=c2, forward=f2, band=band)
#     t = GTerm(pos1, pos2, coeff=1)
#     rs = zeros(eltype(lattice), length(A))
#     Zs = zeros(eltype(lattice), length(A))
#     a, b = max(i, j), min(i, j)
#     center = div(a+b, 2)
#     for l in 1:length(A)
#         Aj = A[l]
#         Aj2 = t * Aj
#         As = [Aj2, B...]
#         left = _bmps_integrate_left(lattice, As, a; kwargs...)
#         right = _bmps_integrate_right(lattice, As, b; kwargs...)
#         rs[l] = _integrate_center(lattice, As, left, a, right, b, center; kwargs...)
#         As = [Aj, B...]
#         Zs[l] = _integrate_center(lattice, As, left, a, right, b, center; kwargs...)
#     end
#     return sum(rs)/sum(Zs)
# end

# function bmps_gf2(lattice::RealGrassmannLattice, i::Int, j::Int, A::Vector{<:GrassmannMPS}, B::GrassmannMPS...; 
#                 f1::Bool, f2::Bool, c1::Bool=true, c2::Bool=false, band::Int=1,
#                 trunc::TruncationScheme=DefaultIntegrationTruncation, Z::Number=bmps_integrate(lattice, A, B..., trunc=trunc))
#     pos1, pos2 = index(lattice, i, conj=c1, forward=f1, band=band), index(lattice, j, conj=c2, forward=f2, band=band)
#     t = GTerm(pos1, pos2, coeff=1)
#     A2 = [t * item for item in A]
#     g = bmps_integrate(lattice, A2, B..., trunc=trunc)/Z
#     return g
# end
