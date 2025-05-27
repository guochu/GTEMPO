# # # the convention here is different from DMRG!!!
TK.dot(psiA::FockMPS, psiB::FockMPS) = _dot(psiA, psiB) * (scaling(psiA) * scaling(psiB))^length(psiA)
function TK.norm(psi::FockMPS) 
	a = real(_dot(psi, psi))
    a = (abs(a) >= 1.0e-14) ? a : zero(a)
	return sqrt(a) * scaling(psi)^(length(psi))
end
DMRG.distance(a::FockMPS, b::FockMPS) = DMRG._distance(a, b)
DMRG.distance2(a::FockMPS, b::FockMPS) = DMRG._distance2(a, b)


function TK.lmul!(f::Number, psi::FockMPS)
    if !isempty(psi)
        psi[1] *= f
    end
    _renormalize!(psi, psi[1], false)
    return psi
end

Base.:*(psi::FockMPS, f::Number) = lmul!(f, copy(psi))
Base.:*(f::Number, psi::FockMPS) = psi * f
Base.:/(psi::FockMPS, f::Number) = psi * (1/f)
Base.:(-)(psi::FockMPS) = (-1) * psi

function _dot(psiA::FockMPS, psiB::FockMPS) 
    (length(psiA) == length(psiB)) || throw(ArgumentError("dimension mismatch"))
    hold = l_LL(psiA, psiB)
    for i in 1:length(psiA)
        hold = updateleft(hold, psiA[i], psiB[i])
    end
    return tr(hold)

end


# the reuslt is also a GrassmannMPS
function Base.:*(x::FockMPS, y::FockMPS)
    (length(x) == length(y)) || throw(DimensionMismatch())
    r = [n_fuse(_mult_site_n(x[i], y[i]), 3) for i in 1:length(x)]
    return FockMPS([tie(rj, (2,1,2)) for rj in r], scaling=scaling(x)*scaling(y))
end

function Base.:+(x::FockMPS, y::FockMPS) 
    (length(x) == length(y)) || throw(DimensionMismatch())
    @assert !isempty(x)
    scaling_x = scaling(x)
    scaling_y = scaling(y)
    (length(x) == 1) && return FockMPS([scaling_x * x[1] + scaling_y * y[1]])

    L = length(x)
    T = promote_type(scalartype(x), scalartype(y))
    r = Vector{Array{T, 3}}(undef, L)
    r[1] = cat(scaling_x*x[1], scaling_y*y[1], dims=3)
    r[L] = cat(scaling_x*x[L], scaling_y*y[L], dims=1)
    for i in 2:L-1
        r[i] = cat(scaling_x*x[i], scaling_y*y[i], dims=(1,3))
    end
    return FockMPS(r)
end
Base.:-(x::FockMPS, y::FockMPS) = x + (-y)


function _permute!(x::AbstractFiniteMPS, perm::Vector{Int}; trunc::TruncationScheme=DefaultIntegrationTruncation)
    @assert length(x) == length(perm)
    if svectors_uninitialized(x)
        canonicalize!(x, alg=Orthogonalize(trunc=trunc, normalize=false))
    end
    p = CoxeterDecomposition(Permutation(perm))
    for i in p.terms
        easy_swap!(x, i, trunc=trunc)
    end
    return x
end
TK.permute!(x::FockMPS, perm::Vector; kwargs...) = _permute!(x, perm; kwargs...)
TK.permute(x::AbstractFiniteMPS, perm::Vector{Int}; kwargs...) = permute!(deepcopy(x), perm; kwargs...)

function naive_permute!(x::AbstractFiniteMPS, perm::Vector{Int}; trunc::TruncationScheme=DefaultIntegrationTruncation)
    @assert length(x) == length(perm)
    p = CoxeterDecomposition(Permutation(perm))
    for i in p.terms
        naive_swap!(x, i, trunc=trunc)
    end
    return x
end
naive_permute(x::AbstractFiniteMPS, perm::Vector{Int}; kwargs...) = naive_permute!(copy(x), perm; kwargs...)

function _mult_site_n(xj, yj)
    @tensor r[1,4,2,5;3,6] := xj[1,2,3] * yj[4,5,6]
    return r
end

