
DMRG.space_l(a::GrassmannMPS) = space_l(a[1])
DMRG.space_r(a::GrassmannMPS) = space_r(a[end])
# the convention here is different from DMRG!!!
DMRG.r_RR(a::GrassmannMPS, b::GrassmannMPS) = DMRG.loose_isometry(Matrix{promote_type(scalartype(a), scalartype(b))}, space_r(a), space_r(b))
DMRG.l_LL(a::GrassmannMPS, b::GrassmannMPS) = DMRG.loose_isometry(Matrix{promote_type(scalartype(a), scalartype(b))}, space_l(a), space_l(b))


TK.dot(psiA::GrassmannMPS, psiB::GrassmannMPS) = _dot(psiA, psiB) * (scaling(psiA) * scaling(psiB))^length(psiA)
function TK.norm(psi::GrassmannMPS) 
	a = real(_dot(psi, psi))
    # println("a is ", a)
	return sqrt(a) * scaling(psi)^(length(psi))
end
DMRG.distance(a::GrassmannMPS, b::GrassmannMPS) = DMRG._distance(a, b)
DMRG.distance2(a::GrassmannMPS, b::GrassmannMPS) = DMRG._distance2(a, b)

function _dot(psiA::GrassmannMPS, psiB::GrassmannMPS)
	(length(psiA) == length(psiB)) || throw(DimensionMismatch())
    hold = l_LL(psiA, psiB)
    for i in 1:length(psiA)
        hold = _updateleft(hold, psiA[i], psiB[i])
    end
    return tr(hold) 
end

# the convention here is different from DMRG!!!
function _updateleft(hold::MPSBondTensor, mpsAj::MPSTensor, mpsBj::MPSTensor) 
    @tensor out[4; 5] := hold[1,2] * conj(mpsAj[1, 3, 4]) * mpsBj[2,3,5]
end

function Base.:*(h::MPO, psi::GrassmannMPS)
    @assert !isempty(h)
    (length(h) == length(psi)) || throw(DimensionMismatch())
    r = [@tensor tmp[-1 -2; -3 -4 -5] := a[-1, -3, -4, 1] * b[-2, 1, -5] for (a, b) in zip(h.data, psi.data)]
    left = isomorphism(fuse(space_l(h), space_l(psi)), space_l(h) ⊗ space_l(psi))
    fusion_ts = [isomorphism(space(item, 4)' ⊗ space(item, 5)', fuse(space(item, 4)', space(item, 5)')) for item in r]
    @tensor tmp[-1 -2; -3] := left[-1, 1, 2] * r[1][1,2,-2,3,4] * fusion_ts[1][3,4,-3]
    mpstensors = Vector{typeof(tmp)}(undef, length(h))
    mpstensors[1] = tmp
    for i in 2:length(h)
        @tensor tmp[-1 -2; -3] := conj(fusion_ts[i-1][1,2,-1]) * r[i][1,2,-2,3,4] * fusion_ts[i][3,4,-3]
        mpstensors[i] = tmp
    end
    return GrassmannMPS(mpstensors, scaling=scaling(psi))
end


# the reuslt is also a GrassmannMPS
function Base.:*(x::GrassmannMPS, y::GrassmannMPS)
    (length(x) == length(y)) || throw(DimensionMismatch())
    out = [_fuse_physical(_mult_site(x[i], y[i])) for i in 1:length(x)]
    fusers = DMRG.PeriodicArray([isomorphism(space(item, 4)' ⊗ space(item, 5)', fuse(space(item, 4), space(item, 5)) ) for item in out])
    return GrassmannMPS([@tensor tmp[3,4;7] := conj(fusers[i-1][1,2,3]) * out[i][1,2,4,5,6] * fusers[i][5,6,7] for i in 1:length(x)], scaling=scaling(x) * scaling(y))
end
# function mult(x::GrassmannMPS, y::GrassmannMPS; trunc::TruncationScheme=DMRG.DefaultTruncation)
#   z = x * y
#   canonicalize!(z, alg=Orthogonalize(SVD(), trunc))
#   return z
# end
# function mult(x::GrassmannMPS, y::GrassmannMPS; trunc::TruncationScheme=DMRG.DefaultTruncation)
#     (length(x) == length(y)) || throw(DimensionMismatch())
#     A = DMRG.mpstensortype(spacetype(x), promote_type(eltype(x), eltype(y)))
#     z = Vector{A}(undef, length(x))
#     left = isomorphism( fuse(space_l(x), space_l(y)), space_l(x) ⊗ space_l(y) )
#     tmp5 = _fuse_physical(_mult_site(x[1], y[1]))
#     @tensor tmp4[1,4;5,6] := left[1,2,3] * tmp5[2,3,4,5,6]
#     for i in 1:length(x)-1
#         q, r = leftorth!(tmp4, alg = QR())
#         z[i] = q
#         tmp5 = _fuse_physical(_mult_site(x[i+1], y[i+1]))
#         @tensor tmp4[1,4;5,6] := r[1,2,3] * tmp5[2,3,4,5,6]
#     end
#     z[length(x)] = @tensor tmp[1,2;5] := tmp4[1,2,3,4] * conj(left[5,3,4])
#     out = GrassmannMPS(z)
#     _rightorth!(out, SVD(), trunc)
#     return _rescaling!(out)
# end

function mult!(x::GrassmannMPS, y::GrassmannMPS; trunc::TruncationScheme=DMRG.DefaultTruncation)
    (length(x) == length(y)) || throw(DimensionMismatch())
    A = mpstensortype(spacetype(x), promote_type(scalartype(x), scalartype(y)))
    left = isomorphism( fuse(space_l(x), space_l(y)), space_l(x) ⊗ space_l(y) )
    tmp5 = _fuse_physical(_mult_site(x[1], y[1]))
    @tensor tmp4[1,4;5,6] := left[1,2,3] * tmp5[2,3,4,5,6]
    for i in 1:length(x)-1
        q, r = leftorth!(tmp4, alg = QR())
        x[i] = q
        tmp5 = _fuse_physical(_mult_site(x[i+1], y[i+1]))
        @tensor tmp4[1,4;5,6] := r[1,2,3] * tmp5[2,3,4,5,6]
    end
    x[end] = @tensor tmp[1,2;5] := tmp4[1,2,3,4] * conj(left[5,3,4])
    _rightorth!(x, SVD(), trunc)
    return _rescaling!(x)
end
mult(x::GrassmannMPS, y::GrassmannMPS; kwargs...) = mult!(copy(x), y; kwargs...)

function Base.:+(x::GrassmannMPS, y::GrassmannMPS) 
    (length(x) == length(y)) || throw(DimensionMismatch())
    @assert !isempty(x)
    (length(x) == 1) && return GrassmannMPS([scaling(x) * x[1] + scaling(y) * y[1]])

    T = promote_type(scalartype(x), scalartype(y))
    A = mpstensortype(spacetype(x), T)
    embedders = [right_embedders(T, space_r(aj)', space_r(bj)') for (aj, bj) in zip(x.data, y.data)]
    r = A[]
    for i in 1:length(x)
        if i == 1
            @tensor m1[-1 -2; -3] := x[i][-1,-2,2] * embedders[i][1][2, -3]
            @tensor m2[-1 -2; -3] := y[i][-1,-2,2] * embedders[i][2][2, -3]
        elseif i == length(x)
            @tensor m1[-1 -2; -3] := (embedders[i-1][1])'[-1, 1] * x[i][1,-2,-3] 
            @tensor m2[-1 -2; -3] := (embedders[i-1][2])'[-1, 1] * y[i][1,-2,-3] 
        else          
            @tensor m1[-1 -2; -3] := (embedders[i-1][1])'[-1, 1] * x[i][1,-2,2] * embedders[i][1][2, -3]
            @tensor m2[-1 -2; -3] := (embedders[i-1][2])'[-1, 1] * y[i][1,-2,2] * embedders[i][2][2, -3]
        end
        push!(r, scaling(x) * m1 + scaling(y) * m2)
    end
    return GrassmannMPS(r)
end

function right_embedders(::Type{T}, a::S...) where {T <: Number, S <: ElementarySpace}
    V = ⊕(a...) 
    ts = [TensorMap(zeros, T, aj, V) for aj in a]
    for c in sectors(V)
        n = 0
        for i in 1:length(ts)
            ni = dim(a[i], c)
            block(ts[i], c)[:, (n+1):(n+ni)] .= Diagonal( ones(ni) )
            n += ni
        end
    end
    return ts
end

function _permute!(x::GrassmannMPS, perm::Vector{Int}; trunc::TruncationScheme=DMRG.DefaultTruncation)
    @assert length(x) == length(perm)
    p = CoxeterDecomposition(Permutation(perm))
    for i in p.terms
        easy_swap!(x, i, trunc=trunc)
    end
    return x
end
_permute(x::GrassmannMPS, perm::Vector{Int}; kwargs...) = _permute!(copy(x), perm; kwargs...)
