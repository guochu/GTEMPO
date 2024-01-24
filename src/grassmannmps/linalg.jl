# # # the convention here is different from DMRG!!!
TK.dot(psiA::GrassmannMPS, psiB::GrassmannMPS) = _dot(psiA, psiB) * (scaling(psiA) * scaling(psiB))^length(psiA)
function TK.norm(psi::GrassmannMPS) 
	a = real(_dot(psi, psi))
    # println("a is ", a)
	return sqrt(a) * scaling(psi)^(length(psi))
end
DMRG.distance(a::AbstractGMPS, b::AbstractGMPS) = DMRG._distance(a, b)
DMRG.distance2(a::AbstractGMPS, b::AbstractGMPS) = DMRG._distance2(a, b)

function TK.lmul!(f::Number, psi::GrassmannMPS)
    if !isempty(psi)
        psi[1] *= f
    end
    _renormalize!(psi, psi[1], false)
    return psi
end

Base.:*(psi::GrassmannMPS, f::Number) = lmul!(f, copy(psi))
Base.:*(f::Number, psi::GrassmannMPS) = psi * f
Base.:/(psi::GrassmannMPS, f::Number) = psi * (1/f)
Base.:(-)(psi::GrassmannMPS) = (-1) * psi

_dot(psiA::GrassmannMPS, psiB::GrassmannMPS) = dot(MPS(psiA.data), MPS(psiB.data))
# 	(length(psiA) == length(psiB)) || throw(DimensionMismatch())
#     hold = l_LL(psiA, psiB)
#     for i in 1:length(psiA)
#         hold = _updateleft(hold, psiA[i], psiB[i])
#     end
#     return tr(hold) 
# end

# # the convention here is different from DMRG!!!
# function _updateleft(hold::MPSBondTensor, mpsAj::MPSTensor, mpsBj::MPSTensor) 
#     @tensor out[4; 5] := hold[1,2] * conj(mpsAj[1, 3, 4]) * mpsBj[2,3,5]
# end

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
    out = [g_fuse(_mult_site(x[i], y[i]), 3) for i in 1:length(x)]
    fusers = PeriodicArray([isomorphism(space(item, 4)' ⊗ space(item, 5)', fuse(space(item, 4), space(item, 5)) ) for item in out])
    return GrassmannMPS([@tensor tmp[3,4;7] := conj(fusers[i-1][1,2,3]) * out[i][1,2,4,5,6] * fusers[i][5,6,7] for i in 1:length(x)], scaling=scaling(x) * scaling(y))
end

function mult2!(x::GrassmannMPS, y::GrassmannMPS; trunc::TruncationScheme=DMRG.DefaultTruncation)
    (length(x) == length(y)) || throw(DimensionMismatch())
    A = mpstensortype(spacetype(x), promote_type(scalartype(x), scalartype(y)))
    left = isomorphism( fuse(space_l(x), space_l(y)), space_l(x) ⊗ space_l(y) )
    tmp5 = g_fuse(_mult_site(x[1], y[1]), 3)
    @tensor tmp4[1,4;5,6] := left[1,2,3] * tmp5[2,3,4,5,6]
    for i in 1:length(x)-1
        q, r = leftorth!(tmp4, alg = QR())
        x[i] = q
        tmp5 = g_fuse(_mult_site(x[i+1], y[i+1]), 3)
        @tensor tmp4[1,4;5,6] := r[1,2,3] * tmp5[2,3,4,5,6]
    end
    x[end] = @tensor tmp[1,2;5] := tmp4[1,2,3,4] * conj(left[5,3,4])
    _rightorth!(x, SVD(), trunc, false)
    return _rescaling!(x)
end
mult2(x::GrassmannMPS, y::GrassmannMPS; kwargs...) = mult2!(copy(x), y; kwargs...)

function mult!(x::GrassmannMPS, y::GrassmannMPS; trunc::TruncationScheme=DMRG.DefaultTruncation)
    (length(x) == length(y)) || throw(DimensionMismatch())
    A = mpstensortype(spacetype(x), promote_type(scalartype(x), scalartype(y)))
    left = isomorphism( fuse(space_l(x), space_l(y)), space_l(x) ⊗ space_l(y) )
    tmp5 = g_fuse(_mult_site(x[1], y[1]), 3)
    @tensor tmp4[1,4;5,6] := left[1,2,3] * tmp5[2,3,4,5,6]
    for i in 1:length(x)-1
        q, r = leftorth!(tmp4, alg = QR())
        x[i] = q
        _renormalize!(x, r, false)
        @tensor tmp1[1,5,4;2] := r[1,2,3] * y[i+1][3,4,5]
        for (f1, f2) in fusiontrees(tmp1)
            coef1 = (isodd(f1.uncoupled[2].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
            coef2 = (isodd(f1.uncoupled[3].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
            coef3 = (isodd(f1.uncoupled[3].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
            # println(coef1, " ", coef2, " ", coef3, " ", coef4, " ", coef5)
            coef = coef1 * coef2 * coef3
            if coef != 1
                lmul!(coef, tmp1[f1, f2])
            end
        end
        @tensor tmp2[1,3,5;6,2] := tmp1[1,2,3,4] * x[i+1][4,5,6]
        for (f1, f2) in fusiontrees(tmp2)
            coef1 = (isodd(f2.uncoupled[2].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
            coef2 = (isodd(f2.uncoupled[2].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
            coef3 = (isodd(f2.uncoupled[2].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
            # println(coef1, " ", coef2, " ", coef3, " ", coef4, " ", coef5)
            coef = coef1 * coef2 * coef3
            if coef != 1
                lmul!(coef, tmp2[f1, f2])
            end
        end
        tmp4 = g_fuse(tmp2, 2)

        # tmp5 = g_fuse(_mult_site(x[i+1], y[i+1]), 3)
        # @tensor tmp4[1,4;5,6] := r[1,2,3] * tmp5[2,3,4,5,6]
    end
    x[end] = @tensor tmp[1,2;5] := tmp4[1,2,3,4] * conj(left[5,3,4])
    _rightorth!(x, SVD(), trunc, false)
    return _rescaling!(x)
end
mult(x::GrassmannMPS, y::GrassmannMPS; kwargs...) = mult!(copy(x), y; kwargs...)

function Base.:+(x::GrassmannMPS, y::GrassmannMPS) 
    (length(x) == length(y)) || throw(DimensionMismatch())
    @assert !isempty(x)
    scaling_x = scaling(x)
    scaling_y = scaling(y)
    (length(x) == 1) && return GrassmannMPS([scaling_x * x[1] + scaling_y * y[1]])

    T = promote_type(scalartype(x), scalartype(y))
    A = mpstensortype(spacetype(x), T)
    embedders = [right_embedders(T, space_r(aj)', space_r(bj)') for (aj, bj) in zip(x.data, y.data)]
    r = A[]
    for i in 1:length(x)
        if i == 1
            @tensor m1[-1 -2; -3] := scaling_x * x[i][-1,-2,2] * embedders[i][1][2, -3]
            @tensor m1[-1 -2; -3] += scaling_y * y[i][-1,-2,2] * embedders[i][2][2, -3]
        elseif i == length(x)
            @tensor m1[-1 -2; -3] := scaling_x * (embedders[i-1][1])'[-1, 1] * x[i][1,-2,-3] 
            @tensor m1[-1 -2; -3] += scaling_y * (embedders[i-1][2])'[-1, 1] * y[i][1,-2,-3] 
        else          
            @tensor m1[-1 -2; -3] := scaling_x * (embedders[i-1][1])'[-1, 1] * x[i][1,-2,2] * embedders[i][1][2, -3]
            @tensor m1[-1 -2; -3] += scaling_y * (embedders[i-1][2])'[-1, 1] * y[i][1,-2,2] * embedders[i][2][2, -3]
        end
        push!(r, m1)
    end
    return GrassmannMPS(r)
end
Base.:-(x::GrassmannMPS, y::GrassmannMPS) = x + (-y)

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

function _permute!(x::AbstractGMPS, perm::Vector{Int}; trunc::TruncationScheme=DMRG.DefaultTruncation)
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
TK.permute!(x::GrassmannMPS, perm::Vector; kwargs...) = _permute!(x, perm; kwargs...)
TK.permute(x::AbstractGMPS, perm::Vector{Int}; kwargs...) = permute!(deepcopy(x), perm; kwargs...)

function _mult_site(xj, yj)
    @tensor r[1,4,2,5;3,6] := xj[1,2,3] * yj[4,5,6]
    for (f1, f2) in fusiontrees(r)
        if isodd(f1.uncoupled[4].n) && isodd(f2.uncoupled[1].n)
            r[f1, f2] .*= -1
        end
    end
    return r
end

