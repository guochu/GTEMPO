# standard multiplication and truncation using SVD
function mult2!(x::GrassmannMPS, y::GrassmannMPS; trunc::TruncationScheme=DefaultTruncation, verbosity::Int=0)
    (length(x) == length(y)) || throw(DimensionMismatch())
    left = GrassmannTensorMap(isomorphism(scalartype(x), fuse(space_l(x), space_l(y)), space_l(x) ⊗ space_l(y) ))
    tmp5 = g_fuse(_mult_site(x[1], y[1]), 3)
    @tensor tmp4[1,4;5,6] := left[1,2,3] * tmp5[2,3,4,5,6]
    for i in 1:length(x)-1
        q, r = leftorth!(tmp4, alg = QR())
        x[i] = get_data(q)
        tmp5 = g_fuse(_mult_site(x[i+1], y[i+1]), 3)
        @tensor tmp4[1,4;5,6] := r[1,2,3] * tmp5[2,3,4,5,6]
    end
    @tensor tmp[1,2;5] := tmp4[1,2,3,4] * conj(left[5,3,4])
    x[end] = get_data(tmp)
    _rightorth!(x, SVD(), trunc, false, verbosity)
    setscaling!(x, scaling(x) * scaling(y))
    return x
end
mult2(x::GrassmannMPS, y::GrassmannMPS; kwargs...) = mult2!(copy(x), y; kwargs...)

"""
    mult!(x::GrassmannMPS, y::GrassmannMPS; trunc::TruncationScheme)

Multiplication of two GMPS x and y, and the result is stored in x
"""
function mult!(x::GrassmannMPS, y::GrassmannMPS; trunc::TruncationScheme=DefaultTruncation, verbosity::Int=0, useHCache::Bool=DefaultUseCache)
    (length(x) == length(y)) || throw(DimensionMismatch())
    left = GrassmannTensorMap(isomorphism(scalartype(x), fuse(space_l(x), space_l(y)), space_l(x) ⊗ space_l(y) ))
    tmp5 = g_fuse(_mult_site(x[1], y[1]), 3)
    @tensor tmp4[1,4;5,6] := left[1,2,3] * tmp5[2,3,4,5,6]
    
    stype = scalartype(tmp5)
    ttype = TensorMap{stype, Z2Space, 2, 1, Vector{stype}}
    res = useHCache ? CachedVector{ttype}(undef, length(x)) : Vector{ttype}(undef, length(x))
    empty_tensor = empty(ttype)

    for i in 1:length(x)-1
        q, r = leftorth!(tmp4, alg = QR())
        res[i] = get_data(q)
        x[i] = empty_tensor
        _renormalize!(x, get_data(r), false)

        @tensor tmp1[1,5,4;2] := r[1,2,3] * GrassmannTensorMap(y[i+1])[3,4,5]
        @tensor tmp2[1,3,5;6,2] := tmp1[1,2,3,4] * GrassmannTensorMap(x[i+1])[4,5,6]
        tmp4 = g_fuse(tmp2, 2)

        # tmp5 = g_fuse(_mult_site(x[i+1], y[i+1]), 3)
        # @tensor tmp4′[1,4;5,6] := r[1,2,3] * tmp5[2,3,4,5,6]
        # println("i=", i, " err=", norm(tmp4.data - tmp4′.data))
    end
    @tensor tmp[1,2;5] := tmp4[1,2,3,4] * conj(left[5,3,4])
    res[end] = get_data(tmp)
    x[end] = empty_tensor
    x′ = GrassmannMPS(res, copy(x.svectors), scaling(x))

    _rightorth!(x′, SVD(), trunc, false, verbosity)
    destory_copy!(x.data, x′.data)
    copy!(x.svectors, x′.svectors)
    setscaling!(x, scaling(x′) * scaling(y))
    return x
end
# function mult!(x::GrassmannMPS, y::GrassmannMPS; trunc::TruncationScheme=DefaultTruncation, verbosity::Int=0)
#     (length(x) == length(y)) || throw(DimensionMismatch())
#     left = GrassmannTensorMap(isomorphism(scalartype(x), fuse(space_l(x), space_l(y)), space_l(x) ⊗ space_l(y) ))
#     tmp5 = g_fuse(_mult_site(x[1], y[1]), 3)
#     @tensor tmp4[1,4;5,6] := left[1,2,3] * tmp5[2,3,4,5,6]
#     for i in 1:length(x)-1
#         q, r = leftorth!(tmp4, alg = QR())
#         x[i] = get_data(q)
#         _renormalize!(x, get_data(r), false)

#         @tensor tmp1[1,5,4;2] := r[1,2,3] * GrassmannTensorMap(y[i+1])[3,4,5]
#         @tensor tmp2[1,3,5;6,2] := tmp1[1,2,3,4] * GrassmannTensorMap(x[i+1])[4,5,6]
#         tmp4 = g_fuse(tmp2, 2)

#         # tmp5 = g_fuse(_mult_site(x[i+1], y[i+1]), 3)
#         # @tensor tmp4′[1,4;5,6] := r[1,2,3] * tmp5[2,3,4,5,6]
#         # println("i=", i, " err=", norm(tmp4.data - tmp4′.data))
#     end
#     @tensor tmp[1,2;5] := tmp4[1,2,3,4] * conj(left[5,3,4])
#     x[end] = get_data(tmp)
#     _rightorth!(x, SVD(), trunc, false, verbosity)
#     setscaling!(x, scaling(x) * scaling(y))
#     return x
# end
mult(x::GrassmannMPS, y::GrassmannMPS; kwargs...) = mult!(copy(x), y; kwargs...)
