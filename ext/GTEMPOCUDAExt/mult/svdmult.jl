
GTEMPO.mult(x::GrassmannMPS, y::GrassmannMPS, alg::CuSVDCompression) = cumult(x, y, trunc=alg.trunc, verbosity=alg.verbosity)
GTEMPO.mult(x::GrassmannMPS, y::GrassmannMPS, alg::CuDMRGMultAlgorithm) = cu_iterativemult(x, y, alg)

GTEMPO.mult!(x::GrassmannMPS, y::GrassmannMPS, alg::CuSVDCompression) = cumult!(x, y, trunc=alg.trunc, verbosity=alg.verbosity)
GTEMPO.mult!(x::GrassmannMPS, y::GrassmannMPS, alg::CuDMRGMultAlgorithm) = cu_iterativemult(x, y, alg)




function cumult!(x::GrassmannMPS, y::GrassmannMPS; trunc::TruncationScheme=DefaultTruncation, verbosity::Int=0, useHCache::Bool=DefaultUseCache)
    (length(x) == length(y)) || throw(DimensionMismatch())
    left = GrassmannTensorMap(isomorphism(scalartype(x), fuse(space_l(x), space_l(y)), space_l(x) ⊗ space_l(y) ))
    tmp5 = g_fuse(_mult_site(x[1], y[1]), 3)
    @tensor tmp4[1,4;5,6] := left[1,2,3] * tmp5[2,3,4,5,6]
    
    stype = scalartype(tmp5)
    ttype = TensorMap{stype, Z2Space, 2, 1, Vector{stype}}
    res = useHCache ? CachedVector{ttype}(undef, length(x)) : Vector{ttype}(undef, length(x))
    empty_tensor = empty(ttype)

	tmp4 = tocu(tmp4)
    for i in 1:length(x)-1
        q, r = leftorth!(tmp4, alg = QR())
        res[i] = fromcu(get_data(q))
        x[i] = empty_tensor
        _renormalize!(x, get_data(r), false)

        @tensor tmp1[1,5,4;2] := r[1,2,3] * tocu(GrassmannTensorMap(y[i+1]))[3,4,5]
        @tensor tmp2[1,3,5;6,2] := tmp1[1,2,3,4] * tocu(GrassmannTensorMap(x[i+1]))[4,5,6]
        tmp4 = g_fuse(tmp2, 2)

        # tmp5 = g_fuse(_mult_site(x[i+1], y[i+1]), 3)
        # @tensor tmp4′[1,4;5,6] := r[1,2,3] * tmp5[2,3,4,5,6]
        # println("i=", i, " err=", norm(tmp4.data - tmp4′.data))
    end
    tmp4 = fromcu(tmp4)
    @tensor tmp[1,2;5] := tmp4[1,2,3,4] * conj(left[5,3,4])
    res[end] = get_data(tmp)
    x[end] = empty_tensor
    x′ = GrassmannMPS(res, copy(x.svectors), scaling(x))

    _cu_rightorth!(x′, SVD(), trunc, false, verbosity)
    destory_copy!(x.data, x′.data)
    copy!(x.svectors, x′.svectors)
    setscaling!(x, scaling(x′) * scaling(y))
    return x
end
# function cumult!(x::GrassmannMPS, y::GrassmannMPS; trunc::TruncationScheme=DefaultTruncation, verbosity::Int=0)
#     (length(x) == length(y)) || throw(DimensionMismatch())
#     left = GrassmannTensorMap(isomorphism(scalartype(x), fuse(space_l(x), space_l(y)), space_l(x) ⊗ space_l(y) ))
#     tmp5 = g_fuse(_mult_site(x[1], y[1]), 3)
#     @tensor tmp4[1,4;5,6] := left[1,2,3] * tmp5[2,3,4,5,6]
# 	tmp4 = tocu(tmp4)
#     for i in 1:length(x)-1
#         q, r = leftorth!(tmp4, alg = QR())
#         x[i] = fromcu(get_data(q))
#         # _renormalize!(x, get_data(r), false)

#         @tensor tmp1[1,5,4;2] := r[1,2,3] * tocu(GrassmannTensorMap(y[i+1]))[3,4,5]
#         @tensor tmp2[1,3,5;6,2] := tmp1[1,2,3,4] * tocu(GrassmannTensorMap(x[i+1]))[4,5,6]
#         tmp4 = g_fuse(tmp2, 2)

#         # tmp5 = g_fuse(_mult_site(x[i+1], y[i+1]), 3)
#         # @tensor tmp4′[1,4;5,6] := r[1,2,3] * tmp5[2,3,4,5,6]
#         # println("i=", i, " err=", norm(tmp4.data - tmp4′.data))
#     end
# 	tmp4 = fromcu(tmp4)
#     @tensor tmp[1,2;5] := tmp4[1,2,3,4] * conj(left[5,3,4])
#     x[end] = get_data(tmp)
#     _cu_rightorth!(x, SVD(), trunc, false, verbosity)
#     setscaling!(x, scaling(x) * scaling(y))
#     return x
# end
cumult(x::GrassmannMPS, y::GrassmannMPS; kwargs...) = cumult!(copy(x), y; kwargs...)

