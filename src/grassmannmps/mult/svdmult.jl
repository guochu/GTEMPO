# standard multiplication and truncation using SVD
"""
    mult!(x::GrassmannMPS, y::GrassmannMPS; trunc::TruncationScheme)

Multiplication of two GMPS x and y, and the result is stored in x
"""
function mult!(x::GrassmannMPS, y::GrassmannMPS; trunc::TruncationScheme=DefaultTruncation, verbosity::Int=0)
    (length(x) == length(y)) || throw(DimensionMismatch())
    left = isomorphism(scalartype(x), fuse(space_l(x), space_l(y)), space_l(x) ⊗ space_l(y) )
    tmp5 = g_fuse(_mult_site(x[1], y[1]), 3)
    @grassmann tmp4[1,4;5,6] := left[1,2,3] * tmp5[2,3,4,5,6]
    for i in 1:length(x)-1
        q, r = leftorth!(tmp4, alg = QR())
        x[i] = q
        tmp5 = g_fuse(_mult_site(x[i+1], y[i+1]), 3)
        @grassmann tmp4[1,4;5,6] := r[1,2,3] * tmp5[2,3,4,5,6]
    end
    @tensor tmp[1,2;5] := tmp4[1,2,3,4] * conj(left[5,3,4])
    x[end] = tmp
    _rightorth!(x, SVD(), trunc, false, verbosity)
    setscaling!(x, scaling(x) * scaling(y))
    return x
end
mult(x::GrassmannMPS, y::GrassmannMPS; kwargs...) = mult!(copy(x), y; kwargs...)
