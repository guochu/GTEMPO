function mult!(x::FockMPS, y::FockMPS; trunc::TruncationScheme=DMRG.DefaultTruncation, verbosity::Int=0)
    (length(x) == length(y)) || throw(DimensionMismatch())
    T = promote_type(scalartype(x), scalartype(y))
    L = length(x)
    z = Vector{Array{T, 3}}(undef, L)
    tmp = tie(n_fuse(_mult_site_n(x[1], y[1]), 3), (2,1,1,1))
    q, r = tqr!(tmp, (1, 2), (3,4))
    x[1] = q
    for i in 2:L-1
        @tensor tmp[1,4,6,5,7] := r[1,2,3] * x[i][2,4,5] * y[i][3,6,7]
        tmp2 = n_fuse(tmp, 2)
        q, r = tqr!(tmp2, (1, 2), (3,4))
        x[i] = q
        _renormalize!(x, r, false)
    end
    @tensor tmp[1,4,6,5,7] := r[1,2,3] * x[L][2,4,5] * y[L][3,6,7]
    x[end] = tie(n_fuse(tmp, 2), (1,2,1))
    _rightorth!(x, SVD(), trunc, false, verbosity)
    setscaling!(x, scaling(x) * scaling(y))
    return x
end
mult(x::FockMPS, y::FockMPS; kwargs...) = mult!(copy(x), y; kwargs...)