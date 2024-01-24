# standard multiplication and truncation using SVD
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
