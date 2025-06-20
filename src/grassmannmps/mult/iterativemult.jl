# abstract type IterativeMultInitialGuess end
# struct SVDGuess <: IterativeMultInitialGuess end

abstract type DMRGMultAlgorithm <: DMRGAlgorithm end

const AllowedInitGuesses = (:svd, :pre, :rand)

struct DMRGMult1 <: DMRGMultAlgorithm
    trunc::TruncationDimCutoff 
    maxiter::Int
    tol::Float64 
    initguess::Symbol
    verbosity::Int 
    callback::Function
end
function DMRGMult1(trunc::TruncationDimCutoff; maxiter::Int=5, tol::Float64=1.0e-12, initguess::Symbol=:svd, verbosity::Int=0, callback::Function=Returns(nothing))
    (initguess in AllowedInitGuesses) || throw(ArgumentError("initguess must be one of $(AllowedInitGuesses)"))
    return DMRGMult1(trunc, maxiter, tol, initguess, verbosity, callback)
end 
DMRGMult1(; trunc::TruncationDimCutoff=DefaultITruncation, kwargs...) = DMRGMult1(trunc; kwargs...)
Base.similar(x::DMRGMult1; trunc::TruncationDimCutoff=x.trunc, maxiter::Int=x.maxiter, tol::Float64=x.tol, initguess::Symbol=x.initguess, verbosity::Int=x.verbosity, callback=x.callback) = DMRGMult1(
            trunc=trunc, maxiter=maxiter, tol=tol, initguess=initguess, verbosity=verbosity, callback=callback)


struct DMRGMult2 <: DMRGMultAlgorithm
    trunc::TruncationDimCutoff 
    maxiter::Int
    tol::Float64 
    initguess::Symbol
    verbosity::Int 
end

function DMRGMult2(trunc::TruncationDimCutoff; maxiter::Int=5, tol::Float64=1.0e-12, initguess::Symbol=:svd, verbosity::Int=0)
    (initguess in AllowedInitGuesses) || throw(ArgumentError("initguess must be one of $(AllowedInitGuesses)"))
    return DMRGMult2(trunc, maxiter, tol, initguess, verbosity)
end 
DMRGMult2(; trunc::TruncationDimCutoff=DefaultITruncation, kwargs...) = DMRGMult2(trunc; kwargs...)
Base.similar(x::DMRGMult2; trunc::TruncationDimCutoff=x.trunc, maxiter::Int=x.maxiter, tol::Float64=x.tol, initguess::Symbol=x.initguess, verbosity::Int=x.verbosity) = DMRGMult2(
            trunc=trunc, maxiter=maxiter, tol=tol, initguess=initguess, verbosity=verbosity)


function Base.getproperty(x::DMRGMultAlgorithm, s::Symbol)
    if s == :D
        return x.trunc.D
    elseif s == :ϵ
        return x.trunc.ϵ
    else
        getfield(x, s)
    end
end


# z is the output GMPS
struct GMPSIterativeMultCache{_O, _A, _B, _H} 
    z::_O
    x::_A
    y::_B
    hstorage::_H
end

function mult_cache(z::GrassmannMPS, x::GrassmannMPS, y::GrassmannMPS)
    @assert length(z) == length(x) == length(y)
    # initialize Hstorage
    L = length(z)
    right = ones(scalartype(z), space_r(y)' ⊗ space_r(x)', space_r(z)')
    hstorage = Vector{typeof(right)}(undef, L+1)
    hstorage[L+1] = right
    hstorage[1] = ones( scalartype(z), space_l(z) ⊗ space_l(x)', space_l(y) )
    for i in L:-1:2
        hstorage[i] = updatemultright(hstorage[i+1], z[i], x[i], y[i])
    end
    # i = 1
    # tmp = updatemultright(hstorage[i+1], z[i], x[i], y[i])
    # println(scaling(z), " ", scaling(x), " ", scaling(y))
    # println("dot is ", scalar(tmp))
    # println("true dot is ", _dot(z, x * y))

    return GMPSIterativeMultCache(z, x, y, hstorage)
end

function iterativemult(x::GrassmannMPS, y::GrassmannMPS, alg::DMRGMultAlgorithm)
    if alg.initguess == :svd
        z = _svd_guess(x, y, alg.D)
    elseif alg.initguess == :rand
        z = randomgmps(promote_type(scalartype(x), scalartype(y)), length(x), D=alg.D)
    elseif alg.initguess == :pre
        z = increase_bond!(copy(x), alg.D)
    else
        error("unsupported initguess $(alg.initguess)")
    end
    cache = mult_cache(z, x, y)
    deltas = compute!(cache, alg)
    z = cache.z
    setscaling!(z, scaling(x) * scaling(y))
    _rescaling!(z)
    return z
end

DMRG.compute!(env::GMPSIterativeMultCache, alg::DMRGMultAlgorithm) = iterative_compute!(env, alg)


function iterative_compute!(m, alg)
    kvals = Float64[]
    iter = 0
    delta = 2 * alg.tol
    while (iter < alg.maxiter) && (delta >= alg.tol)
        _kvals = sweep!(m, alg)
        delta = iterative_error_2(_kvals)
        push!(kvals, delta)
        iter += 1
        (alg.verbosity > 1) && println("finish the $iter-th sweep with error $delta", "\n")
    end
    if (alg.verbosity >= 2) && (iter < alg.maxiter)
        println("early converge in $iter-th sweeps with error $delta")
    end
    if (alg.verbosity > 0) && (delta >= alg.tol)
        println("fail to converge, required precision: $(alg.tol), actual precision $delta in $iter sweeps")
    end
    finalize!(m ,alg)
    return kvals
end
iterative_error_2(m::AbstractVector) = std(m) / abs(mean(m))

DMRG.sweep!(m::GMPSIterativeMultCache, alg::DMRGMultAlgorithm) = vcat(leftsweep!(m, alg), rightsweep!(m, alg))

function finalize!(m::GMPSIterativeMultCache, alg::DMRGMultAlgorithm) end
function finalize!(m::GMPSIterativeMultCache, alg::DMRGMult1)
    leftsweep!(m, alg)
    rightsweep_final!(m, alg)
end
function DMRG.leftsweep!(m::GMPSIterativeMultCache, alg::DMRGMult1)
    z, x, y = m.z, m.x, m.y
    hstorage = m.hstorage
    L = length(z)
    kvals = Float64[]
    for site in 1:L-1
        (alg.verbosity > 3) && println("Sweeping from left to right at site: $site")
        
        # mpsj = g_ac_prime(x[site], y[site], hstorage[site], hstorage[site+1])
        left_xy = get_left_xy(hstorage[site], x[site], y[site])
        @tensor mpsj[1,2;5] := left_xy[1,2,3,4] * hstorage[site+1][4,3,5]
        
        push!(kvals, norm(mpsj))
        (alg.verbosity > 1) && println("residual is $(kvals[end])...")
        z[site], r = leftorth!(mpsj, alg = QR())
        
        # hstorage[site+1] = updatemultleft(hstorage[site], z[site], x[site], y[site])
        @tensor tmp[5,3;4] := left_xy[1,2,3,4] * conj(z[site][1,2,5])
        hstorage[site+1] = tmp 
    end
    return kvals    
end

function DMRG.rightsweep!(m::GMPSIterativeMultCache, alg::DMRGMult1)
    z, x, y = m.z, m.x, m.y
    hstorage = m.hstorage
    L = length(z)
    kvals = Float64[]
    local r
    for site in L:-1:2
        (alg.verbosity > 3) && println("Sweeping from right to left at site: $site")
        
        # mpsj = g_ac_prime(x[site], y[site], hstorage[site], hstorage[site+1])
        xy_right = get_xy_right(hstorage[site+1], x[site], y[site])
        @tensor mpsj[1,4;5] := hstorage[site][1,2,3] * xy_right[3,2,4,5]

        push!(kvals, norm(mpsj))
        (alg.verbosity > 1) && println("residual is $(kvals[end])...")

        r, zj = rightorth(mpsj, (1,), (2,3), alg=LQ())
        z[site] = permute(zj, (1,2), (3,))
        
        # hstorage[site] = updatemultright(hstorage[site+1], z[site], x[site], y[site])
        @tensor tmp[4,5;1] := conj(z[site][1,2,3]) * xy_right[4,5,2,3]
        hstorage[site] = tmp
    end
    # println("norm of r is $(norm(r))")
    z[1] = @tensor tmp[1,2;4] := z[1][1,2,3] * r[3,4]
    return kvals    
end

function rightsweep_final!(m::GMPSIterativeMultCache, alg::DMRGMult1)
    z, x, y = m.z, m.x, m.y
    hstorage = m.hstorage
    L = length(z)
    kvals = Float64[]
    trunc = alg.trunc
    for site in L:-1:2
        (alg.verbosity > 3) && println("Sweeping from right to left at site: $site")
        
        # mpsj = g_ac_prime(x[site], y[site], hstorage[site], hstorage[site+1])
        xy_right = get_xy_right(hstorage[site+1], x[site], y[site])
        @tensor mpsj[1,4;5] := hstorage[site][1,2,3] * xy_right[3,2,4,5]

        push!(kvals, norm(mpsj))
        (alg.verbosity > 1) && println("residual is $(kvals[end])...")

        u, s, v = stable_tsvd(mpsj, (1,), (2,3), trunc=trunc)
        z[site] = permute(v, (1,2), (3,))
        if site == 2
            r = u * s
            z[1] = @tensor tmp[1,2;4] := z[1][1,2,3] * r[3,4]
        end
        z.s[site] = normalize!(s)
        
        # hstorage[site] = updatemultright(hstorage[site+1], z[site], x[site], y[site])
        @tensor tmp[4,5;1] := conj(z[site][1,2,3]) * xy_right[4,5,2,3]
        hstorage[site] = tmp
    end
    # println("norm of r is $(norm(r))")
    return kvals    
end

function DMRG.leftsweep!(m::GMPSIterativeMultCache, alg::DMRGMult2)
    z, x, y = m.z, m.x, m.y
    hstorage = m.hstorage
    L = length(z)
    kvals = Float64[]
    trunc = alg.trunc
    for site in 1:L-1
        (alg.verbosity > 3) && println("Sweeping from left to right at bond: $site")
        twositemps = g_ac_prime2(x[site], x[site+1], y[site], y[site+1], hstorage[site], hstorage[site+2])
        push!(kvals, norm(twositemps))
        (alg.verbosity > 1) && println("residual is $(kvals[end])...")
        u, s, v = stable_tsvd!(twositemps, trunc=trunc)
        z[site] = u
        hstorage[site+1] = updatemultleft(hstorage[site], z[site], x[site], y[site])
    end
    return kvals    
end

function DMRG.rightsweep!(m::GMPSIterativeMultCache, alg::DMRGMult2)
    z, x, y = m.z, m.x, m.y
    hstorage = m.hstorage
    L = length(z)
    kvals = Float64[]
    trunc = alg.trunc
    for site in L-1:-1:1
        (alg.verbosity > 3) && println("Sweeping from right to left at bond: $site")
        twositemps = g_ac_prime2(x[site], x[site+1], y[site], y[site+1], hstorage[site], hstorage[site+2])
        push!(kvals, norm(twositemps))
        (alg.verbosity > 1) && println("residual is $(kvals[end])...")
        u, s, v = stable_tsvd!(twositemps, trunc=trunc)
        z[site+1] = permute(v, (1,2), (3,))
        if site == 1
            z[1] = u * s
        end
        z.s[site+1] = normalize!(s)
        hstorage[site+1] = updatemultright(hstorage[site+2], z[site+1], x[site+1], y[site+1])
    end
    # println("norm of r is $(norm(r))")
    return kvals    
end

# provide the initial guess
_svd_guess(x::GrassmannMPS, y::GrassmannMPS, D::Int) = _svd_guess!(copy(x), y, D)
function _svd_guess!(x::GrassmannMPS, y::GrassmannMPS, D::Int)
    (length(x) == length(y)) || throw(DimensionMismatch())
    left = GrassmannTensorMap(isomorphism(scalartype(x), fuse(space_l(x), space_l(y)), space_l(x) ⊗ space_l(y) ))
    tmp5 = g_fuse(_mult_site(x[1], y[1]), 3)
    @tensor tmp4[1,4;5,6] := left[1,2,3] * tmp5[2,3,4,5,6]
    trunc = truncdim(D)
    for i in 1:length(x)-1
        u, s, v = stable_tsvd!(tmp4, trunc=trunc)
        x[i] = get_data(u)
        _renormalize!(x, get_data(s), false)
        r = s * v
        @tensor tmp1[1,5,4;2] := r[1,2,3] * GrassmannTensorMap(y[i+1])[3,4,5]
        @tensor tmp2[1,3,5;6,2] := tmp1[1,2,3,4] * GrassmannTensorMap(x[i+1])[4,5,6]
        tmp4 = g_fuse(tmp2, 2)

    end
    @tensor tmp[1,2;5] := tmp4[1,2,3,4] * conj(left[5,3,4])
    x[end] = get_data(tmp)
    _rightorth!(x, SVD(), trunc, false, 0)
    setscaling!(x, 1)
    return x
end

function g_ac_prime(xj::MPSTensor, yj::MPSTensor, left::MPSTensor, right::MPSTensor)
    @tensor tmp1[1,5,4;2] := left[1,2,3] * yj[3,4,5]
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
    @tensor tmp2[1,3,5;6,2] := tmp1[1,2,3,4] * xj[4,5,6]
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
    tmp3 = g_fuse(tmp2, 2)
    @tensor tmp2[1,2;5] := tmp3[1,2,3,4] * right[4,3,5]
    return tmp2
end

function g_c_prime(left::MPSTensor, right::MPSTensor)
    @tensor tmp[1; 4] := left[1,2,3] * right[3,2,4]
end

function g_ac_prime2(xj1::MPSTensor, xj2::MPSTensor, yj1::MPSTensor, yj2::MPSTensor, left::MPSTensor, right::MPSTensor)
    @tensor tmp1[1,5,4;2] := left[1,2,3] * yj1[3,4,5]
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
    @tensor tmp2[1,3,5;6,2] := tmp1[1,2,3,4] * xj1[4,5,6]
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

    tmp3 = g_fuse(tmp2, 2)


    @tensor tmp4[4; 1 2 5] := yj2[1,2,3] * right[3,4,5]
    for (f1, f2) in fusiontrees(tmp4)
        coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
        coef2 = (isodd(f2.uncoupled[2].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1

        coef = coef1 * coef2 
        if coef != 1
            lmul!(coef, tmp4[f1, f2])
        end
    end 
    @tensor tmp5[4 1 2 5; 6] := xj2[1,2,3] * tmp4[3,4,5,6]
    for (f1, f2) in fusiontrees(tmp5)
        coef1 = (isodd(f1.uncoupled[2].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
        coef2 = (isodd(f1.uncoupled[3].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1

        coef = coef1 * coef2 
        if coef != 1
            lmul!(coef, tmp5[f1, f2])
        end
    end 
    tmp4 = g_fuse(tmp5, 3)


    @tensor tmp2[1,2;5,6] := tmp3[1,2,3,4] * tmp4[4,3,5,6]
    return tmp2
end

function updatemultleft(left::MPSTensor, zj::MPSTensor, xj::MPSTensor, yj::MPSTensor)
    @tensor tmp1[1,5,4;2] := left[1,2,3] * yj[3,4,5]
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
    @tensor tmp2[1,3,5;6,2] := tmp1[1,2,3,4] * xj[4,5,6]
    for (f1, f2) in fusiontrees(tmp2)
        coef1 = (isodd(f2.uncoupled[2].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
        coef2 = (isodd(f2.uncoupled[2].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
        coef3 = (isodd(f2.uncoupled[2].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1

        coef = coef1 * coef2 * coef3 
        if coef != 1
            lmul!(coef, tmp2[f1, f2])
        end
    end
    tmp3 = g_fuse(tmp2, 2)
    @tensor tmp2[5,3;4] := tmp3[1,2,3,4] * conj(zj[1,2,5])
    # for (f1, f2) in fusiontrees(tmp2)
    #     coef1 = (isodd(f1.uncoupled[2].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
    #     coef2 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
        
    #     coef = coef1 * coef2 
    #     if coef != 1
    #         lmul!(coef, tmp2[f1, f2])
    #     end
    # end
    return tmp2
end

function updatemultright(right::MPSTensor, zj::MPSTensor, xj::MPSTensor, yj::MPSTensor)
    @tensor tmp1[4; 1 2 5] := yj[1,2,3] * right[3,4,5]
    for (f1, f2) in fusiontrees(tmp1)
        coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
        coef2 = (isodd(f2.uncoupled[2].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1

        coef = coef1 * coef2 
        if coef != 1
            lmul!(coef, tmp1[f1, f2])
        end
    end 
    @tensor tmp2[4 1 2 5; 6] := xj[1,2,3] * tmp1[3,4,5,6]
    for (f1, f2) in fusiontrees(tmp2)
        coef1 = (isodd(f1.uncoupled[2].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
        coef2 = (isodd(f1.uncoupled[3].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1

        coef = coef1 * coef2 
        if coef != 1
            lmul!(coef, tmp2[f1, f2])
        end
    end 
    tmp3 = g_fuse(tmp2, 3)
    @tensor tmp2[4,5;1] := conj(zj[1,2,3]) * tmp3[4,5,2,3]
    
    return tmp2
end
