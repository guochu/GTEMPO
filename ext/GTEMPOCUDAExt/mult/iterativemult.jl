

function cu_mult_cache(z::GrassmannMPS, x::GrassmannMPS, y::GrassmannMPS; useHCache::Bool=DefaultUseCache)
    @assert length(z) == length(x) == length(y)
    # initialize Hstorage
    L = length(z)
    right = ones(scalartype(z), space_r(y)' ⊗ space_r(x)', space_r(z)')
    hstorage = useHCache ? CachedVector{typeof(right)}(undef, L+1) : Vector{typeof(right)}(undef, L+1)

    hstorage[1] = ones( scalartype(z), space_l(z) ⊗ space_l(x)', space_l(y) )
    hstorage[L+1] = right
    hip1 = tocu(hstorage[L+1])
    for i in L:-1:2
        xy_right = get_xy_right(hip1, tocu(x[i]), tocu(y[i]))
        @tensor hip1[4,5;1] := conj(tocu(z[i])[1,2,3]) * xy_right[4,5,2,3]
        hstorage[i] = fromcu(hip1)
    end

    return GMPSIterativeMultCache(z, x, y, hstorage)
end

function cu_iterativemult(x::GrassmannMPS, y::GrassmannMPS, alg::CuDMRGMultAlgorithm)
    if alg.initguess == :svd
        z = _cu_svd_guess(x, y, alg.D)
    elseif alg.initguess == :rand
        z = randomgmps(promote_type(scalartype(x), scalartype(y)), length(x), D=alg.D)
    elseif alg.initguess == :pre
        z = increase_bond!(copy(x), alg.D)
    else
        error("unsupported initguess $(alg.initguess)")
    end
    cache = cu_mult_cache(z, x, y)
    deltas = cu_compute!(cache, alg)
    z = cache.z
    setscaling!(z, scaling(x) * scaling(y))
    _rescaling!(z)
    return z
end

cu_compute!(env::GMPSIterativeMultCache, alg::CuDMRGMultAlgorithm) = iterative_compute!(env, alg)


GTEMPO.sweep!(m::GMPSIterativeMultCache, alg::CuDMRGMultAlgorithm) = vcat(leftsweep!(m, alg), rightsweep!(m, alg))

function GTEMPO.finalize!(m::GMPSIterativeMultCache, alg::CuDMRGMult1)
    leftsweep!(m, alg)
    rightsweep_final!(m, alg)
    # rightsweep!(m, alg)
end
function leftsweep!(m::GMPSIterativeMultCache, alg::CuDMRGMult1)
    z, x, y = m.z, m.x, m.y
    hstorage = m.hstorage
    L = length(z)
    kvals = Float64[]
    hi = tocu(hstorage[1])
    for site in 1:L-1
        (alg.verbosity >= 4) && println("Sweeping from left to right at site: $site")

        # mpsj = g_ac_prime(x[site], y[site], hstorage[site], hstorage[site+1])
        left_xy = get_left_xy(hi, tocu(x[site]), tocu(y[site]))
        @tensor mpsj[1,2;5] := left_xy[1,2,3,4] * tocu(hstorage[site+1])[4,3,5]
        
        push!(kvals, norm(mpsj))
        (alg.verbosity >= 3) && println("residual is $(kvals[end])...")
        q, r = leftorth!(mpsj, alg = QR())
        z[site] = fromcu(q)
        
        # hstorage[site+1] = updatemultleft(hstorage[site], z[site], x[site], y[site])
        @tensor hi[5,3;4] := left_xy[1,2,3,4] * conj(q[1,2,5])
        hstorage[site+1] = fromcu(hi) 
    end
    return kvals    
end

function rightsweep!(m::GMPSIterativeMultCache, alg::CuDMRGMult1)
    z, x, y = m.z, m.x, m.y
    hstorage = m.hstorage
    L = length(z)
    kvals = Float64[]
    local r
    hi = tocu(hstorage[L+1])
    for site in L:-1:2
        (alg.verbosity >= 4) && println("Sweeping from right to left at site: $site")
        
        # mpsj = g_ac_prime(x[site], y[site], hstorage[site], hstorage[site+1])
        xy_right = get_xy_right(hi, tocu(x[site]), tocu(y[site]))
        @tensor mpsj[1,4;5] := tocu(hstorage[site])[1,2,3] * xy_right[3,2,4,5]

        push!(kvals, norm(mpsj))
        (alg.verbosity >= 3) && println("residual is $(kvals[end])...")

        r, zj = rightorth(mpsj, (1,), (2,3), alg=LQ())
        zj = permute(zj, (1,2), (3,))
        z[site] = fromcu(zj)
        
        # hstorage[site] = updatemultright(hstorage[site+1], z[site], x[site], y[site])
        @tensor hi[4,5;1] := conj(zj[1,2,3]) * xy_right[4,5,2,3]
        hstorage[site] = fromcu(hi)
    end
    # println("norm of r is $(norm(r))")
    z[1] = @tensor tmp[1,2;4] := z[1][1,2,3] * fromcu(r)[3,4]
    return kvals    
end

function rightsweep_final!(m::GMPSIterativeMultCache, alg::CuDMRGMult1)
    z, x, y = m.z, m.x, m.y
    hstorage = m.hstorage
    L = length(z)
    kvals = Float64[]
    trunc = alg.trunc
    hi = tocu(hstorage[L+1])
    for site in L:-1:2
        (alg.verbosity >= 4) && println("Sweeping from right to left at site: $site")
        
        xy_right = get_xy_right(hi, tocu(x[site]), tocu(y[site]))
        @tensor mpsj[1,4;5] := tocu(hstorage[site])[1,2,3] * xy_right[3,2,4,5]

        push!(kvals, norm(mpsj))
        (alg.verbosity >= 3) && println("residual is $(kvals[end])...")

        u, s, v = stable_tsvd(mpsj, (1,), (2,3), trunc=trunc)
        v = permute(v, (1,2), (3,))
        z[site] = fromcu(v)
        if site == 2
            r = fromcu(u * s)
            z[1] = @tensor tmp[1,2;4] := z[1][1,2,3] * r[3,4]
        end
        z.s[site] = normalize!(fromcu(s))
        
        # hstorage[site] = updatemultright(hstorage[site+1], z[site], x[site], y[site])
        @tensor hi[4,5;1] := conj(v[1,2,3]) * xy_right[4,5,2,3]
        hstorage[site] = fromcu(hi)
    end
    # println("norm of r is $(norm(r))")
    return kvals    
end





# provide the initial guess
_cu_svd_guess(x::GrassmannMPS, y::GrassmannMPS, D::Int) = _cu_svd_guess!(copy(x), y, D)
function _cu_svd_guess!(x::GrassmannMPS, y::GrassmannMPS, D::Int)
    (length(x) == length(y)) || throw(DimensionMismatch())
    left = GrassmannTensorMap(isomorphism(scalartype(x), fuse(space_l(x), space_l(y)), space_l(x) ⊗ space_l(y) ))
    tmp5 = g_fuse(_mult_site(x[1], y[1]), 3)
    @tensor tmp4[1,4;5,6] := left[1,2,3] * tmp5[2,3,4,5,6]
    trunc = truncdim(D)
    tmp4 = tocu(tmp4)
    for i in 1:length(x)-1
        u, s, v = stable_tsvd!(tmp4, trunc=trunc)
        x[i] = get_data(fromcu(u))
        _renormalize!(x, get_data(s), false)
        r = s * v
        @tensor tmp1[1,5,4;2] := r[1,2,3] * tocu(GrassmannTensorMap(y[i+1]))[3,4,5]
        @tensor tmp2[1,3,5;6,2] := tmp1[1,2,3,4] * tocu(GrassmannTensorMap(x[i+1]))[4,5,6]
        tmp4 = g_fuse(tmp2, 2)

    end
    @tensor tmp[1,2;5] := fromcu(tmp4)[1,2,3,4] * conj(left[5,3,4])
    x[end] = get_data(tmp)
    _cu_rightorth!(x, SVD(), trunc, false, 0)
    setscaling!(x, 1)
    return x
end


