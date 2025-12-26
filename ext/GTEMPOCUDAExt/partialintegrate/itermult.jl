



function cu_parint_cache(z::GrassmannMPS, xs::GrassmannMPS...; cidx::Vector{Int}, verbosity::Int=0, useHCache::Bool=DefaultUseCache)
    Lxs = length(xs[1])
    (unique(length.(xs)) == [Lxs,]) || throw(DimensionMismatch("unique($(length.(xs))) != [$Lxs,]"))
    @assert Lxs == length(z) + 2*length(cidx)
    check_contract_idx(Lxs, cidx)

    Lz = length(z)
    left = isomorphism(scalartype(xs[1]), space_l(z), ⊗(space_l.(xs)...) )
    right = isomorphism(scalartype(z), ⊗(space_r.(xs)...)', space_r(z)')
    hstorage = useHCache ? CachedVector{Union{typeof(left),typeof(right)}}(undef, Lz+1) : Vector{Union{typeof(left),typeof(right)}}(undef, Lz+1)
    hstorage[1] = left
    hstorage[Lz+1] = right

    ixs = Lxs
    iz = Lz
    hip1 = tocu(hstorage[iz+1])
	while true
		rt = @elapsed if insorted(ixs-1, cidx)
            j = ixs ÷ 2
            @assert 2*j == ixs
			hip1 = (tocu(GrassmannTransferMatrix(j, xs...)) * GrassmannTensorMap(hip1)).data
            # hstorage[iz+1] = hip1
            normalize!(hip1)
            ixs -= 2
		else
            hstorage[iz+1] = fromcu(hip1)
			(iz == 1) && break
            below_right = get_below_right(hip1, tocu.(getindex.(xs, ixs))...)
            hip1 = above_below_right(tocu(z[iz]), below_right)
            normalize!(hip1)
			iz -= 1
            ixs -= 1
		end
        (verbosity >= 2) && println("$ixs / $Lxs cost $rt Seconds")
	end
    for j in 2:2:ixs
        hstorage[1] = left_m(hstorage[1], GrassmannTransferMatrix(j÷2, xs...))
        normalize!(hstorage[1])
    end

    return PartialIntegrateIterativeMultCache(z, xs, cidx, hstorage)
end

function cu_parint_iterativemult(xs::GrassmannMPS...; cidx::Vector{Int}, alg::CuDMRGMultAlgorithm)
    rt = @elapsed if alg.initguess == :svd
        z = _cu_parint_svd_guess(xs...; cidx=cidx, trunc=alg.trunc, verbosity=alg.verbosity)
    else
        error("unsupported initguess $(alg.initguess)")
    end
    (alg.verbosity >= 1) && println("initial guess cost $rt Seconds")

	rt = @elapsed cache = cu_parint_cache(z, xs..., cidx=cidx, verbosity=alg.verbosity)
    (alg.verbosity >= 1) && println("build cache cost $rt Seconds")

    deltas = compute!(cache, alg)
    z = cache.o
    _rescaling!(z)
    setscaling!(z, 1)
    return z
end

compute!(env::PartialIntegrateIterativeMultCache, alg::CuDMRGMult1) = iterative_compute!(env, alg)
GTEMPO.sweep!(m::PartialIntegrateIterativeMultCache, alg::CuDMRGMult1) = vcat(leftsweep!(m, alg), rightsweep!(m, alg))
function GTEMPO.finalize!(m::PartialIntegrateIterativeMultCache, alg::CuDMRGMult1)
    leftsweep!(m, alg)
    rightsweep_final!(m, alg)
end

function leftsweep!(m::PartialIntegrateIterativeMultCache, alg::CuDMRGMult1)
    z, xs, cidx, hstorage = m.o, m.xs, m.cidx, m.hstorage

    Lxs = length(xs[1])
    ixs = findfirst(x->!insorted(x, cidx), 1:2:Lxs) * 2 - 1
    Lz, iz = length(z), 1
    kvals = Float64[]
    hi = tocu(hstorage[iz])
    while true
        rt = @elapsed if insorted(ixs, cidx)
            j = (ixs+1) ÷ 2
            hi = left_m(hi, tocu(GrassmannTransferMatrix(j, xs...)))
            normalize!(hi)
            ixs += 2
        else
            hstorage[iz] = fromcu(hi)
            (iz == Lz) && break
            (alg.verbosity >= 4) && println("Sweeping from left to right at site: $iz")
            tmp = get_left_below(hi, tocu.(getindex.(xs, ixs))...)
            mpsj = left_below__right(tmp, tocu(hstorage[iz+1]))
    
            push!(kvals, norm(mpsj))
            (alg.verbosity >= 3) && println("residual is $(kvals[end])...")
            q, r = leftorth!(mpsj, alg = QR())
            z[iz] = fromcu(q)
            
            hi = left_below_above(tmp, q)
            normalize!(hi)
            iz += 1
            ixs += 1
        end
        (alg.verbosity >= 2) && println("$ixs / $Lxs cost $rt Seconds")
    end
    
    (alg.verbosity >= 2) && println("z of bond dimension: ", bond_dimension(z))
	# println(kvals)
    return kvals    
end

function rightsweep!(m::PartialIntegrateIterativeMultCache, alg::CuDMRGMult1)
    z, xs, cidx, hstorage = m.o, m.xs, m.cidx, m.hstorage
    Lxs = length(xs[1])
    ixs = findlast(x->!insorted(x, cidx), 1:2:Lxs) * 2
    Lz, iz = length(z), length(z)
    kvals = Float64[]
    local l
    hip1 = tocu(hstorage[iz+1])
    while true
        rt = @elapsed if insorted(ixs-1, cidx)
            j = (ixs+1) ÷ 2
            hip1 = m_right(tocu(GrassmannTransferMatrix(j, xs...)), hip1)
            normalize!(hip1)
            ixs -= 2
        else
            hstorage[iz+1] = fromcu(hip1)
            (iz == 1) && break
            (alg.verbosity >= 4) && println("Sweeping from left to right at site: $iz")
            tmp = get_below_right(hip1, tocu.(getindex.(xs, ixs))...)
            mpsj = left__below_right(tocu(hstorage[iz]), tmp)
    
            push!(kvals, norm(mpsj))
            (alg.verbosity >= 3) && println("residual is $(kvals[end])...")
            l, q = rightorth(mpsj, (1,), (2,3), alg=LQ())
            zi = permute(q, (1,2), (3,))
            z[iz] = fromcu(zi)
    
            hip1 = above_below_right(zi, tmp)
            normalize!(hip1)
            iz -= 1
            ixs -= 1
        end
        (alg.verbosity >= 2) && println("$ixs / $Lxs cost $rt Seconds")
    end
    z[1] = @tensor tmp[1,2;4] := z[1][1,2,3] * fromcu(l)[3,4]
    (alg.verbosity >= 2) && println("z of bond dimension: ", bond_dimension(z))
    return kvals
end

# TODO; check svectors
function rightsweep_final!(m::PartialIntegrateIterativeMultCache, alg::CuDMRGMult1)
    z, xs, cidx, hstorage = m.o, m.xs, m.cidx, m.hstorage
    trunc = alg.trunc

    Lxs = length(xs[1])
    ixs = findlast(x->!insorted(x, cidx), 1:2:Lxs) * 2
    Lz, iz = length(z), length(z)
    kvals = Float64[]
    local u, s
    hip1 = tocu(hstorage[iz+1])
    while true
        rt = @elapsed if insorted(ixs-1, cidx)
            j = (ixs+1) ÷ 2
            hip1 = m_right(tocu(GrassmannTransferMatrix(j, xs...)), hip1)
            normalize!(hip1)
            ixs -= 2
        else
            hstorage[iz+1] = fromcu(hip1)
            (iz == 1) && break
            (alg.verbosity >= 4) && println("Sweeping from left to right at site: $iz")
            tmp = get_below_right(hip1, tocu.(getindex.(xs, ixs))...)
            mpsj = left__below_right(tocu(hstorage[iz]), tmp)
    
            push!(kvals, norm(mpsj))
            (alg.verbosity >= 3) && println("residual is $(kvals[end])...")
            u, s, v = stable_tsvd(mpsj, (1,), (2,3), trunc=trunc)
            zi = permute(v, (1,2), (3,))
            z[iz] = fromcu(zi)
            z.s[iz] = fromcu(normalize!(s))

            hip1 = above_below_right(zi, tmp)
            normalize!(hip1)
            iz -= 1
            ixs -= 1
        end
        (alg.verbosity >= 2) && println("$ixs / $Lxs cost $rt Seconds")
    end
    r = fromcu(u * s)
    z[1] = @tensor tmp[1,2;4] := z[1][1,2,3] * r[3,4]

    # z.svectors[1] = Diagonal(id(space_l(z[1])))
	# z.svectors[end] = Diagonal(id(space_r(z[end])'))
    z.svectors[1] = DiagonalTensorMap{Float64}(ones, space_l(z[1]) )
	z.svectors[end] = DiagonalTensorMap{Float64}(ones, space_r(z[end])' )

    (alg.verbosity >= 2) && println("z of bond dimension: ", bond_dimension(z))
    return kvals
end




function _cu_parint_svd_guess(xs::GrassmannMPS...; cidx::Vector{Int}, trunc::TruncationScheme=DefaultTruncation, verbosity::Int=0)
    Lxs = length(xs[1])
    (unique(length.(xs)) == [Lxs,]) || throw(DimensionMismatch("unique($(length.(xs))) != [$Lxs,]"))
    check_contract_idx(Lxs, cidx)

    Lz = Lxs-2*length(cidx)
    A = xs[1][1]
	B = bondtensortype(spacetype(A), real(scalartype(A)))
	z = GrassmannMPS(similar(xs[1].data, Lz), Vector{Union{Missing, B}}(missing, Lz), Ref(1.0))
    left = isomorphism(scalartype(xs[1]), fuse(space_l.(xs)...), ⊗(space_l.(xs)...) )
    right = isomorphism(scalartype(xs[1]), ⊗(space_r.(xs)...)', fuse(space_r.(xs)...))

    ixs = 1
    iz = 1
    left = tocu(left)
	while true
		rt = @elapsed if insorted(ixs, cidx)
            left = left_m(left, tocu(GrassmannTransferMatrix((ixs+1)÷2, xs...)))
            ixs += 2
		else
            tmp = get_left_below(left, tocu.(getindex.(xs, ixs))...)
            # z[iz], left = leftorth!(tmp, alg=QR())
            u, s, v = stable_tsvd!(tmp, trunc=trunc)
            z[iz] = fromcu(u)
            left = s * v
            iz += 1
            ixs += 1
		end
        (verbosity >= 2) && println("$ixs / $Lxs cost $rt Seconds")
        _renormalize!(z, left, false)
        (ixs > Lxs) && break
	end
    left = fromcu(left)
    @tensor tmp[1 2;4] := z[end][1 2 3] * left_right(left, right)[3 4]
    z[end] = tmp
    _rightorth!(z, SVD(), trunc, false, verbosity)
    setscaling!(z, 1)

    (verbosity >= 2) && println("initial z of bond dimension: ", bond_dimension(z))
    return z
end
