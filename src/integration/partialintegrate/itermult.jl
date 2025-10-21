
##########################################################################
# The struct of the hstorage                                             #
# left enviroment ▶ of z[iz] is the tensor network on the left of z[iz]  #
# left enviroment ◀ of z[iz] is the tensor network on the right of z[iz] #
#                                                                        #
#              ▶ ◀   hstorage[4], hstorage[5]                            #
#            ▶ ◀     hstorage[3], hstorage[4]                            #
#      ▶ ◀           hstorage[2], hstorage[3]                            #
#    ▶ ◀             hstorage[1], hstorage[2]                            #
# ___ ⊤ ⊤ ___ ⊤ ⊤    z                                                   #
# ⊥ ⊥ ⊥ ⊥ ⊥ ⊥ ⊥ ⊥    xs[1]                                               #
# ⊥ ⊥ ⊥ ⊥ ⊥ ⊥ ⊥ ⊥    xs[2]                                               #
#     1 2     3 4    iz                                                  #
# 1 2 3 4 5 6 7 8    ixs                                                 #
#                                                                        #
##########################################################################

# Using iterative multiply, the scaling will lost



# multiply xs, integrate out cidx, result in o
struct PartialIntegrateIterativeMultCache{M<:GrassmannMPS, G<:Tuple, H} 
    o::M
    xs::G
    cidx::Vector{Int}
    hstorage::H
end
function parint_cache(z::GrassmannMPS, xs::GrassmannMPS...; cidx::Vector{Int})
    Lxs = length(xs[1])
    (unique(length.(xs)) == [Lxs,]) || throw(DimensionMismatch())
    @assert Lxs == length(z) + 2*length(cidx)
    chack_contract_idx(Lxs, cidx)

    Lz = length(z)
    left = isomorphism(scalartype(xs[1]), space_l(z), ⊗(space_l.(xs)...) )
    right = isomorphism(scalartype(z), ⊗(space_r.(xs)...)', space_r(z)')
    hstorage = Vector{Union{typeof(left),typeof(right)}}(undef, Lz+1)
    hstorage[Lz+1] = right
    hstorage[1] = left

    ixs = Lxs
    iz = Lz
	while true
		if insorted(ixs-1, cidx)
            j = ixs ÷ 2
            @assert 2*j == ixs
			hstorage[iz+1] = (GrassmannTransferMatrix(j, xs...) * GrassmannTensorMap(hstorage[iz+1])).data
            normalize!(hstorage[iz+1])
            ixs -= 2
		else
			(iz == 1) && break
            below_right = get_below_right(hstorage[iz+1], getindex.(xs, ixs)...)
            hstorage[iz] = above_below_right(z[iz], below_right)
            normalize!(hstorage[iz])
			iz -= 1
            ixs -= 1
		end
	end
    for j in 2:2:ixs
        hstorage[1] = left_m(hstorage[1], GrassmannTransferMatrix(j÷2, xs...))
        normalize!(hstorage[1])
    end

    return PartialIntegrateIterativeMultCache(z, xs, cidx, hstorage)
end

function parint_iterativemult(xs::GrassmannMPS...; cidx::Vector{Int}, alg::DMRGMultAlgorithm)
    if alg.initguess == :svd
        z = _parint_svd_guess(xs...; cidx=cidx, trunc=alg.trunc)
    else
        error("unsupported initguess $(alg.initguess)")
    end

	cache = parint_cache(z, xs..., cidx=cidx)

    deltas = compute!(cache, alg)
    z = cache.o
    _rescaling!(z)
    setscaling!(z, 1)
    return z
end

DMRG.compute!(env::PartialIntegrateIterativeMultCache, alg::DMRGMult1) = iterative_compute!(env, alg)
DMRG.sweep!(m::PartialIntegrateIterativeMultCache, alg::DMRGMult1) = vcat(leftsweep!(m, alg), rightsweep!(m, alg))
function finalize!(m::PartialIntegrateIterativeMultCache, alg::DMRGMult1)
    leftsweep!(m, alg)
    rightsweep_final!(m, alg)
end

function DMRG.leftsweep!(m::PartialIntegrateIterativeMultCache, alg::DMRGMult1)
    z, xs, cidx, hstorage = m.o, m.xs, m.cidx, m.hstorage

    Lxs = length(xs[1])
    ixs = findfirst(x->!insorted(x, cidx), 1:2:Lxs) * 2 - 1
    Lz, iz = length(z), 1
    kvals = Float64[]
    while true
        if insorted(ixs, cidx)
            j = (ixs+1) ÷ 2
            hstorage[iz] = left_m(hstorage[iz], GrassmannTransferMatrix(j, xs...))
            normalize!(hstorage[iz])
            ixs += 2
        else
            (iz == Lz) && break
            (alg.verbosity > 3) && println("Sweeping from left to right at site: $iz")
            tmp = get_left_below(hstorage[iz], getindex.(xs, ixs)...)
            mpsj = left_below__right(tmp, hstorage[iz+1])
    
            push!(kvals, norm(mpsj))
            (alg.verbosity > 1) && println("residual is $(kvals[end])...")
            z[iz], r = leftorth!(mpsj, alg = QR())
            
            tmp = left_below_above(tmp, z[iz])
            hstorage[iz+1] = tmp
            normalize!(hstorage[iz+1])
            iz += 1
            ixs += 1
        end
    end
    
	# println(kvals)
    return kvals    
end

function DMRG.rightsweep!(m::PartialIntegrateIterativeMultCache, alg::DMRGMult1)
    z, xs, cidx, hstorage = m.o, m.xs, m.cidx, m.hstorage
    Lxs = length(xs[1])
    ixs = findlast(x->!insorted(x, cidx), 1:2:Lxs) * 2
    Lz, iz = length(z), length(z)
    kvals = Float64[]
    local l
    while true
        if insorted(ixs-1, cidx)
            j = (ixs+1) ÷ 2
            hstorage[iz+1] = m_right(GrassmannTransferMatrix(j, xs...), hstorage[iz+1])
            normalize!(hstorage[iz+1])
            ixs -= 2
        else
            (iz == 1) && break
            (alg.verbosity > 3) && println("Sweeping from left to right at site: $iz")
            tmp = get_below_right(hstorage[iz+1], getindex.(xs, ixs)...)
            mpsj = left__below_right(hstorage[iz], tmp)
    
            push!(kvals, norm(mpsj))
            (alg.verbosity > 1) && println("residual is $(kvals[end])...")
            l, zj = rightorth(mpsj, (1,), (2,3), alg=LQ())
            z[iz] = permute(zj, (1,2), (3,))
    
            tmp = above_below_right(z[iz], tmp)
            hstorage[iz] = tmp
            normalize!(hstorage[iz])
            iz -= 1
            ixs -= 1
        end
    end
    z[1] = @tensor tmp[1,2;4] := z[1][1,2,3] * l[3,4]
    return kvals
end

# TODO; check svectors
function rightsweep_final!(m::PartialIntegrateIterativeMultCache, alg::DMRGMult1)
    z, xs, cidx, hstorage = m.o, m.xs, m.cidx, m.hstorage
    trunc = alg.trunc

    Lxs = length(xs[1])
    ixs = findlast(x->!insorted(x, cidx), 1:2:Lxs) * 2
    Lz, iz = length(z), length(z)
    kvals = Float64[]
    local u, s
    while true
        if insorted(ixs-1, cidx)
            j = (ixs+1) ÷ 2
            hstorage[iz+1] = m_right(GrassmannTransferMatrix(j, xs...), hstorage[iz+1])
            normalize!(hstorage[iz+1])
            ixs -= 2
        else
            (iz == 1) && break
            (alg.verbosity > 3) && println("Sweeping from left to right at site: $iz")
            tmp = get_below_right(hstorage[iz+1], getindex.(xs, ixs)...)
            mpsj = left__below_right(hstorage[iz], tmp)
    
            push!(kvals, norm(mpsj))
            (alg.verbosity > 1) && println("residual is $(kvals[end])...")
            u, s, v = stable_tsvd(mpsj, (1,), (2,3), trunc=trunc)
            z[iz] = permute(v, (1,2), (3,))
            z.s[iz] = normalize!(s)

            tmp = above_below_right(z[iz], tmp)
            hstorage[iz] = tmp
            normalize!(hstorage[iz])
            iz -= 1
            ixs -= 1
        end
    end
    r = u * s
    z[1] = @tensor tmp[1,2;4] := z[1][1,2,3] * r[3,4]
    return kvals
end




function _parint_svd_guess(xs::GrassmannMPS...; cidx::Vector{Int}, trunc::TruncationScheme=DMRG.DefaultTruncation, verbosity::Int=0)
    L = length(xs[1])
    (unique(length.(xs)) == [L,]) || throw(DimensionMismatch())
    chack_contract_idx(L, cidx)

    Lz = L-2*length(cidx)
    A = xs[1][1]
	B = bondtensortype(spacetype(A), real(scalartype(A)))
	z = GrassmannMPS(similar(xs[1].data, Lz), Vector{Union{Missing, B}}(missing, Lz), Ref(1.0))
    left = isomorphism(scalartype(xs[1]), fuse(space_l.(xs)...), ⊗(space_l.(xs)...) )
    right = isomorphism(scalartype(xs[1]), ⊗(space_r.(xs)...)', fuse(space_r.(xs)...))

    ixs = 1
    iz = 1
	while true
		if insorted(ixs, cidx)
            left = left_m(left, GrassmannTransferMatrix((ixs+1)÷2, xs...))
            ixs += 2
		else
            tmp = get_left_below(left, getindex.(xs, ixs)...)
            # z[iz], left = leftorth!(tmp, alg=QR())
            z[iz], s, v = stable_tsvd!(tmp, trunc=trunc)
            left = s * v
            iz += 1
            ixs += 1
		end
        _renormalize!(z, left, false)
        (ixs > L) && break
	end
    @tensor tmp[1 2;4] := z[end][1 2 3] * left_right(left, right)[3 4]
    z[end] = tmp
    _rightorth!(z, SVD(), trunc, false, verbosity)
    setscaling!(z, 1)
    return z
end

# function _parint_svd_guess(xs::GrassmannMPS...; cidx::Vector{Int}, trunc::TruncationScheme=DMRG.DefaultTruncation, verbosity::Int=0)
#     L = length(xs[1])
#     (unique(length.(xs)) == [L,]) || throw(DimensionMismatch())
#     chack_contract_idx(L, cidx)

# 	data = similar(xs[1].data, L-2*length(cidx))
#     left = isomorphism(scalartype(xs[1]), fuse(space_l.(xs)...), ⊗(space_l.(xs)...) )
#     right = isomorphism(scalartype(xs[1]), ⊗(space_r.(xs)...)', fuse(space_r.(xs)...))

#     ixs = 1
#     iz = 1
# 	while true
# 		if insorted(ixs, cidx)
#             left = left_m(left, GrassmannTransferMatrix((ixs+1)÷2, xs...))
#             ixs += 2
# 		else
#             tmp = get_left_below(left, getindex.(xs, ixs)...)
#             data[iz], s, v = stable_tsvd!(tmp, trunc=trunc)
#             left = s * v
#             iz += 1
#             ixs += 1
# 		end
#         (ixs > L) && break
# 	end
#     @tensor tmp[1 2;4] := data[end][1 2 3] * left_right(left, right)[3 4]
#     data[end] = tmp
#     return GrassmannMPS(data)
# end


