

# multiply x and y, integrate out band, result in z
struct IntegrateBandIterativeMultCache{_O, _A, _B, L, _H} 
    z::_O
    x::_A
    y::_B
    lattice::L
    band::Int
    hstorage::_H
end


function mult_cache(z::GrassmannMPS, x::GrassmannMPS, y::GrassmannMPS, lattice::AbstractGrassmannLattice; band::Int=1)
    @assert length(x) == length(y) == length(lattice)
    lattice2 = similar(lattice, bands=lattice.bands-1)
    @assert length(z) == length(lattice2)

    r2 = indexmappings(lattice2)
	r1 = indexmappings(lattice)
	mm = Dict(r1[(j, c, b, ifelse(bj<band, bj, bj+1))]=>pos for ((j, c, b, bj), pos) in r2)

    L = length(z)
    right = ones(scalartype(z), space_r(y)' ⊗ space_r(x)', space_r(z)')
    hstorage = Vector{typeof(right)}(undef, L+1)
    hstorage[L+1] = right
    hstorage[1] = ones( scalartype(z), space_l(z) ⊗ space_l(x)', space_l(y) )

    ixy = length(lattice)
    iz = L+1
	while true
		pos2 = get(mm, ixy, nothing)
		if isnothing(pos2)
            j = ixy ÷ 2
            @assert 2*j == ixy
			hstorage[iz] = (GrassmannTransferMatrix(j, x, y) * GrassmannTensorMap(hstorage[iz])).data
            ixy -= 2
		else
			iz -= 1
			(iz == 1) && break
            hstorage[iz] = updatemultright(hstorage[iz+1], z[iz], x[ixy], y[ixy])
            ixy -= 1
		end
	end
    for j in 1:2:(ixy-1)
		# tmp = GrassmannTensorMap(permute(hstorage[iz], ((1,), (2,3)))) * GrassmannTransferMatrix(j, x, y)
        # hstorage[1] = permute(tmp.data, ((1,2), (3,)))
		tmp = permute(GrassmannTensorMap(hstorage[iz]), ((1,), (2,3))) * GrassmannTransferMatrix(j, x, y)
        hstorage[1] = permute(tmp, ((1,2), (3,))).data
    end

    return IntegrateBandIterativeMultCache(z, x, y, lattice, band, hstorage)
end

# integrateband(lattice::AbstractGrassmannLattice, x::GrassmannMPS, y::GrassmannMPS; alg::DMRGMultAlgorithm=DefaultMultAlg, band::Int=1) = integrateband(lattice, x, y, alg, band=band)
integrateband(lattice::AbstractGrassmannLattice, x::GrassmannMPS, y::GrassmannMPS, alg::DMRGMultAlgorithm; band::Int=1) = iterativemult(x, y, lattice, alg, band=band)


# multiply x and y, integrate out band
function iterativemult(x::GrassmannMPS, y::GrassmannMPS, lattice::AbstractGrassmannLattice, alg::DMRGMultAlgorithm; band::Int=1)
    if alg.initguess == :svd
        z = _integrateband_svd_guess(lattice, x, y, alg.D; band=band)
    else
        error("unsupported initguess $(alg.initguess)")
    end

	cache = mult_cache(z, x, y, lattice, band=band)

    deltas = compute!(cache, alg)
    z = cache.z
    setscaling!(z, scaling(x) * scaling(y))
    _rescaling!(z)
    return z
end


compute!(env::IntegrateBandIterativeMultCache, alg::DMRGMult1) = iterative_compute!(env, alg)
sweep!(m::IntegrateBandIterativeMultCache, alg::DMRGMult1) = vcat(leftsweep!(m, alg), rightsweep!(m, alg))
function finalize!(m::IntegrateBandIterativeMultCache, alg::DMRGMult1)
    leftsweep!(m, alg)
    rightsweep_final!(m, alg)
end




function leftsweep!(m::IntegrateBandIterativeMultCache, alg::DMRGMult1)
    z, x, y = m.z, m.x, m.y
    hstorage = m.hstorage

	lattice2 = similar(m.lattice, bands=m.lattice.bands-1)
	mm = _lattice_mapping(m.lattice, lattice2, band=m.band)

    L = length(z)
    kvals = Float64[]
    for site in 1:L-1
        (alg.verbosity > 3) && println("Sweeping from left to right at site: $site")
        # mpsj = g_ac_prime(x[mm[site]], y[mm[site]], hstorage[site], hstorage[site+1])
        left_xy = get_left_xy(hstorage[site], x[mm[site]], y[mm[site]])
        @tensor mpsj[1,2;5] := left_xy[1,2,3,4] * hstorage[site+1][4,3,5]

        push!(kvals, norm(mpsj))
        (alg.verbosity > 1) && println("residual is $(kvals[end])...")
        z[site], r = leftorth!(mpsj, alg = QR())
        
        # hstorage[site+1] = updatemultleft(hstorage[site], z[site], x[mm[site]], y[mm[site]])
        @tensor tmp[5,3;4] := left_xy[1,2,3,4] * conj(z[site][1,2,5])
        hstorage[site+1] = tmp

        for i in (mm[site]+1):2:(mm[site+1]-1)
			j = (i+1) ÷ 2
			# tmp = GrassmannTensorMap(permute(hstorage[site+1], ((1,), (2,3)))) * GrassmannTransferMatrix(j, x, y)
			# hstorage[site+1] = permute(tmp.data, ((1,2), (3,)))
			tmp = permute(GrassmannTensorMap(hstorage[site+1]), ((1,), (2,3))) * GrassmannTransferMatrix(j, x, y)
			hstorage[site+1] = permute(tmp, ((1,2), (3,))).data

		end
    end
	# println(kvals)
    return kvals    
end

function rightsweep!(m::IntegrateBandIterativeMultCache, alg::DMRGMult1)
    z, x, y = m.z, m.x, m.y
    hstorage = m.hstorage

	lattice2 = similar(m.lattice, bands=m.lattice.bands-1)
	mm = _lattice_mapping(m.lattice, lattice2, band=m.band)

	L = length(z)
    kvals = Float64[]
    local r
    for site in L:-1:2
        (alg.verbosity > 3) && println("Sweeping from right to left at site: $site")
        # mpsj = g_ac_prime(x[mm[site]], y[mm[site]], hstorage[site], hstorage[site+1])
        xy_right = get_xy_right(hstorage[site+1], x[mm[site]], y[mm[site]])
        @tensor mpsj[1,4;5] := hstorage[site][1,2,3] * xy_right[3,2,4,5]
        
        push!(kvals, norm(mpsj))
        (alg.verbosity > 1) && println("residual is $(kvals[end])...")

        r, zj = rightorth(mpsj, (1,), (2,3), alg=LQ())
        z[site] = permute(zj, (1,2), (3,))
        # hstorage[site] = updatemultright(hstorage[site+1], z[site], x[mm[site]], y[mm[site]])
        @tensor tmp[4,5;1] := conj(z[site][1,2,3]) * xy_right[4,5,2,3]
        hstorage[site] = tmp
        for i in (mm[site]-1):-2:(mm[site-1]+1)
			j = i ÷ 2
			hstorage[site] = (GrassmannTransferMatrix(j, x, y) * GrassmannTensorMap(hstorage[site])).data
		end
    end
    # println("norm of r is $(norm(r))")
    z[1] = @tensor tmp[1,2;4] := z[1][1,2,3] * r[3,4]
	# println(kvals)
    return kvals    
end

function rightsweep_final!(m::IntegrateBandIterativeMultCache, alg::DMRGMult1)
    z, x, y = m.z, m.x, m.y
    hstorage = m.hstorage

	lattice2 = similar(m.lattice, bands=m.lattice.bands-1)
	mm = _lattice_mapping(m.lattice, lattice2, band=m.band)

	L = length(z)
    kvals = Float64[]
    trunc = alg.trunc
    for site in L:-1:2
        (alg.verbosity > 3) && println("Sweeping from right to left at site: $site")
        # mpsj = g_ac_prime(x[mm[site]], y[mm[site]], hstorage[site], hstorage[site+1])
        xy_right = get_xy_right(hstorage[site+1], x[mm[site]], y[mm[site]])
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
        # hstorage[site] = updatemultright(hstorage[site+1], z[site], x[mm[site]], y[mm[site]])
        @tensor tmp[4,5;1] := conj(z[site][1,2,3]) * xy_right[4,5,2,3]
        hstorage[site] = tmp
		for i in (mm[site]-1):-2:(mm[site-1]+1)
			j = i ÷ 2
			hstorage[site] = (GrassmannTransferMatrix(j, x, y) * GrassmannTensorMap(hstorage[site])).data
		end
    end
    # println("norm of r is $(norm(r))")
    return kvals    
end




function _integrateband_svd_guess(lattice::AbstractGrassmannLattice, x::GrassmannMPS, y::GrassmannMPS, D::Int; band::Int=1)
    (ConjugationStyle(lattice) isa AdjacentConjugation) || throw(ArgumentError("integrateband only supports AdjacentConjugation style"))
	(1 <= band <= lattice.bands) || throw(BoundsError(1:lattice.bands, band))
	(length(x) == length(lattice)) || throw(DimensionMismatch())
    (length(x) == length(y)) || throw(DimensionMismatch())

	lattice2 = similar(lattice, bands=lattice.bands-1)
	r2 = indexmappings(lattice2)
	r1 = indexmappings(lattice)
	mm = Dict(r1[(j, c, b, ifelse(bj<band, bj, bj+1))]=>pos for ((j, c, b, bj), pos) in r2)

	data = similar(x.data, length(lattice2))
    trunc = truncdim(D)
    fuser = GrassmannTensorMap(isomorphism(scalartype(x), fuse(space_l(x), space_l(y)), space_l(x) ⊗ space_l(y) ))
    left = fuser

    i = 1
	idx = 1
	while i <= length(lattice)
		pos2 = get(mm, i, nothing)
		if isnothing(pos2)
            j = (i+1)÷2
            @assert 2*j == i+1
			left = left * GrassmannTransferMatrix(j, x, y)
            i += 2
		else
            @tensor tmp1[1,5,4;2] := left[1,2,3] * GrassmannTensorMap(y[i])[3,4,5]
            @tensor tmp2[1,3,5;6,2] := tmp1[1,2,3,4] * GrassmannTensorMap(x[i])[4,5,6]
            tmp3 = g_fuse(tmp2, 2)
            
            u, s, v = stable_tsvd!(tmp3, trunc=trunc)
			data[idx] = u.data
			idx += 1
            left = s * v
    
            i += 1
		end		
	end
    @tensor tmp[1,2;6] := GrassmannTensorMap(data[end])[1,2,3] * left[3,4,5] * conj(fuser[6,4,5])
    data[end] = tmp.data

    return GrassmannMPS(data)
end


function _lattice_mapping(lattice, lattice2; band::Int=1)
	r2 = indexmappings(lattice2)
	r1 = indexmappings(lattice)
	mm = Dict(r1[(j, c, b, ifelse(bj<band, bj, bj+1))]=>pos for ((j, c, b, bj), pos) in r2)
	rmm = Dict(v=>k for (k,v) in mm)
	return [rmm[i] for i in 1:length(lattice2)]
end














# integrateband(lattice::AbstractGrassmannLattice, x::GrassmannMPS, y::GrassmannMPS; alg::SVDCompression, band::Int=1) = mult(lattice, x, y, alg.trunc, band=band)
integrateband(lattice::AbstractGrassmannLattice, x::GrassmannMPS, y::GrassmannMPS, alg::SVDCompression; band::Int=1) = mult(lattice, x, y, alg.trunc, band=band, verbosity=alg.verbosity)
integrateband(lattice::AbstractGrassmannLattice, x::GrassmannMPS, y::GrassmannMPS; trunc::TruncationScheme, band::Int=1, verbosity::Int=0) = mult(lattice, x, y, trunc, band=band, verbosity=verbosity)

function mult(lattice::AbstractGrassmannLattice, x::GrassmannMPS, y::GrassmannMPS, trunc::TruncationScheme; band::Int=1, verbosity::Int=0)
    (ConjugationStyle(lattice) isa AdjacentConjugation) || throw(ArgumentError("integrateband only supports AdjacentConjugation style"))
	(1 <= band <= lattice.bands) || throw(BoundsError(1:lattice.bands, band))
	(length(x) == length(lattice)) || throw(DimensionMismatch())
    (length(x) == length(y)) || throw(DimensionMismatch())

	lattice2 = similar(lattice, bands=lattice.bands-1)
	r2 = indexmappings(lattice2)
	r1 = indexmappings(lattice)
	mm = Dict(r1[(j, c, b, ifelse(bj<band, bj, bj+1))]=>pos for ((j, c, b, bj), pos) in r2)

	data = similar(x.data, length(lattice2))
    fuser = GrassmannTensorMap(isomorphism(scalartype(x), fuse(space_l(x), space_l(y)), space_l(x) ⊗ space_l(y) ))
    left = fuser

    i = 1
	idx = 1
	while i <= length(lattice)
		pos2 = get(mm, i, nothing)
		if isnothing(pos2)
            j = (i+1)÷2
            @assert 2*j == i+1
			left = left * GrassmannTransferMatrix(j, x, y)
            i += 2
		else
            @tensor tmp1[1,5,4;2] := left[1,2,3] * GrassmannTensorMap(y[i])[3,4,5]
            @tensor tmp2[1,3,5;6,2] := tmp1[1,2,3,4] * GrassmannTensorMap(x[i])[4,5,6]
            tmp3 = g_fuse(tmp2, 2)
            
            q, r = leftorth!(tmp3, alg=QR())
			data[idx] = q.data
			idx += 1
            left = r
    
            i += 1
		end		
	end
    @tensor tmp[1,2;6] := GrassmannTensorMap(data[end])[1,2,3] * left[3,4,5] * conj(fuser[6,4,5])
    data[end] = tmp.data

    z = GrassmannMPS(data)
    (verbosity >= 2) && println("bond dimension of intermediate GMPS: ", bond_dimension(z))
    _rightorth!(z, SVD(), trunc, false, verbosity)
    return z
end
