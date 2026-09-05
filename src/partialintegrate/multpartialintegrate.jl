

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
			hstorage[iz] = GrassmannTransferMatrix(j, x, y) * hstorage[iz]
            ixy -= 2
		else
			iz -= 1
			(iz == 1) && break
            hstorage[iz] = updatemultright(hstorage[iz+1], z[iz], x[ixy], y[ixy])
            ixy -= 1
		end
	end
    for j in 1:2:(ixy-1)
		# single fermionic permutes bracketing the transfer-matrix product (via @grassmann)
		@grassmann tmp0[1; 2 3] := hstorage[iz][1,2,3]
		tmp = tmp0 * GrassmannTransferMatrix(j, x, y)
		@grassmann tmp1[1 2; 3] := tmp[1,2,3]
        hstorage[1] = tmp1
    end

    return IntegrateBandIterativeMultCache(z, x, y, lattice, band, hstorage)
end

# multintegrateband(lattice::AbstractGrassmannLattice, x::GrassmannMPS, y::GrassmannMPS; alg::DMRGMultAlgorithm=DefaultMultAlg, band::Int=1) = multintegrateband(lattice, x, y, alg, band=band)
# multintegrateband(lattice::AbstractGrassmannLattice, x::GrassmannMPS, y::GrassmannMPS, alg::DMRGMultAlgorithm; band::Int=1) = iterativemultintegrate(x, y, lattice, alg, band=band)


# multiply x and y, integrate out band
function multintegrateband(lattice::AbstractGrassmannLattice, x::GrassmannMPS, y::GrassmannMPS, alg::DMRGMultAlgorithm; band::Int=1)
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


compute!(env::IntegrateBandIterativeMultCache, alg::DMRG1) = iterative_compute!(env, alg)
sweep!(m::IntegrateBandIterativeMultCache, alg::DMRG1) = vcat(leftsweep!(m, alg), rightsweep!(m, alg))
function finalize!(m::IntegrateBandIterativeMultCache, alg::DMRG1)
    leftsweep!(m, alg)
    rightsweep_final!(m, alg)
end




function leftsweep!(m::IntegrateBandIterativeMultCache, alg::DMRG1)
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
			# single fermionic permutes bracketing the transfer-matrix product (via @grassmann)
			@grassmann tmp0[1; 2 3] := hstorage[site+1][1,2,3]
			tmp = tmp0 * GrassmannTransferMatrix(j, x, y)
			@grassmann tmp1[1 2; 3] := tmp[1,2,3]
			hstorage[site+1] = tmp1
		end
    end
	# println(kvals)
    return kvals    
end

function rightsweep!(m::IntegrateBandIterativeMultCache, alg::DMRG1)
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
			hstorage[site] = GrassmannTransferMatrix(j, x, y) * hstorage[site]
		end
    end
    # println("norm of r is $(norm(r))")
    z[1] = @tensor tmp[1,2;4] := z[1][1,2,3] * r[3,4]
	# println(kvals)
    return kvals
end

function rightsweep_final!(m::IntegrateBandIterativeMultCache, alg::DMRG1)
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

        u, s, v, _ = tsvd(mpsj, (1,), (2,3); alg=SDD(), trunc=trunc)
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
			hstorage[site] = GrassmannTransferMatrix(j, x, y) * hstorage[site]
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
    fuser = isomorphism(scalartype(x), fuse(space_l(x), space_l(y)), space_l(x) ⊗ space_l(y) )
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
            @grassmann tmp1[1,5,4;2] := left[1,2,3] * y[i][3,4,5]
            @grassmann tmp2[1,3,5;6,2] := tmp1[1,2,3,4] * x[i][4,5,6]
            tmp3 = g_fuse(tmp2, 2)

            u, s, v, _ = tsvd(tmp3; alg=SDD(), trunc=trunc)
			data[idx] = u
			idx += 1
            left = s * v

            i += 1
		end
	end
    @tensor tmp[1,2;6] := data[end][1,2,3] * left[3,4,5] * conj(fuser[6,4,5])
    data[end] = tmp

    return GrassmannMPS(data)
end


function _lattice_mapping(lattice, lattice2; band::Int=1)
	r2 = indexmappings(lattice2)
	r1 = indexmappings(lattice)
	mm = Dict(r1[(j, c, b, ifelse(bj<band, bj, bj+1))]=>pos for ((j, c, b, bj), pos) in r2)
	rmm = Dict(v=>k for (k,v) in mm)
	return [rmm[i] for i in 1:length(lattice2)]
end














# multintegrateband(lattice::AbstractGrassmannLattice, x::GrassmannMPS, y::GrassmannMPS; alg::SVDCompression, band::Int=1) = mult(lattice, x, y, alg.trunc, band=band)
multintegrateband(lattice::AbstractGrassmannLattice, x::GrassmannMPS, y::GrassmannMPS, alg::SVDCompression; band::Int=1) = multintegrateband(lattice, x, y, alg.trunc, band=band, verbosity=alg.verbosity)
multintegrateband(lattice::AbstractGrassmannLattice, x::GrassmannMPS, y::GrassmannMPS; trunc::TruncationScheme, band::Int=1, verbosity::Int=0) = multintegrateband(lattice, x, y, trunc, band=band, verbosity=verbosity)

function multintegrateband(lattice::AbstractGrassmannLattice, x::GrassmannMPS, y::GrassmannMPS, trunc::TruncationScheme; band::Int=1, verbosity::Int=0)
    (ConjugationStyle(lattice) isa AdjacentConjugation) || throw(ArgumentError("multintegrateband only supports AdjacentConjugation style"))
	(1 <= band <= lattice.bands) || throw(BoundsError(1:lattice.bands, band))
	(length(x) == length(lattice)) || throw(DimensionMismatch())
    (length(x) == length(y)) || throw(DimensionMismatch())

	lattice2 = similar(lattice, bands=lattice.bands-1)
	r2 = indexmappings(lattice2)
	r1 = indexmappings(lattice)
	mm = Dict(r1[(j, c, b, ifelse(bj<band, bj, bj+1))]=>pos for ((j, c, b, bj), pos) in r2)

	data = similar(x.data, length(lattice2))
    fuser = isomorphism(scalartype(x), fuse(space_l(x), space_l(y)), space_l(x) ⊗ space_l(y) )
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
            @grassmann tmp1[1,5,4;2] := left[1,2,3] * y[i][3,4,5]
            @grassmann tmp2[1,3,5;6,2] := tmp1[1,2,3,4] * x[i][4,5,6]
            tmp3 = g_fuse(tmp2, 2)

            q, r = leftorth!(tmp3, alg=QR())
			data[idx] = q
			idx += 1
            left = r

            i += 1
		end
	end
    @tensor tmp[1,2;6] := data[end][1,2,3] * left[3,4,5] * conj(fuser[6,4,5])
    data[end] = tmp

    z = GrassmannMPS(data)
    (verbosity >= 2) && println("bond dimension of intermediate GMPS: ", bond_dimension(z))
    _rightorth!(z, SVD(), trunc, false, verbosity)
    return z
end
