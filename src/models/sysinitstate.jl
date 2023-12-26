

function systhermalstate_iterative!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, model::AbstractImpurityModel; trunc::TruncationScheme=DefaultKTruncation, δτ::Real=0.01)
	δτ = convert(Float64, δτ)
	β = model.bath.β
	n = round(Int, β/δτ)
	@assert n >= 1
	δτ = β / n
	# println("n=", n)

	lattice_i = GrassmannLattice(δτ=β, N=1, contour=:imag, order=1, ordering=A1B1B1A1(), bands=lattice.bands)
	gmps_i = acc_sysdynamics2(lattice_i, model, scaling=n, trunc=trunc)
	# lattice_i = GrassmannLattice(δτ=δτ, N=n, contour=:imag, order=1, ordering=A1B1B1A1(), bands=lattice.bands)
	# gmps_i = acc_sysdynamics2(lattice_i, model, scaling=10, trunc=trunc)


	# for band in 1:lattice_i.bands
	# 	gmps_i = boundarycondition(gmps_i, lattice_i, band=band, trunc=trunc)
	# end
	@assert length(gmps_i) == length(lattice_i)

	# right = _contract_band(gmps_i, lattice_i, 1)
	# for i in 2:lattice_i.k
	# 	right = _contract_band(gmps_i, lattice_i, i) * right
	# end
	# posa, posb = band_boundary(lattice_i, 0)
	# data = [gmps_i[i] for i in posa:posb]
	# data[end] = data[end] * right

	left = _contract_band(gmps_i, lattice_i, 0) 
	posa, posb = band_boundary(lattice_i, lattice_i.k)
	_x, posb = band_boundary(lattice_i, 1)
	# println("posa=", posa, " posb=", posb)
	data = [gmps_i[i] for i in posa:posb]
	data[1] = @tensor tmp[1,3;4] := left[1,2] * data[1][2,3,4]
	# for item in data
	# 	println(space_l(item), " ", space_r(item))
	# end


	# # band, conj, forward
	inds1 = [index(lattice, 1, band=band, conj=false, forward=true) for band in 1:lattice.bands]
	inds2 = [index(lattice, 1, band=band, conj=true, forward=true) for band in lattice.bands:-1:1]
	inds3 = [index(lattice, 1, band=band, conj=false, forward=false) for band in 1:lattice.bands]
	inds4 = [index(lattice, 1, band=band, conj=true, forward=false) for band in lattice.bands:-1:1]

	inds = vcat(inds1, inds2, inds3, inds4)
	# println(inds)

	gmps2 = _creategmps(length(lattice), inds, data, trunc=trunc)
	_rescaling!(gmps2, scale(gmps_i)^(length(gmps_i)))
	# println(bond_dimensions(gmps2))
	# return gmps2
	return mult!(gmps, gmps2, trunc=trunc)
end


function _creategmps(L::Int, inds::Vector{Int}, data::Vector{<:MPSTensor}; trunc)
	@assert length(inds) == length(data)
	if !issorted(inds)
		# println("here...")
		perm = sortperm(inds)
		# println(perm)
		inds = inds[perm]
		gmps = _permute(GrassmannMPS(data), perm, trunc=trunc)
		# gmps = GrassmannMPS(data)
		# easy_swap!(gmps, 1, trunc=trunc)
		data = gmps.data
	end
	# println(inds)
	# for item in data
	# 	println(space_l(item), " ", space_r(item))
	# end
	new_data = similar(data, L)
	leftind = oneunit(_ph)
	for i in 1:L
		# println("i=", i, " ", leftind)
		pos = findfirst(x -> x==i, inds)
		if isnothing(pos)
			new_data[i] = isometry(leftind ⊗ _ph ← leftind)
			# new_data[i] = _trivial_site_tensor(leftind)
		else
			new_data[i] = copy(data[pos])
		end
		leftind = space_r(new_data[i])'
	end
	return GrassmannMPS(new_data)
end


function _trivial_site_tensor(leftind)
	m = isometry(leftind ⊗ _ph ← leftind)
	for (k, v) in blocks(m)
		# println(k)
		if isodd(k.n)
			v .*= (-1)
		end
	end
	return m
end


