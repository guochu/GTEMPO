"""
	focktograssmann(lattice::AbstractFockLattice, mps::FockMPS; trunc::TruncationScheme=DefaultKTruncation)
"""
function focktograssmann(lattice::AbstractFockLattice, mps::FockMPS)
	lattice2 = similargrassmannlattice(lattice)
	mps2 = vacuumstate(lattice2)
	for band in 1:lattice.bands, b in branches(lattice)
		if b == :τ
			k = lattice.Nτ
		else
			k = lattice.Nt
		end
		for j in 1:k
			focktograssmann_one_step!(lattice2, mps2, lattice, mps, j, band, b)
		end
	end
	new_scale = scaling(mps)^(length(mps) / length(mps2))
	setscaling!(mps2, new_scale)
	return lattice2, mps2
end


function focktograssmann(::Type{O}, lattice::AbstractFockLattice, mps::FockMPS; kwargs...) where {O<:GrassmannOrdering}
	lattice2, mps2 = focktograssmann(lattice, mps)
	if !(lattice2.ordering isa O)
		lattice2, mps2 = changeordering(O, lattice2, mps2; kwargs...)
	end
	return lattice2, mps2
end
focktograssmann(o::GrassmannOrdering, lattice::AbstractFockLattice, mps::FockMPS; kwargs...) = focktograssmann(typeof(o), lattice, mps; kwargs...)

function focktograssmann_one_step!(lattice2::AbstractGrassmannLattice, mps2::GrassmannMPS, lattice::AbstractFockLattice, mps::FockMPS, j::Int, band::Int, b::Symbol)
	ph = grassmannpspace()
	pos = index(lattice, j, band=band, branch=b)
	t = mps[pos]
	if b == :-
		pos1, pos2 = index(lattice2, j, band=band, branch=b, conj=true), index(lattice2, j+1, band=band, branch=b, conj=false)
	else
		pos1, pos2 = index(lattice2, j+1, band=band, branch=b, conj=true), index(lattice2, j, band=band, branch=b, conj=false)
	end
	# println("j=", j, ",band=", band, ",branch=", b, ",pos1=", pos1, ",pos2=", pos2)
	(pos1 + 1 == pos2) || error("something wrong")
	left = Rep[ℤ₂](0=>space_l(t))
	right = Rep[ℤ₂](0=>space_r(t))
	t2 = zeros(scalartype(mps), left ⊗ ph ← ph' ⊗ right)
	copy!(block(t2, Z2Irrep(0)), t[:, 1, :])
	copy!(block(t2, Z2Irrep(1)), t[:, 2, :])
	x, y = rightorth!(t2, alg=LQ())
	mps2[pos1] = x
	mps2[pos2] = permute(y, (1,2), (3,))
end

reweighting(lattice2::AbstractGrassmannLattice, mps2::GrassmannMPS, lattice::AbstractFockLattice, mps::FockMPS; kwargs...) = reweighting!(
			lattice2, copy(mps2), lattice, mps; kwargs...)

function reweighting!(lattice2::AbstractGrassmannLattice{O}, mps2::GrassmannMPS, lattice::AbstractFockLattice, mps::FockMPS; kwargs...) where O
	fo = similargrassmannordering(OrderingStyle(lattice))
	if fo isa O
		return _reweighting!(lattice2, mps2, lattice, mps)
	else
		lattice3, mps3 = changeordering(fo, lattice2, mps2; kwargs...)
		mps3 = _reweighting!(lattice3, mps3, lattice, mps)
		lattice2′, mps3′ = changeordering(O, lattice3, mps3; kwargs...)
		mps2.data[:] = mps3′.data[:]
		mps2.svectors[:] = mps3′.svectors[:]
		setscaling!(mps2, scaling(mps3′))
		return mps2
	end
end

function _reweighting!(lattice2::AbstractGrassmannLattice{O}, mps2::GrassmannMPS, lattice::AbstractFockLattice, mps::FockMPS) where O
	(similargrassmannordering(OrderingStyle(lattice)) isa O) || throw(ArgumentError("reweighting requires ordering matching"))
	for band in 1:lattice.bands, b in branches(lattice)
		if b == :τ
			k = lattice.Nτ
		else
			k = lattice.Nt
		end
		for j in 1:k
			reweighting_one_step!(lattice2, mps2, lattice, mps, j, band, b)
		end
	end
	new_scale = scaling(mps)^(length(mps) / length(mps2))
	setscaling!(mps2, new_scale)
	unset_svectors!(mps2)
	return mps2	
end


function reweighting_one_step!(lattice2::AbstractGrassmannLattice, mps2::GrassmannMPS, lattice::AbstractFockLattice, mps::FockMPS, j::Int, band::Int, b::Symbol)
	# (lattice.N == lattice2.N) || throw(ArgumentError("lattice size mismatch"))
	pos = index(lattice, j, band=band, branch=b)
	t = mps[pos]
	if b == :-
		pos1, pos2 = index(lattice2, j, band=band, branch=b, conj=true), index(lattice2, j+1, band=band, branch=b, conj=false)
	else
		pos1, pos2 = index(lattice2, j+1, band=band, branch=b, conj=true), index(lattice2, j, band=band, branch=b, conj=false)
	end
	(pos1 + 1 == pos2) || error("something wrong")

	s3 = space_r(t)
	t2 = zeros(scalartype(t), s3, 2, s3)
	for i in 1:s3
		t2[i, :, i] .= 1
	end
	mps2[pos1] = mult_f_g_tensor(mps2[pos1], t)
	mps2[pos2] = mult_f_g_tensor(mps2[pos2], t2)
end

function mult_f_g_tensor(g::AbstractParityTensorMap, f::AbstractArray{<:Number, 3})
	@assert size(f, 2) == 2
	al, a2, ar = size(f)
	sl, sr = space_l(g), space_r(g)'
	sl′ = spacetype(sl)(0=>dim(sl, Z2Irrep(0))*al, 1=>dim(sl, Z2Irrep(1))*al)
	sr′ = spacetype(sr)(0=>dim(sr, Z2Irrep(0))*ar, 1=>dim(sr, Z2Irrep(1))*ar)
	g′ = zeros(scalartype(g), sl′ ⊗ grassmannpspace() ← sr′)
	for (fl, fr) in fusiontrees(g′)
		if isodd(fl.uncoupled[2].n)
			g′[fl, fr] = kron(g[fl, fr], f[:, 1:1, :])
		else
			g′[fl, fr] = kron(g[fl, fr], f[:, 2:2, :])
		end
	end
	return g′
end
