# # move left and right most physical indices after the i-th position
# function pif_util(data::Vector{M}, i::Int; trunc) where {M}
# 	L = length(data) 
# 	ndata = Vector{M}(undef, L+1)
# 	util = TensorMap(ds->ones(Float64, ds), oneunit(_ph) ← one(_ph))
# 	@tensor left[1,3;4,5,2] := util[1] * data[1][2,3,4,5]
# 	for j in 1:i-1
# 		ndata[j], left = move_to_right(left, data[j+1], trunc=trunc)
# 	end
# 	if i == L
# 		u, s, v, _ = tsvd(left, (1,2,4), (3,5); alg=SDD(), trunc=trunc)
# 		ndata[i] = permute(u, (1,2), (4,3))
# 		ndata[i+1] = @tensor tmp[1,4;5,3] := s[1,2] * v[2,3,4] * conj(util[5])
# 	else
# 		@tensor right[3,1,2;5,4] := data[L][1,2,3,4] * conj(util[5])
# 		for j in (L+1):-1:(i+3)
# 			right, ndata[j] = move_to_left(data[j-2], right, trunc=trunc)
# 		end
# 		if i == 0
# 			u, s, v, _ = tsvd(right, (1,2), (3,4,5); alg=SDD(), trunc=trunc)
# 			u = u * s
# 			ndata[2] = permute(v, (1,2), (3,4))
# 			ndata[1] = @tensor tmp[1,3;4,2] := util[1] * u[2,3,4]
# 		else
# 			ndata[i], ndata[i+1], ndata[i+2] = central(left, right, trunc=trunc)	
# 		end			
# 	end
# 	return ndata
# end

# # t1 is rank-4, t2 is rank-5
# function move_to_left(t1, t2; trunc)
# 	@tensor tmp[5,1,2,4; 6,7,8] := t1[1,2,3,4] * t2[5,3,6,7,8]
# 	u, s, v, err = tsvd(tmp; alg=SDD(), trunc=trunc)
# 	return permute(u*s, (1,2,3), (5,4)), permute(v, (1,2), (3,4))
# end

# # t1 is rank-5, t2 is rank-4
# function move_to_right(t1, t2; trunc)
# 	@tensor tmp[1,2,4; 6,7,8,5] := t1[1,2,3,4,5] * t2[3,6,7,8]
# 	u, s, v, err = tsvd(tmp; alg=SDD(), trunc=trunc)
# 	return permute(u, (1,2), (4,3)), permute(s * v, (1,2), (3,4,5))
# end

# # t1 is rank-5, t2 is rank-5
# function central(t1, t2; trunc)
# 	u1, s1, v1, err = tsvd(t1, (1,2,4), (3,5); alg=SDD(), trunc=trunc)
# 	v1 = s1 * v1
# 	u2, s2, v2, err = tsvd(t2, (1,2), (3,4,5); alg=SDD(), trunc=trunc)
# 	u2 = u2 * s2
# 	@tensor center[1,3;5,4] := v1[1,2,3] * u2[4,2,5]
# 	return permute(u1, (1,2), (4,3)), center, permute(v2, (1,2), (3,4))
# end


"""
	_latticeindex(x::GrassmannContourIndex; conj=x.conj)

Map a correlation-function node index to the corresponding lattice node
index. The imaginary-time branch of a Grassmann lattice carries one more
node than the correlation function (`kτ = Nτ + 1`), so its time index is
shifted by one; the real-time branches are unchanged.
"""
_latticeindex(x::GrassmannContourIndex; conj::Bool=x.conj) =
	GrassmannContourIndex(x.branch === :τ ? x.j + 1 : x.j, x.branch, x.band, conj)

# The influence-functional exponent pairs the row leg āᵢ with the column legs
# aⱼ, while the QuAPI correlation function is read with the a-ā convention
# (x[a, ā] = Δ); conj flips the conjugation flags when looking up the contraction.

"""
	partialif_hybrid(lattice::AbstractGrassmannLattice, i::GrassmannContourIndex,
					 cols::AbstractVector{<:GrassmannContourIndex}, corr)

Build the partial IF contribution of the `i`-th row of the hybridization
function, `exp(āᵢ Σⱼ Δᵢⱼ aⱼ)`, analytically as an MPO of bond dimension 2.

`i` is the row leg and `cols` are the column legs it couples to. For a normal
fermionic bath the row leg is the conjugate field `ā` (`conj=true`) and the
columns are `a` (`conj=false`); BCS baths additionally use the anomalous
combinations, where either leg can carry either conjugation. Each leg carries
its own branch/time/band, so imaginary, real and mixed contours are all covered
by this one function. `Δᵢⱼ` is read from `corr` as the contraction between the
row leg and each column leg.
"""
function partialif_hybrid(lattice::AbstractGrassmannLattice, i::GrassmannContourIndex,
						  cols::AbstractVector{<:GrassmannContourIndex}, corr::AbstractCorrelationFunction)
	row = lattice[_latticeindex(i)]
	col_pos = [lattice[_latticeindex(c)] for c in cols]
	# QuAPI 统一实时约定：跨支块值内嵌 Keldysh 负号，IF 指数用物理块，读数补回
	#（两支 Keldysh 轮廓才嵌入；由 lattice 的分支判断，虚时 (:τ,) / 混合 (+,−,τ) 不补）
	keldysh = branches(lattice) == (:+, :-)
	coefs = [corr[conj(i), conj(c)] * (keldysh && branch(i) !== branch(c) ? -1 : 1) for c in cols]
	mpo = partialmpo(row, col_pos, coefs)
	return mpo * vacuumstate(promote_type(scalartype(lattice), eltype(coefs)), lattice)
end

# BCS 的指数系数读取在 QuAPI 的 `getindex(::BCSCorrelationFunction, i, j)` 里
# 完成（与 Dense 的 `_delta_terms` 同一分解：四类单项式 (ā↑,a↑)/(ā↓,a↓)/(ā↑,ā↓)/
# (a↓,a↑) 全 + 系数，其余腿对为 0；有向的指数系数语义，无需普通 bath 的 conj
# 翻转技巧）。旧装配在读取处补 a-ā 符号、且不区分指数项与误配对，迫使 ↓↓ 正常
# 项写成 (a↓,ā↓)·(−Δ)——Grassmann 单项式意义上 = ↓↓ 核转置：U=0 时 band-2 只
# 贡献标量、在 Z 归约里相消，U≠0 起产生 δt 无关、随时间累积的系统偏差
#（2026-11-03 firstorderif bcs_realtime 基准定位，bcs_ed ED 判决）。
function partialif_hybrid(lattice::AbstractGrassmannLattice, i::GrassmannContourIndex,
						  cols::AbstractVector{<:GrassmannContourIndex}, corr::BCSCorrelationFunction)
	keldysh = branches(lattice) == (:+, :-)
	row = lattice[_latticeindex(i)]
	# 行腿自身不进列集合：(band2,·) 既是行又是列，同节点项 ξ² = 0 且位置不可重复
	col_pos = [lattice[_latticeindex(c)] for c in cols if c != i]
	coefs = [corr[i, c] * (keldysh && branch(i) !== branch(c) ? -1 : 1)
			 for c in cols if c != i]
	mpo = partialmpo(row, col_pos, coefs)
	return mpo * vacuumstate(promote_type(scalartype(lattice), eltype(coefs)), lattice)
end

"""
	partialif_hybrid(lattice::AbstractGrassmannLattice,
					 rows::AbstractVector{<:GrassmannContourIndex},
					 col::GrassmannContourIndex, corr)

Transposed counterpart of the row-based method: build
`exp(Σᵢ Δᵢ,col āᵢ a_col)`, i.e. one column leg coupling to many row legs.
Used by the stepper to add the terms of the newest time step.
"""
function partialif_hybrid(lattice::AbstractGrassmannLattice, rows::AbstractVector{<:GrassmannContourIndex},
						  col::GrassmannContourIndex, corr::BCSCorrelationFunction)
	# 转置（列驱动）装配的单项式次序与行驱动不同，旧的 (a↓,ā↓) 结构在 U≠0 有
	# 系统偏差；BCS 尚无步进式调用方，先显式拒绝而非静默产出错误 IF。
	error("partialif_hybrid(rows, col, corr::BCSCorrelationFunction): 列驱动的 BCS 装配未实现，请用行驱动 hybriddynamics")
end

function partialif_hybrid(lattice::AbstractGrassmannLattice, rows::AbstractVector{<:GrassmannContourIndex},
						  col::GrassmannContourIndex, corr::AbstractCorrelationFunction)
	col_pos = lattice[_latticeindex(col)]
	row_pos = [lattice[_latticeindex(r)] for r in rows]
	keldysh = branches(lattice) == (:+, :-)   # 跨支负号补偿，同 partialif_hybrid
	coefs = [corr[conj(r), conj(col)] * (keldysh && branch(r) !== branch(col) ? -1 : 1) for r in rows]
	mpo = partialmpo(col_pos, row_pos, -coefs)
	return mpo * vacuumstate(promote_type(scalartype(lattice), eltype(coefs)), lattice)
end

"""
	partialif_hybrid_naive(lattice::AbstractGrassmannLattice, i::GrassmannContourIndex,
						   cols::AbstractVector{<:GrassmannContourIndex}, corr;
						   trunc::TruncationScheme=Defaults.alg_trunc_if())

Same partial IF as [`partialif_hybrid`](@ref), but built by applying the
exponentials `exp(Δᵢⱼ āᵢ aⱼ)` one by one and canonicalizing after each step.
"""
function partialif_hybrid_naive(lattice::AbstractGrassmannLattice, i::GrassmannContourIndex,
						  cols::AbstractVector{<:GrassmannContourIndex}, corr::AbstractCorrelationFunction;
						  trunc::TruncationScheme=Defaults.alg_trunc_if())
	alg = Orthogonalize(orth=TK.SVD(), trunc=trunc)
	gmps = vacuumstate(promote_type(scalartype(lattice), scalartype(corr)), lattice)
	pos1 = lattice[_latticeindex(i)]
	for c in cols
		pos2 = lattice[_latticeindex(c)]
		v = corr[conj(i), conj(c)] * (branches(lattice) == (:+, :-) && branch(i) !== branch(c) ? -1 : 1)
		t = exp(GTerm(pos1, pos2, coeff=v))
		apply!(t, gmps)
		canonicalize!(gmps, alg=alg)
	end
	return gmps
end

function partialif_hybrid_naive(lattice::AbstractGrassmannLattice, i::GrassmannContourIndex,
						  cols::AbstractVector{<:GrassmannContourIndex}, corr::BCSCorrelationFunction;
						  trunc::TruncationScheme=Defaults.alg_trunc_if())
	alg = Orthogonalize(orth=TK.SVD(), trunc=trunc)
	gmps = vacuumstate(promote_type(scalartype(lattice), scalartype(corr)), lattice)
	pos1 = lattice[_latticeindex(i)]
	for c in cols
		(c == i || iszero(corr[i, c])) && continue
		pos2 = lattice[_latticeindex(c)]
		v = corr[i, c] * (branches(lattice) == (:+, :-) && branch(i) !== branch(c) ? -1 : 1)
		t = exp(GTerm(pos1, pos2, coeff=v))
		apply!(t, gmps)
		canonicalize!(gmps, alg=alg)
	end
	return gmps
end

# column legs shared by all rows of the partial IF: the imaginary-time legs
# first, then the two real-time branches interleaved per time step (the terms
# within one row commute, so only the set of legs matters)
function _columns(lattice::AbstractGrassmannLattice, band::Int)
	cols = GrassmannContourIndex[]
	if :τ in branches(lattice)
		for j in 1:lattice.Nτ
			push!(cols, GrassmannContourIndex(j, :τ, band, false))
		end
	end
	if :+ in branches(lattice)
		for j in 1:lattice.kt
			push!(cols, GrassmannContourIndex(j, :+, band, false))
			push!(cols, GrassmannContourIndex(j, :-, band, false))
		end
	end
	return cols
end


function partialmpo(row::Int, cols::Vector{Int}, coefs::Vector{<:Number})
	# println("row=", row, " cols ", cols)
	@assert length(cols) == length(coefs) >= 1
	# @assert issorted(cols)
	p = sortperm(cols)
	cols = cols[p]
	coefs = coefs[p]
	# I2 = one(JW)

    virtual = isomorphism(eltype(coefs), Z2Space(1=>1), Z2Space(1=>1))
    pspace = z2space()
    T = scalartype(virtual)
    @tensor m22I[1,3;2,4] := virtual[1,2] * I2[3,4] 
	if row < cols[1]
	    tmp = Matrix{Any}(undef, 1, 2)
	    tmp[1, 1] = I2
	    tmp[1, 2] = one(eltype(coefs)) * σ₊
	    mpoj = SparseMPOTensor(tmp, T, pspace)
	    data = [mpoj]
	    for i in 1:length(cols)-1
	        tmp = Matrix{Any}(undef, 2, 2)
	        tmp[1,1] = I2
	        tmp[2,2] = m22I
	        tmp[1,2] = 0.
	        tmp[2,1] = adjoint(σ₋) * coefs[i]
	        push!(data, SparseMPOTensor(tmp, T, pspace))
	    end
	    tmp = Matrix{Any}(undef, 2, 1)
	    tmp[1,1] = I2
	    tmp[2,1] = adjoint(σ₋) * coefs[end]
	    mpof = SparseMPOTensor(tmp, T, pspace)
	    push!(data, mpof)
	    positions = [row; cols]
	elseif row > cols[end]
	    tmp = Matrix{Any}(undef, 1, 2)
	    tmp[1, 1] = I2
	    tmp[1, 2] = (-coefs[1]) * σ₊
	    mpoj = SparseMPOTensor(tmp, T, pspace)
	    data = [mpoj]
	    for i in 2:length(cols)
	        tmp = Matrix{Any}(undef, 2, 2)
	        tmp[1,1] = I2
	        tmp[2,2] = m22I
	        tmp[2,1] = 0.
	        tmp[1,2] = (-coefs[i]) * σ₊
	        push!(data, SparseMPOTensor(tmp, T, pspace))
	    end	
	    tmp = Matrix{Any}(undef, 2, 1)
	    tmp[1,1] = I2
	    tmp[2,1] = adjoint(σ₋) * one(eltype(coefs))
	    mpof = SparseMPOTensor(tmp, T, pspace)
	    push!(data, mpof)
	    positions = [cols; row]
	else
	    tmp = Matrix{Any}(undef, 1, 2)
	    tmp[1, 1] = I2
	    tmp[1, 2] = (-coefs[1]) * σ₊
	    mpoj = SparseMPOTensor(tmp, T, pspace)
	    data = [mpoj]

		row_pos = findfirst(x -> row < x, cols)
		# println("row pos is ", row_pos)
		for i in 2:row_pos-1
	        tmp = Matrix{Any}(undef, 2, 2)
	        tmp[1,1] = I2
	        tmp[2,2] = m22I
	        tmp[2,1] = 0.
	        tmp[1,2] = (-coefs[i]) * σ₊
	        push!(data, SparseMPOTensor(tmp, T, pspace))			
		end
		tmp = Matrix{Any}(undef, 2, 2)
		tmp[1,1] = I2
		tmp[2,2] = 0.
		tmp[1,2] = σ₊
		tmp[2,1] = adjoint(σ₋) * one(eltype(coefs))
		push!(data, SparseMPOTensor(tmp, T, pspace))

	    for i in row_pos:length(cols)-1
	        tmp = Matrix{Any}(undef, 2, 2)
	        tmp[1,1] = I2
	        tmp[2,2] = m22I
	        tmp[1,2] = 0.
	        tmp[2,1] = adjoint(σ₋) * coefs[i]
	        push!(data, SparseMPOTensor(tmp, T, pspace))
	    end
	    tmp = Matrix{Any}(undef, 2, 1)
	    tmp[1,1] = I2
	    tmp[2,1] = adjoint(σ₋) * coefs[end]
	    mpof = SparseMPOTensor(tmp, T, pspace)
	    push!(data, mpof)
	    positions = insert!(copy(cols), row_pos, row)
	end
	mpo = MPO(MPOHamiltonian(data))
	return PartialMPO(mpo.data, positions)
end



function split_mpotensor(mpoj::MPOTensor, trunc)
	ph = z2space()
	f = isomorphism(scalartype(mpoj), fuse(ph, ph), ph ⊗ ph)
	@tensor mpoj6[1,5,7;6,3,8] := mpoj[1,2,3,4] * conj(f[2,5,6]) * f[4,7,8]
	u, s, v = tsvd!(mpoj6, trunc=trunc)
	ss = sqrt(s)
	return permute(u * ss, (1,2), (4,3)), permute(ss * v, (1,2), (3,4))
end


