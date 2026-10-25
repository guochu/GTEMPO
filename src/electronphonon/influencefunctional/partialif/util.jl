# 单条腿（FockContourIndex = 时间 + 分支 + band）的密度型 partial IF 构造
# + Fock lattice 上的统一入口。
#
# `partialif_hybrid` 是 electron-phonon 路径的公共入口（与费米子的同名函数
# 对应）：它按 `FockContourIndex` 接收一行/一列，内部再调用下面的
# `_partialif_densemps` 做**单个格点**的 density-MPS 构造（每条腿只构造一次）。

"""
	_partialif_densemps(lat_size::Int, row::Int, cols::Vector{Int}, coefs::Vector{<:Number})

单个 Fock 格点 `row` 的 partial IF：`exp(Σ_{c∈cols} coefs[c] n̂_row n̂_c)`
（Strathearn-Lovett 2018 的 density-MPS 构造，键维 2）。

`coefs` 是**未指数化**的 Δ；内部取 `exp.(coefs)` 作为
`Λ = e^{Δ}` 的 two-body 转移矩阵元。`cols` 必须包含 `row`。
"""
function _partialif_densemps(lat_size::Int, row::Int, cols::Vector{Int}, coefs::Vector{<:Number})
	# println("row=", row, " cols ", cols)
	@assert length(cols) == length(coefs) > 1
	# @assert (row in cols)
	# @assert issorted(cols)
	p = sortperm(cols)
	cols = cols[p]
	coefs = exp.(coefs[p])
	T = scalartype(coefs)

	function onebody(m)
		tmp = ones(T,2)
		tmp[1] = m
		return tmp
	end
	function twobody(m)
		tmp = ones(T, 2,2)
		tmp[1,1] = m
		return tmp
	end
	L = length(cols)
	mpsdata = Vector{Array{T, 3}}(undef, L)
	pos = findfirst(x->row==x, cols)
	isnothing(pos) && throw(ArgumentError("$(row) is not a member of $cols"))
	if pos == 1
		m = onebody(coefs[1])
		tmp = zeros(T, 2, 2)
		for i in 1:2
			tmp[i, i] = m[i]
		end
		mpsdata[1] = reshape(tmp, 1,2,2)
		for j in 2:L-1
			m = twobody(coefs[j])
			tmp = zeros(T, 2,2,2)
			for i in 1:2
				tmp[i,:,i] = m[i,:]
			end
			mpsdata[j] = tmp
		end
		m = twobody(coefs[L])
		mpsdata[L] = reshape(m,2,2,1)
	elseif pos == L
		m = twobody(coefs[1])
		mpsdata[1] = reshape(m, 1,2,2)
		for j in 2:L-1
			m = twobody(coefs[j])
			tmp = zeros(T, 2,2,2)
			for i in 1:2
				tmp[i,:,i] = m[i,:]
			end
			mpsdata[j] = tmp			
		end
		m = onebody(coefs[L])
		tmp = zeros(T, 2, 2)
		for i in 1:2
			tmp[i, i] = m[i]
		end
		mpsdata[L] = reshape(tmp,2,2,1)
	else
		m = twobody(coefs[1])
		mpsdata[1] = reshape(m, 1,2,2)
		for j in 2:pos-1
			m = twobody(coefs[j])
			tmp = zeros(T, 2,2,2)
			for i in 1:2
				tmp[i,:,i] = m[i,:]
			end
			mpsdata[j] = tmp
		end
		m = onebody(coefs[pos])
		tmp = zeros(T,2,2,2)
		for i in 1:2
			tmp[i,i,i] = m[i]
		end
		mpsdata[pos] = tmp
		for j in pos+1:L-1
			m = twobody(coefs[j])
			tmp = zeros(T, 2,2,2)
			for i in 1:2
				tmp[i,:,i] = m[i,:]
			end
			mpsdata[j] = tmp			
		end
		m = twobody(coefs[L])
		mpsdata[L] = reshape(m,2,2,1)
	end
	return _fit_to_full(lat_size, cols, mpsdata)
end

function _fit_to_full(L::Int, pos, mpsdata)
	r = similar(mpsdata, L)
	for j in 1:pos[1]-1
		r[j] = ones(1,2,1)
	end
	for j in pos[end]+1:L
		r[j] = ones(1,2,1)
	end
	tmp = zeros(2,2,2)
	for i in 1:2
		tmp[:,i,:] = FMA.isometry(2)
	end
	for j in pos[1]:pos[end]
		posj = findfirst(x->x==j, pos)
		r[j] = isnothing(posj) ? tmp : mpsdata[posj]
	end
	return FockMPS(r)
end

# ---------------------------------------------------------------- 统一入口 --

"""
	_focknumnodes(lattice::AbstractFockLattice, branch::Symbol) -> Int

某一 contour 分支的节点数：`:τ` 用 `Nτ`，`:+`/`:-` 用 `Nt`
（三种 Fock lattice 的 `Nτ`/`Nt` 都指向各自的步数）。
"""
_focknumnodes(lattice::AbstractFockLattice, branch::Symbol) = branch === :τ ? lattice.Nτ : lattice.Nt

"""
	_fockcolumns(lattice::AbstractFockLattice) -> Vector{FockContourIndex}

Fock lattice 的全部 contour 腿（分支 → 时间 → band），即 partial IF 的
列腿集合。
"""
function _fockcolumns(lattice::AbstractFockLattice)
	cols = FockContourIndex[]
	for branch in branches(lattice), j in 1:_focknumnodes(lattice, branch), band in 1:lattice.bands
		push!(cols, FockContourIndex(j, branch, band))
	end
	return cols
end

"""
	partialif_hybrid(lattice::AbstractFockLattice, rowind::FockContourIndex,
	                 colinds::AbstractVector{<:FockContourIndex},
	                 corr::AbstractCorrelationFunction) -> FockMPS

Fock lattice 上**一条腿**（一个 `(时间, 分支, band)`）的 partial IF：
`exp(Σ_{c∈colinds} Δ_{row,c} n̂_row n̂_c)`。由
[`_partialif_densemps`](@ref) 构造，每个 `(节点, band)` 只做一次
density-MPS 构造（键维 2），与旧的逐行实现开销相同。

`Δ_{row,c}` 由 `corr[rowind, colind]` 读取（`FockContourIndex` 的 `band`
不参与取值，普通声子 bath 与 band 无关）。
"""
function partialif_hybrid(lattice::AbstractFockLattice, rowind::FockContourIndex,
						  colinds::AbstractVector{<:FockContourIndex},
						  corr::AbstractCorrelationFunction)
	row = index(lattice, rowind.j, branch=rowind.branch, band=rowind.band)
	colpos = [index(lattice, c.j, branch=c.branch, band=c.band) for c in colinds]
	coefs = [corr[rowind, c] for c in colinds]
	return _partialif_densemps(length(lattice), row, colpos, coefs)
end
