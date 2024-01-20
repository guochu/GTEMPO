# function contract_center(left::AbstractTensorMap{S, 1, 1}, right::AbstractTensorMap) where S 
# 	@tensor r = left[1,2] * right[2,1]
# 	return r
# end

function contract_center(left::AbstractTensorMap{S, 1, N}, right::AbstractTensorMap{S, N, 1}) where {S, N}
	@assert space(left, 1) == space(right, N+1)'
	r = TensorMap(zeros, scalartype(left), space(left, 1), space(right, N+1)')
	cindA = ntuple(x->x+1, N)
	cindB = ntuple(x->N-x+1, N)
	TK.contract!(r, left, ((1,), cindA), right, (cindB, (N+1,)), ((1,), (2,)), true, false)
	return tr(r)
end 


# DMRG.l_LL(x::GrassmannMPS, y::GrassmannMPS, z::GrassmannMPS...) = error("l_LL not implemented for $(2+length(z)) mps")
# DMRG.r_RR(x::GrassmannMPS, y::GrassmannMPS, z::GrassmannMPS...) = error("r_RR not implemented for $(2+length(z)) mps")

# function DMRG.l_LL(x::GrassmannMPS)
# 	vacuum = space_l(x)
# 	left = isomorphism(Matrix{scalartype(x)}, one(vacuum), vacuum)
# 	return left
# end

# function DMRG.l_LL(x::GrassmannMPS, y::GrassmannMPS)
# 	vacuum = space_l(x)
# 	left = isomorphism(Matrix{scalartype(x)}, vacuum', vacuum)
# 	return left
# end

# function DMRG.r_RR(x::GrassmannMPS)
# 	vacuum = space_l(x)
# 	left = isomorphism(Matrix{scalartype(x)}, vacuum, one(vacuum))
# 	return left
# end

# function DMRG.r_RR(x::GrassmannMPS, y::GrassmannMPS)
# 	vacuum = space_l(x)
# 	right = isomorphism(Matrix{scalartype(x)}, vacuum, vacuum')
# 	return right
# end

# function DMRG.l_LL(x::GrassmannMPS, y::GrassmannMPS, z::GrassmannMPS)
# 	vacuum = space_l(x)
# 	left = isomorphism(Matrix{scalartype(x)}, vacuum' ⊗ vacuum',  vacuum )
# 	return left
# end

# function DMRG.r_RR(x::GrassmannMPS, y::GrassmannMPS, z::GrassmannMPS)
# 	vacuum = space_l(x)
# 	right = isomorphism(Matrix{scalartype(x)}, vacuum ⊗ vacuum,  vacuum')
# 	return right
# end
# function contract_center(left::AbstractTensorMap{S, 2, 1}, right::AbstractTensorMap) where S 
# 	@tensor r = left[1,2,3] * right[3,2,1]
# 	return r
# end

# function DMRG.l_LL(x::GrassmannMPS, y::GrassmannMPS, z::GrassmannMPS, u::GrassmannMPS)
# 	vacuum = space_l(x)
# 	left = isomorphism(Matrix{scalartype(x)}, vacuum' ⊗ vacuum' ,  vacuum ⊗ vacuum )
# 	return left
# end

# function DMRG.l_LL(x::GrassmannMPS, y::GrassmannMPS, z::GrassmannMPS, u::GrassmannMPS, v::GrassmannMPS, w::GrassmannMPS)
# 	vacuum = space_l(x)
# 	left = isomorphism(Matrix{scalartype(x)}, vacuum' ⊗ vacuum' ⊗ vacuum',  vacuum ⊗ vacuum ⊗ vacuum )
# 	return left
# end

# function DMRG.l_LL(x::GrassmannMPS, y::GrassmannMPS, z::GrassmannMPS, u::GrassmannMPS, v::GrassmannMPS, w::GrassmannMPS, r::GrassmannMPS)
# 	vacuum = space_l(x)
# 	left = isomorphism(Matrix{scalartype(x)}, vacuum' ⊗ vacuum' ⊗ vacuum' ⊗ vacuum',  vacuum ⊗ vacuum ⊗ vacuum )
# 	return left
# end

_mult_A(t, A::GrassmannMPS) = t * A
_mult_A(t, A::Vector{<:GrassmannMPS}) = [t * Aj for Aj in A]


# the 2 and 3-th indices are physical ones to be traced out
# nt = true means normal trace aa*
function _trace_physical(a::AbstractTensorMap{S, 3, 1}; nt::Bool=true) where {S}
	r = TensorMap(ds->zeros(scalartype(a), ds), space(a, 1) ← space(a, 4)' )
	I = sectortype(a)
	for (f1, f2) in fusiontrees(a)
		if f1.uncoupled[2] == f1.uncoupled[3]
			f0 = FusionTree{I}((f1.uncoupled[1],), f1.coupled, (f1.isdual[1],), ())
			coef = (isodd(f1.uncoupled[2].n) && (!nt)) ? -1 : 1
			@tensor r[f0, f2][1,3] += coef * a[f1, f2][1,2,2,3]
			# axpy!(1, tmp, r[f0, f2]) 
		end
	end
	return r
end


function swap12!(twositemps::AbstractTensorMap{S, 3, 1}) where S
	for (f1, f2) in fusiontrees(twositemps)
		if isodd(f1.uncoupled[1].n) && isodd(f1.uncoupled[2].n)
			twositemps[f1, f2] .*= -1
		end
	end	
	return twositemps
end


function _fuse_physical(m::AbstractTensorMap{S, 3, 1}) where S
	cod = ProductSpace{S}((space(m, 1), space(m, 2)))
	dom = ProductSpace{S}((space(m, 4)',))
	r = TensorMap(ds->zeros(scalartype(m), ds), cod ← dom)
	for (f1, f2) in fusiontrees(m)
		if f1.uncoupled[2] == f1.uncoupled[3]
			if iseven(f1.uncoupled[2].n)
				f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2]), f1.coupled, (f1.isdual[1],f1.isdual[2]))
				r[f1′, f2] .+= m[f1, f2][:, :, 1, :]
			end
		else
			isdual2 = isodd(f1.uncoupled[2].n) ? f1.isdual[2] : f1.isdual[3]
			f1′ = FusionTree((f1.uncoupled[1], Z2Irrep(1)), f1.coupled, (f1.isdual[1], isdual2))
			r[f1′, f2] .+= m[f1, f2][:, :, 1, :]
		end
	end
	return r
end


function swap34!(twositemps::AbstractTensorMap{S, 1, 3}) where S
	for (f1, f2) in fusiontrees(twositemps)
		if isodd(f2.uncoupled[2].n) && isodd(f2.uncoupled[3].n)
			twositemps[f1, f2] .*= -1
		end
	end	
	return twositemps
end

function _fuse_physical(m::AbstractTensorMap{S, 1, 3}) where S
	cod = ProductSpace{S}((space(m, 1),))
	dom = ProductSpace{S}((space(m, 3)', space(m, 4)'))
	r = TensorMap(ds->zeros(scalartype(m), ds), cod ← dom)
	for (f1, f2) in fusiontrees(m)
		if f2.uncoupled[1] == f2.uncoupled[2]
			if iseven(f2.uncoupled[1].n)
				f2′ = FusionTree((f2.uncoupled[2],f2.uncoupled[3]), f2.coupled, (f2.isdual[2],f2.isdual[3]))
				r[f1, f2′] .+= m[f1, f2][:, 1, :, :]
			end
		else
			isdual2 = isodd(f2.uncoupled[2].n) ? f2.isdual[2] : f2.isdual[1]
			f2′ = FusionTree((Z2Irrep(1), f2.uncoupled[3]), f2.coupled, (isdual2, f2.isdual[3]))
			r[f1, f2′] .+= m[f1, f2][:, 1, :, :]
		end
	end
	return r
end

# fuse i and i+1 into a single index
function g_fuse(m::AbstractTensorMap{S, M, N}, i::Int) where {S, M, N}
	@assert (i != M) && (i < M+N)
	@assert space(m, i) == space(m, i+1)
	# @assert (i < M) || (M <= i < N)
	local tmp
	if i < M
		idx = ntuple(x -> (x <= i) ? x : x+1, M-1)
		cod = ProductSpace{S,M-1}(map(n -> space(m, n), idx))
		dom = domain(m)

		tmp = TensorMap(zeros, scalartype(m), cod ← dom) 
		for (f1, f2) in fusiontrees(m)
			n = f1.uncoupled[i].n + f1.uncoupled[i+1].n
			uncoupled = map(n->f1.uncoupled[n], idx)
			isdual = map(n->f1.isdual[n], idx)
			if n == 0
				f1′ = FusionTree(uncoupled, f1.coupled, isdual)
				tmp[f1′, f2] .+= StridedView(dropdims(m[f1, f2], dims=i+1) )
			elseif n == 1
				isdual2 = ifelse(isodd(f1.uncoupled[i].n), f1.isdual[i], f1.isdual[i+1])
				uncoupled = TupleTools.setindex(uncoupled, Z2Irrep(1), i)
				isdual = TupleTools.setindex(isdual, isdual2, i)
				f1′ = FusionTree(uncoupled, f1.coupled, isdual)
				tmp[f1′, f2] .+= StridedView(dropdims(m[f1, f2], dims=i+1))
			end
		end
	else
		i2 = i - M
		idx = ntuple(x -> (x <= i2) ? x : x+1, N-1)
		cod = codomain(m)
		dom = ProductSpace{S,N-1}(map(n->domain(m)[n], idx))

		tmp = TensorMap(zeros, scalartype(m), cod ← dom) 
		for (f1, f2) in fusiontrees(m)
			n = f2.uncoupled[i2].n + f2.uncoupled[i2+1].n
			uncoupled = map(n->f2.uncoupled[n], idx)
			isdual = map(n->f2.isdual[n], idx)
			if n == 0
				f2′ = FusionTree(uncoupled, f2.coupled, isdual)
				tmp[f1, f2′] .+= StridedView(dropdims(m[f1, f2], dims=i+1))
			elseif n == 1
				isdual2 = ifelse(isodd(f2.uncoupled[i2].n), f2.isdual[i2], f2.isdual[i2+1])
				uncoupled = TupleTools.setindex(uncoupled, Z2Irrep(1), i2)
				isdual = TupleTools.setindex(isdual, isdual2, i2)
				f2′ = FusionTree(uncoupled, f2.coupled, isdual)
				tmp[f1, f2′] .+= StridedView(dropdims(m[f1, f2], dims=i+1))
			end
		end	
	end
	return tmp
end

function g_trace(m::AbstractTensorMap{S, M, N}, i::Int) where {S, M, N}
	@assert (i != M) && (i < M+N)
	@assert space(m, i) == space(m, i+1)
	local tmp
	if i < M
		idx = ntuple(x -> (x < i) ? x : x+2, M-2)
		cod = ProductSpace{S,M-2}(map(n -> space(m, n), idx))
		dom = domain(m)

		tmp = TensorMap(zeros, scalartype(m), cod ← dom) 

		for (f1, f2) in fusiontrees(m)
			if f1.uncoupled[i] == f1.uncoupled[i+1]
				uncoupled = map(n->f1.uncoupled[n], idx)
				isdual = map(n->f1.isdual[n], idx)

				f0 = FusionTree(uncoupled, f1.coupled, isdual)
				tmp[f0, f2] += StridedView(dropdims(m[f1, f2], dims=(i, i+1)))
			end
		end	
	else
		i2 = i - M
		idx = ntuple(x -> (x < i2) ? x : x+2, N-2)
		cod = codomain(m)
		dom = ProductSpace{S,N-2}(map(n->domain(m)[n], idx))

		tmp = TensorMap(zeros, scalartype(m), cod ← dom) 

		for (f1, f2) in fusiontrees(m)
			if f2.uncoupled[i2] == f2.uncoupled[i2+1]
				uncoupled = map(n->f2.uncoupled[n], idx)
				isdual = map(n->f2.isdual[n], idx)

				f0 = FusionTree(uncoupled, f2.coupled, isdual)

				tmp[f1, f0] += StridedView(dropdims(m[f1, f2], dims=(i, i+1)))

			end
		end	

	end

	return tmp
end

