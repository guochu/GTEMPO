
g_fuse(t::GrassmannTensorMap, i::Int) = GrassmannTensorMap(g_fuse(t.data, i))

# fuse i and i+1 into a single index
function g_fuse(m::AbstractTensorMap{T, S, M, N}, i::Int) where {T<:Number, S, M, N}
	@assert (i != M) && (i < M+N)
	@assert space(m, i) == space(m, i+1)
	# @assert (i < M) || (M <= i < N)
	local tmp
	if i < M
		idx = ntuple(x -> (x <= i) ? x : x+1, M-1)
		cod = ProductSpace{S,M-1}(map(n -> space(m, n), idx))
		dom = domain(m)

		# tmp = zeros(scalartype(m), cod ← dom) 
		tmp = fill!(similar(m, cod ← dom), zero(T)) 
		for (f1, f2) in fusiontrees(m)
			n = f1.uncoupled[i].n + f1.uncoupled[i+1].n
			uncoupled = map(n->f1.uncoupled[n], idx)
			# isdual = map(n->f1.isdual[n], idx)
			if n == 0
				f1′ = FusionTree(uncoupled, f1.coupled)
				# tmp[f1′, f2] .+= StridedView(dropdims(m[f1, f2], dims=i+1) )
				axpy!(true, StridedView(dropdims(m[f1, f2], dims=i+1) ), tmp[f1′, f2])
			elseif n == 1
				# isdual2 = ifelse(isodd(f1.uncoupled[i].n), f1.isdual[i], f1.isdual[i+1])
				uncoupled = TupleTools.setindex(uncoupled, Z2Irrep(1), i)
				# isdual = TupleTools.setindex(isdual, isdual2, i)
				f1′ = FusionTree(uncoupled, f1.coupled)
				# tmp[f1′, f2] .+= StridedView(dropdims(m[f1, f2], dims=i+1))
				axpy!(true, StridedView(dropdims(m[f1, f2], dims=i+1)), tmp[f1′, f2])
			end
		end
	else
		i2 = i - M
		idx = ntuple(x -> (x <= i2) ? x : x+1, N-1)
		cod = codomain(m)
		dom = ProductSpace{S,N-1}(map(n->domain(m)[n], idx))

		# tmp = zeros(scalartype(m), cod ← dom) 
		tmp = fill!(similar(m, cod ← dom), zero(T)) 
		for (f1, f2) in fusiontrees(m)
			n = f2.uncoupled[i2].n + f2.uncoupled[i2+1].n
			uncoupled = map(n->f2.uncoupled[n], idx)
			# isdual = map(n->f2.isdual[n], idx)
			if n == 0
				f2′ = FusionTree(uncoupled, f2.coupled)
				# tmp[f1, f2′] .+= StridedView(dropdims(m[f1, f2], dims=i+1))
				axpy!(true, StridedView(dropdims(m[f1, f2], dims=i+1)), tmp[f1, f2′])
			elseif n == 1
				# isdual2 = ifelse(isodd(f2.uncoupled[i2].n), f2.isdual[i2], f2.isdual[i2+1])
				uncoupled = TupleTools.setindex(uncoupled, Z2Irrep(1), i2)
				# isdual = TupleTools.setindex(isdual, isdual2, i2)
				f2′ = FusionTree(uncoupled, f2.coupled)
				# tmp[f1, f2′] .+= StridedView(dropdims(m[f1, f2], dims=i+1))
				axpy!(true, StridedView(dropdims(m[f1, f2], dims=i+1)), tmp[f1, f2′])
			end
		end	
	end
	return tmp
end

g_trace(t::GrassmannTensorMap, i::Int) = GrassmannTensorMap(_g_trace(t.data, i))

function _g_trace(m::AbstractTensorMap{T, S, M, N}, i::Int) where {T<:Number, S, M, N}
	@assert (i != M) && (i < M+N)
	@assert space(m, i) == space(m, i+1)
	local tmp
	if i < M
		idx = ntuple(x -> (x < i) ? x : x+2, M-2)
		cod = ProductSpace{S,M-2}(map(n -> space(m, n), idx))
		dom = domain(m)

		# tmp = zeros(scalartype(m), cod ← dom) 
		tmp = fill!(similar(m, cod ← dom), zero(T)) 

		for (f1, f2) in fusiontrees(m)
			if f1.uncoupled[i] == f1.uncoupled[i+1]
				uncoupled = map(n->f1.uncoupled[n], idx)
				# isdual = map(n->f1.isdual[n], idx)

				f0 = FusionTree(uncoupled, f1.coupled)
				# tmp[f0, f2] += StridedView(dropdims(m[f1, f2], dims=(i, i+1)))
				axpy!(true, StridedView(dropdims(m[f1, f2], dims=(i, i+1))), tmp[f0, f2])
			end
		end	
	else
		i2 = i - M
		idx = ntuple(x -> (x < i2) ? x : x+2, N-2)
		cod = codomain(m)
		dom = ProductSpace{S,N-2}(map(n->domain(m)[n], idx))

		# tmp = zeros(scalartype(m), cod ← dom) 
		tmp = fill!(similar(m, cod ← dom), zero(T)) 

		for (f1, f2) in fusiontrees(m)
			if f2.uncoupled[i2] == f2.uncoupled[i2+1]
				uncoupled = map(n->f2.uncoupled[n], idx)
				# isdual = map(n->f2.isdual[n], idx)

				f0 = FusionTree(uncoupled, f2.coupled)

				# tmp[f1, f0] += StridedView(dropdims(m[f1, f2], dims=(i, i+1)))
				coeff = isodd(f2.uncoupled[i2].n) ? -1 : 1
				axpy!(coeff, StridedView(dropdims(m[f1, f2], dims=(i, i+1))), tmp[f1, f0])

			end
		end	

	end

	return tmp
end
