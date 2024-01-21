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


# # the 2 and 3-th indices are physical ones to be traced out
# # nt = true means normal trace aa*
# function _trace_physical(a::AbstractTensorMap{S, 3, 1}; nt::Bool=true) where {S}
# 	r = TensorMap(ds->zeros(scalartype(a), ds), space(a, 1) ← space(a, 4)' )
# 	I = sectortype(a)
# 	for (f1, f2) in fusiontrees(a)
# 		if f1.uncoupled[2] == f1.uncoupled[3]
# 			f0 = FusionTree{I}((f1.uncoupled[1],), f1.coupled, (f1.isdual[1],), ())
# 			coef = (isodd(f1.uncoupled[2].n) && (!nt)) ? -1 : 1
# 			@tensor r[f0, f2][1,3] += coef * a[f1, f2][1,2,2,3]
# 			# axpy!(1, tmp, r[f0, f2]) 
# 		end
# 	end
# 	return r
# end


# function swap12!(twositemps::AbstractTensorMap{S, 3, 1}) where S
# 	for (f1, f2) in fusiontrees(twositemps)
# 		if isodd(f1.uncoupled[1].n) && isodd(f1.uncoupled[2].n)
# 			twositemps[f1, f2] .*= -1
# 		end
# 	end	
# 	return twositemps
# end


# function _fuse_physical(m::AbstractTensorMap{S, 3, 1}) where S
# 	cod = ProductSpace{S}((space(m, 1), space(m, 2)))
# 	dom = ProductSpace{S}((space(m, 4)',))
# 	r = TensorMap(ds->zeros(scalartype(m), ds), cod ← dom)
# 	for (f1, f2) in fusiontrees(m)
# 		if f1.uncoupled[2] == f1.uncoupled[3]
# 			if iseven(f1.uncoupled[2].n)
# 				f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2]), f1.coupled, (f1.isdual[1],f1.isdual[2]))
# 				r[f1′, f2] .+= m[f1, f2][:, :, 1, :]
# 			end
# 		else
# 			isdual2 = isodd(f1.uncoupled[2].n) ? f1.isdual[2] : f1.isdual[3]
# 			f1′ = FusionTree((f1.uncoupled[1], Z2Irrep(1)), f1.coupled, (f1.isdual[1], isdual2))
# 			r[f1′, f2] .+= m[f1, f2][:, :, 1, :]
# 		end
# 	end
# 	return r
# end


# function swap34!(twositemps::AbstractTensorMap{S, 1, 3}) where S
# 	for (f1, f2) in fusiontrees(twositemps)
# 		if isodd(f2.uncoupled[2].n) && isodd(f2.uncoupled[3].n)
# 			twositemps[f1, f2] .*= -1
# 		end
# 	end	
# 	return twositemps
# end

# function _fuse_physical(m::AbstractTensorMap{S, 1, 3}) where S
# 	cod = ProductSpace{S}((space(m, 1),))
# 	dom = ProductSpace{S}((space(m, 3)', space(m, 4)'))
# 	r = TensorMap(ds->zeros(scalartype(m), ds), cod ← dom)
# 	for (f1, f2) in fusiontrees(m)
# 		if f2.uncoupled[1] == f2.uncoupled[2]
# 			if iseven(f2.uncoupled[1].n)
# 				f2′ = FusionTree((f2.uncoupled[2],f2.uncoupled[3]), f2.coupled, (f2.isdual[2],f2.isdual[3]))
# 				r[f1, f2′] .+= m[f1, f2][:, 1, :, :]
# 			end
# 		else
# 			isdual2 = isodd(f2.uncoupled[2].n) ? f2.isdual[2] : f2.isdual[1]
# 			f2′ = FusionTree((Z2Irrep(1), f2.uncoupled[3]), f2.coupled, (isdual2, f2.isdual[3]))
# 			r[f1, f2′] .+= m[f1, f2][:, 1, :, :]
# 		end
# 	end
# 	return r
# end

