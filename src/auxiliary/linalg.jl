

function _distance2(x, y)
	sA = real(dot(x, x))
	sB = real(dot(y, y))
	c = dot(x, y)
	r = sA+sB-2*real(c)
	return abs(r)
end

_distance(x, y) = sqrt(_distance2(x, y))


function stable_tsvd(m::AbstractTensorMap, args...; trunc::TruncationScheme=NoTruncation())
	try
		return tsvd(m, args...; trunc=trunc, alg=TK.SDD())
	catch
		return tsvd(m, args...; trunc=trunc, alg=TK.SVD())
	end
end

function stable_tsvd!(m::AbstractTensorMap; trunc::TruncationScheme=NoTruncation())
	try
		return tsvd!(copy(m), trunc=trunc, alg=TK.SDD())
	catch
		return tsvd!(m, trunc=trunc, alg=TK.SVD())
	end
end


# loose_isometry(cod::TensorSpace, dom::TensorSpace) =
#     loose_isometry(Matrix{Float64}, cod, dom)
# loose_isometry(P::TensorSpace) = loose_isometry(codomain(P), domain(P))
# loose_isometry(A::Type{<:DenseMatrix}, P::TensorSpace) =
#     loose_isometry(A, codomain(P), domain(P))
# function loose_isometry(::Type{A},
#                     cod::TensorSpace,
#                     dom::TensorSpace) where {A<:DenseMatrix}
#     t = TensorMap{scalartype(A)}(undef, cod, dom)
#     for (c, b) in blocks(t)
#         TK.MatrixAlgebra.one!(b)
#     end
#     return t
# end

function right_embedders(::Type{T}, a::S...) where {T <: Number, S <: ElementarySpace}
    V = ⊕(a...) 
    ts = [zeros(T, aj, V) for aj in a]
    for c in sectors(V)
    	n = 0
    	for i in 1:length(ts)
    		ni = dim(a[i], c)
    		block(ts[i], c)[:, (n+1):(n+ni)] .= Diagonal( ones(ni) )
    		n += ni
    	end
    end
    return ts
end

function left_embedders(::Type{T}, a::S...) where {T <: Number, S <: ElementarySpace}
    V = ⊕(a...) 
    ts = [zeros(T, V, aj) for aj in a]
    for c in sectors(V)
    	n = 0
    	for i in 1:length(ts)
    		ni = dim(a[i], c)
    		block(ts[i], c)[(n+1):(n+ni), :] .= Diagonal( ones(ni) )
    		n += ni
    	end
    end
    return ts	
end



# TK.:⊗(::Type{I}) where {I<:Sector} = (one(I),)
# TK.:⊗(::Type{I}, a::I, rest::Vararg{I}) where {I<:Sector} = ⊗(a, rest...)

