"""
	struct MPOHamiltonian{M <: AbstractSparseMPOTensor}

A generic MPO which stores a chain of AbstractSparseMPOTensor (Matrix of MPOTensors)

For finite system, the first site tensor is understood as the first row of the 
first AbstractSparseMPOTensor, and the last site tensor is understood as the last 
column of the last AbstractSparseMPOTensor
"""
struct MPOHamiltonian{M <: AbstractSparseMPOTensor}
	data::Vector{M}

function MPOHamiltonian{M}(data::AbstractVector) where {M <: AbstractSparseMPOTensor}
	@assert !isempty(data) 
	(data[1].leftspaces == data[end].rightspaces) || throw(SpaceMismatch())
	for i in 1:length(data)-1
		(size(data[i], 2) == size(data[i+1], 1)) || throw(DimensionMismatch())

	end
	new{M}(convert(Vector{M}, data))
end

end

Base.length(x::MPOHamiltonian) = length(x.data)
Base.getindex(x::MPOHamiltonian, i::Int) = getindex(x.data, i)
Base.firstindex(x::MPOHamiltonian) = firstindex(x.data)
Base.lastindex(x::MPOHamiltonian) = lastindex(x.data)
Base.copy(x::MPOHamiltonian) = MPOHamiltonian(copy(x.data))

TK.spacetype(::Type{MPOHamiltonian{M}}) where M = spacetype(M)
TK.spacetype(x::MPOHamiltonian) = spacetype(typeof(x))
TK.scalartype(::Type{MPOHamiltonian{M}}) where {M} = scalartype(M)

MPOHamiltonian(data::AbstractVector{M}) where {M <: AbstractSparseMPOTensor} = MPOHamiltonian{M}(data)
function MPOHamiltonian(data::Vector{Matrix{Any}}, virtualspaces::Vector{Vector{S}}, pspaces::Vector{S}) where {S<:ElementarySpace}
	@assert length(data) == length(pspaces) == length(virtualspaces)-1
	return MPOHamiltonian([SparseMPOTensor(data[i], virtualspaces[i], virtualspaces[i+1], pspaces[i]) for i in 1:length(data)])
end
MPOHamiltonian(data::Vector{Matrix{Any}}) = MPOHamiltonian([SparseMPOTensor(item) for item in data])


function Base.getindex(x::MPOHamiltonian, i::Int, j::Int, k::Int)
	x[i][j, k]
end 

Base.getindex(x::MPOHamiltonian, i::Colon, j::Int, k::Int) = [getindex(x, i,j,k) for i in 1:length(x)]


"""
	MPO(h::MPOHamiltonian, L::Int) 
	
Conversion of an MPOHamiltonian into a finite dense MPO
"""
MPO(h::MPOHamiltonian{<:SchurMPOTensor}) = _tompo(h, 1, size(h[end], 2))
MPO(h::MPOHamiltonian{<:SparseMPOTensor}; rowl::Int=1, colr::Int=1) = _tompo(h, rowl, colr)

function _tompo(h::MPOHamiltonian, leftrow::Int, rightcol::Int) 
	L = length(h)
	(L >= 2) || throw(ArgumentError("size of MPO must at least be 2"))
	# isstrict(h) || throw(ArgumentError("only strict MPOHamiltonian is allowed"))
	T = scalartype(h)
	S = spacetype(h)

	mpotensors = Vector{mpotensortype(S, T)}(undef, L)
	embedders = [right_embedders(T, h[i].rightspaces...) for i in 1:length(h)-1]

	tmp = zeros(T, oneunit(S) ⊗ h[1].pspace ← space(embedders[1][1], 2)' ⊗ h[1].pspace )
	for i in 1:length(embedders[1])
		@tensor tmp[-1, -2, -3, -4] += h[1, leftrow, i][-1,-2,1,-4] * embedders[1][i][1, -3]
	end
	mpotensors[1] = tmp
	for n in 2:L-1
		tmp = zeros(T, space(mpotensors[n-1], 3)' ⊗ h[n].pspace ← space(embedders[n][1], 2)' ⊗ h[n].pspace )
		for (i, j) in opkeys(h[n])
			@tensor tmp[-1, -2, -3, -4] += conj(embedders[n-1][i][1, -1]) * h[n, i, j][1,-2,2,-4] * embedders[n][j][2, -3]
		end
		for (i, j) in scalkeys(h[n])
			# iden = h[n].Os[i, j] * isomorphism(T, h[n].pspace, h[n].pspace)
			# @tensor tmp[-1, -2, -3, -4] += conj(embedders[n-1][i][1, -1]) * embedders[n][j][1, -3] * iden[-2, -4] 
			@tensor tmp[-1, -2, -3, -4] += conj(embedders[n-1][i][1, -1]) * h[n, i, j][1,-2,2,-4] * embedders[n][j][2, -3]
		end
		mpotensors[n] = tmp
	end
	tmp = zeros(T, space(embedders[L-1][1], 2)' ⊗ h[L].pspace, h[L].rightspaces[rightcol] ⊗ h[L].pspace )
	# _a = size(h[L], 2)
	for i in 1:size(h[L], 1)
		@tensor tmp[-1, -2, -3, -4] += conj(embedders[L-1][i][1, -1]) * h[L, i, rightcol][1,-2,-3,-4]
	end
	mpotensors[L] = tmp
	return MPO(mpotensors)
end
