abstract type AbstractNTerm end

struct ExpNTerm{N, T <: Number} <: AbstractNTerm
	positions::NTuple{N, Int}
	coeff::T

function ExpNTerm(positions::NTuple{N, Int}; coeff::Number) where {N}
	(length(Set(positions)) == N) || throw(ArgumentError("multiple n̂ on the same position not allowed"))
	p = sortperm([positions...])
	positions = positions[p]
	new{N, typeof(coeff)}(positions, coeff)
end
end

TK.scalartype(::Type{ExpNTerm{N, T}}) where {N, T} = T
DMRG.positions(x::ExpNTerm) = x.positions
ExpNTerm(a::Vararg{Int}; kwargs...) = ExpNTerm(a; kwargs...)


# function Base.convert(::Type{<:PartialDenseMPO}, x::NTerm)
# 	nop2 = Matrix{scalartype(x)}([0 0; 1 1])
# 	nop = reshape(nop2, 1,2,1,2)

# 	return PartialDenseMPO([nop for i in 1:length(positions(x))], [positions(x)...]) * x.coeff
# end

function Base.convert(::Type{<:PartialDenseMPO}, x::ExpNTerm{N, T}) where {N, T}
	Λ = exp(x.coeff)
	if N == 1
		nop = reshape(Matrix{T}([Λ 0; 0 1]), 1,2,1,2)
		return PartialDenseMPO([nop], [positions(x)[1]])
	else
		I2 = Matrix{T}([1 0; 0 1])
		ml = Matrix{T}([Λ 1; 1 1])
		opl = zeros(T, 1, 2, 2, 2)
		for i in 1:2
			opl[1,i,:,i] = ml[i,:]
		end
		@tensor opm[1,3,2,4] := I2[1,2] * I2[3,4] 
		opr = zeros(T, 2,2,1,2)
		for i in 1:2
			opr[:,i,1,i] = I2[:,i]
		end
		oplist = [opl, ntuple(v->opm, N-2)..., opr]
		return PartialDenseMPO(oplist, [positions(x)...])
	end
end

function DMRG.apply!(x::ExpNTerm, mps::FockMPS) 
    all(s -> 1 <= s <= length(mps), positions(x)) || throw(BoundsError())
    t = convert(PartialDenseMPO, x)
    apply!(t, mps)
    return mps
end

