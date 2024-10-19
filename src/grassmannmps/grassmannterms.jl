abstract type AbstractGTerm end


function spin_site_ops_z2()
    ph = _ph
    vacuum = oneunit(ph)
    σ₊ = TensorMap(zeros, vacuum ⊗ ph ← Rep[ℤ₂](1=>1) ⊗ ph)
    blocks(σ₊)[Irrep[ℤ₂](1)] = ones(1, 1)
    σ₋ = TensorMap(zeros, vacuum ⊗ ph ← Rep[ℤ₂](1=>1) ⊗ ph)
    blocks(σ₋)[Irrep[ℤ₂](0)] = ones(1, 1)
    σz = TensorMap(ones, ph ← ph)
    blocks(σz)[Irrep[ℤ₂](0)] = -ones(1, 1)
    JW = -σz
    n = TensorMap(zeros, ph ← ph)
    blocks(n)[Irrep[ℤ₂](1)] = ones(1, 1)
    return σ₊, σ₋,σz, JW, one(JW), n
end

const σ₊, σ₋, σz, JW, I2, n̂ = spin_site_ops_z2()

struct GTerm{N, T <: Number} <: AbstractGTerm
	positions::NTuple{N, Int}
	coeff::T

function GTerm(positions::NTuple{N, Int}; coeff::Number) where {N}
	(length(Set(positions)) == N) || throw(ArgumentError("multiple grassmann on the same position not allowed"))
	(N % 2 == 0) || throw(ArgumentError("only support even number of grassmann variables"))
	p = sortperm([positions...])
	positions = positions[p]
	p = Permutation(p)
	coeff *= sign(p)
	new{N, typeof(coeff)}(positions, coeff)
end
end

TK.scalartype(::Type{GTerm{N, T}}) where {N, T} = T
DMRG.positions(x::GTerm) = x.positions
GTerm(a::Vararg{Int}; kwargs...) = GTerm(a; kwargs...)

function Base.convert(::Type{<:PartialMPO}, x::GTerm)
	n = x.positions[end] - x.positions[1] + 1
	pos = collect(x.positions[1]:x.positions[end])
	A = mpotensortype(grassmannpspacetype(), scalartype(x))
	ops = Vector{A}(undef, n)
	ops[1] = σ₊
	ops[end] = adjoint(σ₋)
	rest_pos = x.positions[2:end-1]
	_s = 1
	leftspace = space_r(ops[1])'
	for (i, pj) in enumerate(x.positions[1]+1:x.positions[end]-1)
		if pj in rest_pos
			_s *= -1
			ops[i+1] = ifelse(_s == -1, adjoint(σ₋), σ₊)
			leftspace = space_r(ops[i+1])'
		else
			mj = ifelse(_s == -1, I2, JW)
			id2 = id(storagetype(A), leftspace)
			ops[i+1] = @tensor tmp[1,3;2,4] := id2[1,2] * mj[3,4]
		end
	end
	return PartialMPO(ops, pos) * x.coeff
end


struct ExpGTerm{N, T <:Number} <: AbstractGTerm
	data::GTerm{N, T}
end
TK.scalartype(::Type{ExpGTerm{N, T}}) where {N, T} = T
DMRG.positions(x::ExpGTerm) = x.data.positions
Base.exp(x::GTerm{N}) where N = ExpGTerm(x)
function Base.convert(::Type{<:PartialMPO}, x::ExpGTerm)
	t1 = convert(PartialMPO, x.data)
	return t1 + id(t1)
end
