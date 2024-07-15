struct MixedCorrelationFunction2{M<:AbstractMatrix} <: AbstractMixedCorrelationFunction
    # real time
	data::Matrix{M}


	function MixedCorrelationFunction2(data::Matrix{M}) where {M <: AbstractMatrix}
		(size(data) == (3,3)) || throw(ArgumentError("input data should be a 3×3 matrix"))
		new{M}(data)
	end
end

Cm(data::Matrix{<:AbstractMatrix}) = MixedCorrelationFunction2(data)

# function Base.show(io::IO, ::MIME"text/plain", A::MixedCorrelationFunction2)
#     print(io, "Mixed Correlation Function [$(isize(A))+$(rsize(A))]")
# end

Base.:+(A::MixedCorrelationFunction2, B::MixedCorrelationFunction2) = MixedCorrelationFunction2(A.data .+ B.data)


branch(x::MixedCorrelationFunction2; b1::Symbol, b2::Symbol) = x.data[branch_to_index(b1), branch_to_index(b2)]

function branch_to_index(b::Symbol)
	if b == :+
		return 1
	elseif b == :-
		return 2
	elseif b == :τ
		return 3
	else
		throw(ArgumentError("branch must be one of :+, :- or :τ"))
	end
end

index(A::MixedCorrelationFunction2, i::Int, j::Int; b1::Symbol, b2::Symbol) = branch(A, b1=b1, b2=b2)[i, j]

isize(x::MixedCorrelationFunction2) = size(branch(x, b1=:τ, b2=:τ), 1)
rsize(x::MixedCorrelationFunction2) = size(branch(x, b1=:+, b2=:+), 1)


