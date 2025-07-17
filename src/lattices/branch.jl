# abstract type Branch end
# abstract type RealTimeBranch <: Branch end
# abstract type ImaginaryTimeBranch <: Branch end

# struct ForwardBranch <: RealTimeBranch end
# struct BackwardBranch <: RealTimeBranch end
# struct ImaginaryBranch <: ImaginaryTimeBranch end

# branch(b::Symbol) = _branch(b)
# function _branch(b::Symbol)
# 	if b == :+
# 		return ForwardBranch()
# 	elseif b ==:-
# 		return BackwardBranch()
# 	elseif b == :τ
# 		return ImaginaryBranch()
# 	else
# 		throw(ArgumentError("branch must be one of :+, :- or :τ"))
# 	end
# end

# struct TimeStep{B <: Branch}
# 	j::Int
# 	branch::B
# end
# TimeStep(j::Int; branch::Branch) = TimeStep(j, branch)


abstract type AbstractLatticeIndex end
# band(x::AbstractLatticeIndex) = x.band

struct ContourIndex <: AbstractLatticeIndex
	j::Int
	band::Int
	conj::Bool
	branch::Symbol


function ContourIndex(j::Int, band::Int, conj::Bool, branch::Symbol)
	(branch in (:+, :-, :τ)) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
	new(j, band, conj, branch)
end

end
ContourIndex(j::Int; conj::Bool, branch::Symbol=:τ, band::Int=1) = ContourIndex(j, band, conj, branch)
branch(x::ContourIndex) = x.branch

Base.:(==)(a::ContourIndex, b::ContourIndex) = (a.j == b.j) && (a.band == b.band) && (a.conj == b.conj) && (branch(a) == branch(b))

function Base.:<(a::ContourIndex, b::ContourIndex)
	((a.j == 0) || (b.j == 0)) && throw(ArgumentError("boundary index are not contour-ordered"))
	# (a.band == b.band) || throw(ArgumentError("band mismatch"))
	if (branch(a) == branch(b)) && (a.j == b.j)
		(a.conj == b.conj) && throw(ArgumentError("the two indices are the same"))
		return a.conj
	end

	if branch(a) == :+
		if branch(b) == :+
			return a.j < b.j
		else
			return true
		end
	elseif branch(a) == :-
		if branch(b) == :+
			return false
		elseif branch(b) == :-
			return a.j > b.j
		else
			return true
		end
	else
		if branch(b) == :τ
			return a.j < b.j
		else
			return false
		end
	end
end

# struct BoundaryIndex <: AbstractLatticeIndex
# 	band::Int
# 	conj::Bool
# end
# BoundaryIndex(; band::Int, conj::Bool) = BoundaryIndex(band, conj)
# timestep(x::BoundaryIndex) = 0