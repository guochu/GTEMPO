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


# `ContourIndex` has been replaced by QuAPI's `GrassmannContourIndex`
# (`j, branch, band, conj`). Lattice indexing is done through
# `lattice[GrassmannContourIndex(...)]`; the contour ordering of these
# indices is defined here because it depends on GTEMPO's contour layout.

# struct BoundaryIndex
# 	band::Int
# 	conj::Bool
# end
# BoundaryIndex(; band::Int, conj::Bool) = BoundaryIndex(band, conj)
# timestep(x::BoundaryIndex) = 0

function Base.isless(a::GrassmannContourIndex, b::GrassmannContourIndex)
	((a.j == 0) || (b.j == 0)) && throw(ArgumentError("boundary index are not contour-ordered"))
	# (a.band == b.band) || throw(ArgumentError("band mismatch"))
	if (a.branch == b.branch) && (a.j == b.j)
		(a.conj == b.conj) && throw(ArgumentError("the two indices are the same"))
		return a.conj
	end

	if a.branch == :+
		if b.branch == :+
			return a.j < b.j
		else
			return true
		end
	elseif a.branch == :-
		if b.branch == :+
			return false
		elseif b.branch == :-
			return a.j > b.j
		else
			return true
		end
	else
		if b.branch == :τ
			return a.j < b.j
		else
			return false
		end
	end
end
Base.:<(a::GrassmannContourIndex, b::GrassmannContourIndex) = isless(a, b)