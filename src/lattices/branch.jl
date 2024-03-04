abstract type Branch end
abstract type RealTimeBranch <: Branch end
abstract type ImaginaryTimeBranch <: Branch end

struct ForwardBranch <: RealTimeBranch end
struct BackwardBranch <: RealTimeBranch end
struct ImaginaryBranch <: ImaginaryTimeBranch end

branch(b::Symbol) = _branch(b)
function _branch(b::Symbol)
	if b == :+
		return ForwardBranch()
	elseif b ==:-
		return BackwardBranch()
	elseif b == :τ
		return ImaginaryBranch()
	else
		throw(ArgumentError("branch must be one of :+, :- or :τ"))
	end
end

# struct TimeStep{B <: Branch}
# 	j::Int
# 	branch::B
# end
# TimeStep(j::Int; branch::Branch) = TimeStep(j, branch)


abstract type AbstractLatticeIndex end
# band(x::AbstractLatticeIndex) = x.band

struct TimeStepIndex{B <: Branch} <: AbstractLatticeIndex
	j::Int
	band::Int
	conj::Bool
	branch::B
end
TimeStepIndex(j::Int; band::Int, conj::Bool, branch::Branch) = TimeStepIndex(j, band, conj, Branch)
timestep(x::TimeStepIndex) = x.j


struct BoundaryIndex <: AbstractLatticeIndex
	band::Int
	conj::Bool
end
BoundaryIndex(; band::Int, conj::Bool) = BoundaryIndex(band, conj)
timestep(x::BoundaryIndex) = 0