abstract type FockOrdering end
abstract type ImagFockOrdering <: FockOrdering end
abstract type RealFockOrdering <: FockOrdering end
abstract type MixedFockOrdering <: FockOrdering end


LayoutStyle(x::FockOrdering) = LayoutStyle(typeof(x))
similargrassmannordering(::Type{T}) where {T<:FockOrdering} = error("similargrassmannordering not implemented for FockOrdering type $T")
similargrassmannordering(x::FockOrdering) = similargrassmannordering(typeof(x))

ImaginaryTimeOrderingStyle(x::FockOrdering) = ImaginaryTimeOrderingStyle(typeof(x))
RealTimeOrderingStyle(x::FockOrdering) = RealTimeOrderingStyle(typeof(x))
TimeOrderingStyle(x::FockOrdering) = TimeOrderingStyle(typeof(x))

"""
	struct M1N1 <: ImagFockOrdering
"""
struct M1N1 <: ImagFockOrdering end
LayoutStyle(::Type{M1N1}) = TimeLocalLayout()
ImaginaryTimeOrderingStyle(::Type{<:ImagFockOrdering}) = TimeDscending()
TimeOrderingStyle(::Type{O}) where {O<:ImagFockOrdering} = ImaginaryTimeOrderingStyle(O)
const MN = M1N1
similargrassmannordering(::Type{M1N1}) = Ā2A1B̄2B1()

"""
	struct M1m1N1n1 <: RealFockOrdering 
"""
struct M1m1N1n1 <: RealFockOrdering end
LayoutStyle(::Type{M1m1N1n1}) = TimeLocalLayout()
RealTimeOrderingStyle(::Type{<:RealFockOrdering}) = TimeDscending()
TimeOrderingStyle(::Type{O}) where {O<:RealFockOrdering} = RealTimeOrderingStyle(O)
const MmNn = M1m1N1n1
similargrassmannordering(::Type{M1m1N1n1}) = Ā2A1ā1a2B̄2B1b̄1b̄2()

"""
	struct M1N1_m1M1n1N1m2M2n2N2 <: MixedFockOrdering
"""
struct M1N1_m1M1n1N1m2M2n2N2 <: MixedFockOrdering end
LayoutStyle(::Type{M1N1_m1M1n1N1m2M2n2N2}) = TimeLocalLayout()
RealTimeOrderingStyle(::Type{M1N1_m1M1n1N1m2M2n2N2}) = TimeAscending()
ImaginaryTimeOrderingStyle(::Type{<:MixedFockOrdering}) = TimeDscending()
const MN_MmNn = M1N1_m1M1n1N1m2M2n2N2
similargrassmannordering(::Type{M1N1_m1M1n1N1m2M2n2N2}) = Ā3A2B̄3B2Ā2A1B̄2B1_ā1a2Ā2A1b̄1b2B̄2B1ā2a3Ā3A2b̄2b3B̄3B2()
