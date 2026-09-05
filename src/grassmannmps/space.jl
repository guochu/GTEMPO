# to do, generalization to SU₂ symmetry
const _ph = Z2Space(0=>1, 1=>1)

z2space() = _ph
z2spacetype() = spacetype(_ph)
# grassmannvacuum() =  Z2Space(0=>1)
z2mpstensortype(::Type{T}) where {T<:Number} = mpstensortype(z2spacetype(), T)
z2mpotensortype(::Type{T}) where {T<:Number} = mpotensortype(z2spacetype(), T)
