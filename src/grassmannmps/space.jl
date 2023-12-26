# to do, generalization to SU₂ symmetry
const _ph = Rep[ℤ₂](0=>1, 1=>1)

grassmannpspace() = _ph
grassmannpspacetype() = spacetype(_ph)
# grassmannvacuum() =  Rep[ℤ₂](0=>1)
grassmannmpstensortype(::Type{T}) where {T<:Number} = mpstensortype(grassmannpspacetype(), T)
grassmannmpotensortype(::Type{T}) where {T<:Number} = mpotensortype(grassmannpspacetype(), T)
