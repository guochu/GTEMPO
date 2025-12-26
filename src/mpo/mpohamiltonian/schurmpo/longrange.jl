abstract type AbstractLongRangeTerm end

space_l(x::AbstractLongRangeTerm) = isa(x.a, MPOTensor) ? space_l(x.a) : oneunit(spacetype(x))
space_r(x::AbstractLongRangeTerm) = isa(x.b, MPOTensor) ? space_r(x.b) : oneunit(spacetype(x))'
TK.spacetype(x::AbstractLongRangeTerm) = spacetype(typeof(x))
coeff(x::AbstractLongRangeTerm) = x.coeff
