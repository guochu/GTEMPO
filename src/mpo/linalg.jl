

# get identity operator

"""
    TK.id(m::MPO)
Retuen an identity MPO from a given MPO
"""
TK.id(m::MPO) = MPO([id(scalartype(m), oneunit(spacetype(m)) ⊗ space(item, 2) ) for item in m.data])

TK.id(m::PartialMPO) = PartialMPO([id(scalartype(m), oneunit(spacetype(m)) ⊗ space(item, 2) ) for item in m.data], positions(m))


function LinearAlgebra.lmul!(f::Number, h::Union{MPO, PartialMPO})
    if !isempty(h)
        h[1] *= f
    end
    _renormalize!(h, h[1], false)
    return h
end
function _renormalize!(psi::Union{MPO, PartialMPO}, r, normalize::Bool)
	normalize && LinearAlgebra.normalize!(r)
end
LinearAlgebra.rmul!(h::Union{MPO, PartialMPO}, f::Number) = lmul!(f, h)

Base.:*(h::Union{MPO, PartialMPO}, f::Number) = lmul!(f, copy(h))
Base.:*(f::Number, h::Union{MPO, PartialMPO}) = h * f
Base.:/(h::Union{MPO, PartialMPO}, f::Number) = h * (1/f)


function Base.:+(hA::PartialMPO, hB::PartialMPO)
    @assert !isempty(hA)
    (positions(hA) == positions(hB)) || throw(ArgumentError("only PartialMPOs with same positions are allowed to add"))
    T = promote_type(scalartype(hA), scalartype(hB))
    S = spacetype(hA)
    M = mpotensortype(S, T)
    scale_a = scale_b = one(T) 
    if length(hA) == 1
        return PartialMPO([scale_a * hA[1] + scale_b * hB[1]], positions(hA))
    end
    embedders = [left_embedders(T, space_l(hA[i]), space_l(hB[i])) for i in 2:length(hA)]
    r = Vector{M}(undef, length(hA))
    for i in 1:length(hA)
        if i == 1
            @tensor m1[-1 -2; -3 -4] := hA[i][-1,-2,1,-4] * (embedders[i][1])'[1, -3]
            @tensor m1[-1 -2; -3 -4] += scale_b * hB[i][-1,-2,1,-4] * (embedders[i][2])'[1, -3]
        elseif i == length(hA)
            @tensor m1[-1 -2; -3 -4] := scale_a * embedders[i-1][1][-1, 1] * hA[i][1,-2,-3,-4] 
            @tensor m1[-1 -2; -3 -4] += scale_b * embedders[i-1][2][-1, 1] * hB[i][1,-2,-3,-4] 
        else          
            @tensor m1[-1 -2; -3 -4] := scale_a * embedders[i-1][1][-1, 1] * hA[i][1,-2,2,-4] * (embedders[i][1])'[2, -3]
            @tensor m1[-1 -2; -3 -4] += scale_b * embedders[i-1][2][-1, 1] * hB[i][1,-2,2,-4] * (embedders[i][2])'[2, -3]
        end
        r[i] = m1 
    end 
    return PartialMPO(r, positions(hA))       
end


