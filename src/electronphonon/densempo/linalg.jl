function TK.lmul!(f::Number, h::Union{DenseMPO, PartialDenseMPO})
    if !isempty(h)
        h[1] *= f
    end
    return h
end
Base.:*(h::Union{DenseMPO, PartialDenseMPO}, f::Number) = lmul!(f, copy(h))
Base.:*(f::Number, h::Union{DenseMPO, PartialDenseMPO}) = h * f
Base.:/(h::Union{DenseMPO, PartialDenseMPO}, f::Number) = h * (1/f)


function Base.:*(h::DenseMPO, psi::FockMPS)
    @assert !isempty(h)
    (length(h) == length(psi)) || throw(DimensionMismatch("mpo mps size mismatch"))
    r = [@tensor tmp[-1 -2; -3 -4 -5] := a[-1, -3, -4, 1] * b[-2, 1, -5] for (a, b) in zip(h.data, psi.data)]
    return FockMPS([tie(item,(2,1,2)) for item in r], scaling=scaling(psi))
end

Base.:*(h::PartialDenseMPO, psi::FockMPS) = apply!(h, copy(psi))
function DMRG.apply!(h::PartialDenseMPO, psi::FockMPS)
    @assert positions(h)[end] <= length(psi)
    T = promote_type(scalartype(h), scalartype(psi))
    _start, _end = positions(h)[1], positions(h)[end]
    leftspace = 1

    for (i, pos) in enumerate(_start:_end)
        _pos = findfirst(x->x==pos, positions(h))
        if !isnothing(_pos)
            @tensor tmp[-1 -2; -3 -4 -5] := h[_pos][-1, -3, -4, 1] * psi[pos][-2, 1, -5]
            leftspace = space_r(h[_pos])
        else
            hj = _eye(T, leftspace, leftspace)
           @tensor tmp[-1 -2; -3 -4 -5] := hj[-1, -4] * psi[pos][-2, -3, -5]
        end
        psi[pos] = tie(tmp, (2, 1, 2))
    end
    unset_svectors!(psi)
    return psi
end

function Base.:+(hA::PartialDenseMPO, hB::PartialDenseMPO)
    @assert !isempty(hA)
    (positions(hA) == positions(hB)) || throw(ArgumentError("only PartialMPOs with same positions are allowed to add"))
    T = promote_type(scalartype(hA), scalartype(hB))
    if length(hA) == 1
        return PartialDenseMPO([hA[1] + hB[1]], positions(hA))
    end
    L = length(hA)
    r = Vector{Array{T, 4}}(undef, L)
    r[1] = cat(hA[1], hB[1], dims=3)
    r[L] = cat(hA[L], hB[L], dims=1)
    for i in 2:L-1
        r[i] = cat(hA[i], hB[i], dims=(1,3))
    end
    return PartialDenseMPO(r, positions(hA))       
end