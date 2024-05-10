# Implement full TensorOperations.jl interface
#----------------------------------------------
TO.tensorstructure(t::GrassmannTensorMap) = space(t.data)
function TO.tensorstructure(t::GrassmannTensorMap, iA::Int, conjA::Symbol)
    return conjA == :N ? space(t.data, iA) : conj(space(t.data, iA))
end

TK.storagetype(::Type{GrassmannTensorMap{P}}) where P = storagetype(P)
TK.spacetype(::Type{GrassmannTensorMap{P}}) where P = spacetype(P)
TK.spacetype(t::GrassmannTensorMap) = spacetype(typeof(t))
TK.scalartype(::Type{GrassmannTensorMap{P}}) where P = scalartype(P)
TK.has_shared_permute(t::GrassmannTensorMap, p::Index2Tuple) = TK.has_shared_permute(t.data, p)

function TO.tensoralloc(::Type{TT}, structure, istemp=false,
                        backend::Backend...) where {TT<:GrassmannTensorMap}
    function blockallocator(d)
        return TO.tensoralloc(storagetype(TT), d, istemp, backend...)
    end
    return GrassmannTensorMap(TensorMap(blockallocator, structure))
end

function TO.tensorfree!(t::GrassmannTensorMap, backend::Backend...)
    for (c, b) in blocks(t.data)
        TO.tensorfree!(b, backend...)
    end
    return nothing
end

TO.tensorscalar(t::GrassmannTensorMap) = scalar(t.data)

function _canonicalize(p::Index2Tuple{N₁,N₂},
                       ::GrassmannTensorMap) where {N₁,N₂}
    return p
end
function _canonicalize(p::Index2Tuple, t::GrassmannTensorMap)
    p′ = linearize(p)
    p₁ = TupleTools.getindices(p′, codomainind(t.data))
    p₂ = TupleTools.getindices(p′, domainind(t.data))
    return (p₁, p₂)
end

@propagate_inbounds function add_permute!(tdst::GrassmannTensorMap,
                                          tsrc::GrassmannTensorMap,
                                          p::Index2Tuple{N₁,N₂},   
                                          α::Number,
                                          β::Number,
                                          backend::Backend...) where {N₁,N₂}
    add_f_permute!(tdst.data, tsrc.data, p, α, β, backend...)
end

# tensoradd!
function TO.tensoradd!(C::GrassmannTensorMap, pC::Index2Tuple,
                       A::GrassmannTensorMap, conjA::Symbol,
                       α::Number, β::Number, backend::Backend...) 
    if conjA == :N
        A′ = A
        pC′ = _canonicalize(pC, C)
    elseif conjA == :C
        A′ = adjoint(A)
        pC′ = TK.adjointtensorindices(A, _canonicalize(pA, C))
    else
        throw(ArgumentError("unknown conjugation flag $conjA"))
    end
    add_permute!(C, A′, pC′, α, β, backend...)
    return C
end

function TO.tensoradd_type(TC, ::Index2Tuple{N₁,N₂}, A::GrassmannTensorMap,
                           ::Symbol) where {N₁,N₂}
    M = TK.similarstoragetype(A.data, TC)
    return GrassmannTensorMap(tensormaptype(spacetype(A.data), N₁, N₂, M))
end

function TO.tensoradd_structure(pC::Index2Tuple{N₁,N₂},
                                A::GrassmannTensorMap, conjA::Symbol) where {N₁,N₂}
    S = spacetype(A.data)
    if conjA == :N
        cod = ProductSpace{S,N₁}(space.(Ref(A), pC[1]))
        dom = ProductSpace{S,N₂}(dual.(space.(Ref(A), pC[2])))
        return dom → cod
    else
        return TO.tensoradd_structure(TK.adjointtensorindices(A, pC), adjoint(A), :N)
    end
end

# tensortrace!
function TO.tensortrace!(C::GrassmannTensorMap, p::Index2Tuple,
                         A::GrassmannTensorMap, q::Index2Tuple, conjA::Symbol,
                         α::Number, β::Number, backend::Backend...) 
    if conjA == :N
        A′ = A
        p′ = _canonicalize(p, C)
        q′ = q
    elseif conjA == :C
        A′ = adjoint(A)
        p′ = TK.adjointtensorindices(A.data, _canonicalize(p, C))
        q′ = TK.adjointtensorindices(A.data, q)
    else
        throw(ArgumentError("unknown conjugation flag $conjA"))
    end
    # TODO: novel syntax for tensortrace?
    # tensortrace!(C, pC′, A′, qA′, α, β, backend...)
    trace_permute!(C, A′, p′, q′, α, β, backend...)
    return C
end

# tensorcontract!
function TO.tensorcontract!(C::GrassmannTensorMap, pAB::Index2Tuple,
                            A::GrassmannTensorMap, pA::Index2Tuple, conjA::Symbol,
                            B::GrassmannTensorMap, pB::Index2Tuple, conjB::Symbol,
                            α::Number, β::Number, backend::Backend...) 
    pAB′ = _canonicalize(pAB, C)
    if conjA == :N
        A′ = A
        pA′ = pA
    elseif conjA == :C
        A′ = A'
        pA′ = TK.adjointtensorindices(A.data, pA)
    else
        throw(ArgumentError("unknown conjugation flag $conjA"))
    end
    if conjB == :N
        B′ = B
        pB′ = pB
    elseif conjB == :C
        B′ = B'
        pB′ = TK.adjointtensorindices(B.data, pB)
    else
        throw(ArgumentError("unknown conjugation flag $conjB"))
    end
    contract!(C, A′, pA′, B′, pB′, pAB′, α, β, backend...)
    return C
end

function TO.tensorcontract_type(TC, ::Index2Tuple{N₁,N₂},
                                A::GrassmannTensorMap, pA, conjA,
                                B::GrassmannTensorMap, pB, conjB) where {N₁,N₂}
    M = TK.similarstoragetype(A.data, TC)
    M == TK.similarstoragetype(B.data, TC) || throw(ArgumentError("incompatible storage types"))
    return GrassmannTensorMap{tensormaptype(spacetype(A), N₁, N₂, M)}
end

function TO.tensorcontract_structure(pC::Index2Tuple{N₁,N₂},
                                     A::GrassmannTensorMap, pA::Index2Tuple, conjA,
                                     B::GrassmannTensorMap, pB::Index2Tuple,
                                     conjB) where {N₁,N₂}
    S = spacetype(A)
    spaces1 = TO.flag2op(conjA).(space.(Ref(A.data), pA[1]))
    spaces2 = TO.flag2op(conjB).(space.(Ref(B.data), pB[2]))
    spaces = (spaces1..., spaces2...)
    cod = ProductSpace{S,N₁}(getindex.(Ref(spaces), pC[1]))
    dom = ProductSpace{S,N₂}(dual.(getindex.(Ref(spaces), pC[2])))
    return dom → cod
end

function TO.checkcontractible(tA::GrassmannTensorMap, iA::Int, conjA::Symbol,
                              tB::GrassmannTensorMap, iB::Int, conjB::Symbol,
                              label) 
    sA = TO.tensorstructure(tA, iA, conjA)'
    sB = TO.tensorstructure(tB, iB, conjB)
    sA == sB ||
        throw(SpaceMismatch("incompatible spaces for $label: $sA ≠ $sB"))
    return nothing
end

TO.tensorcost(t::GrassmannTensorMap, i::Int) = dim(space(t.data, i))

#----------------
# IMPLEMENTATONS
#----------------

# Trace implementation
#----------------------
function trace_permute!(tdst::GrassmannTensorMap,
                        tsrc::GrassmannTensorMap,
                        (p₁, p₂)::Index2Tuple{N₁,N₂},
                        (q₁, q₂)::Index2Tuple{N₃,N₃},
                        α::Number,
                        β::Number,
                        backend::Backend...) where {N₁,N₂,N₃}
    @boundscheck begin
        all(i -> space(tsrc.data, p₁[i]) == space(tdst.data, i), 1:N₁) ||
            throw(SpaceMismatch("trace: tsrc = $(codomain(tsrc.data))←$(domain(tsrc.data)),
                    tdst = $(codomain(tdst.data))←$(domain(tdst.data)), p₁ = $(p₁), p₂ = $(p₂)"))
        all(i -> space(tsrc.data, p₂[i]) == space(tdst.data, N₁ + i), 1:N₂) ||
            throw(SpaceMismatch("trace: tsrc = $(codomain(tsrc.data))←$(domain(tsrc.data)),
                    tdst = $(codomain(tdst.data))←$(domain(tdst.data)), p₁ = $(p₁), p₂ = $(p₂)"))
        all(i -> space(tsrc.data, q₁[i]) == dual(space(tsrc.data, q₂[i])), 1:N₃) ||
            throw(SpaceMismatch("trace: tsrc = $(codomain(tsrc.data))←$(domain(tsrc.data)),
                    q₁ = $(q₁), q₂ = $(q₂)"))
    end
    S = spacetype(tdst)
    I = sectortype(S)
    # TODO: is it worth treating UniqueFusion separately? Is it worth to add multithreading support?
    cod = codomain(tsrc.data)
    dom = domain(tsrc.data)
    n = length(cod)
    if iszero(β)
        fill!(tdst.data, β)
    elseif β != 1
        mul!(tdst.data, β, tdst.data)
    end
    r₁ = (p₁..., q₁...)
    r₂ = (p₂..., q₂...)
    for (f₁, f₂) in fusiontrees(tsrc.data)
        for ((f₁′, f₂′), coeff) in f_permute(f₁, f₂, r₁, r₂)
            f₁′′, g₁ = split(f₁′, N₁)
            f₂′′, g₂ = split(f₂′, N₂)
            g₁ == g₂ || continue
            coeff *= dim(g₁.coupled) / dim(g₁.uncoupled[1])
            # for i in 2:length(g₁.uncoupled)
            #     if !(g₁.isdual[i])
            #         coeff *= twist(g₁.uncoupled[i])
            #     end
            # end
            C = tdst.data[f₁′′, f₂′′]
            A = tsrc.data[f₁, f₂]
            α′ = α * coeff
            TO.tensortrace!(C, (p₁, p₂), A, (q₁, q₂), :N, α′, true, backend...)
        end
    end
    return tdst
end

# Contract implementation
#-------------------------
# TODO: contraction with either A or B a rank (1, 1) tensor does not require to
# permute the fusion tree and should therefore be special cased. This will speed
# up MPS algorithms
function contract!(C::GrassmannTensorMap,
                   A::GrassmannTensorMap,
                   (oindA, cindA)::Index2Tuple{N₁,N₃},
                   B::GrassmannTensorMap,
                   (cindB, oindB)::Index2Tuple{N₃,N₂},
                   (p₁, p₂)::Index2Tuple,
                   α::Number,
                   β::Number,
                   backend::Backend...) where {N₁,N₂,N₃}
    
    # find optimal contraction scheme
    hsp = TK.has_shared_permute
    ipC = TupleTools.invperm((p₁..., p₂...))
    oindAinC = TupleTools.getindices(ipC, ntuple(n -> n, N₁))
    oindBinC = TupleTools.getindices(ipC, ntuple(n -> n + N₁, N₂))

    qA = TupleTools.sortperm(cindA)
    cindA′ = TupleTools.getindices(cindA, qA)
    cindB′ = TupleTools.getindices(cindB, qA)

    qB = TupleTools.sortperm(cindB)
    cindA′′ = TupleTools.getindices(cindA, qB)
    cindB′′ = TupleTools.getindices(cindB, qB)

    dA, dB, dC = dim(A.data), dim(B.data), dim(C.data)

    # keep order A en B, check possibilities for cind
    memcost1 = memcost2 = dC * (!hsp(C, (oindAinC, oindBinC)))
    memcost1 += dA * (!hsp(A, (oindA, cindA′))) +
                dB * (!hsp(B, (cindB′, oindB)))
    memcost2 += dA * (!hsp(A, (oindA, cindA′′))) +
                dB * (!hsp(B, (cindB′′, oindB)))

    # reverse order A en B, check possibilities for cind
    memcost3 = memcost4 = dC * (!hsp(C, (oindBinC, oindAinC)))
    memcost3 += dB * (!hsp(B, (oindB, cindB′))) +
                dA * (!hsp(A, (cindA′, oindA)))
    memcost4 += dB * (!hsp(B, (oindB, cindB′′))) +
                dA * (!hsp(A, (cindA′′, oindA)))

    if min(memcost1, memcost2) <= min(memcost3, memcost4)
        if memcost1 <= memcost2
            return _contract!(α, A, B, β, C, oindA, cindA′, oindB, cindB′, p₁, p₂)
        else
            return _contract!(α, A, B, β, C, oindA, cindA′′, oindB, cindB′′, p₁, p₂)
        end
    else
        p1′ = map(n -> ifelse(n > N₁, n - N₁, n + N₂), p₁)
        p2′ = map(n -> ifelse(n > N₁, n - N₁, n + N₂), p₂)
        if memcost3 <= memcost4
            return _contract!(α, B, A, β, C, oindB, cindB′, oindA, cindA′, p1′, p2′)
        else
            return _contract!(α, B, A, β, C, oindB, cindB′′, oindA, cindA′′, p1′, p2′)
        end
    end
end

# TODO: also transform _contract! into new interface, and add backend support
function _contract!(α, A::GrassmannTensorMap, B::GrassmannTensorMap,
                    β, C::GrassmannTensorMap,
                    oindA::IndexTuple{N₁}, cindA::IndexTuple,
                    oindB::IndexTuple{N₂}, cindB::IndexTuple,
                    p₁::IndexTuple, p₂::IndexTuple) where {N₁,N₂}
    copyA = false
    A′ = permute(A, (oindA, cindA); copy=copyA)
    B′ = permute(B, (cindB, oindB))
    ipC = TupleTools.invperm((p₁..., p₂...))
    oindAinC = TupleTools.getindices(ipC, ntuple(n -> n, N₁))
    oindBinC = TupleTools.getindices(ipC, ntuple(n -> n + N₁, N₂))
    if TK.has_shared_permute(C, (oindAinC, oindBinC))
        C′ = permute(C, (oindAinC, oindBinC))
        mul!(C′.data, A′.data, B′.data, α, β)
    else
        C′ = GrassmannTensorMap(A′.data * B′.data)
        add_permute!(C, C′, (p₁, p₂), α, β)
    end
    return C
end

# Scalar implementation
#-----------------------
scalar(t::GrassmannTensorMap) = scalar(t.data)

