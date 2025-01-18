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

function TO.tensoralloc(::Type{TT},
                        structure::TensorMapSpace{S,N₁,N₂},
                        istemp::Val,
                        allocator=TO.DefaultAllocator()) where
         {S,N₁,N₂,TT<:GrassmannTensorMap}
    A = storagetype(TT)
    dim = TK.fusionblockstructure(structure).totaldim
    data = TO.tensoralloc(A, dim, istemp, allocator)
    # return TT(data, structure)
    return GrassmannTensorMap(TensorMap{scalartype(TT)}(data, structure))
end

function TO.tensorfree!(t::GrassmannTensorMap, allocator=TO.DefaultAllocator())
    TO.tensorfree!(t.data.data, allocator)
    return nothing
end

TO.tensorscalar(t::GrassmannTensorMap) = scalar(t.data)

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
                                          backend...) where {N₁,N₂}
    add_f_permute!(tdst.data, tsrc.data, p, α, β, backend...)
end

# tensoradd!
function TO.tensoradd!(C::GrassmannTensorMap,
                       A::GrassmannTensorMap, pA::Index2Tuple, conjA::Bool,
                       α::Number, β::Number,
                       backend, allocator)
    if conjA
        A′ = adjoint(A)
        pA′ = TK.adjointtensorindices(A.data, _canonicalize(pA, C))
        add_permute!(C, A′, pA′, α, β, backend)
    else
        add_permute!(C, A, _canonicalize(pA, C), α, β, backend)
    end
    return C
end

function TO.tensoradd_type(TC, A::GrassmannTensorMap, ::Index2Tuple{N₁,N₂},
                           ::Bool) where {N₁,N₂}
    M = TK.similarstoragetype(A.data, TC)
    return GrassmannTensorMap(tensormaptype(spacetype(A.data), N₁, N₂, M))
end

function TO.tensoradd_structure(A::GrassmannTensorMap, pA::Index2Tuple{N₁,N₂},
                                conjA::Bool) where {N₁,N₂}
    if !conjA
        # don't use `permute` as this is also used when indices are traced
        return TK.select(space(A.data), pA)
    else
        return TO.tensoradd_structure(adjoint(A), TK.adjointtensorindices(A.data, pA), false)
    end
end

# tensortrace!
function TO.tensortrace!(C::GrassmannTensorMap,
                         A::GrassmannTensorMap, p::Index2Tuple, q::Index2Tuple,
                         conjA::Bool,
                         α::Number, β::Number, backend, allocator) 
    if conjA
        A′ = adjoint(A)
        p′ = TK.adjointtensorindices(A.data, _canonicalize(p, C))
        q′ = TK.adjointtensorindices(A.data, q)
        trace_permute!(C, A′, p′, q′, α, β, backend)
    else
        trace_permute!(C, A, _canonicalize(p, C), q, α, β, backend)
    end
    return C
end

# tensorcontract!
function TO.tensorcontract!(C::GrassmannTensorMap,
                            A::GrassmannTensorMap, pA::Index2Tuple, conjA::Bool,
                            B::GrassmannTensorMap, pB::Index2Tuple, conjB::Bool,
                            pAB::Index2Tuple, α::Number, β::Number,
                            backend, allocator)
    pAB′ = _canonicalize(pAB, C)
    if conjA && conjB
        A′ = A'
        pA′ = TK.adjointtensorindices(A.data, pA)
        B′ = B'
        pB′ = TK.adjointtensorindices(B.data, pB)
        contract!(C, A′, pA′, B′, pB′, pAB′, α, β, backend, allocator)
    elseif conjA
        A′ = A'
        pA′ = TK.adjointtensorindices(A.data, pA)
        contract!(C, A′, pA′, B, pB, pAB′, α, β, backend, allocator)
    elseif conjB
        B′ = B'
        pB′ = TK.adjointtensorindices(B.data, pB)
        contract!(C, A, pA, B′, pB′, pAB′, α, β, backend, allocator)
    else
        contract!(C, A, pA, B, pB, pAB′, α, β, backend, allocator)
    end
    return C
end

function TO.tensorcontract_type(TC,
                                A::GrassmannTensorMap, ::Index2Tuple, ::Bool,
                                B::GrassmannTensorMap, ::Index2Tuple, ::Bool,
                                ::Index2Tuple{N₁,N₂}) where {N₁,N₂}
    M = TK.similarstoragetype(A.data, TC)
    M == TK.similarstoragetype(B.data, TC) ||
        throw(ArgumentError("incompatible storage types:\n$(M) ≠ $(TK.similarstoragetype(B, TC))"))
    spacetype(A) == spacetype(B) || throw(SpaceMismatch("incompatible space types"))
    return GrassmannTensorMap{tensormaptype(spacetype(A), N₁, N₂, M)}
end

function TO.tensorcontract_structure(A::GrassmannTensorMap, pA::Index2Tuple, conjA::Bool,
                                     B::GrassmannTensorMap, pB::Index2Tuple, conjB::Bool,
                                     pAB::Index2Tuple{N₁,N₂}) where {N₁,N₂}
    sA = TO.tensoradd_structure(A, pA, conjA)
    sB = TO.tensoradd_structure(B, pB, conjB)
    return permute(TK.compose(sA, sB), pAB)
end

function TO.checkcontractible(tA::GrassmannTensorMap, iA::Int, conjA::Bool,
                              tB::GrassmannTensorMap, iB::Int, conjB::Bool,
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
                        backend=TO.DefaultBackend()) where {N₁,N₂,N₃}
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
            TO.tensortrace!(C, (p₁, p₂), A, (q₁, q₂), :N, α′, true, backend)
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
                   backend, allocator) where {N₁,N₂,N₃}
    
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


    if memcost1 <= memcost2
        return _contract!(α, A, B, β, C, oindA, cindA′, oindB, cindB′, p₁, p₂)
    else
        return _contract!(α, A, B, β, C, oindA, cindA′′, oindB, cindB′′, p₁, p₂)
    end

    # # reverse order A en B, check possibilities for cind
    # memcost3 = memcost4 = dC * (!hsp(C, (oindBinC, oindAinC)))
    # memcost3 += dB * (!hsp(B, (oindB, cindB′))) +
    #             dA * (!hsp(A, (cindA′, oindA)))
    # memcost4 += dB * (!hsp(B, (oindB, cindB′′))) +
    #             dA * (!hsp(A, (cindA′′, oindA)))

    # if min(memcost1, memcost2) <= min(memcost3, memcost4)
    #     if memcost1 <= memcost2
    #         return _contract!(α, A, B, β, C, oindA, cindA′, oindB, cindB′, p₁, p₂)
    #     else
    #         return _contract!(α, A, B, β, C, oindA, cindA′′, oindB, cindB′′, p₁, p₂)
    #     end
    # else
    #     p1′ = map(n -> ifelse(n > N₁, n - N₁, n + N₂), p₁)
    #     p2′ = map(n -> ifelse(n > N₁, n - N₁, n + N₂), p₂)
    #     if memcost3 <= memcost4
    #         return _contract!(α, B, A, β, C, oindB, cindB′, oindA, cindA′, p1′, p2′)
    #     else
    #         return _contract!(α, B, A, β, C, oindB, cindB′′, oindA, cindA′′, p1′, p2′)
    #     end
    # end
end

# function contract!(C::GrassmannTensorMap,
#                    A::GrassmannTensorMap,
#                    (oindA, cindA)::Index2Tuple{N₁,N₃},
#                    B::GrassmannTensorMap,
#                    (cindB, oindB)::Index2Tuple{N₃,N₂},
#                    (p₁, p₂)::Index2Tuple,
#                    α::Number,
#                    β::Number,
#                    backend::Backend...) where {N₁,N₂,N₃}
    
#     # find optimal contraction scheme
#     hsp = TK.has_shared_permute
#     ipC = TupleTools.invperm((p₁..., p₂...))
#     oindAinC = TupleTools.getindices(ipC, ntuple(n -> n, N₁))
#     oindBinC = TupleTools.getindices(ipC, ntuple(n -> n + N₁, N₂))

#     qA = TupleTools.sortperm(cindA)
#     cindA′ = TupleTools.getindices(cindA, qA)
#     cindB′ = TupleTools.getindices(cindB, qA)

#     qB = TupleTools.sortperm(cindB)
#     cindA′′ = TupleTools.getindices(cindA, qB)
#     cindB′′ = TupleTools.getindices(cindB, qB)

#     dA, dB, dC = dim(A.data), dim(B.data), dim(C.data)

#     # keep order A en B, check possibilities for cind
#     memcost1 = memcost2 = dC * (!hsp(C, (oindAinC, oindBinC)))
#     memcost1 += dA * (!hsp(A, (oindA, cindA′))) +
#                 dB * (!hsp(B, (cindB′, oindB)))
#     memcost2 += dA * (!hsp(A, (oindA, cindA′′))) +
#                 dB * (!hsp(B, (cindB′′, oindB)))

#     # reverse order A en B, check possibilities for cind
#     memcost3 = memcost4 = dC * (!hsp(C, (oindBinC, oindAinC)))
#     memcost3 += dB * (!hsp(B, (oindB, cindB′))) +
#                 dA * (!hsp(A, (cindA′, oindA)))
#     memcost4 += dB * (!hsp(B, (oindB, cindB′′))) +
#                 dA * (!hsp(A, (cindA′′, oindA)))

#     if min(memcost1, memcost2) <= min(memcost3, memcost4)
#         if memcost1 <= memcost2
#             return _contract!(α, A, B, β, C, oindA, cindA′, oindB, cindB′, p₁, p₂)
#         else
#             return _contract!(α, A, B, β, C, oindA, cindA′′, oindB, cindB′′, p₁, p₂)
#         end
#     else
#         p1′ = map(n -> ifelse(n > N₁, n - N₁, n + N₂), p₁)
#         p2′ = map(n -> ifelse(n > N₁, n - N₁, n + N₂), p₂)
#         if memcost3 <= memcost4
#             return _contract!(α, B, A, β, C, oindB, cindB′, oindA, cindA′, p1′, p2′)
#         else
#             return _contract!(α, B, A, β, C, oindB, cindB′′, oindA, cindA′′, p1′, p2′)
#         end
#     end
# end

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
TK.scalar(t::GrassmannTensorMap) = scalar(t.data)

