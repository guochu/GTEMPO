#---------------------------------------------------------------
# Fermionic (Grassmann) execution methods of the TensorOperations
# interface, keyed on `GrassmannBackend`.
#
# All type-level methods of the TensorOperations interface (structures,
# types, allocation, costs, contractibility checks) are sign-independent and
# are inherited from the Z2Tensors implementations for plain
# `AbstractTensorMap`s. Only the three execution methods `tensoradd!`,
# `tensortrace!` and `tensorcontract!` carry fermionic signs, entering
# through `f_permute` (see grassmanntensor.jl).
#
# Plain `@tensor` expressions (without `GrassmannBackend`) are unaffected and
# keep the bosonic (sign-free) Z2Tensors semantics.
#---------------------------------------------------------------

function _canonicalize(p::Index2Tuple, t::AbstractTensorMap)
    p′ = linearize(p)
    p₁ = TupleTools.getindices(p′, codomainind(t))
    p₂ = TupleTools.getindices(p′, domainind(t))
    return (p₁, p₂)
end

# NOTE: the tensor arguments are annotated as `AbstractTensorMap{<:Number}`
# (a partial application constraining only the scalar type, which suffices
# since Z2Space is the only space type) rather than via the
# `AbstractParityTensorMap` alias: the alias produces a `where`-clause with a
# different variable order than the Z2Tensors methods, which would make these
# definitions ambiguous with theirs instead of overriding them.

# tensoradd!
function TO.tensoradd!(C::AbstractTensorMap{<:Number},
                       A::AbstractTensorMap{<:Number},
                       pA::Index2Tuple, conjA::Bool,
                       α::Number, β::Number,
                       ::GrassmannBackend, allocator)
    if conjA
        A′ = adjoint(A)
        pA′ = TK.adjointtensorindices(A, _canonicalize(pA, C))
        add_f_permute!(C, A′, pA′, α, β)
    else
        add_f_permute!(C, A, _canonicalize(pA, C), α, β)
    end
    return C
end

# tensortrace!
function TO.tensortrace!(C::AbstractTensorMap{<:Number},
                         A::AbstractTensorMap{<:Number},
                         p::Index2Tuple, q::Index2Tuple, conjA::Bool,
                         α::Number, β::Number,
                         ::GrassmannBackend, allocator)
    if conjA
        A′ = adjoint(A)
        p′ = TK.adjointtensorindices(A, _canonicalize(p, C))
        q′ = TK.adjointtensorindices(A, q)
        trace_permute!(C, A′, p′, q′, α, β)
    else
        trace_permute!(C, A, _canonicalize(p, C), q, α, β)
    end
    return C
end

# tensorcontract!
function TO.tensorcontract!(C::AbstractTensorMap{<:Number},
                            A::AbstractTensorMap{<:Number},
                            pA::Index2Tuple, conjA::Bool,
                            B::AbstractTensorMap{<:Number},
                            pB::Index2Tuple, conjB::Bool,
                            pAB::Index2Tuple, α::Number, β::Number,
                            ::GrassmannBackend, allocator)
    pAB′ = _canonicalize(pAB, C)
    if conjA && conjB
        A′ = adjoint(A)
        pA′ = TK.adjointtensorindices(A, pA)
        B′ = adjoint(B)
        pB′ = TK.adjointtensorindices(B, pB)
        contract!(C, A′, pA′, B′, pB′, pAB′, α, β)
    elseif conjA
        A′ = adjoint(A)
        pA′ = TK.adjointtensorindices(A, pA)
        contract!(C, A′, pA′, B, pB, pAB′, α, β)
    elseif conjB
        B′ = adjoint(B)
        pB′ = TK.adjointtensorindices(B, pB)
        contract!(C, A, pA, B′, pB′, pAB′, α, β)
    else
        contract!(C, A, pA, B, pB, pAB′, α, β)
    end
    return C
end

#----------------
# IMPLEMENTATONS
#----------------

# Trace implementation
#----------------------
function trace_permute!(tdst::AbstractParityTensorMap,
                        tsrc::AbstractParityTensorMap,
                        (p₁, p₂)::Index2Tuple{N₁,N₂},
                        (q₁, q₂)::Index2Tuple{N₃,N₃},
                        α::Number,
                        β::Number) where {N₁,N₂,N₃}
    @boundscheck begin
        all(i -> space(tsrc, p₁[i]) == space(tdst, i), 1:N₁) ||
            throw(SpaceMismatch("trace: tsrc = $(codomain(tsrc))←$(domain(tsrc)),
                    tdst = $(codomain(tdst))←$(domain(tdst)), p₁ = $(p₁), p₂ = $(p₂)"))
        all(i -> space(tsrc, p₂[i]) == space(tdst, N₁ + i), 1:N₂) ||
            throw(SpaceMismatch("trace: tsrc = $(codomain(tsrc))←$(domain(tsrc)),
                    tdst = $(codomain(tdst))←$(domain(tdst)), p₁ = $(p₁), p₂ = $(p₂)"))
        all(i -> space(tsrc, q₁[i]) == dual(space(tsrc, q₂[i])), 1:N₃) ||
            throw(SpaceMismatch("trace: tsrc = $(codomain(tsrc))←$(domain(tsrc)),
                    q₁ = $(q₁), q₂ = $(q₂)"))
    end
    if iszero(β)
        fill!(tdst, β)
    elseif β != 1
        mul!(tdst, β, tdst)
    end
    r₁ = (p₁..., q₁...)
    r₂ = (p₂..., q₂...)
    for (f₁, f₂) in fusiontrees(tsrc)
        for ((f₁′, f₂′), coeff) in f_permute(f₁, f₂, r₁, r₂)
            f₁′′, g₁ = split(f₁′, N₁)
            f₂′′, g₂ = split(f₂′, N₂)
            g₁ == g₂ || continue
            coeff *= dim(g₁.coupled) / dim(g₁.uncoupled[1])
            C = tdst[f₁′′, f₂′′]
            A = tsrc[f₁, f₂]
            α′ = α * coeff
            TO.tensortrace!(C, (p₁, p₂), A, (q₁, q₂), false, α′, true)
        end
    end
    return tdst
end

# Contract implementation
#-------------------------
# TODO: contraction with either A or B a rank (1, 1) tensor does not require to
# permute the fusion tree and should therefore be special cased. This will speed
# up MPS algorithms
function contract!(C::AbstractParityTensorMap,
                   A::AbstractParityTensorMap,
                   (oindA, cindA)::Index2Tuple{N₁,N₃},
                   B::AbstractParityTensorMap,
                   (cindB, oindB)::Index2Tuple{N₃,N₂},
                   (p₁, p₂)::Index2Tuple,
                   α::Number,
                   β::Number) where {N₁,N₂,N₃}

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

    dA, dB, dC = dim(A), dim(B), dim(C)

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
end

function _contract!(α, A::AbstractParityTensorMap, B::AbstractParityTensorMap,
                    β, C::AbstractParityTensorMap,
                    oindA::IndexTuple{N₁}, cindA::IndexTuple,
                    oindB::IndexTuple{N₂}, cindB::IndexTuple,
                    p₁::IndexTuple, p₂::IndexTuple) where {N₁,N₂}
    A′ = f_permute(A, (oindA, cindA))
    B′ = f_permute(B, (cindB, oindB))
    ipC = TupleTools.invperm((p₁..., p₂...))
    oindAinC = TupleTools.getindices(ipC, ntuple(n -> n, N₁))
    oindBinC = TupleTools.getindices(ipC, ntuple(n -> n + N₁, N₂))
    if TK.has_shared_permute(C, (oindAinC, oindBinC))
        C′ = f_permute(C, (oindAinC, oindBinC))
        mul!(C′, A′, B′, α, β)
    else
        C′ = A′ * B′
        add_f_permute!(C, C′, (p₁, p₂), α, β)
    end
    return C
end
