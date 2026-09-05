const AbstractParityTensorMap{T, N₁, N₂} = AbstractTensorMap{T, N₁, N₂}

"""
    GrassmannBackend()

TensorOperations backend that implements `tensoradd!`, `tensortrace!` and
`tensorcontract!` for Z2-graded (parity) `TensorMap`s with the fermionic
(Grassmann) sign convention: every internal index permutation is carried out
by the fermionic `gpermute` instead of the bosonic `permute`. It is inserted
automatically by the [`@grassmann`](@ref) macro, and can equally be passed
explicitly as `backend = GrassmannBackend()` to the function-based
TensorOperations API.

The backend instance is also threaded through the internal helpers
`trace_permute!`, `contract!` and `_contract!` as their last argument, both
to supply the sign information and to mark these functions as carrying the
fermionic convention, in contrast to the sign-free Z2Tensors operations with
otherwise similar signatures.

In addition to the fermionic `gpermute` reordering signs, the backend applies
the fermionic twists that put every contracted pair into the canonical
`(a, ā)` order before contracting, in the same way as the GrassmannTensors /
TensorKit fermionic conventions:
* `contract!`: every contracted pair whose A-side leg carries a non-dual
  space (an `a`, contracted against an `ā` on the B side) contributes a twist
  of `-1` per odd sector, applied on the A side — independent of which
  region the B-side legs live in.
* `trace_permute!`: closing the trace loop, every traced codomain leg with a
  non-dual space and odd sector contributes a twist of `-1`.
"""
struct GrassmannBackend <: AbstractBackend end

function gpermute(t::AbstractParityTensorMap, (p₁, p₂)::Index2Tuple{N₁,N₂};
                  copy::Bool=false) where {N₁,N₂}
    cod = ProductSpace{N₁}(map(n -> space(t, n), p₁))
    dom = ProductSpace{N₂}(map(n -> dual(space(t, n)), p₂))
    # share data if possible
    if (!copy) && (p₁ === codomainind(t) && p₂ === domainind(t)) 
        return t
    end
    # general case
    @inbounds begin
        return gpermute!(similar(t, cod ← dom), t, (p₁, p₂))
    end
end
function gpermute(t::AdjointTensorMap, (p₁, p₂)::Index2Tuple; copy::Bool=false) 
    p₁′ = TK.adjointtensorindices(t, p₂)
    p₂′ = TK.adjointtensorindices(t, p₁)
    return adjoint(gpermute(adjoint(t), (p₁′, p₂′); copy=copy))
end
# convenience: two separate leg-index tuples
gpermute(t::AbstractParityTensorMap, p1::IndexTuple, p2::IndexTuple; kwargs...) = gpermute(t, (p1, p2); kwargs...)


@propagate_inbounds function gpermute!(tdst::AbstractParityTensorMap{<:Number, N₁, N₂},
                                       tsrc::AbstractParityTensorMap,
                                       p::Index2Tuple{N₁,N₂}) where {N₁,N₂}
    return add_gpermute!(tdst, tsrc, p, true, false)
end

"""
    g_twist!(t, inds)

Apply the fermion-parity (Z2) twist to the legs of the parity tensor `t` at
linearized positions `inds`: every fusion-tree block is multiplied by
``(-1)^{\\text{number of odd sectors among the twisted legs}}``, i.e. by
`-1` for each fermion (odd sector) on the twisted legs.
"""
function g_twist!(t::AbstractParityTensorMap, inds)
    isempty(inds) && return t
    N₁ = numout(t)
    for (f₁, f₂) in fusiontrees(t)
        θ = 1
        @inbounds for i in inds
            sect = i <= N₁ ? f₁.uncoupled[i] : f₂.uncoupled[i - N₁]
            isodd(sect.n) && (θ = -θ)
        end
        θ == 1 || lmul!(θ, t[f₁, f₂])
    end
    return t
end

"""
    compensate_twists!(t, (i₁, j₁), (i₂, j₂), ...)

Apply one or several fermion-pair twists to the tensor `t` in a single
pass over its fusion trees: every block is multiplied by the product of
the factors `(-1)^{p_i p_j}` of all given index pairs, with `p_i, p_j ∈
{0,1}` the parities of the two legs of each pair. Legs are numbered
linearly, codomain legs `1:N₁` first and then domain legs `N₁+1:N₁+N₂`.

This replaces the repetitive hand-written loops

    for (f₁, f₂) in fusiontrees(t)
        coef = (isodd(f₁.uncoupled[i].n) && isodd(f₂.uncoupled[j].n)) ? -1 : 1
        coef != 1 && lmul!(coef, t[f₁, f₂])
    end

(equivalently with both legs on `f₁` or both on `f₂`) that manually
compensate a missing fermionic sign after a bosonic `@tensor` contraction.
"""
function compensate_twists!(t::AbstractParityTensorMap, pairs::Vararg{Tuple{Int,Int}})
    isempty(pairs) && return t
    N₁ = numout(t)
    for (f₁, f₂) in fusiontrees(t)
        coef = 1
        for (i, j) in pairs
            pᵢ = i <= N₁ ? f₁.uncoupled[i].n : f₂.uncoupled[i - N₁].n
            pⱼ = j <= N₁ ? f₁.uncoupled[j].n : f₂.uncoupled[j - N₁].n
            (isodd(pᵢ) && isodd(pⱼ)) && (coef = -coef)
        end
        coef == 1 || lmul!(coef, t[f₁, f₂])
    end
    return t
end


@propagate_inbounds function add_gpermute!(tdst::AbstractParityTensorMap{<:Number, N₁, N₂},
                                          	tsrc::AbstractParityTensorMap,
                                         	p::Index2Tuple{N₁,N₂},
                                         	α::Number,
                                         	β::Number,
                                         	backend::AbstractBackend...) where {N₁,N₂}
    treepermuter(f₁, f₂) = gpermute(f₁, f₂, p[1], p[2])
    return TK.add_transform!(tdst, tsrc, p, treepermuter, α, β, backend...)
end


function gpermute(f1::FusionTree, f2::FusionTree,
                            p1::IndexTuple{N₁}, p2::IndexTuple{N₂}) where {N₁, N₂}
    uncoupled = (f1.uncoupled..., dual.(f2.uncoupled)...)
    uncoupled1′, uncoupled2′ = TupleTools.getindices(uncoupled, p1), TupleTools.getindices(uncoupled, p2)
    uncoupled2′ = ntuple(i->dual(uncoupled2′[i]), Val(N₂))
    coupled1′ = TK.couple(uncoupled1′)
    coupled2′ = TK.couple(uncoupled2′)
    f1′ = FusionTree(uncoupled1′, coupled1′)
    f2′ = FusionTree(uncoupled2′, coupled2′)
    # compute the sign
    if coupled1′ == coupled2′
    	p = TK.linearizepermutation(p1, p2, length(f1), length(f2))
    	swaps = TK.permutation2swaps(p)
    	coeff = 1
        uncoupled = (f1.uncoupled..., dual.(reverse(f2.uncoupled))...)
    	for s in swaps
    		v = uncoupled[s]
    		coeff = (isodd(uncoupled[s].n) && isodd(uncoupled[s+1].n)) ? -coeff : coeff
    		uncoupled = TupleTools.setindex(uncoupled, uncoupled[s+1], s)
    		uncoupled = TupleTools.setindex(uncoupled, v, s+1)
    	end
        tmp = (uncoupled1′..., reverse(uncoupled2′)...)
        # println(uncoupled, " ", tmp)
        (uncoupled == tmp) || error("something wrong")
    else
   		coeff = 0
    end

    return TK.SingletonDict((f1′, f2′)=>coeff)
end
