const AbstractParityTensorMap{T, N₁, N₂} = AbstractTensorMap{T, S, N₁, N₂} where {S <: GradedSpace{ZNIrrep{2}}}


struct GrassmannTensorMap{P<:AbstractParityTensorMap}
	data::P
end

TK.permute(t::GrassmannTensorMap, p1::IndexTuple, p2::IndexTuple; copy::Bool=false) = permute(t, (p1, p2), copy=copy)
TK.permute(t::GrassmannTensorMap, p::Index2Tuple; copy::Bool=false) = GrassmannTensorMap(f_permute(t.data, p, copy=copy))
Base.adjoint(t::GrassmannTensorMap) = GrassmannTensorMap(adjoint(t.data))
get_data(t::GrassmannTensorMap) = t.data

function f_permute(t::AbstractParityTensorMap, (p₁, p₂)::Index2Tuple{N₁,N₂};
                   copy::Bool=false) where {N₁,N₂}
	S = spacetype(t)
    cod = ProductSpace{S,N₁}(map(n -> space(t, n), p₁))
    dom = ProductSpace{S,N₂}(map(n -> dual(space(t, n)), p₂))
    # share data if possible
    if (!copy) && (p₁ === codomainind(t) && p₂ === domainind(t)) 
        return t
    end
    # general case
    @inbounds begin
        return f_permute!(similar(t, cod ← dom), t, (p₁, p₂))
    end
end
function f_permute(t::AdjointTensorMap, (p₁, p₂)::Index2Tuple; copy::Bool=false) 
    p₁′ = TK.adjointtensorindices(t, p₂)
    p₂′ = TK.adjointtensorindices(t, p₁)
    return adjoint(f_permute(adjoint(t), (p₁′, p₂′); copy=copy))
end


@propagate_inbounds function f_permute!(tdst::AbstractParityTensorMap{<:Number, N₁, N₂},
                                        tsrc::AbstractParityTensorMap,
                                        p::Index2Tuple{N₁,N₂}) where {N₁,N₂}
    return add_f_permute!(tdst, tsrc, p, true, false)
end


@propagate_inbounds function add_f_permute!(tdst::AbstractParityTensorMap{<:Number, N₁, N₂},
                                          	tsrc::AbstractParityTensorMap,
                                         	p::Index2Tuple{N₁,N₂},
                                         	α::Number,
                                         	β::Number,
                                         	backend::AbstractBackend...) where {N₁,N₂}
    treepermuter(f₁, f₂) = f_permute(f₁, f₂, p[1], p[2])
    return TK.add_transform!(tdst, tsrc, p, treepermuter, α, β, backend...)
end


function f_permute(f1::FusionTree{ZNIrrep{2}}, f2::FusionTree{ZNIrrep{2}},
                            p1::IndexTuple{N₁}, p2::IndexTuple{N₂}) where {N₁, N₂}
    isdual = (f1.isdual..., (!).(f2.isdual)...)
    uncoupled = (f1.uncoupled..., dual.(f2.uncoupled)...)
    isdual1′, isdual2′ = TupleTools.getindices(isdual, p1), TupleTools.getindices(isdual, p2)
    uncoupled1′, uncoupled2′ = TupleTools.getindices(uncoupled, p1), TupleTools.getindices(uncoupled, p2)
    uncoupled2′ = ntuple(i->dual(uncoupled2′[i]), Val(N₂))
    isdual2′ = (!).(isdual2′)
    coupled1′ = first(⊗(ZNIrrep{2}, uncoupled1′...))
    coupled2′ = first(⊗(ZNIrrep{2}, uncoupled2′...))
    f1′ = FusionTree(uncoupled1′, coupled1′, isdual1′)
    f2′ = FusionTree(uncoupled2′, coupled2′, isdual2′)
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
