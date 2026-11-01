# BCS 快速诊断：partialif (GTEMPO) vs QuAPI.Dense（τ / 实时 / 混合逐轮廓）
using Random
using GTEMPO
using Z2Tensors
using ImpurityModelBase
using QuAPI
using QuAPI.Dense
using Z2Tensors: fusiontrees, dim, sectors
using GTEMPO: space_l, space_r

Random.seed!(12354)

# ---- 复制自 test/api/quapi_crosscheck.jl 的辅助函数 ----
function _sector_offset(V, c)
    o = 0
    for d in sectors(V)
        d.n < c.n || break
        o += dim(V, d)
    end
    return o
end
function _site_matrix(A, σ::Int)
    Vl, Vr = space_l(A), space_r(A)
    M = zeros(GTEMPO.scalartype(A), dim(Vl), dim(Vr))
    for (fl, fr) in fusiontrees(A)
        sl, sp = fl.uncoupled
        sp.n == σ || continue
        fl.coupled.n == fr.coupled.n || continue
        sr = fr.uncoupled[1]
        ol, o = _sector_offset(Vl, sl), _sector_offset(Vr, sr)
        dl, dr = dim(Vl, sl), dim(Vr, sr)
        M[ol+1:ol+dl, o+1:o+dr] .= reshape(collect(A[fl, fr]), dl, dr)
    end
    return M
end
function _gmps_to_dense(psi)
    L = length(psi)
    Ms = [[_site_matrix(psi[p], σ) for σ in 0:1] for p in 1:L]
    out = zeros(GTEMPO.scalartype(psi), ntuple(_ -> 2, L))
    for I in CartesianIndices(out)
        v = Ms[1][I[1]]
        for p in 2:L
            v = v * Ms[p][I[p]]
        end
        out[I] = v[1, 1]
    end
    return out .* (GTEMPO.scaling(psi)^L)
end
function _gt_vs_quapi(dense, lat, corr)
    qif = GrassmannIF(GrassmannHybridizationMatrix(corr))
    Lq = length(qif.basis)
    sites = [lat[GrassmannContourIndex((leg.branch === :τ) ? leg.j + 1 : leg.j,
                                       leg.branch, leg.band, leg.conj)]
             for leg in qif.basis]
    out = zeros(eltype(dense), 1 << Lq)
    for m in 0:(1 << Lq)-1
        σ = falses(length(dense))
        occ = Int[]
        for q in 1:Lq
            iszero((m >> (q - 1)) & 1) && continue
            σ[sites[q]] = true
            push!(occ, sites[q])
        end
        sgn = 1
        for a in 1:length(occ)-1, b in a+1:length(occ)
            occ[a] > occ[b] && (sgn = -sgn)
        end
        out[m+1] = sgn * dense[ntuple(p -> σ[p] ? 2 : 1, ndims(dense))...]
    end
    return maximum(abs.(out - qif.data)) / norm(qif.data)
end

# ---- BCS 三轮廓（小规模：dense IF 2^腿数，N=1 时 τ/real 仅 2^8）----
μ = 0.7
bath = fermionicbath(DiracDelta(ω=1.0, α=0.5), β=0.6, μ=μ)
bcs = bcsbath(bath, Δ=0.6)

for (label, lat) in (("τ",   GrassmannLattice(N=1, δτ=0.3, contour=:imag, bands=2)),
                     ("real", GrassmannLattice(N=1, δt=0.1, contour=:real, bands=2)),
                     ("mixed", GrassmannLattice(Nt=1, δt=0.1, Nτ=1, δτ=0.3, contour=:mixed, bands=2)))
    corr = correlationfunction(bcs, lat)
    IF = hybriddynamics(lat, corr, orbital=1, trunc=NoTruncation())
    e = _gt_vs_quapi(_gmps_to_dense(IF), lat, corr)
    println("BCS ", label, ": relerr=", e)
end
