# Cross-check of QuAPI.Dense.GrassmannADT against GTEMPO's untruncated
# Grassmann-MPS influence functionals: normal fermionic bath on the imaginary,
# real and mixed contours, plus the BCS bath (two bands). electron-phonon
# (Fock lattice) is intentionally not covered.

using QuAPI.Dense
using Z2Tensors: fusiontrees, dim, sectors
using GTEMPO: space_l, space_r

function _sector_offset(V, c)
    o = 0
    for d in sectors(V)
        d.n < c.n || break
        o += dim(V, d)
    end
    return o
end

# dense bond matrix of a site tensor for the Grassmann occupation σ ∈ {0,1}
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

# dense coefficient tensor of a GrassmannMPS (axis p = lattice site p,
# index 1 = empty, index 2 = occupied); includes scaling^L
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

# map GTEMPO's dense coefficient tensor onto QuAPI's leg basis and return the
# relative maximum deviation from GrassmannIF(corr)
function _gt_vs_quapi(dense, lat, corr)
    qif = GrassmannIF(corr)
    Lq = length(qif.basis)
    out = zeros(eltype(dense), 1 << Lq)
    for m in 0:(1 << Lq)-1
        σ = falses(length(dense))
        for q in 1:Lq
            iszero((m >> (q - 1)) & 1) && continue
            leg = qif.basis[q]
            i = (leg.branch === :τ) ? leg.j + 1 : leg.j
            σ[index(lat, i, conj=leg.conj, branch=leg.branch, band=leg.band)] = true
        end
        out[m+1] = dense[ntuple(p -> σ[p] ? 2 : 1, ndims(dense))...]
    end
    return maximum(abs.(out - qif.data)) / norm(qif.data)
end

@testset "QuAPI dense IF cross-check" begin
    μ, ω, α = 0.7, 1.0, 0.5
    spec = DiracDelta(ω=ω, α=α)
    β = 0.6
    bath = fermionicbath(spec, β=β, μ=μ)
    tol = 1.0e-10

    # ---------------- normal fermionic bath ----------------
    lat = GrassmannLattice(N=2, δτ=0.3, contour=:imag)
    corr = correlationfunction(bath, lat)
    IF = hybriddynamics(lat, corr, trunc=NoTruncation())
    @test _gt_vs_quapi(_gmps_to_dense(IF), lat, corr) < tol

    latr = GrassmannLattice(N=2, δt=0.1, contour=:real)
    corrr = correlationfunction(bath, latr)
    IFr = hybriddynamics(latr, corrr, trunc=NoTruncation())
    @test _gt_vs_quapi(_gmps_to_dense(IFr), latr, corrr) < tol

    latm = GrassmannLattice(Nt=1, δt=0.1, Nτ=2, δτ=0.3, contour=:mixed)
    corrm = correlationfunction(bath, latm)
    IFm = hybriddynamics(latm, corrm, trunc=NoTruncation())
    @test _gt_vs_quapi(_gmps_to_dense(IFm), latm, corrm) < tol

    # ---------------- BCS bath (two bands) ----------------
    latb = GrassmannLattice(N=2, δτ=0.3, contour=:imag, bands=2)
    corrb = correlationfunction(bcsbath(bath, Δ=0.6), latb)
    IFb = hybriddynamics(latb, corrb, orbital=1, trunc=NoTruncation())
    @test _gt_vs_quapi(_gmps_to_dense(IFb), latb, corrb) < tol
end
