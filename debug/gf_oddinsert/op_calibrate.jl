# 标定:golden greater/lesser 数组 vs 算符路线关联(全部索引组合)
using Printf
import GTEMPO; using GTEMPO; using ImpurityModelBase
import GrassmannMPS as GM
using LinearAlgebra: Diagonal, tr
μ, ω, α, β = 0.7, 1.0, 0.5, 1.0
bath = discretefermionicbath([1.0], [0.5]; β=β)
golden_gt = [complex(-0.06945314205699218, -0.6608068990426897),
             complex(-0.13549191017047374, -0.6390680792460552),
             complex(-0.1948904437761827, -0.6041604002540943)]
golden_lt = [complex(0.03448940964763923, 0.328146994708573),
             complex(0.06728329137328251, 0.3173518162757938),
             complex(0.09677973022863828, 0.30001717589891896)]
lat = GrassmannLattice(N=3, δt=0.15, bands=1, contour=:real)
ρ0 = Diagonal(ComplexF64[1.0, exp(-β * μ)]); ρ0 /= tr(ρ0)
cn = α
Is = (hybriddynamics(lat, bath, PartialIF(); band=1, coupling_norm=cn),)
K = sysdynamics_weyl2order(lat, ToulouseIM(ϵ_d=μ); coupling_norm=cn, initial_state=ρ0)
cache = environments_cache(K, Is...)
idx(i, b, c) = GrassmannContourIndex(i, b, 1, c)
# greater 候选:-i·⟨a(+,p) ā(+,q)⟩ 对 (p,q) ∈ {(k+1,1),(k,1),(k+1,0→无)…}
for k in 1:3
    ref = golden_gt[k]
    for (p, q) in ((k + 1, 1), (k, 1), (k + 1, 0), (k, 0))
        (1 <= p <= 3 && 1 <= q <= 3) || continue
        v = -1.0im * cached_gf(lat, (idx(p, :+, false), idx(q, :+, true)), K, Is...)
        println("gt[$k] op(p=$p,q=$q): ", round(real(v), digits=6), "+", round(imag(v), digits=6), "im  diff=", round(abs(v - ref), sigdigits=3))
    end
    println("  golden_gt[$k] = ", round(real(ref), digits=6), " + ", round(imag(ref), digits=6), "im")
end
for k in 1:3
    ref = golden_lt[k]
    for (p, q) in ((1, k + 1), (1, k), (0, k))
        (1 <= p <= 3 && 1 <= q <= 3) || continue
        v = 1.0im * cached_gf(lat, (idx(p, :-, true), idx(q, :+, false)), K, Is...)
        @printf("lt[%d] op(p=%d,q=%d): %.6f%+.6fi  diff=%.2e\n", k, p, q, real(v), imag(v), abs(v - ref))
    end
    println("  golden_lt[$k] = ", round(real(ref), digits=6), " + ", round(imag(ref), digits=6), "im")
end
