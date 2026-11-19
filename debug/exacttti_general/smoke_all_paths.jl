using GTEMPO
using LinearAlgebra: diagind
import GrassmannMPS as GM
using Random
gd(x, y) = maximum(abs.(GM.Dense.todense(x).data - GM.Dense.todense(y).data))

# ---- BCS 虚时 ----
latb = GTEMPO.GrassmannLattice(N=4, δτ=0.25, bands=2, contour=:imag, symmetry=Z2Symmetry())
bathb = bcsbath(fermionicbath(DiracDelta(ω=1.0, α=0.5), β=latb.β), Δ=0.4)
kb = correlationfunction(bathb, latb)
Ib, db = hybriddynamics(latb, kb, ExactTTI())
Fb, _ = hybriddynamics(latb, kb, PartialIF())
println("BCS imag swap dev = ", gd(Ib, Fb))

# ---- 跨带 normal 虚时 ----
function ring_block(n)
    lower = [0.3 * 0.7^(m-1) * cis(0.4m) for m in 1:n-1]
    M = zeros(ComplexF64, n, n)
    for d in 1:n-1, k in 1:n-d
        M[k+d, k] = lower[d]
        M[k, k+d] = -lower[n-d]
    end
    M[diagind(M)] .= 0.05
    return M
end
n = 4
blocks = Dict{Tuple{ContourBlockIndex,ContourBlockIndex},Union{Matrix{ComplexF64},CorrelationMatrix{ComplexF64}}}()
for (ba, bb) in ((1, 1), (1, 2), (2, 1), (2, 2))
    blocks[(ContourBlockIndex(:τ, ba, true), ContourBlockIndex(:τ, bb, false))] = ring_block(n)
end
kc = HybridizationKernel(blocks, 0.1, [1, 2])
latc = GTEMPO.GrassmannLattice(N=n, δτ=0.25, bands=2, contour=:imag)
Ic, dc = hybriddynamics(latc, kc, ExactTTI())
Fc, _ = hybriddynamics(latc, kc, PartialIF())
println("跨带 normal imag swap dev = ", gd(Ic, Fc))

# ---- 同味（T,T）跨带 ----
blocks2 = Dict{Tuple{ContourBlockIndex,ContourBlockIndex},Union{Matrix{ComplexF64},CorrelationMatrix{ComplexF64}}}()
for (ba, bb) in ((1, 2), (2, 1))
    blocks2[(ContourBlockIndex(:τ, ba, true), ContourBlockIndex(:τ, bb, true))] = ring_block(n)
end
ks = HybridizationKernel(blocks2, 0.1, [1, 2])
Is, ds = hybriddynamics(latc, ks, ExactTTI())
Fs, _ = hybriddynamics(latc, ks, PartialIF())
println("同味 (T,T) imag swap dev = ", gd(Is, Fs))

# ---- 带对角 2-band（tuple 契约）----
latd = GTEMPO.GrassmannLattice(N=4, δτ=0.25, bands=2, contour=:imag)
kd = correlationfunction(discretefermionicbath([1.0], [0.5]; β=1.0), latd)
Id, ddg = hybriddynamics(latd, kd, ExactTTI())
Ifd, _ = hybriddynamics(latd, kd, PartialIF())
println("带对角 2-band 返回 = ", typeof(Id), "  integrate dev = ", abs(integrate(Id...) - integrate(Ifd)))

# ---- 跨带 normal 实时（手造 TTI + Keldysh 结构核——用 normal bath 多带）----
latd_r = GTEMPO.GrassmannLattice(N=2, δt=0.15, bands=2, contour=:real)
kr2 = correlationfunction(discretefermionicbath([1.0], [0.5]; β=1.0), latd_r)
I2, d2 = hybriddynamics(latd_r, kr2, ExactTTI())
F2, _ = hybriddynamics(latd_r, kr2, PartialIF())
println("带对角 2-band real: ", typeof(I2), "  integrate dev = ", abs(integrate(I2...) - integrate(F2)))

# ---- 跨带实时（changebands 扩展 + 手工跨带：取单带核重排出双带 + 加跨带块）----
k1 = correlationfunction(discretefermionicbath([1.0], [0.5]; β=1.0), latd_r)
kc_r = changebands(k1, [1 => 1, 1 => 2])
blocks3 = Dict{Tuple{ContourBlockIndex,ContourBlockIndex},Union{Matrix{ComplexF64},CorrelationMatrix{ComplexF64}}}()
merge!(blocks3, Dict(kc_r.blocks))
# 加跨带键：带间 Toeplitz 块（拷带 1 块、缩幅值——保持 Toeplitz/TTI）
rng = Random.MersenneTwister(3)
for b1 in (:+, :-), b2 in (:+, :-)
    base = blocks3[(ContourBlockIndex(b1, 1, true), ContourBlockIndex(b2, 1, false))]
    blocks3[(ContourBlockIndex(b1, 1, true), ContourBlockIndex(b2, 2, false))] = base
    blocks3[(ContourBlockIndex(b1, 2, true), ContourBlockIndex(b2, 1, false))] = base
end
kx = HybridizationKernel(blocks3, 0.2, [1, 2])
I3, d3 = hybriddynamics(latd_r, kx, ExactTTI())
F3, _ = hybriddynamics(latd_r, kx, PartialIF())
println("跨带 normal real swap dev = ", gd(I3, F3))

# ---- 实时单带 + BCS 实时 ----
latr = GTEMPO.GrassmannLattice(N=3, δt=0.15, contour=:real)
kr = correlationfunction(fermionicbath(DiracDelta(ω=1.0, α=0.5), β=1.0), latr)
Ir, dr = hybriddynamics(latr, kr, ExactTTI())
Fr, _ = hybriddynamics(latr, kr, PartialIF())
println("real single-band swap dev = ", gd(Ir, Fr))
latbr = GTEMPO.GrassmannLattice(N=2, δt=0.15, bands=2, contour=:real, symmetry=Z2Symmetry())
kbr = correlationfunction(bathb, latbr)
Ibr, dbr = hybriddynamics(latbr, kbr, ExactTTI())
Fbr, _ = hybriddynamics(latbr, kbr, PartialIF())
println("BCS real swap dev = ", gd(Ibr, Fbr), "  modes = ", [e.mode for e in dbr[:real_expansions]])

# ---- 混合格拒绝 ----
latm = GTEMPO.GrassmannLattice(Nt=2, δt=0.15, Nτ=3, δτ=0.25, contour=:mixed)
spec_c = spectrum(ϵ -> sqrt(1 - (ϵ / 10)^2) / π, lb=-10, ub=10)
km = correlationfunction(fermionicbath(spec_c, β=latm.β), latm)
try
    hybriddynamics(latm, km, ExactTTI())
    println("混合格: 未拒绝（错误！）")
catch e
    println("混合格拒绝: ", e isa ArgumentError)
end
