include("/home/guochu/Documents/Missile/GTEMPO/test/util.jl")
import GTEMPO; using GTEMPO, ImpurityModelBase
using LinearAlgebra
μ = 0.5; β = 0.8; Nτ = 6; δτ = β / Nτ
# 原子极限 ED（α=0）：能级 μ=+0.5，4 带
Z_atomic = (1 + exp(-β * μ))^4
println("原子极限 Z（4 带）= ", Z_atomic)
for sym in (Z2Symmetry(), U1Symmetry(), U1U1Symmetry(), Z2Z2Symmetry())
    lat = GrassmannLattice(N=Nτ, δτ=δτ, bands=4, contour=:imag, symmetry=sym)
    K = sysdynamics(lat, KanamoriIM(U=0.0, J=0.0, norb=2, μ=μ))
    zk = integrate(K)
    println("$sym: Z(K) = ", real(zk))
end
