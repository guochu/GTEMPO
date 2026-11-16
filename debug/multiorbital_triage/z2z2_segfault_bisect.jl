include("/home/guochu/Documents/Missile/GTEMPO/test/util.jl")
import GTEMPO; using GTEMPO, ImpurityModelBase
import GrassmannMPS as GM
using LinearAlgebra

N = parse(Int, get(ENV, "NN", "6"))
μ = 0.5
println("== stage 1: lattice =="); flush(stdout)
lat = GrassmannLattice(N=N, δτ=0.8/N, bands=4, contour=:imag, symmetry=Z2Z2Symmetry())
println("   ok; k=", lat.k, " len=", length(lat)); flush(stdout)
println("== stage 2: vacuumstate =="); flush(stdout)
v = vacuumstate(lat)
println("   ok"); flush(stdout)
println("== stage 3: sysdynamics (K build) =="); flush(stdout)
K = sysdynamics(lat, KanamoriIM(U=0.0, J=0.0, norb=2, μ=μ))
println("   ok; len=", length(K), " bonddims=", GM.bonddims(K)); flush(stdout)
println("== stage 4: integrate =="); flush(stdout)
z = integrate(K)
println("   ok; Z = ", real(z)); flush(stdout)
