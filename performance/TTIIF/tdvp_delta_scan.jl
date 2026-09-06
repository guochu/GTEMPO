using GTEMPO, ImpurityModelBase, LinearAlgebra, Printf

f(D, ϵ) = sqrt(1 - (ϵ / D)^2) / π

# ---------------- 虚时 ----------------
N = 25; δτ = 0.01; β = N * δτ; ϵ_d = 1.25π
trunc = truncdimcutoff(D = 100, ϵ = 1.0e-8, add_back = 0)
bath = fermionicbath(spectrum(ϵ -> f(10, ϵ), lb = -10, ub = 10), β = β, μ = 0)
lattice = GrassmannLattice(N = N, δτ = β / N, contour = :imag, ordering = A1Ā1B1B̄1())
corr = correlationfunction(bath, lattice)
b2 = discretebath(bath, δw = 0.2)
exactGτ = toulouse_Gτ(Toulouse(b2, ϵ_d = ϵ_d), collect(0:δτ:β))

println("=" ^ 70)
println("虚时 TDVPIF δ 扫描")
println("=" ^ 70)
for δ in (0.1, 0.05, 0.02)
	t = @elapsed mpsI = hybriddynamics(lattice, corr, TDVPIF(trunc = trunc, δ = δ))
	mpsI = boundarycondition(mpsI, lattice)
	tobs = @elapsed begin
		mpsK = sysdynamics(lattice, AndersonIM(μ = ϵ_d, U = 0), trunc = trunc)
		g = Gτ(lattice, mpsK, mpsI)
	end
	err = norm(g - exactGτ) / norm(exactGτ)
	println(@sprintf("δ=%.2f  IF构建: %6.1f s, 观测量: %5.1f s, 最大键维: %3d, 相对误差: %.3e",
		δ, t, tobs, maximum(bond_dimensions(mpsI)), err))
end

# ---------------- 实时 ----------------
Nt = 10; δt = 0.02; β = 1.0; ϵ_d = -1.0
trunc2 = truncdimcutoff(D = 100, ϵ = 1.0e-9, add_back = 0)
bath2 = fermionicbath(spectrum(ϵ -> f(1, ϵ), lb = -1, ub = 1), β = β, μ = 0)
lattice2 = GrassmannLattice(N = Nt, δt = δt, contour = :real, order = 1, ordering = A1Ā1a1ā1B1B̄1b1b̄1())
corr2 = correlationfunction(bath2, lattice2)
b22 = discretebath(bath2, δw = 0.02)
model = Toulouse(b22, ϵ_d = ϵ_d)
ts = collect(δt:δt:Nt*δt)
exactGt, _ = toulouse_neq_greater_lesser(model, ts, nsys = 0)

println("=" ^ 70)
println("实时 TDVPIF δ 扫描")
println("=" ^ 70)
for δ in (0.1, 0.05, 0.02)
	t = @elapsed mpsI = hybriddynamics(lattice2, corr2, TDVPIF(trunc = trunc2, δ = δ))
	mpsI = boundarycondition(mpsI, lattice2)
	tobs = @elapsed begin
		mpsK = sysdynamics(lattice2, AndersonIM(μ = ϵ_d, U = 0), trunc = trunc2)
		g = [greater(lattice2, i, 1, mpsK, mpsI) for i in 2:Nt]
	end
	ref = im .* exactGt[1:Nt-1]
	err = norm(g - ref) / norm(ref)
	println(@sprintf("δ=%.2f  IF构建: %6.1f s, 观测量: %5.1f s, 最大键维: %3d, 相对误差: %.3e",
		δ, t, tobs, maximum(bond_dimensions(mpsI)), err))
end
