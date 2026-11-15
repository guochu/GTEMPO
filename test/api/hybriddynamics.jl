β0 = 1.0
bath0 = discretefermionicbath([1.0], [0.5]; β=β0)

@testset "API: hybridization kernel" begin
	lat = GrassmannLattice(N=4, δt=0.1, contour=:real)
	kernel = hybridization(lat, bath0; coupling_norm=0.5)
	@test kernel isa HybridizationKernel
	for b1 in (:+, :-), b2 in (:+, :-)
		@test size(kernel.blocks[(b1, b2)]) == (4, 4)
	end
	@test kernel.coupling_norm == 0.5
	@test_throws ArgumentError hybridization(lat, bath0; coupling_norm=-1.0)
	# 形状守卫：块尺寸须与 lattice 一致
	lat2 = GrassmannLattice(N=6, δt=0.1, contour=:real)
	blocks = Dict((b1, b2) => zeros(ComplexF64, 4, 4) for b1 in (:+, :-), b2 in (:+, :-))
	@test_throws DimensionMismatch HybridizationKernel(lat2, blocks)

	# 虚时只有 τ 块
	latτ = GrassmannLattice(N=4, δτ=0.25, contour=:imag)
	kernelτ = hybridization(latτ, bath0)
	@test size(kernelτ.blocks[(:τ, :τ)]) == (4, 4)

	# 混合九块
	latm = GrassmannLattice(Nt=3, δt=0.1, Nτ=4, δτ=0.25, contour=:mixed)
	kernelm = hybridization(latm, bath0)
	for b1 in (:+, :-, :τ), b2 in (:+, :-, :τ)
		@test haskey(kernelm.blocks, (b1, b2))
	end
end

@testset "API: hybriddynamics（IF 链）" begin
	lat = GrassmannLattice(N=4, δt=0.1, bands=2, contour=:real)
	for b in 1:2
		I = hybriddynamics(lat, bath0, PartialIF(); band=b, coupling_norm=1.0)
		@test I isa GM.GMPS && length(I) == length(lat) ÷ 2
		@test check_charge(I)
	end
	# band 越界显式拒绝
	@test_throws BoundsError hybriddynamics(lat, bath0, PartialIF(); band=3)

	# 虚时：PartialIF 与 ExactTTI 同一泛函的两种物化 → 同 Z
	latτ = GrassmannLattice(N=4, δτ=0.25, contour=:imag)
	Iτ = hybriddynamics(latτ, bath0, PartialIF(); band=1)
	@test length(Iτ) == length(latτ) ÷ 2
	Ie = hybriddynamics(latτ, bath0, ExactTTI(); band=1)
	@test length(Ie) == length(latτ) ÷ 2
	@test integrate(Iτ) ≈ integrate(Ie) rtol = 1.0e-8
end
