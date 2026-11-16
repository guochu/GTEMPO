using LinearAlgebra: I
using TensorKit: scalartype

@testset "API: GrassmannLattice 几何与对称性" begin
	# 虚时
	lat = GrassmannLattice(N=6, δτ=0.1, contour=:imag)
	@test lat isa ImagGrassmannLattice
	@test lat.k == lat.kτ == 6 && lat.N == lat.Nτ == 6
	@test lat.β ≈ 0.6 && lat.bands == 1
	@test length(lat) == 2 * lat.bands * lat.k
	@test branches(lat) == (:τ,)
	@test symmetry(lat) === :z2
	@test scalartype(lat) == ComplexF64
	lat2 = similar(lat, bands=2)
	@test lat2.bands == 2 && lat2.δτ == lat.δτ && lat2.N == lat.N
	lat3 = similar(lat, N=8, δτ=0.05)
	@test lat3.N == 8 && lat3.δτ == 0.05

	# 实时
	lat = GrassmannLattice(N=5, δt=0.1, contour=:real)
	@test lat isa RealGrassmannLattice
	@test lat.k == lat.kt == lat.Nt == 5
	@test lat.t == 0.5 && scalartype(lat) == ComplexF64
	@test branches(lat) == (:+, :-)

	# 混合（Kadanoff）
	lat = GrassmannLattice(Nt=3, δt=0.1, Nτ=4, δτ=0.25, contour=:mixed)
	@test lat isa MixedGrassmannLattice
	@test lat.kt == lat.Nt == 3 && lat.kτ == lat.Nτ == 4 && lat.β ≈ 1.0

	# 对称性手选（D1）：四形 + 荷分辨 bands 守卫
	for (sym, G) in ((:z2, GTEMPO.Z2Symmetry()), (:u1, GTEMPO.U1Symmetry()),
	                 (:u1u1, GTEMPO.U1U1Symmetry()), (:z2z2, GTEMPO.Z2Z2Symmetry()))
		l = GrassmannLattice(N=4, δt=0.1, bands=2, contour=:real, symmetry=sym)
		@test symmetry(l) === sym
		@test GTEMPO.symmetry_spec(sym) isa typeof(G)
	end
	@test_throws ArgumentError GTEMPO.symmetry_spec(:u1x)
	# 荷分辨形要求偶数带（轨道交错：奇带 up / 偶带 down；多轨道 bands = 2·轨道数）
	@test_throws ArgumentError GrassmannLattice(N=4, δt=0.1, bands=3, contour=:real, symmetry=:u1u1)
	@test_throws ArgumentError GrassmannLattice(N=4, δt=0.1, bands=3, contour=:real, symmetry=:z2z2)
	# 多轨道（bands=4 = 2 轨道）合法构造（D1 修订）
	for sym in (:u1u1, :z2z2)
		l4 = GrassmannLattice(N=4, δt=0.1, bands=4, contour=:real, symmetry=sym)
		@test l4.bands == 4 && length(l4) == 4 * 4 * l4.k   # 实时 = 2k 块 × 2M 位
	end
	# BCS 守卫（D7）：仅 :z2
	for sym in (:u1, :u1u1, :z2z2)
		@test_throws ArgumentError GTEMPO.bcs_symmetry_guard(sym)
	end
	@test GTEMPO.bcs_symmetry_guard(:z2)

	# 站空间序列与变量荷
	l = GrassmannLattice(N=4, δt=0.1, bands=2, contour=:real, symmetry=:u1u1)
	ss = site_spaces(GTEMPO.symmetry_spec(:u1u1), 4, 2)
	@test length(ss) == 4 && ss[1] != ss[2]   # up/down 交替站形
	cc = GTEMPO.cell_charges(GTEMPO.symmetry_spec(:u1u1), 2)
	@test cc == [(1, 0), (-1, 0), (0, 1), (0, -1)]
	# 时序格序 / 存储格序 / 切位索引
	l = GrassmannLattice(N=4, δt=0.1, contour=:real)
	cells = GTEMPO.chronological_cells(l)
	@test length(cells) == 2 * l.Nt
	@test GTEMPO.insertion_cut_index(l, OperatorInsertion(:+, 2, Matrix{ComplexF64}(I, 2, 2))) == 2
	@test GTEMPO.insertion_cut_index(l, OperatorInsertion(:-, 3, Matrix{ComplexF64}(I, 2, 2))) == 2 * l.Nt - 3
end
