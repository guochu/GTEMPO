# Phase 1 烟测（plan §Phase 1 验收）：
#   ① 包可加载（runtests 顶层面已完成）
#   ② 真空态/randomgmps 烟测（四形各一）
#   ③ 对称性校验拒绝路径（未知 Symbol、荷分辨形 bands>2、BCS 守卫）
#   ④ lattice _index 单元测试（位置映射健全性：双射 + 结构公式）
#   ⑤ vacuumstate/coefficient_scaled 往返
using Test, Random
using GTEMPO
using GTEMPO: branches, sectortype, cell_charges, monomial, coefficient
import TensorKit as TK
import GrassmannMPS as GM

Random.seed!(12354)

@testset "Phase 1 smoke: lattices × 4 symmetries" begin
    # ②③ 真空/随机 + 拒绝路径
    # 单荷形任意带数（bands=1）；荷分辨形偶数带（轨道交错）
    for (sym, bands) in ((Z2(), 1), (U1(), 1), (Z2Z2(), 2), (U1U1(), 2))
        lat = GrassmannLattice(N=3, δτ=0.1, bands=bands, contour=:imag, symmetry=sym)
        v = vacuumstate(lat)
        @test length(v) == length(lat) ÷ 2
        @test GM.grassmann_basis(GM.physpace(v, 1)) !== nothing
        r = randomgmps(Float64, site_spaces(lat), fill(2, length(v) - 1))
        @test length(r) == length(v)
    end
    # :z2/:u1 bands=3 合法（荷合并）
    for sym in (Z2(), U1())
        lat = GrassmannLattice(N=2, δτ=0.1, bands=3, contour=:imag, symmetry=sym)
        @test length(vacuumstate(lat)) == 3 * lat.kτ   # 站数 = 变量数/2
    end
    # ③ 拒绝路径（Symbol 手选已退役：非法输入在 kwarg 类型层即拒）
    @test_throws TypeError GrassmannLattice(N=2, δτ=0.1, contour=:imag, symmetry=:su2)
    @test_throws ArgumentError GrassmannLattice(N=2, δτ=0.1, bands=3, contour=:imag, symmetry=Z2Z2())
    @test_throws ArgumentError GrassmannLattice(N=2, δτ=0.1, bands=3, contour=:imag, symmetry=U1U1())
end

@testset "Phase 1 smoke: lattice _index 位置映射" begin
    # imag Chrono：格 i ⟷ 块 k−i+1；块内 2M 位 (a,ā) 相邻
    lat = GrassmannLattice(N=3, δτ=0.1, bands=2, contour=:imag)
    @test length(lat) == 2 * 2 * 3
    m = indexmappings(lat)
    @test length(m) == length(lat)
    @test length(unique(values(m))) == length(m)               # 双射
    # 格 1 = 末块（k−1+1 = k）
    @test m[(1, false, :τ, 1)] == 2 * 2 * 3 - 2 * 2 + 1
    @test m[(1, true, :τ, 2)] == 2 * 2 * 3 - 2 * 2 + 4
    # 实时 Chrono：块 (+,i) → i；(−,i) → 2k−i+1
    latr = GrassmannLattice(N=2, δt=0.1, bands=1, contour=:real)
    mr = indexmappings(latr)
    @test mr[(1, false, :+, 1)] == 1
    @test mr[(1, true, :+, 1)] == 2
    @test mr[(1, false, :-, 1)] == 3 * 2 + 1   # block 2k−i+1 = 4
    @test length(unique(values(mr))) == length(mr) == length(latr)
    # 混合 Chrono：τ j → 块 j；+ j → Nτ+j；− j → Nτ+Nt+(Nt−j+1)
    latm = GrassmannLattice(Nτ=2, δτ=0.1, Nt=1, δt=0.1, bands=1, contour=:mixed)
    mm = indexmappings(latm)
    @test mm[(1, false, :τ, 1)] == 1
    @test mm[(2, false, :τ, 1)] == 3   # 块 2 起点 = (2−1)*2M+1
    @test mm[(1, false, :+, 1)] == 2 * 2 + 1
    @test mm[(1, false, :-, 1)] == 2 * 3 + 1
    @test length(unique(values(mm))) == length(mm) == length(latm)
end

@testset "Phase 1 smoke: vacuumstate/coefficient_scaled 往返" begin
    for sym in (Z2(), U1())
        lat = GrassmannLattice(N=3, δτ=0.1, bands=1, contour=:imag, symmetry=sym)
        v = vacuumstate(lat)
        # 全空占据的系数 = 1（站级 digit 全 0 = vac；恒等真空链）
        z = coefficient_scaled(v, zeros(Int, length(v)))
        @test abs(value(z) - 1.0) < 1e-14
    end
    # U1U1 交替站形真空往返
    lat = GrassmannLattice(N=2, δτ=0.1, bands=2, contour=:imag, symmetry=U1U1())
    v = vacuumstate(lat)
    z = coefficient_scaled(v, zeros(Int, length(v)))
    @test abs(value(z) - 1.0) < 1e-14
    # 单项式系数：monomial(带荷 positions) 的对应组态 = coeff
    spec = GTEMPO.sectortype(Z2())
    charges = cell_charges(spec, 2)  # bands=2
    mono = monomial(spec, repeat(charges, 3), [1, 6]; coeff=1.0)
    occ = zeros(Int, length(mono))
    occ[1] = 2   # 站 1 = η（digit 2：vac0/pair1/eta2/etabar3）
    occ[3] = 3   # 站 3 = η̄（digit 3；变量 6 = 全链第 6 变量 = 站 3 的第二变量）
    c = coefficient(mono, occ)
    @test abs(c) ≈ 1.0
end
