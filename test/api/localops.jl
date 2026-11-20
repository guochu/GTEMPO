import LinearAlgebra
import LinearAlgebra: I
# localops.jl：符号算符 LocalOp/FockOp 单测（plan/fockmpo_dress_selfread_plan.md
# Phase 1：乘法表 36 项 vs 2×2 矩阵乘逐位 + 结合律 + 物化逐位锚 + 荷派生锚）
@testset "LocalOp/FockOp: 符号算符" begin
    syms = (:I, :a, :adag, :n, :aadg, :empty)
    d2 = [0.0im 1.0; 0.0 0.0]   # 手工 2×2 湮灭算符（|0⟩,|1⟩ 基）
    mats = Dict(:a => d2, :adag => d2', :n => d2'd2, :aadg => d2 * d2',
                :I => LinearAlgebra.Matrix{ComplexF64}(LinearAlgebra.I, 2, 2),
                :empty => zeros(ComplexF64, 2, 2))

    # 乘法表 36 项 vs 2×2 矩阵乘逐位（代数像唯一，only 校验）
    for x in syms, y in syms
        m = mats[x] * mats[y]
        hit = only([s for s in syms if mats[s] == m])
        @test (GTEMPO.LocalOp(x) * GTEMPO.LocalOp(y)).s === hit
    end

    # 结合律（穷举 6³）
    for x in syms, y in syms, z in syms
        @test (GTEMPO.LocalOp(x) * GTEMPO.LocalOp(y)) * GTEMPO.LocalOp(z) ===
              GTEMPO.LocalOp(x) * (GTEMPO.LocalOp(y) * GTEMPO.LocalOp(z))
    end

    # 非对称注记锚（plan §2.1：a·n = a 而 n·a = empty）
    @test (GTEMPO.LocalOp(:a) * GTEMPO.LocalOp(:n)).s === :a
    @test (GTEMPO.LocalOp(:n) * GTEMPO.LocalOp(:a)).s === :empty
    @test (GTEMPO.LocalOp(:adag) * GTEMPO.LocalOp(:adag)).s === :empty

    # 物化：单带因子 vs 手工矩阵逐位 ==（bands{1,2} × 六符号 × 每带位置）
    for bands in (1, 2), s in syms, band in 1:bands
        d = GTEMPO._fock_annihilators(bands)
        manual = s === :a ? d[band] : s === :adag ? d[band]' :
                 s === :n ? d[band]' * d[band] : s === :aadg ? d[band] * d[band]' :
                 s === :I ? Matrix{ComplexF64}(I, 2^bands, 2^bands) :
                 zeros(ComplexF64, 2^bands, 2^bands)
        @test GTEMPO._fockop_matrix(GTEMPO._single_band_op(s, band, bands), bands) == manual
    end

    # fold 序锚：[:a, :adag] ≡ d₁·d₂†（band 1 最左）+ push 序 ≡ 拼接
    # （同 side 同切位列积 foldr 首 push 最左 ≡ 单 FockOp 物化）
    d = GTEMPO._fock_annihilators(2)
    op_fold = GTEMPO.FockOp([GTEMPO.LocalOp(:a), GTEMPO.LocalOp(:adag)])
    @test GTEMPO._fockop_matrix(op_fold, 2) == d[1] * d[2]'
    @test GTEMPO._fockop_matrix(op_fold, 2) ==
          foldr((m, A) -> m * A, [d[1], d[2]'];
                init=Matrix{ComplexF64}(LinearAlgebra.I, 4, 4))

    # JW 实证（d₁·d₂ = −d₂·d₁；band 1 = MSB、JW 穿越更低索引带）
    @test d[1] * d[2] == -(d[2] * d[1])

    # :empty 短路：整积零矩阵 + 荷 0
    op_empty = GTEMPO.FockOp([GTEMPO.LocalOp(:a), GTEMPO.LocalOp(:empty)])
    @test GTEMPO._fockop_matrix(op_empty, 2) == zeros(ComplexF64, 4, 4)
    @test GTEMPO._fockop_charge(op_empty, GTEMPO.sectortype(U1Symmetry()), 2) ==
          one(GTEMPO.sectortype(U1Symmetry()))

    # 荷派生 vs var_charge（spec{:z2,:u1,:u1u1} × bands{1,2} × 每带位置）
    specs = ((Z2Symmetry(), GTEMPO.sectortype(Z2Symmetry())),
             (U1Symmetry(), GTEMPO.sectortype(U1Symmetry())),
             (U1U1Symmetry(), GTEMPO.sectortype(U1U1Symmetry())))
    for bands in (1, 2), (_, G) in specs, band in 1:bands
        @test GTEMPO._fockop_charge(GTEMPO._single_band_op(:a, band, bands), G, bands) ==
              GTEMPO.qneg(GTEMPO.var_charge(G, band, bands))
        @test GTEMPO._fockop_charge(GTEMPO._single_band_op(:adag, band, bands), G, bands) ==
              GTEMPO.var_charge(G, band, bands)
        @test GTEMPO._fockop_charge(GTEMPO._single_band_op(:n, band, bands), G, bands) == one(G)
    end
    # 双带积荷 = 逐带和（:a@band1 + :adag@band2）
    G2 = GTEMPO.sectortype(U1U1Symmetry())
    @test GTEMPO._fockop_charge(op_fold, G2, 2) ==
          GTEMPO.qadd(GTEMPO.qneg(GTEMPO.var_charge(G2, 1, 2)), GTEMPO.var_charge(G2, 2, 2))

    # 构造守卫：未知符号 / band 越界 / FockOp 带数失配（物化点校验）
    @test_throws ArgumentError GTEMPO.LocalOp(:x)
    @test_throws BoundsError GTEMPO._single_band_op(:n, 3, 2)
    @test_throws ArgumentError GTEMPO._fockop_matrix(GTEMPO.FockOp([GTEMPO.LocalOp(:a)]), 2)

    # 结构化 ==/hash（插入缓存键先决；Vector 逐元素）
    @test GTEMPO.FockOp([GTEMPO.LocalOp(:a), GTEMPO.LocalOp(:I)]) ==
          GTEMPO.FockOp([GTEMPO.LocalOp(:a), GTEMPO.LocalOp(:I)])
    @test hash(GTEMPO.FockOp([GTEMPO.LocalOp(:a), GTEMPO.LocalOp(:I)])) ==
          hash(GTEMPO.FockOp([GTEMPO.LocalOp(:a), GTEMPO.LocalOp(:I)]))
    @test GTEMPO.FockOp([GTEMPO.LocalOp(:a)]) != GTEMPO.FockOp([GTEMPO.LocalOp(:adag)])
    @test GTEMPO.FockOp([GTEMPO.LocalOp(:n)]) != GTEMPO.FockOp([GTEMPO.LocalOp(:n), GTEMPO.LocalOp(:I)])
end
