import TensorKit as TK
import TensorKit: ⊗, ←
import Random
import LinearAlgebra
# fockmpo.jl：FockMPO 直积 MPO + fockmpo2gmps 闸门（plan V1 Phase 1；
# 2026-11-20 起穿戴内联 fock_sysdynamics + converter 荷账自读——
# plan/fockmpo_dress_selfread_plan.md；绝对数值锚由 observables golden
# （U1GTEMPO2 机器精度）与 ED 电池承接）
#
# 对拍 = 公开 sysdynamics（穿戴内联 + 自读）vs 手工闭环（fock_sysdynamics 裸
# 调用（真空穿戴）→ 手工稀疏乘法（公开 sparsefockmpo *）→ 手工 ρ₀ 三明治折
# 首站 → 手工 oddlater 号 → fockmpo2gmps）todense 机器精度——穿戴与自读两条
# 独立簿记路径互验。配置矩阵：
#   {imag, real, mixed} × bands{1,2} × spec{:z2,:u1,:u1u1}（:z2z2 抽查）
#   × 插入{无、偶（n̂ cut0/内/末）、奇（greater/lesser 式双插）、多插}
#   × order{1,2}
# 端到端：integrate(K_new) vs integrate(K_base)（无插入 + 单偶插 + 奇双插）。
# 荷账自读锚（Phase 3 先立，改造前后都绿）：未 dressing 链 δ ≡ 0；单带荷插入
# δ = 符号荷、他站零；prefix 恒等式（多插偶+奇：键 k = Σ_{j≤k} δ_j）。
@testset "FockMPO: fockmpo2gmps 对拍矩阵" begin

    function _wf_build(lattice, model, order, coupling_norm, insertions, initial_state)
        # 手工闭环：**原始装配复刻**（fock_sysdynamics 步骤 ①——裸调用已按 P3
        # 恒折边界态，不能作手工基座）→ 手工稀疏乘法（公开 sparsefockmpo *）→
        # 手工 ρ₀ 三明治折首站（**在稀疏乘法之后**——红线 1）→ 手工 oddlater 号
        # （cut ≥ k）→ fockmpo2gmps（荷账自读，无插入元数据）——与公开入口
        # （穿戴内联）两条独立簿记路径互验
        bands = lattice.bands
        dim = 2^bands
        G = TK.sectortype(lattice)
        z = one(G)
        weights = GTEMPO.local_weights(G, bands, Int(order))
        Kf = GTEMPO.FockMPO([begin
            S = GTEMPO._fock_propagator_at(lattice, model, b, j, Int(order))
            mats = order == 2 ? [S * W * S for W in weights] :
                   b === :- ? [W * S for W in weights] : [S * W for W in weights]
            t = TK.zeros(ComplexF64, GTEMPO._fock_bondspace(G, dim) ⊗ GTEMPO._fock_bondspace(G, dim^2) ←
                                 GTEMPO._fock_bondspace(G, dim^2) ⊗ GTEMPO._fock_bondspace(G, dim))
            blk = t[(z, z, TK.dual(z), TK.dual(z))]
            for m in 0:dim^2-1
                blk[:, m+1, m+1, :] .= mats[m+1]
            end
            t
        end for (b, j) in GTEMPO.chronological_cells(lattice)])
        cut0 = filter(i -> GTEMPO.insertion_cut_index(lattice, i) == 0, insertions)
        ge1 = filter(i -> GTEMPO.insertion_cut_index(lattice, i) >= 1, insertions)
        Wl = sparsefockmpo(lattice, filter(i -> i.side === :left, ge1))
        Wr = sparsefockmpo(lattice, filter(i -> i.side === :right, ge1))
        Kf = Wr * (Wl * Kf)
        # ρ₀ 三明治（initial_state 仅实时覆盖——公开入口同门控；push 序 = 0⁺ 序；
        # 边界态构造内联复刻——fock_vacuumstate 已删，公开路径同式在
        # fock_sysdynamics 内联）
        rho0 = lattice isa GTEMPO.AbstractRealGrassmannLattice ? initial_state : nothing
        rho = if lattice isa GTEMPO.AbstractRealGrassmannLattice
            rho0 === nothing ? Matrix{ComplexF64}(LinearAlgebra.Diagonal(vcat([1.0], zeros(dim - 1)))) :
            (rho0 isa GTEMPO.FockMatrix ? Matrix{ComplexF64}(rho0.data) : Matrix{ComplexF64}(rho0))
        else
            Matrix{ComplexF64}(LinearAlgebra.I, dim, dim)
        end
        rho = foldr((i, A) -> GTEMPO._fockop_matrix(i.op, bands) * A,
                    filter(i -> i.side === :left, cut0); init=rho)
        rho = foldl((A, i) -> A * GTEMPO._fockop_matrix(i.op, bands),
                    filter(i -> i.side === :right, cut0); init=rho)
        blk = Kf.data[1][(z, z, TK.dual(z), TK.dual(z))]
        for m in 0:dim^2-1
            blk[:, m+1, m+1, :] .= blk[:, m+1, m+1, :] * rho
        end
        # oddlater（谓词 cut ≥ k 逐字；奇偶源 = 符号荷）
        bycut = Dict{Int,Vector{eltype(insertions)}}()
        for ins in insertions
            push!(get!(bycut, GTEMPO.insertion_cut_index(lattice, ins), eltype(insertions)[]), ins)
        end
        for (k, t) in enumerate(Kf.data)
            oddlater = 0
            for (cut, items) in bycut, ins in items
                cut >= k || continue
                oddlater += GM.fermionparity(GTEMPO._fockop_charge(ins.op, G, bands)) ? 1 : 0
            end
            oddlater == 0 && continue
            blk = t[(z, z, TK.dual(z), TK.dual(z))]
            for m in 0:dim^2-1
                blk[:, m+1, m+1, :] .*= (-1)^(oddlater * count_ones(m))
            end
        end
        out = fockmpo2gmps(Kf, lattice)
        # Weyl 标量由调用方后置累加补上（⚠️ logscale += lscale——压缩归一化因子
        # 已在 logscale 内，覆盖会抹掉）
        lscale = order == 2 ? real(coupling_norm / 4 *
                                   sum(GTEMPO._cell_step2(lattice, b) for (b, _) in GTEMPO.chronological_cells(lattice))) : 0.0
        GM.setlogscale!(out, GM.logscale(out) + lscale)
        return out
    end

    function _wf_check(lattice, model, order, coupling_norm, insertions; initial_state=nothing)
        # 基线 = 公开入口（穿戴内联实现体）——对拍手工闭环 workflow（_wf_build）：
        # 穿戴装配合同（切位过滤/三明治序/oddlater 谓词/轮廓门控）与荷账自读一致
        K_base = sysdynamics(lattice, model; order=order, coupling_norm=coupling_norm,
                             insertions=insertions, initial_state=initial_state)
        K_new = _wf_build(lattice, model, order, coupling_norm, insertions, initial_state)
        nst = length(K_base)
        if 4^nst <= 2^20
            tb = GM.Dense.todense(K_base).data
            tn = GM.Dense.todense(K_new).data
            @test length(tn) == length(tb)
            Δ = maximum(abs.(tn .- tb))
            @test Δ <= 1e-10
        else
            # 长链：coefficient_scaled 抽样（全 vac、全 pair、随机 3 组）
            occs = [fill(0, nst), fill(1, nst)]
            Random.shuffle!(vcat(occs, [rand(0:3, nst) for _ in 1:3]))
            local Δ = 0.0
            for occ in occs
                Δ = max(Δ, abs(value(GM.coefficient_scaled(K_base, occ)) - value(GM.coefficient_scaled(K_new, occ))))
            end
            @test Δ <= 1e-10
        end
        return Δ
    end

    "链 δ 读数（荷账自读锚用）"
    function _deltas(lattice, model, order, insertions=OperatorInsertion[], initial_state=nothing)
        G = TK.sectortype(lattice)
        bands = lattice.bands
        z = one(G)
        fqs1 = GTEMPO.fock_charges(G, bands)
        qms = GTEMPO._slot_charges(G, bands)
        Kf = fock_sysdynamics(lattice, model; order=order, insertions=insertions,
                              initial_state=initial_state)
        return [GTEMPO._site_charge_shift(Kf.data[k][(z, z, TK.dual(z), TK.dual(z))],
                                          fqs1, qms, bands) for k in 1:length(Kf.data)], z
    end

    model1 = ImpurityHamiltonian([tunneling(1, 1, coeff=0.5)]; bands=1)
    model2 = KanamoriIM(U=1.0, J=0.2, norb=1)

    # ---- 荷账自读锚（Phase 3 先立）：未 dressing 链 δ ≡ 0 ----
    # （全轮廓 × bands × spec × order——槽荷表符号约定与 local_weights 逐位核对
    # 的机器可检不变量）
    @testset "荷账自读锚：未 dressing 链 δ ≡ 0" begin
        for order in (1, 2)
            for (bands, model, syms) in ((1, model1, (Z2Symmetry(), U1Symmetry())),
                                         (2, model2, (Z2Symmetry(), U1Symmetry(), U1U1Symmetry())))
                for sym in syms
                    lat_i = GrassmannLattice(N=2, δτ=0.25, bands=bands, contour=:imag, symmetry=sym)
                    lat_r = GrassmannLattice(N=2, δt=0.1, bands=bands, contour=:real, symmetry=sym)
                    lat_m = GrassmannLattice(Nt=2, Nτ=2, δt=0.1, δτ=0.25, bands=bands, contour=:mixed, symmetry=sym)
                    for (tag, lat) in (("imag", lat_i), ("real", lat_r), ("mixed", lat_m))
                        δs, z = _deltas(lat, model, order)
                        @test all(==(z), δs)
                    end
                end
            end
        end
    end

    @testset "荷账自读锚：单带荷插入 δ = 符号荷 + prefix 恒等式" begin
        # U1 bands=1 虚时 N=2：d@cut1（−q@站1）、ddag@cut2（+q@站2）他站零
        lat = GrassmannLattice(N=2, δτ=0.25, bands=1, contour=:imag, symmetry=U1Symmetry())
        G = TK.sectortype(lat)
        q = GTEMPO.var_charge(G, 1, 1)
        z = one(G)
        δs, _ = _deltas(lat, model1, 1, [OperatorInsertion(:τ, 1, GTEMPO._single_band_op(:a, 1, 1))])
        @test δs[1] == GTEMPO.qneg(q) && δs[2] == z
        δs, _ = _deltas(lat, model1, 1, [OperatorInsertion(:τ, 2, GTEMPO._single_band_op(:adag, 1, 1))])
        @test δs[1] == z && δs[2] == q
        # prefix 恒等式（plan §2.2 验证例：多插偶+奇 n̂@cut0 + d@cut1(−q) +
        # ddag@cut2(+q)）——键 1 = δ₁ = −q、键 2 = δ₁+δ₂ = 0（≡ 旧 cut 累计）
        δs, _ = _deltas(lat, model1, 1,
                        [OperatorInsertion(:τ, 0, GTEMPO._single_band_op(:n, 1, 1)),
                         OperatorInsertion(:τ, 1, GTEMPO._single_band_op(:a, 1, 1)),
                         OperatorInsertion(:τ, 2, GTEMPO._single_band_op(:adag, 1, 1))])
        @test δs[1] == GTEMPO.qneg(q)
        @test GTEMPO.qadd(δs[1], δs[2]) == z   # prefix(键2) = δ₁+δ₂ = −q+q = 0（≡ 旧 cut 累计）
        # 奇单插链 Σδ ≠ 0：fock_sysdynamics 可建（P5 逐站自检过），converter Σδ
        # 守卫拒绝
        Kf = fock_sysdynamics(lat, model1; insertions=[OperatorInsertion(:τ, 1, GTEMPO._single_band_op(:a, 1, 1))])
        @test_throws ArgumentError fockmpo2gmps(Kf, lat)
    end

    @testset "荷账自读锚：电荷非对角初态拒绝（A1′）" begin
        lat = GrassmannLattice(N=2, δt=0.1, bands=2, contour=:real, symmetry=U1Symmetry())
        rho_bad = Matrix{ComplexF64}(LinearAlgebra.I, 4, 4); rho_bad[1, 2] = 0.3
        @test_throws ArgumentError fock_sysdynamics(lat, model2; initial_state=rho_bad)
        @test_throws ArgumentError sysdynamics(lat, model2; initial_state=rho_bad)
        # 对角初态（含真空投影缺省）不受影响
        δs, z = _deltas(lat, model2, 1)
        @test all(==(z), δs)
    end

    # ---- 对拍矩阵 ----
    for order in (1, 2)
        cn = order == 2 ? 0.7 : 0.0
        # ---- imag：bands{1,2} × spec{:z2,:u1,:u1u1}（:z2z2 抽查）----
        for (bands, model, syms) in ((1, model1, (Z2Symmetry(), U1Symmetry())),
                                     (2, model2, (Z2Symmetry(), U1Symmetry(), U1U1Symmetry(), Z2Z2Symmetry())))
            for sym in syms
                lat = GrassmannLattice(N=2, δτ=0.25, bands=bands, contour=:imag, symmetry=sym)
                @testset "imag bands=$bands $sym order=$order 无插入" begin
                    _wf_check(lat, model, order, cn, OperatorInsertion[])
                end
                @testset "imag bands=$bands $sym order=$order 偶插入三切位" begin
                    for band in 1:bands, jj in (0, 1, 2)
                        _wf_check(lat, model, order, cn,
                                  [OperatorInsertion(:τ, jj, GTEMPO._single_band_op(:n, band, bands))])
                    end
                end
            end
        end
        # ---- imag：奇插入（符号双插）+ 多插 ----
        for (bands, model, sym) in ((1, model1, U1Symmetry()), (2, model2, U1U1Symmetry()))
            lat = GrassmannLattice(N=2, δτ=0.25, bands=bands, contour=:imag, symmetry=sym)
            a = GTEMPO._single_band_op(:a, 1, bands)
            adag = GTEMPO._single_band_op(:adag, 1, bands)
            @testset "imag bands=$bands $sym order=$order 奇双插(:τ)" begin
                ins = [OperatorInsertion(:τ, 2, a), OperatorInsertion(:τ, 1, adag)]
                _wf_check(lat, model, order, cn, ins)
            end
            @testset "imag bands=$bands $sym order=$order 多插(偶+奇)" begin
                ins = [OperatorInsertion(:τ, 0, GTEMPO._single_band_op(:n, 1, bands)),
                       OperatorInsertion(:τ, 1, a),
                       OperatorInsertion(:τ, 2, adag)]
                _wf_check(lat, model, order, cn, ins)
            end
        end
    end

    # ---- real：真空 ρ₀ + 偶/奇插入（bands=2 覆盖 U1U1 与单 U1——单 U1×2 带
    # 曾漏网：:- 支 j = 0 = 链末切位，lesser(j=0) 路径）----
    for (bands, model, sym) in ((1, model1, U1Symmetry()), (2, model2, U1U1Symmetry()),
                                (2, model2, U1Symmetry()))
        lat = GrassmannLattice(N=2, δt=0.1, bands=bands, contour=:real, symmetry=sym)
        a = GTEMPO._single_band_op(:a, 1, bands)
        adag = GTEMPO._single_band_op(:adag, 1, bands)
        for order in (1, 2)
            cn = order == 2 ? 0.7 : 0.0
            @testset "real bands=$bands $sym order=$order 无插入" begin
                _wf_check(lat, model, order, cn, OperatorInsertion[])
            end
            @testset "real bands=$bands $sym order=$order 偶插入(:+,1)" begin
                _wf_check(lat, model, order, cn,
                          [OperatorInsertion(:+, 1, GTEMPO._single_band_op(:n, 1, bands))])
            end
            @testset "real bands=$bands $sym order=$order greater 式奇双插" begin
                ins = [OperatorInsertion(:+, 2, a), OperatorInsertion(:+, 1, adag)]
                _wf_check(lat, model, order, cn, ins)
            end
            @testset "real bands=$bands $sym order=$order lesser 式奇双插" begin
                ins = [OperatorInsertion(:-, 1, adag), OperatorInsertion(:+, 2, a)]
                _wf_check(lat, model, order, cn, ins)
            end
            @testset "real bands=$bands $sym order=$order lesser 式奇双插(j=0 链末切位)" begin
                ins = [OperatorInsertion(:-, 0, adag), OperatorInsertion(:+, 2, a)]
                _wf_check(lat, model, order, cn, ins)
            end
            # 实时专属：热化对角初态 + 奇双插（initial_state 覆盖路径）
            q = exp(-1.0 * (-0.5))
            ρ0 = LinearAlgebra.Diagonal(vcat([1.0],
                    [q^count_ones(st) for st in 1:(2^bands - 1)]) ./ (1 + 2q + q * q))
            @testset "real bands=$bands $sym order=$order 热化初态 + 奇双插" begin
                ins = [OperatorInsertion(:+, 2, a), OperatorInsertion(:+, 1, adag)]
                _wf_check(lat, model, order, cn, ins; initial_state=ρ0)
            end
        end
    end

    # ---- mixed：Kadanoff 轮廓抽查（对称性覆盖同 real 段）----
    for (bands, model, sym) in ((1, model1, U1Symmetry()), (2, model2, U1U1Symmetry()),
                                (2, model2, U1Symmetry()))
        lat = GrassmannLattice(Nt=2, Nτ=2, δt=0.1, δτ=0.25, bands=bands, contour=:mixed, symmetry=sym)
        a = GTEMPO._single_band_op(:a, 1, bands)
        adag = GTEMPO._single_band_op(:adag, 1, bands)
        for order in (1, 2)
            cn = order == 2 ? 0.7 : 0.0
            @testset "mixed bands=$bands $sym order=$order 无插入" begin
                _wf_check(lat, model, order, cn, OperatorInsertion[])
            end
            @testset "mixed bands=$bands $sym order=$order τ 奇双插 + rt 偶插" begin
                ins = [OperatorInsertion(:τ, 2, a),
                       OperatorInsertion(:τ, 1, adag),
                       OperatorInsertion(:+, 1, GTEMPO._single_band_op(:n, 1, bands))]
                _wf_check(lat, model, order, cn, ins)
            end
            @testset "mixed bands=$bands $sym order=$order rt 链末奇双插(:-,0)" begin
                ins = [OperatorInsertion(:-, 0, adag), OperatorInsertion(:+, 1, a)]
                _wf_check(lat, model, order, cn, ins)
            end
        end
    end

    # ---- 端到端：integrate(K_new) vs integrate(K_base) ----
    for (bands, model, sym) in ((1, model1, U1Symmetry()), (2, model2, U1U1Symmetry()))
        lat = GrassmannLattice(N=2, δτ=0.25, bands=bands, contour=:imag, symmetry=sym)
        for insertions in (OperatorInsertion[],
                           [OperatorInsertion(:τ, 1, GTEMPO._single_band_op(:n, 1, bands))],
                           [OperatorInsertion(:τ, 2, GTEMPO._single_band_op(:a, 1, bands)),
                            OperatorInsertion(:τ, 1, GTEMPO._single_band_op(:adag, 1, bands))])
            K_base = sysdynamics(lat, model; order=1, insertions=insertions)
            K_new = _wf_build(lat, model, 1, 0.0, insertions, nothing)
            @testset "integrate 端到端 bands=$bands ins=$(length(insertions))" begin
                Zb = integrate(K_base)
                Zn = integrate(K_new)
                @test abs(Zn - Zb) <= 1e-10 * (1 + abs(Zb))
            end
        end
    end
end
