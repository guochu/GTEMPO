import TensorKit as TK
import Random
import LinearAlgebra
# fockmpo.jl：FockMPO 直积 MPO + fockmpo2gmps 闸门（plan V1 Phase 1；
# 2026-11-18 起 = 公开 sysdynamics vs 手工闭环 workflow 接线对拍——迁移期的
# 旧径锚定（_system_factor_reference，203 项全绿）已随 V2 删除退役，绝对
# 数值锚由 observables golden（U1GTEMPO2 机器精度）与 ED 电池承接）
#
# 对拍 = 公开 sysdynamics（纯委托实现体）vs 手工闭环（fock_sysdynamics →
# sparse * → fockmpo2gmps kwargs）todense 机器精度。配置矩阵：
#   {imag, real, mixed} × bands{1,2} × spec{:z2,:u1,:u1u1}（:z2z2 抽查）
#   × 插入{无、偶（n̂ cut0/内/末）、奇（greater/lesser 式双插带荷）、多插}
#   × order{1,2}
# 端到端：integrate(K_new) vs integrate(K_base)（无插入 + 单偶插 + 奇双插）。
@testset "FockMPO: fockmpo2gmps 对拍矩阵" begin

    function _wf_build(lattice, model, order, coupling_norm, insertions, initial_state)
        # 新路径：切位（cut index）≥ 1 走 sparsefockmpo 乘法（左右各一）；切位 0
        # 三明治走 fockmpo2gmps 的 rho0/insertions kwargs（fold 在转换内做，不预折）。
        # ⚠️ 过滤按 insertion_cut_index 而非 boundary 字段——:- 支 boundary 0 =
        # 链末切位（cut = ni+2Nt）非 ρ₀ 切位（与 src/fockmpo.jl 同口径）
        cut0 = filter(i -> GTEMPO.insertion_cut_index(lattice, i) == 0, insertions)
        ge1 = filter(i -> GTEMPO.insertion_cut_index(lattice, i) >= 1, insertions)
        Kf = fock_sysdynamics(lattice, model; order=order)
        Wl = sparsefockmpo(lattice, filter(i -> i.side === :left, ge1))
        Wr = sparsefockmpo(lattice, filter(i -> i.side === :right, ge1))
        Kf = Wr * (Wl * Kf)
        rho0 = initial_state === nothing && isempty(cut0) ? nothing :
               (initial_state === nothing ? nothing :
                (initial_state isa FockMatrix ? Matrix{ComplexF64}(initial_state.data) :
                 Matrix{ComplexF64}(initial_state)))
        lscale = order == 2 ? real(coupling_norm / 4 *
                                   sum(GTEMPO._cell_step2(lattice, b) for (b, _) in GTEMPO.chronological_cells(lattice))) : 0.0
        # Weyl 标量由调用方后置累加补上（2026-11-18 起转换不再收 weyl_lscale；
        # ⚠️ logscale += lscale——压缩归一化因子已在 logscale 内，覆盖会抹掉）
        out = fockmpo2gmps(Kf, lattice; rho0=rho0, insertions=insertions)
        GM.setlogscale!(out, GM.logscale(out) + lscale)
        return out
    end

    function _wf_check(lattice, model, order, coupling_norm, insertions; initial_state=nothing)
        # 基线 = 公开入口（纯委托 _system_factor_fock）——对拍手工闭环 workflow
        # （_wf_build）：kwargs 接线合同（rho0 轮廓门控/切位过滤/prefix/oddlater/
        # lscale 路由）与规范装配一致
        K_base = sysdynamics(lattice, model; order=order, coupling_norm=coupling_norm,
                             insertions=insertions, initial_state=initial_state)
        K_new = _wf_build(lattice, model, order, coupling_norm, insertions, initial_state)
        nst = length(K_base)
        if 4^nst <= 2^20
            tb = GM.todense(K_base)
            tn = GM.todense(K_new)
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

    model1 = ImpurityHamiltonian([tunneling(1, 1, coeff=0.5)]; bands=1)
    model2 = KanamoriIM(U=1.0, J=0.2, norb=1)

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
                    n = GTEMPO._number_op(bands, nothing)
                    for jj in (0, 1, 2)
                        _wf_check(lat, model, order, cn, [OperatorInsertion(:τ, jj, n)])
                    end
                end
            end
        end
        # ---- imag：奇插入（带荷双插）+ 多插 ----
        for (bands, model, sym) in ((1, model1, U1Symmetry()), (2, model2, U1U1Symmetry()))
            lat = GrassmannLattice(N=2, δτ=0.25, bands=bands, contour=:imag, symmetry=sym)
            qc = var_charge(TK.sectortype(lat), 1, bands)
            d = GTEMPO._fock_ops(bands)[1]; ddag = d'
            @testset "imag bands=$bands $sym order=$order 奇双插(:τ)" begin
                ins = [OperatorInsertion(:τ, 2, d; charge=GTEMPO.qneg(qc)),
                       OperatorInsertion(:τ, 1, ddag; charge=qc)]
                _wf_check(lat, model, order, cn, ins)
            end
            @testset "imag bands=$bands $sym order=$order 多插(偶+奇)" begin
                ins = [OperatorInsertion(:τ, 0, GTEMPO._number_op(bands, nothing)),
                       OperatorInsertion(:τ, 1, d; charge=GTEMPO.qneg(qc)),
                       OperatorInsertion(:τ, 2, ddag; charge=qc)]
                _wf_check(lat, model, order, cn, ins)
            end
        end
    end

    # ---- real：真空 ρ₀ + 偶/奇插入（bands=2 覆盖 U1U1 与单 U1——单 U1×2 带
    # 曾漏网：:- 支 boundary 0 = 链末切位，lesser(j=0) 路径）----
    for (bands, model, sym) in ((1, model1, U1Symmetry()), (2, model2, U1U1Symmetry()),
                                (2, model2, U1Symmetry()))
        lat = GrassmannLattice(N=2, δt=0.1, bands=bands, contour=:real, symmetry=sym)
        qc = var_charge(TK.sectortype(lat), 1, bands)
        d = GTEMPO._fock_ops(bands)[1]; ddag = d'
        for order in (1, 2)
            cn = order == 2 ? 0.7 : 0.0
            @testset "real bands=$bands $sym order=$order 无插入" begin
                _wf_check(lat, model, order, cn, OperatorInsertion[])
            end
            @testset "real bands=$bands $sym order=$order 偶插入(:+,1)" begin
                _wf_check(lat, model, order, cn, [OperatorInsertion(:+, 1, GTEMPO._number_op(bands, nothing))])
            end
            @testset "real bands=$bands $sym order=$order greater 式奇双插" begin
                ins = [OperatorInsertion(:+, 2, d; charge=GTEMPO.qneg(qc)),
                       OperatorInsertion(:+, 1, ddag; charge=qc)]
                _wf_check(lat, model, order, cn, ins)
            end
            @testset "real bands=$bands $sym order=$order lesser 式奇双插" begin
                ins = [OperatorInsertion(:-, 1, ddag; charge=qc),
                       OperatorInsertion(:+, 2, d; charge=GTEMPO.qneg(qc))]
                _wf_check(lat, model, order, cn, ins)
            end
            @testset "real bands=$bands $sym order=$order lesser 式奇双插(j=0 链末切位)" begin
                ins = [OperatorInsertion(:-, 0, ddag; charge=qc),
                       OperatorInsertion(:+, 2, d; charge=GTEMPO.qneg(qc))]
                _wf_check(lat, model, order, cn, ins)
            end
        end
    end

    # ---- mixed：Kadanoff 轮廓抽查（对称性覆盖同 real 段）----
    for (bands, model, sym) in ((1, model1, U1Symmetry()), (2, model2, U1U1Symmetry()),
                                (2, model2, U1Symmetry()))
        lat = GrassmannLattice(Nt=2, Nτ=2, δt=0.1, δτ=0.25, bands=bands, contour=:mixed, symmetry=sym)
        d = GTEMPO._fock_ops(bands)[1]; ddag = d'
        qc = var_charge(TK.sectortype(lat), 1, bands)
        for order in (1, 2)
            cn = order == 2 ? 0.7 : 0.0
            @testset "mixed bands=$bands $sym order=$order 无插入" begin
                _wf_check(lat, model, order, cn, OperatorInsertion[])
            end
            @testset "mixed bands=$bands $sym order=$order τ 奇双插 + rt 偶插" begin
                ins = [OperatorInsertion(:τ, 2, d; charge=GTEMPO.qneg(qc)),
                       OperatorInsertion(:τ, 1, ddag; charge=qc),
                       OperatorInsertion(:+, 1, GTEMPO._number_op(bands, nothing))]
                _wf_check(lat, model, order, cn, ins)
            end
            @testset "mixed bands=$bands $sym order=$order rt 链末奇双插(:-,0)" begin
                ins = [OperatorInsertion(:-, 0, ddag; charge=qc),
                       OperatorInsertion(:+, 1, d; charge=GTEMPO.qneg(qc))]
                _wf_check(lat, model, order, cn, ins)
            end
        end
    end

    # ---- fock_initstate! 与转换 kwargs 双路一致 ----
    lat = GrassmannLattice(N=2, δτ=0.25, bands=2, contour=:imag, symmetry=U1U1Symmetry())
    Kf = fock_sysdynamics(lat, model2)
    rho_fold = Matrix{ComplexF64}(LinearAlgebra.I, 4, 4)
    rho_fold[3, 3] = 0.5
    K_a = fockmpo2gmps(fock_initstate!(lat, copy(Kf), rho_fold), lat)
    K_b = fockmpo2gmps(copy(Kf), lat; rho0=rho_fold)
    @testset "fock_initstate! 与 kwargs 双路一致" begin
        @test GM.todense(K_a) ≈ GM.todense(K_b) atol = 1e-12
    end

    # ---- 端到端：integrate(K_new) vs integrate(K_base) ----
    for (bands, model, sym) in ((1, model1, U1Symmetry()), (2, model2, U1U1Symmetry()))
        lat = GrassmannLattice(N=2, δτ=0.25, bands=bands, contour=:imag, symmetry=sym)
        qc = var_charge(TK.sectortype(lat), 1, bands)
        d = GTEMPO._fock_ops(bands)[1]; ddag = d'
        for insertions in (OperatorInsertion[],
                           [OperatorInsertion(:τ, 1, GTEMPO._number_op(bands, nothing))],
                           [OperatorInsertion(:τ, 2, d; charge=GTEMPO.qneg(qc)),
                            OperatorInsertion(:τ, 1, ddag; charge=qc)])
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
