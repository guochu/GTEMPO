# =============================================================================
# validation/formb_legacy_ifgtempo.jl —— 形态 B 验证电池的 IFGTEMPO 依赖部分
#
# 从 formb_dense.jl 拆出的两个需要 IFGTEMPO 在作用域内的函数：
#   - fock_chain_two_spin：两自旋离散链参照（FiniteBath/TwoSpinImpurityModel 类型注解）
#   - step_d_two_spin：Step D 两自旋消元式 vs Fock vs IFGTEMPO GMPS 三方对齐
#     （timelocallattice / build_impurity_path / gtempo_density_matrix_direct）
#
# **先 `using IFGTEMPO` 再 include 本文件**（与 formb_dense.jl 同序即可）。
# 消费者：validation/verify_formb.jl（legacy 验证记录）。
# =============================================================================
using IFGTEMPO

"""两自旋离散链（mode 序：imp_up=0, imp_dn=1, up-bath, dn-bath）。"""
function fock_chain_two_spin(up::FiniteBath, dn::FiniteBath,
                             imp::TwoSpinImpurityModel, sched::TrotterSchedule)
    mu, md = bath_modes(up), bath_modes(dn)
    nt = 2 + mu + md
    dim = 1 << nt
    ops = [annih_op(dim, k) for k in 0:(nt - 1)]
    Hq = zeros(ComplexF64, dim, dim)
    Hb = zeros(ComplexF64, dim, dim)
    for (model, off, impmode) in ((up, 2, 0), (dn, 2 + mu, 1))
        for i in 1:bath_modes(model), j in 1:bath_modes(model)
            Hq .+= model.bath_hamiltonian[i, j] .* (ops[off + i]' * ops[off + j])
            Hb .+= model.bath_hamiltonian[i, j] .* (ops[off + i]' * ops[off + j])
        end
        for i in 1:bath_modes(model)
            Hq .+= conj(model.hybridization[i]) .* (ops[impmode + 1]' * ops[off + i])
            Hq .+= model.hybridization[i] .* (ops[off + i]' * ops[impmode + 1])
        end
    end
    Him = imp.epsilon_up .* (ops[1]' * ops[1]) +
          imp.epsilon_down .* (ops[2]' * ops[2]) +
          imp.interaction_u .* (ops[1]' * ops[1]) * (ops[2]' * ops[2])
    # ρ_b 嵌入全空间，按 bath 配分函数归一（杂质维度 4）
    Zb = tr(exp(-up.beta .* Hb)) / 4
    rb = exp(-up.beta .* Hb) ./ Zb
    U = chain_operator(Hq, Him, sched.cells, sched.links)
    r0 = zeros(ComplexF64, dim, dim)
    for b in 0:((1 << (nt - 2)) - 1), a in 0:3, a2 in 0:3
        r0[(b << 2 | a) + 1, (b << 2 | a2) + 1] = imp.initial_density[a + 1, a2 + 1]
    end
    rtot = U * (r0 * rb) * adjoint(U)
    rho = zeros(ComplexF64, 4, 4)
    for b in 0:((1 << (nt - 2)) - 1), a in 0:3, a2 in 0:3
        rho[a + 1, a2 + 1] += rtot[(b << 2 | a) + 1, (b << 2 | a2) + 1]
    end
    return rho
end

"""Step D：两自旋消元式收缩 vs Fock vs IFGTEMPO GMPS（权威对齐）。"""
function step_d_two_spin(elim_sign, pf_conv, int_conv)
    up = FiniteBath([0.4 0.1; 0.1 0.3], [0.23 + 0.11im, 0.15 - 0.05im], BETA)
    dn = FiniteBath([0.38 0.08; 0.08 0.28], [0.19 - 0.07im, 0.12 + 0.04im], BETA)
    imp = TwoSpinImpurityModel(0.35, -0.22, 1.0)
    println("═"^72)
    println("Step D：两自旋消元式收缩 vs Fock vs IFGTEMPO GMPS（gtempo_density_matrix_direct）")
    println("═"^72)
    # 断言：我们的两 flavor 布局 == timelocallattice
    for N in (1, 2)
        smap, vars, L = build_sites(N, 2)
        lat = timelocallattice(N)
        for v in vars
            g = IFGTEMPO.index(lat, GrassmannVariable(v.typ, v.t, v.br, v.stat, v.f))
            g == smap[v] || error("lattice 布局与 timelocallattice 不一致：$v")
        end
    end
    println("  lattice 布局 == timelocallattice（N=1,2 全变量逐位断言通过）")
    all_ok = true
    quick = get(ENV, "FORMB_QUICK", "0") == "1"
    # 两自旋消元 L = 12 + 8·nl：order 4（nl=3N）N=1 已是 2^36，dense 不做
    grid = quick ? [(1, 1)] : [(o, n) for o in (1, 2) for n in (1, 2)]
    for (order, N) in grid
        dt = T_TOT / N
        sched = trotter_schedule(order, N, dt)
        rho_ref = fock_chain_two_spin(up, dn, imp, sched)
        path = build_impurity_path(up, dn, N, dt; order = order,
                                   impurity_model = imp, cutoff = 1e-12)
        rho_gmps = gtempo_density_matrix_direct(path.factors...)
        data_up = influence_matrix_G(up, sched.cells)
        data_dn = influence_matrix_G(dn, sched.cells)
        nl4 = length(sched.cells)
        pf = pf_mine(data_up, nl4, 1; step_sign = pf_conv.step_sign,
                     order = pf_conv.order) *
             pf_mine(data_dn, nl4, 2; step_sign = pf_conv.step_sign,
                     order = pf_conv.order)
        res = dense_case(order = order, N = N, nflav = 2,
                         models = [(Ginv = inv(data_up.G),), (Ginv = inv(data_dn.G),)],
                         imp = imp, mode = :elim, conv = int_conv,
                         elim_sign = elim_sign, pf_scale = pf)
        errz = abs(res.Z - 1)
        errf = maxabsdiff(res.rho, rho_ref)
        errg = maxabsdiff(res.rho, rho_gmps)
        errgf = maxabsdiff(rho_gmps, rho_ref)
        ok = errz < 1e-8 && errf < 1e-8 && errg < 1e-7
        all_ok &= ok
        @printf("  order %d N=%d : Z-1 = %.2e | 消元-Fock %.2e | 消元-GMPS %.2e | GMPS-Fock %.2e  %s\n",
                order, N, errz, errf, errg, errgf, passfail(ok))
    end
    println()
    return all_ok
end
