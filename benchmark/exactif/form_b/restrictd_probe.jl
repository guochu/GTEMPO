# =============================================================================
# form_b/restrictd_probe.jl —— Q1 判决实验：消元核 D 的 bar→eta 块是否足够
#
# GTEMPO 的 partialif 行扫（hybriddynamics/partialif_hybrid）只能装载
# exp(Σ_{bar 行, eta 列} C[r,c] θ̄_r θ_c)——即 D 的 (bar,eta) 有向块。
# 由于 D 反对称，(bar,eta) 腿对的单项式在 site 序下与全 D 逐位一致
# （θ_eta θ_bar ≡ −θ_bar θ_eta = −D[bar,eta]·(…)），因此限制后**丢弃的只有
# bar-bar 与 eta-eta 单项式**。本实验在稠密系数级直接对照：
#
#   全 D 消元式收缩   vs   丢弃同 kind 对的消元式收缩
#
# 若 (Z, ρ) 逐位一致 ⇒ 行扫可完整承载消元核 ⇒ partialif_hybrid 路线可行。
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/form_b/restrictd_probe.jl
# =============================================================================
using Printf
using GTEMPO

const IFT_ROOT = joinpath(dirname(dirname(dirname(pathof(GTEMPO)))), "IFGTEMPO")
isdir(IFT_ROOT) || error("IFGTEMPO 仓库未找到：$IFT_ROOT")
push!(LOAD_PATH, IFT_ROOT)
using IFGTEMPO

# 复用稠密机器（formb_dense.jl 仅在作为主脚本时运行 main）
include(joinpath(@__DIR__, "formb_dense.jl"))

"""消元核因子（可按 leg kind 过滤）。

`same_kind::Symbol`：:keep = 全 D；:drop = 丢弃 (bar,bar) 与 (eta,eta) 对
（即只保留 bar→eta 有向块；反对称性保证混合对的 site 序无歧义）。"""
function delta_eff_factors_q1(smap, Ginv::AbstractMatrix, N::Int, f::Int,
                              elim_sign::Int; mode::Symbol=:keep)
    GA = (Ginv .- transpose(Ginv)) ./ 2
    legs = Tuple{Int,Int,Symbol}[]   # (partner site, G 槽位, stat)
    for t in 1:N, br in (:+, :-), stat in (:bar, :eta)
        push!(legs, (smap[partner_var(t, br, stat, f)], gslot(t, br, stat), stat))
    end
    sort!(legs, by = x -> x[1])
    pairs = Tuple{Int,Int,ComplexF64}[]
    ndrop = 0
    for i in eachindex(legs), k in (i + 1):length(legs)
        v = elim_sign * GA[legs[i][2], legs[k][2]]
        iszero(v) && continue
        if mode === :drop && legs[i][3] == legs[k][3]
            ndrop += 1
            continue
        end
        push!(pairs, (legs[i][1], legs[k][1], v))
    end
    return bilinear_row_factors(pairs), ndrop
end

"""消元式稠密收缩（Q1 版）：mode = :keep / :drop。返回 (Z, ρ)。"""
function dense_case_q1(; order::Int, N::Int, models, imp, elim_sign::Int,
                       pf_scale, mode::Symbol, conv)
    sched = trotter_schedule(order, N, T_TOT / N)
    nl = length(sched.cells)
    smap0, vars, Lfull = build_sites(nl, 1)
    cmap, L = compact_map(smap0, vars, v -> v.typ === :imp)
    smap = cmap
    rho0 = imp.rho0
    Him = ComplexF64[0 0; 0 imp.eps]

    dfac, ndrop = delta_eff_factors_q1(smap, models[1].Ginv, nl, 1, elim_sign; mode)
    stack = Monomials[]
    append!(stack, dfac)
    append!(stack, layer0_measure_factors(smap, 1))
    append!(stack, k_kernel_factors(smap, nl, 1, sched.links, Him, rho0))

    dense = embed(L, stack[1])
    dense .*= pf_scale
    buf = similar(dense)
    for k in 2:length(stack)
        mul_monomials!(buf, dense, L, stack[k])
        dense, buf = buf, dense
    end

    openvars = [Var(:imp, nl + 1, :+, :bar, 1), Var(:imp, nl + 1, :-, :eta, 1)]
    keep = [smap[v] for v in openvars]
    intsites = [s for s in 0:(L - 1) if all(s != k for k in keep)]
    poly = integrate_out(dense, L, intsites; step_sign = conv.step_sign,
                         order = conv.order)
    rho = open_rho(poly, smap, nl, 1)

    zclosures = closure_factors(smap, nl, 1, conv.closure_k)
    Zpoly = dense
    for zc in zclosures
        mul_monomials!(buf, Zpoly, L, zc)
        Zpoly, buf = buf, Zpoly
    end
    Z = integrate_out(Zpoly, L, collect(0:(L - 1)); step_sign = conv.step_sign,
                      order = conv.order)[1]
    return (Z = Z, rho = rho, ndrop = ndrop)
end

function main()
    up = single_models().up
    imp = single_imp()
    println("═"^72)
    println("Q1 判决：消元核 D 丢弃 (bar,bar)/(eta,eta) 对后收缩是否逐位不变")
    println("═"^72)
    # 锁定约定（formb_dense 全量运行的结果）：step_sign=1, order=:inc, closure_k=1
    conv = Conv(1, :inc, 1)
    all_ok = true
    grid = vcat([(o, n) for o in (1, 2) for n in (1, 2)], [(4, 1)])
    for (order, N) in grid
        sched = trotter_schedule(order, N, T_TOT / N)
        nl = length(sched.cells)
        data = influence_matrix_G(up, sched.cells)
        Ginv = inv(data.G)
        pf = pf_mine(data, nl, 1; step_sign = conv.step_sign, order = conv.order)
        full = dense_case_q1(order = order, N = N, models = [(Ginv = Ginv,)],
                             imp = imp, elim_sign = -1, pf_scale = pf,
                             mode = :keep, conv = conv)
        drop = dense_case_q1(order = order, N = N, models = [(Ginv = Ginv,)],
                             imp = imp, elim_sign = -1, pf_scale = pf,
                             mode = :drop, conv = conv)
        ref = fock_chain_single(up, EPS_D, sched, imp.rho0)
        dz = abs(full.Z - drop.Z)
        dr = maxabsdiff(full.rho, drop.rho)
        okr = maxabsdiff(drop.rho, ref) < 1e-9
        all_ok &= (dz < 1e-12 && dr < 1e-12 && okr)
        @printf("  order %d N=%d (nl=%2d, 丢弃 %2d 个同 kind 对) : ΔZ = %.2e, Δρ = %.2e, 丢弃版 vs Fock %.2e  %s\n",
                order, N, nl, drop.ndrop, dz, dr,
                maxabsdiff(drop.rho, ref), (dz < 1e-12 && dr < 1e-12) ? "PASS" : "FAIL")
    end
    println("═"^72)
    println("结论：", all_ok ? "bar→eta 块（经反对称性）完整承载消元核 —— partialif 行扫可行" :
            "存在真实结构差异 —— 行扫不充分，需要进一步处理")
    return all_ok
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
