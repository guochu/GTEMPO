# =============================================================================
# src/verify_formb.jl —— FormB（GTEMPO GMPS 消元装配）与 IFGTEMPO 对齐验证
#
# 参照链：同一 schedule 的 Fock 离散链（../form_b/formb_dense.jl 的
# fock_chain_single / fock_chain_two_spin）+ IFGTEMPO 两扇区 GMPS
# （gtempo_density_matrix_direct，权威机器）。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/src/verify_formb.jl
# 快测：FORMB_QUICK=1 同上（只跑 order 1 N=1 单自旋 + 两自旋）
# =============================================================================
using LinearAlgebra, Printf
using GTEMPO

const IFT_ROOT = joinpath(dirname(dirname(dirname(pathof(GTEMPO)))), "IFGTEMPO")
isdir(IFT_ROOT) || error("IFGTEMPO 仓库未找到：$IFT_ROOT")
push!(LOAD_PATH, IFT_ROOT)
using IFGTEMPO

include(joinpath(@__DIR__, "formb_dense.jl"))   # Fock 参照链（含 T_TOT/BETA/EPS_D）
include(joinpath(@__DIR__, "..", "src", "FormB.jl"))
using .FormB

passfail(ok) = ok ? "PASS" : "FAIL"
maxabsdiff(A, B) = maximum(abs.(A .- B))

"""单自旋 Fock 参照链（转发 formb_dense 实现）。"""
fock1(up, eps, sched, rho0) = fock_chain_single(up, eps, sched, rho0)

"""单自旋 FormB 路径 + 收缩。"""
function single_case(up, eps, rho0, order, N; counter_c=1, trunc=nothing)
    dt = T_TOT / N
    path = FormB.formb_path((up,), (eps = eps, rho0 = rho0), N, dt;
                            order = order, counter_c = counter_c,
                            trunc = something(trunc, GTEMPO.notrunc()))
    X = path.factors
    Z = FormB.formb_path_integral(X...)
    rho = FormB.formb_density_matrix_single(X...)
    return (Z = Z, rho = rho, sched = path.sched,
            maxbond = path.diagnostics.max_bond)
end

"""两自旋 FormB 路径 + 收缩。"""
function two_spin_case(up, dn, imp, order, N; counter_c=1, trunc=nothing)
    dt = T_TOT / N
    path = FormB.formb_path((up, dn), imp, N, dt;
                            order = order, counter_c = counter_c,
                            trunc = something(trunc, GTEMPO.notrunc()))
    X = path.factors
    Z = FormB.formb_path_integral(X...)
    rho = FormB.formb_density_matrix(X...)
    return (Z = Z, rho = rho, sched = path.sched,
            maxbond = path.diagnostics.max_bond)
end

"""约定标定：反测度符号 × 因子阶（单自旋 N=1 order 1，Z=1 证书）。"""
function calibrate()
    up = FiniteBath([0.4;;], [0.23 + 0.11im], BETA)
    rho0 = ComplexF64[1 0; 0 0]
    println("═"^72)
    println("标定：反测度符号 counter_c × ε_d=0.7（单自旋 N=1 order 1，Z=1 + ρ=Fock）")
    println("═"^72)
    sched = trotter_schedule(1, 1, T_TOT)
    rho_ref = fock1(up, EPS_D, sched, rho0)
    found = nothing
    for cc in (-1, 1, 0)
        res = single_case(up, EPS_D, rho0, 1, 1; counter_c = cc)
        errz = abs(res.Z - 1)
        errr = maxabsdiff(res.rho, rho_ref)
        @printf("  counter_c = %+d : Z-1 = %.3e, ρ-err = %.3e\n", cc, errz, errr)
        if errz < 1e-10 && errr < 1e-10 && found === nothing
            found = cc
        end
    end
    found === nothing && error("标定失败：没有反测度符号给出 Z=1 且 ρ=Fock")
    @printf("  → counter_c = %+d 锁定。\n\n", found)
    return found
end

"""单自旋电池：orders 1/2 × N=1..8，order 4 × N=1..3。"""
function battery_single(cc)
    up = FiniteBath([0.4;;], [0.23 + 0.11im], BETA)
    rho0 = ComplexF64[1 0; 0 0]
    println("═"^72)
    println("单自旋 FormB GMPS vs Fock 离散链（T=$T_TOT, ε_d=$EPS_D）")
    println("═"^72)
    quick = get(ENV, "FORMB_QUICK", "0") == "1"
    # 精确档（notrunc）：Δ_eff 的 G⁻¹ 非局域 → bond 按 4^N 增长，小 N 逐位
    return quick ? [(1, 1)] :
           vcat([(o, n) for o in (1, 2) for n in 1:4], [(4, 1)])
end

function battery_single_trunc(cc)
    up = FiniteBath([0.4;;], [0.23 + 0.11im], BETA)
    rho0 = ComplexF64[1 0; 0 0]
    println("═"^72)
    println("单自旋 FormB GMPS vs Fock——截断档（D≤512, ϵ=1e-12，同 twosector 生产用法）")
    println("═"^72)
    quick = get(ENV, "FORMB_QUICK", "0") == "1"
    trunc_grid = quick ? [] : [(o, 6) for o in (1, 2)]
    all_ok = true
    for (order, N) in trunc_grid
        sched = trotter_schedule(order, N, T_TOT / N)
        rho_ref = fock1(up, EPS_D, sched, rho0)
        trunc = GTEMPO.truncdimcutoff(D = 512, ϵ = 1e-12, add_back = 0)
        res = single_case(up, EPS_D, rho0, order, N; counter_c = cc, trunc = trunc)
        errz = abs(res.Z - 1)
        errr = maxabsdiff(res.rho, rho_ref)
        ok = errz < 1e-6 && errr < 1e-6
        all_ok &= ok
        @printf("  order %d N=%d [截断] : Z-1 = %.2e, ρ-err = %.2e, bond = %d  %s\n",
                order, N, errz, errr, res.maxbond, passfail(ok))
    end
    println()
    return all_ok
end

function battery_single_exact(cc)
    up = FiniteBath([0.4;;], [0.23 + 0.11im], BETA)
    rho0 = ComplexF64[1 0; 0 0]
    println("═"^72)
    println("单自旋 FormB GMPS vs Fock 离散链——精确档（T=$T_TOT, ε_d=$EPS_D）")
    println("═"^72)
    grid = battery_single(cc)
    all_ok = true
    for (order, N) in grid
        sched = trotter_schedule(order, N, T_TOT / N)
        rho_ref = fock1(up, EPS_D, sched, rho0)
        res = single_case(up, EPS_D, rho0, order, N; counter_c = cc)
        errz = abs(res.Z - 1)
        errr = maxabsdiff(res.rho, rho_ref)
        n_ex = real(rho_ref[2, 2])
        n_pb = real(res.rho[2, 2])
        ok = errz < 1e-7 && errr < 1e-7
        all_ok &= ok
        @printf("  order %d N=%2d : Z-1 = %.2e, ρ-err = %.2e, ⟨n⟩ = %.10f (Fock %.10f), bond = %d  %s\n",
                order, N, errz, errr, n_pb, n_ex, res.maxbond, passfail(ok))
    end
    println()
    return all_ok
end

"""两自旋模型构造。"""
function two_spin_models()
    up = FiniteBath([0.4 0.1; 0.1 0.3], [0.23 + 0.11im, 0.15 - 0.05im], BETA)
    dn = FiniteBath([0.38 0.08; 0.08 0.28], [0.19 - 0.07im, 0.12 + 0.04im], BETA)
    imp = TwoSpinImpurityModel(0.35, -0.22, 1.0)
    return up, dn, imp
end

"""两自旋电池（精确档）：FormB vs Fock vs IFGTEMPO GMPS，N ≤ 2。"""
function battery_two_spin_exact(cc)
    up, dn, imp = two_spin_models()
    println("═"^72)
    println("两自旋 FormB GMPS vs Fock vs IFGTEMPO GMPS——精确档（三方对齐）")
    println("═"^72)
    quick = get(ENV, "FORMB_QUICK", "0") == "1"
    exact_grid = quick ? [(1, 1)] : [(o, n) for o in (1, 2) for n in 1:2]
    all_ok = true
    for (order, N) in exact_grid
        dt = T_TOT / N
        sched = trotter_schedule(order, N, dt)
        rho_ref = fock_chain_two_spin(up, dn, imp, sched)
        path2 = build_impurity_path(up, dn, N, dt; order = order,
                                    impurity_model = imp, cutoff = 1e-12)
        rho_gmps = gtempo_density_matrix_direct(path2.factors...)
        res = two_spin_case(up, dn, imp, order, N; counter_c = cc)
        errz = abs(res.Z - 1)
        errf = maxabsdiff(res.rho, rho_ref)
        errg = maxabsdiff(res.rho, rho_gmps)
        ok = errz < 1e-7 && errf < 1e-7 && errg < 1e-8
        all_ok &= ok
        @printf("  order %d N=%d : Z-1 = %.2e | FormB-Fock %.2e | FormB-GMPS %.2e | bond = %d  %s\n",
                order, N, errz, errf, errg, res.maxbond, passfail(ok))
    end
    println()
    return all_ok
end

"""两自旋电池（截断档，D≤512）：N=3。"""
function battery_two_spin_trunc(cc)
    up, dn, imp = two_spin_models()
    println("═"^72)
    println("两自旋 FormB GMPS vs Fock vs IFGTEMPO GMPS——截断档（D≤512, ϵ=1e-12）")
    println("═"^72)
    quick = get(ENV, "FORMB_QUICK", "0") == "1"
    trunc_grid = quick ? [] : [(o, 3) for o in (1, 2)]
    all_ok = true
    for (order, N) in trunc_grid
        dt = T_TOT / N
        sched = trotter_schedule(order, N, dt)
        rho_ref = fock_chain_two_spin(up, dn, imp, sched)
        path2 = build_impurity_path(up, dn, N, dt; order = order,
                                    impurity_model = imp, cutoff = 1e-12)
        rho_gmps = gtempo_density_matrix_direct(path2.factors...)
        trunc = GTEMPO.truncdimcutoff(D = 512, ϵ = 1e-12, add_back = 0)
        res = two_spin_case(up, dn, imp, order, N; counter_c = cc, trunc = trunc)
        errz = abs(res.Z - 1)
        errf = maxabsdiff(res.rho, rho_ref)
        errg = maxabsdiff(res.rho, rho_gmps)
        ok = errz < 1e-6 && errf < 1e-6 && errg < 1e-6
        all_ok &= ok
        @printf("  order %d N=%d [截断]: Z-1 = %.2e | FormB-Fock %.2e | FormB-GMPS %.2e | bond = %d  %s\n",
                order, N, errz, errf, errg, res.maxbond, passfail(ok))
    end
    println()
    return all_ok
end

"""分块运行（限单进程峰值内存）：ARGS = calib | single_exact |
single_trunc | two_exact | two_trunc（默认全跑，建议分块调用）。"""
function main()
    stage = isempty(ARGS) ? "all" : ARGS[1]
    println("FormB（GTEMPO GMPS 消元装配）与 IFGTEMPO 对齐验证  [stage: ", stage, "]")
    @printf("T = %.2f, β = %.1f, ε_d = %.2f\n\n", T_TOT, BETA, EPS_D)
    cc = -1   # 已标定：counter_c = −1（GTEMPO 隐式 g_trace 与 dense 约定的取向差）
    if stage in ("all", "calib")
        cc = calibrate()
        stage == "calib" && return
    end
    if stage in ("all", "single_exact")
        battery_single_exact(cc)
    end
    if stage in ("all", "single_trunc")
        battery_single_trunc(cc)
    end
    if stage in ("all", "two_exact")
        battery_two_spin_exact(cc)
    end
    if stage in ("all", "two_trunc")
        battery_two_spin_trunc(cc)
    end
end

main()
