include(joinpath("/home/guochu/Documents/Missile/GTEMPO", "benchmark", "exactif", "src", "NativeExact.jl"))
using .NativeExact
using GTEMPO
using Printf
using Z2Tensors: NoTruncation
const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(GT_ROOT, "benchmark", "exactif", "form_b", "native_probe.jl"))

function gmps_marginal(K, I)
    env = NativeExact._right_sweep1((K, I))
    env = NativeExact._absorb_site2(env, K.data[4], I.data[4], 0)
    env = NativeExact._absorb_site2(env, K.data[3], I.data[3], 1)
    env = NativeExact._absorb_site2(env, K.data[2], I.data[2], 2)
    env = NativeExact._absorb_site2(env, K.data[1], I.data[1], 3)
    env = GTEMPO.TK.lmul!(GTEMPO.scaling(K) * GTEMPO.scaling(I), env)
    poly = zeros(ComplexF64, 16)
    for (f1, f2) in GTEMPO.TK.fusiontrees(env)
        unc1 = f1.uncoupled; unc2 = f2.uncoupled
        (unc2[5].n == 0 && unc1[1].n == 0 && unc1[2].n == 0) || continue
        mask = 0
        for (leg, bit) in ((2, 2), (3, 3), (4, 0), (5, 1))
            (unc2[leg].n == 1) && (mask |= 1 << bit)
        end
        poly[mask + 1] += env[f1, f2][1, 1]
    end
    return poly
end

function dense_marginal(sched, up0, imp0, n)
    G = native_G(up0.bath_hamiltonian, up0.hybridization, BETA, collect(Float64.(sched.cells)), n; trace_style = :edge)
    pf = pf_mine_native(G; step_sign = 1, order = :inc)
    pairs = native_D_pairs(n, G)
    stack = Monomials[]
    append!(stack, bilinear_row_factors(pairs))
    # fold=neg：折返因子包含在内（与 GMPS 管线的 boundarycondition 对应）
    append!(stack, native_skeleton_factors(n, collect(Float64.(sched.links)), zeros(ComplexF64,2,2), imp0.rho0; fold = :neg, res = :boundary))
    Ld = n.L
    d = embed(Ld, stack[1]); d .*= pf
    b = similar(d)
    for kk in 2:length(stack)
        mul_monomials!(b, d, Ld, stack[kk]); d, b = b, d
    end
    keep = [8, 9, 10, 11]
    intsites = [s for s in 0:(Ld-1) if all(s != kk for kk in keep)]
    full = integrate_out(d, Ld, intsites; step_sign = 1, order = :inc)
    return [full[(m << 8) + 1] for m in 0:15]
end

function main()
    T_TOT = 0.4; BETA = 1.3
    HB = ComplexF64[0.4;;]
    sched = NativeExact.native_schedule(2, 1, T_TOT)
    model = ToulouseIM(ϵ_d=0.7)
    pipe = exactif_native_pipeline(((hb = HB, v = ComplexF64[0.0+0im]),), model, 1, T_TOT;
                                   order=2, trunc=NoTruncation(), β=BETA, bands=1)
    mg = gmps_marginal(pipe.mpsK, pipe.mpsIs[1])
    n = Nat(length(sched.cells), length(sched.cells) + 2, 4 * (length(sched.cells) + 2))
    up0, imp0 = setup(eps_d = 0.0)
    md = dense_marginal(sched, up0, imp0, n)
    println("bar⁺eta⁺bar⁻eta⁻ :        dense        |         gmps")
    for m in 0:15
        d = md[m + 1]; g = mg[m + 1]
        (abs(d) > 1e-10 || abs(g) > 1e-10) || continue
        s = join(reverse([(m >> i) & 1 for i in 0:3]), "")
        @printf("%s : %+9.5f%+9.5f | %+9.5f%+9.5f\n", s, real(d), imag(d), real(g), imag(g))
    end
end
main()

# =============================================================================
# calibrate_native.jl —— 开放腿读出的约定标定工作台（进行中）
#
# 用法：julia --project=<GTEMPO 根> benchmark/exactif/src/calibrate_native.jl
#
# 当前状态（2026-11-05）：
#   * T1 已过：右扫 update_pair_right 全收缩 == integrate Z（机器精度）；
#   * T3（本脚本，V=0 自由情形，含折返对照）判定：GMPS 复合多项式的边界
#     4 腿边缘与稠密参照存在真实内容差异——稠密 = {0000, 0110(折返双腿),
#     1001, 1111}；GMPS = 6 分量（折返双腿模式缺失 + 单占据伪分量）。
#     这不是全局符号（比例不随 sc/符号变），是约定组合未对齐。
#   * 已排除：GTerm 插入式读出（原生 occupation/gf 路线，值互补错误）；
#     简单旋钮 {counter_c, node1_c, ρ_sign} 网格（22 组合未命中）。
#   * 下一步候选（按优先级）：
#     1. node-1 隐式测度符号：lattice pair = (eta,bar) 序，隐式 g_trace
#        ≡ (1−θηθb)（FormB compact 标定），而稠密 M₀ = (1−θ̄θ) = (1+θηθb)
#        ——修正需 (1−θηθb)(1+c·θηθb) = (1+θηθb) 即 c=2，网格已试未过，
#        需与 fold/吸收序联合重推；
#     2. 折返积分方式：g_trace_phy 迹 vs 稠密"折返因子+素积分"（两种语义
#        在 fold 腿上未必等价——素积分 = 单腿 Berezin 抽取，需实现）；
#     3. 吸收序（T1 先/T2 先）与 g_fuse 的费米符号；
#     4. partialmpo 反测度的净符号（case 处理内建排列符号）。
#   * 判决协议：V=0 → ρ = diag(1,0) 逐位；εd≠0 → ρ₂₂ = 0.003856；
#     参照 = form_b/native_probe.jl 的稠密机器（注意其对 g_trace 语义的
#     稠密仿真有 δt³ 级偏差，结构性判决可用、数值判决以真实机器为准）。
# =============================================================================
