# =============================================================================
# imag_dense.jl —— 严格 IF 虚轴·稠密路线（已全部验证，2026-11-07）
#
# 自洽的虚轴严格 IF 基准：全 kernel 多项式 + 暴力 Berezin 积分，机器精度对齐
# Fock 环 Trotter 链（Z 与 n(β)），orders 1/2/4 调度全支持。此路线同时是
# GMPS 管线（proto_imag.jl）的约定基准。
#
# 约定（全部锁定）：
#   1. 配对 = 每（位置 p，模 s）一对（s=1 杂质、2..1+m bath）；无门拆分。
#   2. 区间算符 W_p = G_p·C_p（门在胞之前、同 index；G₁ = 环合并门
#      links[1]+links[end]）。均匀 schedule 与 C·G 序循环等价，非均匀必须此序。
#   3. W_p 的多模相干核 = kernel_monomials（impurity_kernel_gmps 配方：
#      coeff = (−1)^{C(ket,2)}·K[bra+1,ket+1]，local 序 (bar…,eta…)）。
#   4. 迹 resolution：每模每边界因子 e^{−ξ̄_pξ_{p+1}}（单项式 coeff = −1，
#      缝位同号）。
#   5. 反周期缝：⟨−ξ̄| 的负号 = 任一区间（取末位）kernel 的 bra 侧取负，
#      即奇数 bar 单项式翻号（等价于原 M 装配的 M[N,N] = +U_N）。
#   6. 积分 = integrate_out 全体生成元（step_sign/order 不敏感——结构偶）。
#   7. n(β)：缝位插 n̂·W₁（Tr[n̂·∏W]/Tr[∏W]，环循环任意位等价）。
#
# 验证（β=1，单模 bath ω=1 α=0.5，εd=0.7）：
#   S0  εd=0（W=C）：Z == 多体链迹，N=2,3,4 机器精度（≤4.4e-16）
#   S2  εd≠0：Z 与 n(β) == Fock 环链，N=2,3,4 机器精度（≤5.6e-17）
#   S4  非均匀调度 orders 2/4（native_schedule + 环门合并）：机器精度
#       （order 4 N_int=3: err 2.2e-16）
#   S3  连续极限：n → ED 热态 ⟨n⟩，二阶收敛（N=4→5 误差 ×0.64 ≈ (4/5)²，
#       与"虚时环基线即二阶"先例一致）
#
# 内存保护：L = 2(1+m)N ≤ 22（闸在 dense_ring 内，N≈5@单模 bath）；
# mul_monomials! 的 out 与 dense 不得同数组（fill! 先清——每轮新分配）。
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_dense.jl
# =============================================================================

using LinearAlgebra
using Printf
using GTEMPO                    # formb_dense 依赖

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(GT_ROOT, "benchmark", "exactif", "validation", "formb_dense.jl"))
include(joinpath(GT_ROOT, "benchmark", "exactif", "src", "NativeExact.jl"))
using .NativeExact

# ------------------------------------------------------------------ 模型 --
const BETA = 1.0
const WB, AB = 1.0, 0.5
HB = ComplexF64[WB;;]
VB = ComplexF64[sqrt(AB) + 0.0im]
EPS_D = 0.7

"""耦合单粒子 h（mode 1 = 杂质 ε=0；mode k+1 = bath k）。"""
function hcoupled()
    m = length(VB)
    h = zeros(ComplexF64, m + 1, m + 1)
    for k in 1:m
        h[1, k+1] = conj(VB[k]); h[k+1, 1] = VB[k]; h[k+1, k+1] = HB[k, k]
    end
    return h
end

# ------------------------------------------------- 全 kernel 稠密收缩 --
"""配对索引：位置 p（1..N）× 模 s（1=imp，2..1+m=bath）→ 对索引。"""
pairidx(p, s, nm) = (p - 1) * nm + s

"""多体区间算符 W_p = C_p·G_p（dim=2^(1+m) Fock 矩阵）。

h = 耦合单粒子（杂质 ε=0）；C_p = exp(-cells_p·Ĥ₀)，G_p = exp(-durs_p·Ĥ_imp)。"""
function interval_operators(cells, durs, eps_d; U=0.0)
    nm = 1 + length(VB)
    dim = 1 << nm
    ops = [annih_op(dim, k) for k in 0:(nm-1)]
    Hq = zeros(ComplexF64, dim, dim)
    h = hcoupled()
    for i in 1:nm, j in 1:nm
        iszero(h[i, j]) && continue
        Hq .+= h[i, j] .* (ops[i]' * ops[j])
    end
    Him = eps_d .* (ops[1]' * ops[1])
    U == 0 || (Him .+= U .* (ops[1]' * ops[1] * ops[1]' * ops[1]))
    n1 = ops[1]' * ops[1]
    # 环链序：g₁ C₁ g₂ C₂ … g_N C_N（g₁ = 环合并门；门在胞之前、同 index）。
    # 均匀 schedule 下与 C·G 循环等价；非均匀（order 4）必须此序。
    N = length(cells)
    Ws = [exp(-durs[p] .* Him) * exp(-cells[p] .* Hq) for p in 1:N]
    return Ws, n1
end

"""全 kernel 稠密收缩：∏_p kernel(W_p) × 迹 resolution（全模）→ 暴力积分。

resolution 因子 = e^{−ξ̄_pξ_{p+1}}（每模，缝位含反周期号已并入原 M 验证结构：
M[p,p+1] = +1 于 −āMa 形式 ⇔ 因子多项式 1 − ξ̄ξ）。返回 (Z_raw, L)。"""
function dense_ring(Ws; step_sign::Int=1, order::Symbol=:inc)
    N = length(Ws)
    nm = size(Ws[1], 1) |> x -> round(Int, log2(x))   # 1 + m
    L = 2nm * N
    L > 22 && error("稠密规模 2^$L 过大（内存保护）")
    dense = zeros(ComplexF64, 1 << L)
    dense[1] = 1.0
    bar(p, s) = 2 * (pairidx(p, s, nm) - 1)
    eta(p, s) = bar(p, s) + 1
    for p in 1:N
        # 区间 kernel（全模，impurity_kernel_gmps 配方）；缝位反周期 ⟨−ξ̄|：
        # 末位置 kernel 的 bra 侧取负 = 奇数 bar 单项式翻号（与原 M 的
        # M[N,N] = +U_N 即 exp 系数 −U_N 等价）
        kmonos = kernel_monomials([bar(p, s) for s in 1:nm],
                                  [eta(p, s) for s in 1:nm], Ws[p])
        if p == N
            kmonos = [(loc, isodd(count(b -> b in [bar(N, s2) for s2 in 1:nm], loc)) ?
                       -c : c) for (loc, c) in kmonos]
        end
        monos = Monomials(kmonos)
        dense = mul_monomials!(zeros(ComplexF64, 1 << L), dense, L, monos)
        # 迹 resolution：每模 1 − ξ̄_p ξ_{p+1}（缝位同号，反周期已验证并入）
        nxt = p == N ? 1 : p + 1
        for s in 1:nm
            rmonos = Monomials([(Int[], 1.0 + 0.0im),
                                ([bar(p, s), eta(nxt, s)], -1.0 + 0.0im)])
            dense = mul_monomials!(zeros(ComplexF64, 1 << L), dense, L, rmonos)
        end
    end
    out = integrate_out(dense, L, collect(0:(L-1)); step_sign=step_sign, order=order)
    return out[1], L
end

# ------------------------------------------------------- Fock 链参照 --
"""环 Fock 链：Tr[∏_p(C_p·G_p)] 与 n̂ 插入（dim = 2·2^m 多体）。

门在 cell 后作用（装配序 C·G：G 先作用）。gateεs[p] = 门系数矩阵生成元。"""
function fock_ring(cells, gate_durs, eps_d)
    m = length(VB)
    nt = 1 + m
    dim = 1 << nt
    ops = [annih_op(dim, k) for k in 0:(nt-1)]
    Hq = zeros(ComplexF64, dim, dim)
    h = hcoupled()
    for i in 1:nt, j in 1:nt
        iszero(h[i, j]) && continue
        Hq .+= h[i, j] .* (ops[i]' * ops[j])
    end
    Him = eps_d .* (ops[1]' * ops[1])
    Z = Matrix{ComplexF64}(I, dim, dim)
    Zn = Matrix{ComplexF64}(I, dim, dim)
    n1 = ops[1]' * ops[1]
    for p in length(cells):-1:1
        W = exp(-cells[p] .* Hq) * exp(-gate_durs[p] .* Him)
        Z = W * Z
        Zn = W * Zn
    end
    Zn = n1 * Zn                       # n̂ 插在链末端（等价任意位置，cyclic）
    return Z, Zn, n1
end

"""连续 ED：热态 ⟨n⟩（εd, U）。"""
function ed_n(eps_d; U=0.0)
    m = length(VB)
    nt = 1 + m
    dim = 1 << nt
    ops = [annih_op(dim, k) for k in 0:(nt-1)]
    H = zeros(ComplexF64, dim, dim)
    h = hcoupled(); h[1, 1] += eps_d
    for i in 1:nt, j in 1:nt
        iszero(h[i, j]) && continue
        H .+= h[i, j] .* (ops[i]' * ops[j])
    end
    U == 0 || (H .+= U .* (ops[1]' * ops[1] * ops[1]' * ops[1]))
    n1 = ops[1]' * ops[1]
    R = exp(-BETA .* H)
    return real(tr(n1 * R) / tr(R))
end

# ================================================================== 主程序 --
function main()
    println("="^78)
    println("S0：εd=0（恒等门，W=C）自检——锁积分定向（对照多体链迹）")
    println("="^78)
    for N in (2, 3, 4)
        cells = fill(BETA / N, N)
        Ws, _ = interval_operators(cells, zeros(N), 0.0)
        FZ, _, _ = fock_ring(cells, zeros(N), 0.0)
        @printf("N=%d  多体链迹 = %.12f\n", N, real(tr(FZ)))
        for ss in (1, -1), ord in (:inc, :dec)
            z, L = dense_ring(Ws; step_sign=ss, order=ord)
            @printf("  ss=%+d %-3s : Z_dense = %14.10f   err = %9.2e\n",
                    ss, ord, real(z), abs(real(z) - real(tr(FZ))))
        end
    end
    println()
    println("="^78)
    println("S2：εd≠0 vs Fock 环链（Z + n̂ 插入）")
    println("="^78)
    for N in (2, 3, 4)
        cells = fill(BETA / N, N); durs = fill(BETA / N, N)
        Ws, n1 = interval_operators(cells, durs, EPS_D)
        Wsn = [p == 1 ? n1 * Ws[1] : Ws[p] for p in 1:N]     # n̂ 插在缝（链首前）
        FZ, FZn, _ = fock_ring(cells, durs, EPS_D)
        nref = real(tr(FZn) / tr(FZ))
        ss, ord = 1, :inc
        z, _ = dense_ring(Ws; step_sign=ss, order=ord)
        zn, _ = dense_ring(Wsn; step_sign=ss, order=ord)
        @printf("N=%d : Z = %12.8f（Fock %12.8f, err %9.2e）  n = %10.8f（Fock %10.8f, err %9.2e）\n",
                N, real(z), real(tr(FZ)), abs(real(z) - real(tr(FZ))),
                real(zn / z), nref, abs(real(zn / z) - nref))
    end
    println()
    println("="^78)
    println("S4：非均匀调度（Trotter orders 2/4，native_schedule + 环门合并）vs Fock 链")
    println("="^78)
    for order in (2, 4), N in (order == 4 ? (1,) : (1, 2))
            sched = NativeExact.native_schedule(order, N, BETA / N)
            links = copy(sched.links)
            ring_links = vcat(links[1] + links[end], links[2:end-1])   # 环门合并
            cells = sched.cells
            length(cells) == length(ring_links) || error("调度表长度不配")
            Ws, n1 = interval_operators(cells, ring_links, EPS_D)
            FZ, _, _ = fock_ring(cells, ring_links, EPS_D)
            z, _ = dense_ring(Ws)
            @printf("N=%d order=%d (N_int=%d) : Z = %12.8f（Fock %12.8f, err %9.2e）\n",
                    N, order, length(cells), real(z), real(tr(FZ)),
                    abs(real(z) - real(tr(FZ))))
    end
    println()
    println("="^78)
    println("S3：连续极限（n_dense vs ED 热态 ⟨n⟩ = ", ed_n(EPS_D), "）")
    println("="^78)
    for N in (4, 5)                       # L = 4N ≤ 20（内存保护，勿升）
        cells = fill(BETA / N, N); durs = fill(BETA / N, N)
        Ws, n1 = interval_operators(cells, durs, EPS_D)
        Wsn = [p == 1 ? n1 * Ws[1] : Ws[p] for p in 1:N]
        z, _ = dense_ring(Ws)
        zn, _ = dense_ring(Wsn)
        @printf("N=%2d : n_dense = %12.10f   |Δ−ED| = %9.2e\n",
                N, real(zn / z), abs(real(zn / z) - ed_n(EPS_D)))
    end
end
main()
