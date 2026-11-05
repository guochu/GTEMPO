# =============================================================================
# src/FormB.jl —— 形态 B（消元式单扇区）的 GTEMPO GrassmannMPS 装配
#
# 把 ../form_b/formb_dense.jl 稠密验证过的消元装配移植到 GTEMPO 的 GMPS
# 机制上（稠密 → 稀疏因子，为将来 GTEMPO 库内集成作准备）：
#
#   P_elim = ν·Pf(G) · exp(Σ_{λ<λ'} D[λ,λ'] θ_λθ_λ') · K[x]
#   D[λ,λ'] = −½(G⁻¹[i,k] − G⁻¹[k,i])    （λ,λ' = 接口轨迹腿，i,k = 其 IF 槽位）
#
# lattice = IMP-only 紧致格（timelocallattice 变量序滤掉 :IF，相邻 pair
# 结构保持；GTEMPO 按相邻 site 配 pair 并施加隐式 resolution g_trace）：
#
#   · layer 0 / 开放对：隐式 resolution 恰为 dense 里的 M₀ / 迹闭合
#     （反周期闭合）——不放任何显式因子；
#   · 接口对 (x_t, x̄_t)：该 resolution 已被 Δ_eff 高斯消元吸收，GTEMPO
#     会重复施加 —— 以幂零反测度 (1 + c·x_t·x̄_t) 补偿（(1−x̄x)(1+x̄x)=1
#     精确抵消；c 的符号经 Z=1 证书标定，默认 +1）。
#
# 因子分组 (Δ_eff+反测度, ρ₀+前向链, 后向链) —— 三个非平凡 GMPS，与
# IFGTEMPO 的 (I1, I2, K) 三因子收缩路径完全同构：Z 用 l_LL ×
# GrassmannTransferMatrix × scalar，ρ(T) 直接复用 gtempo_open_kernel_direct
# / density_matrix_direct。与 IFGTEMPO 对齐的数据源：G/ν 来自同一
# influence_matrix_G(model, sched.cells)，links/cells 来自同一
# trotter_schedule(order, N, dt)。
#
# IFGTEMPO 依赖仅限 impurity_kernel_gmps / PathLattice 与收缩函数（纯
# GTEMPO 底层；库内集成时随 pair-g_trace 语义一并移入 GTEMPO）。
# =============================================================================
module FormB

using GTEMPO
using GTEMPO: TK
using LinearAlgebra: isdiag, Diagonal
using TensorOperations: @tensor
using IFGTEMPO: trotter_schedule, influence_matrix_G, impurity_hamiltonian,
                PathLattice, GrassmannVariable, impurity_kernel_gmps, index,
                nvariables, gtempo_open_kernel_direct, gtempo_density_matrix_direct

const GV = GrassmannVariable
IMP(t, br, stat, f) = GV(:imp, t, br, stat, f)
const GMPS = GTEMPO.GrassmannMPS

"""IF 变量的 G 矩阵槽位（1-based）：flat = 4(t−1)+slot+1，slot 0..3 =
(+bar, +eta, −bar, −eta)——influence_matrix_G 的约定。"""
gslot(t::Int, br::Symbol, stat::Symbol) =
    4 * (t - 1) + (br === :+ ? 1 : 3) + (stat === :bar ? 0 : 1)

if_var(f::Int, i::Int) =
    GV(:IF, div(i - 1, 4) + 1, mod(i - 1, 4) < 2 ? :+ : :-,
        iseven(mod(i - 1, 4)) ? :bar : :eta, f)

"""IF slot 的测度 partner（IMP 腿，bar↔eta 交叉）。"""
partner_var(t::Int, br::Symbol, stat::Symbol, f::Int) =
    IMP(t, br, stat === :bar ? :eta : :bar, f)

# ---------------------------------------------------------------------------
# lattice
# ---------------------------------------------------------------------------

"""IMP-only 紧致 lattice：timelocallattice 的变量序滤掉 :IF。

相邻 pair 结构保持：layer 0 的 (x̄₀±, x₀±)、接口的 (x_t±, x̄_t±)、开放的
(x̄⁺_{N+1}, x⁻_{N+1})——与 GTEMPO 的逐 pair 隐式 resolution 一一对应。"""
function formb_lattice(N::Int, nflav::Int)
    vars = GV[]
    for f in 1:nflav
        push!(vars, IMP(0, :+, :bar, f)); push!(vars, IMP(0, :+, :eta, f))
    end
    for f in 1:nflav
        push!(vars, IMP(0, :-, :bar, f)); push!(vars, IMP(0, :-, :eta, f))
    end
    for t in 1:N, f in 1:nflav
        push!(vars, IMP(t, :+, :eta, f)); push!(vars, IMP(t, :+, :bar, f))
        push!(vars, IMP(t, :-, :eta, f)); push!(vars, IMP(t, :-, :bar, f))
    end
    for f in 1:nflav
        push!(vars, IMP(N + 1, :+, :bar, f)); push!(vars, IMP(N + 1, :-, :eta, f))
    end
    position = Dict(v => i - 1 for (i, v) in enumerate(vars))   # 0-based
    return PathLattice(vars, position, N), vars
end

# ---------------------------------------------------------------------------
# 因子构造
# ---------------------------------------------------------------------------

_vacuum(n::Int) = GMPS(ComplexF64, n)

"""kernel 因子。两 flavor：转发 IFGTEMPO.impurity_kernel_gmps（分级构造）。

单 flavor：对角 2×2 kernel（踢门/对角初态）是纯双线性
1 + K[2,2]·x̄x —— 以 partialmpo 直构（递增序 (eta, bar) 上系数
−K[2,2]，local (bar,eta)→全局 (eta,bar) 的排列符号已并入）。"""
function kernel_gmps(lattice, bars::Tuple, etas::Tuple, K::AbstractMatrix;
                     trunc=GTEMPO.notrunc())
    n = nvariables(lattice)
    if length(bars) == 1
        isdiag(K) || error("单 flavor kernel 仅支持对角矩阵（得 1-mode 分级构造）")
        bar = bars[1]; eta = etas[1]
        sb, se = index(lattice, bar), index(lattice, eta)
        sb != se || error("kernel 两个 site 重合")
        # local 序 (bar, eta) → 全局递增序的排列符号
        sgn = sb < se ? 1 : -1
        c = sgn * K[2, 2]
        iszero(c) && return _vacuum(n)
        lo, hi = minmax(sb, se)
        return GTEMPO.partialmpo(lo + 1, [hi + 1], [ComplexF64(c)]) * _vacuum(n)
    elseif length(bars) == 2
        return impurity_kernel_gmps(lattice, bars, etas, K)
    else
        error("不支持 $nflav flavor")
    end
end

"""acc × factor（GTEMPO mult，返回 acc）。"""
function _mul(acc, fac; trunc=GTEMPO.notrunc())
    acc, _ = GTEMPO.mult(acc, fac; trunc = trunc)
    return acc
end

"""把标量 c 折进 GMPS scaling（GTEMPO 收缩按 scaling^nsites 计）。"""
function _fold_scalar!(mps, c::Real)
    GTEMPO.setscaling!(mps, GTEMPO.scaling(mps) * c^(1 / length(mps)))
    return mps
end

"""数值 Pfaffian（Parlett–Reid：斜对称矩阵的 Gauss 消元到三对角，
Pf = 对角 pivots 之积；O(n³)，无 2^n 展开——稠密展开法在 nl 大时是
2^{4nl} 的内存炸弹（N=8 截断档 = 54GB），勿回退）。"""
function pfaffian(A0::Matrix{ComplexF64})
    n = size(A0, 1)
    (n % 2 == 0) || return zero(ComplexF64)
    A = copy(A0)
    pf = one(ComplexF64)
    for k in 1:2:n
        kp = k + 1
        kp > n && break
        # pivot：|A[k+1,k]| 过小则与下方某行同时交换行列（congruence，Pf 不变）
        if abs(A[kp, k]) < 1e-14
            r = 0
            for i in (kp + 1):n
                abs(A[i, k]) > 1e-14 && (r = i; break)
            end
            r == 0 && return zero(ComplexF64)   # 奇异 → Pf = 0
            A[[kp, r], :] = A[[r, kp], :]
            A[:, [kp, r]] = A[:, [r, kp]]
        end
        piv = A[kp, k]
        for i in (kp + 1):n
            c = A[i, k] / piv
            c == 0 && continue
            for j in (kp + 1):n
                A[i, j] -= c * A[kp, j]
            end
            for j in (kp + 1):n
                A[j, i] -= c * A[j, kp]
            end
        end
        pf *= piv
    end
    return pf
end

"""ν·Pf(G)：归一化 × Gaussian 归一化（解析 Pfaffian）。

校验基准：与完美匹配显式求和及 formb_dense 的行因子 Berezin 积分
三方一致（≤1e-14）。"""
function bath_pf_norm(data, nl::Int, f::Int)
    return data.normalization * pfaffian(data.G)
end

"""Δ_eff 因子 + 接口反测度，折进同一个 GMPS。

D[λ,λ'] = elim_sign·(G⁻¹[i,k] − G⁻¹[k,i])/2，行因子化
（1 + θ_head·Σ D θ_partner，get_if_gmps 逐行 partialmpo 同款，精确），
随后逐接口对乘反测度 (1 + counter_c·x_t·x̄_t)。`norm` = ν·Pf(G) 折进
scaling。"""
function delta_eff_gmps(lattice, Ginv::AbstractMatrix, nl::Int, f::Int;
                        elim_sign::Int=-1, counter_c::Int=1, norm::Real=1.0,
                        trunc=GTEMPO.notrunc())
    n = nvariables(lattice)
    GA = (Ginv .- transpose(Ginv)) ./ 2
    legs = Tuple{Int,Int}[]   # (partner site, G 槽位)
    for t in 1:nl, br in (:+, :-), stat in (:bar, :eta)
        push!(legs, (index(lattice, partner_var(t, br, stat, f)),
                     gslot(t, br, stat)))
    end
    sort!(legs, by = x -> x[1])   # site 递增 → monomial 序 / 行分组
    acc = _vacuum(n)
    _fold_scalar!(acc, norm)
    seen = Int[]
    for (site, slot) in legs
        site in seen && continue
        push!(seen, site)
        cols = Int[]; coefs = ComplexF64[]
        for (site2, slot2) in legs
            site2 > site || continue
            v = elim_sign * GA[slot, slot2]
            iszero(v) && continue
            push!(cols, site2 + 1)
            push!(coefs, v)
        end
        isempty(coefs) && continue
        group = GTEMPO.partialmpo(site + 1, cols, coefs) * _vacuum(n)
        acc = _mul(acc, group; trunc = trunc)
    end
    # 接口反测度 (1 + counter_c·x_t·x̄_t)：紧凑格接口对 = (IMP eta, IMP bar)
    for t in 1:nl, br in (:+, :-)
        s_eta = index(lattice, IMP(t, br, :eta, f)) + 1
        s_bar = index(lattice, IMP(t, br, :bar, f)) + 1
        group = GTEMPO.partialmpo(s_eta, [s_bar], [ComplexF64(counter_c)]) *
                _vacuum(n)
        acc = _mul(acc, group; trunc = trunc)
    end
    return acc
end

"""K 因子：初态 ρ₀ + 前向/后向链 kernel（get_k_gmps 的紧致格版）。

返回 (x_init_fwd, x_bwd)：初态+前向链折为一个因子，后向链为另一个
（三因子收缩需要三个非平凡 GMPS；乘积与整链相同）。单 flavor：kernel
矩阵 2×2；两 flavor：4×4（mode 序 f=1,2）。"""
function k_gmps_split(lattice, nl::Int, nflav::Int, links::Vector{Float64},
                      Him::AbstractMatrix, rho0::AbstractMatrix;
                      trunc=GTEMPO.notrunc())
    bars(t, br) = ntuple(f -> IMP(t, br, :bar, f), nflav)
    etas(t, br) = ntuple(f -> IMP(t, br, :eta, f), nflav)
    x_if = _vacuum(nvariables(lattice))
    x_if = _mul(x_if, kernel_gmps(lattice, bars(0, :+), etas(0, :-), rho0);
                trunc)
    for t in 1:(nl + 1)
        x_if = _mul(x_if, kernel_gmps(lattice, bars(t, :+), etas(t - 1, :+),
                                      exp((-1im * links[t]) .* Him)); trunc)
    end
    x_bwd = _vacuum(nvariables(lattice))
    for t in (nl + 1):-1:1
        x_bwd = _mul(x_bwd, kernel_gmps(lattice, bars(t - 1, :-),
                                        etas(t, :-),
                                        exp((+1im * links[t]) .* Him)); trunc)
    end
    return x_if, x_bwd
end

# ---------------------------------------------------------------------------
# 路径组装与收缩
# ---------------------------------------------------------------------------

"""组装形态 B 路径：返回 (Δ_eff+反测度, ρ₀+前向链, 后向链, lattice, sched)。

`models`：单 flavor 传 (bath,)；两 flavor 传 (bath_up, bath_dn)。
`imp`：单 flavor NamedTuple(eps=…, rho0=…)；两 flavor TwoSpinImpurityModel。
G/ν 取自 influence_matrix_G(model, sched.cells)（与 IFGTEMPO 的
build_impurity_path 同一数据源）。"""
function formb_path(models::Tuple, imp, N::Int, dt::Real; order::Integer=2,
                    counter_c::Int=1, trunc=GTEMPO.notrunc())
    sched = trotter_schedule(order, N, dt)
    nflav = length(models)
    nl = length(sched.cells)              # 接口层数（order 4 = 3N）
    lattice, vars = formb_lattice(nl, nflav)
    rho0 = nflav == 1 ? imp.rho0 : imp.initial_density
    Him = nflav == 1 ? ComplexF64[0 0; 0 imp.eps] : impurity_hamiltonian(imp)

    x_delta = _vacuum(nvariables(lattice))
    for f in 1:nflav
        data = influence_matrix_G(models[f], sched.cells)
        Ginv = inv(data.G)
        pf = bath_pf_norm(data, nl, f)    # ν·Pf(G)
        x_delta = _mul(x_delta, delta_eff_gmps(lattice, Ginv, nl, f;
                                               elim_sign = -1,
                                               counter_c = counter_c,
                                               norm = real(pf), trunc = trunc);
                       trunc)
    end
    x_if, x_bwd = k_gmps_split(lattice, nl, nflav, sched.links, Him, rho0;
                               trunc)
    # 因子排列 (Δ, bwd, ρ₀+前向)：graded 读出符号经 Z=1 / ρ=Fock 双证书
    # 标定（bwd 置中；与 timelocallattice 的 (I1,I2,K) 排列的偶宇称扇区
    # 符号差由 formb_density_matrix 的 S 承载）
    return (factors = (x_delta, x_bwd, x_if), lattice = lattice,
            sched = sched,
            diagnostics = (nvariables = nvariables(lattice),
                           max_bond = maximum(GTEMPO.bonddim.(
                               (x_delta, x_bwd, x_if)))))
end

"""全收缩 Z（任意因子数，GTEMPO pair transfer + scalar）。"""
function formb_path_integral(factors::GMPS...)
    tmp = GTEMPO.l_LL(factors...) * GTEMPO.GrassmannTransferMatrix(factors...)
    return TK.scalar(tmp)
end

"""单 flavor 开放核：前 L−2 个 site 逐 pair transfer，最后 2 个 site
（开放腿 x̄⁺_{N+1}, x⁻_{N+1}）逐一吸收——gtempo_open_kernel_direct 的
单 flavor 版（后者硬编码 4 条开放腿）。逐 site 收缩模式与
Contraction.jl 逐块一致：(X1, X2, X3) = (Δ, ρ₀+前向, 后向)。"""
function open_kernel_single(x1::GMPS, x2::GMPS, x3::GMPS)
    L = length(x1)
    (length(x2) == L && length(x3) == L) ||
        throw(DimensionMismatch("因子长度不一致"))
    iseven(L) || throw(ArgumentError("pair 收缩要求偶数 site"))
    env = GTEMPO.l_LL(x1, x2, x3)
    for pair in 1:div(L - 2, 2)
        env = env * GTEMPO.GrassmannTransferMatrix(pair, x1, x2, x3)
    end
    # ---- site L−1（次末 = x̄⁺_{N+1}）----
    Xt, It, Kt = x1[L - 1], x2[L - 1], x3[L - 1]
    @tensor a1[v, b1, pk, rk; b2] := env[v, b1, b2, b3] * Kt[b3, pk, rk]
    @tensor a2[v, pk, pi, rk, ri; b1] := a1[v, b1, pk, rk, b2] * It[b2, pi, ri]
    a3 = GTEMPO.g_fuse(a2, 2)
    @tensor a4[v, pk, px, rk, ri; rx] := a3[v, pk, rk, ri, b1] * Xt[b1, px, rx]
    a5 = GTEMPO.g_fuse(a4, 2)
    env = TK.permute(a5, (1, 2), (5, 4, 3))     # cod(v, P1), dom(rx, ri, rk)
    # ---- site L（末 = x⁻_{N+1}）----
    Xt, It, Kt = x1[L], x2[L], x3[L]
    @tensor b1[v, P1, b1, pk, rk; b2] := env[v, P1, b1, b2, b3] * Kt[b3, pk, rk]
    @tensor b2[v, P1, pk, pi, rk, ri; b1] := b1[v, P1, b1, pk, rk, b2] *
                                             It[b2, pi, ri]
    b3 = GTEMPO.g_fuse(b2, 3)
    @tensor b4[v, P1, Pki, px, rk, ri; rx] := b3[v, P1, Pki, rk, ri, b1] *
                                              Xt[b1, px, rx]
    b5 = GTEMPO.g_fuse(b4, 3)
    env = TK.permute(b5, (1, 2, 3), (6, 5, 4))  # cod(v, P1, P2), dom(rx, ri, rk)
    sc = GTEMPO.scaling(x1) * GTEMPO.scaling(x2) * GTEMPO.scaling(x3)
    return TK.lmul!(sc^2, env)
end

"""开放核 ρ(T)——两 flavor（复用 IFGTEMPO 的收缩与系数提取）。

紧凑格 + 因子排列 (Δ, bwd, ρ₀+前向) 的 graded 读出符号 S：偶宇称扇区
（对角初态类）= diag(1, 1, 1, −1)，即仅 |↑↓⟩ 对角扇区差一个费米子符号
（与单 flavor 的 (−1)^{bra·ket} 同性质），由 Z=1 / ρ=Fock 双证书标定，
并在全部 order/N/模型上一致。奇宇称扇区是 impurity_kernel_gmps 电荷零
分级构造的固有限制（IFGTEMPO 测试均为对角初态），不在支持类内。"""
const TWO_SPIN_EXTRACT_SIGN = Diagonal(ComplexF64[1, 1, 1, -1])

function formb_density_matrix(factors::GMPS...)
    return TWO_SPIN_EXTRACT_SIGN .* gtempo_density_matrix_direct(factors...)
end

"""开放核 ρ(T)——单 flavor（开放腿 2 条，4 系数直接映射 2×2 矩阵）。

env 的 codomain 腿 = (aux, P1, P2)（P1 = x̄⁺_{N+1}，P2 = x⁻_{N+1}）。
开放 monomial (bra on x̄⁺, ket on x⁻) ↔ ρ[bra+1,ket+1]——在
(Δ, bwd, ρ₀+前向) 因子排列下直接映射（graded 符号恰好为 +；该排列由
Z=1 / ρ=Fock 双证书标定）。"""
function formb_density_matrix_single(x1::GMPS, x2::GMPS, x3::GMPS)
    env = open_kernel_single(x1, x2, x3)
    rho = zeros(ComplexF64, 2, 2)
    for (f1, f2) in TK.fusiontrees(env)
        u = f1.uncoupled
        rho[u[2].n + 1, u[3].n + 1] = env[f1, f2][1, 1]
    end
    return rho
end

end # module FormB
