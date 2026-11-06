# =============================================================================
# form_b/native_probe.jl —— Q2 判决实验：native 腿集上的两扇区消元 × 原生骨架
#
# 目标：把 FormB 的消元式构造搬到 GTEMPO 原生实时骨架（每节点 4 腿、ρ₀ 在
# node 1、折返在开放端 node k、门核链连接相邻节点）上：
#
#   两扇区（native 腿集）  P(a,y) = ν·exp(½aᵀGa) · ∏_c M_c(a_c, y_c) · S_native[y]
#   --∫Da-->               P_elim(y) = ν·Pf(G) · exp(D) · S_native[y]
#
# 其中 G = 两扇区 bath Gaussian（bath pair per (branch, node) 布局，QuAPI
# exact-Δt 同构），D[λ,λ'] = −½(G⁻¹[i,k] − G⁻¹[k,i])（λ,λ' = 接口轨迹腿，
# i,k = 其 bar↔eta 交叉 partner 槽位），S_native = 门核链 + ρ₀ + 折返 +
# （旋钮）逐节点分辨率。
#
# lattice 结构（对齐 GTEMPO RealGrassmannLattice + FormB schedule）：
#   Nt = nl+1 切片，k = Nt+1 = nl+2 个节点；node 1 = t=0（ρ₀ 层）、
#   node 2..nl+1 = 接口 1..nl（IF 槽位/测度/D 源所在）、node k = 开放端
#   （物理 T：折返 + ρ(T) 提取）。gap g（node g → g+1）承载门 A(links[g])，
#   cell c 连接 node c → c+1（与 FormB 链 A(l₁)·Q₁·A(l₂)···Q_nl·A(l_{nl+1})
#   逐因子一致）。每 node 腿序 (bar⁺, eta⁺, bar⁻, eta⁻)。
#   bath pair：per (mode, branch, node 1..nl+1)——node 1 挂热边、node nl+1 挂
#   迹边（物理 T 端的 +/− pair），开放端 node k 无 bath pair（无 cell 耦合）。
#
# 判决协议：Z（ε_d=0，Z·ν·Pf == 1）锁定测度/闭合旋钮；ρ（ε_d≠0 vs Fock 链）
# 锁定 D/源结构。约定旋钮：trace_style(:edge/:flip)、fold(:none/:neg/:pos)、
# res(:none/:all)、conv(step_sign × order)。
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/form_b/native_probe.jl
# =============================================================================
using Printf
using LinearAlgebra
using GTEMPO

const IFT_ROOT = joinpath(dirname(dirname(dirname(pathof(GTEMPO)))), "IFGTEMPO")
isdir(IFT_ROOT) || error("IFGTEMPO 仓库未找到：$IFT_ROOT")
push!(LOAD_PATH, IFT_ROOT)
using IFGTEMPO

include(joinpath(@__DIR__, "formb_dense.jl"))   # 稠密机器 + Fock 链参照

# ---- 本地工具（IFGTEMPO 同名函数的移植；生产模块将直接使用这些）----

"""反对称矩阵占位（1-based 槽位）：上三角 += v、下三角 −= v（i≠j）。"""
function put_G!(A::Matrix{ComplexF64}, i::Int, j::Int, value)
    i == j && return A
    if i < j
        A[i, j] += value
        A[j, i] -= value
    else
        A[j, i] -= value
        A[i, j] += value
    end
    return A
end

"""占据态模式列表。"""
_occ_modes(state::Int, m::Int) = [mode for mode in 0:(m - 1) if state & (1 << mode) != 0]

"""单粒子矩阵 S 的二次量子化（公开 Fock 基）。"""
function second_quantized(S::AbstractMatrix)
    m = size(S, 1)
    dim = 1 << m
    out = zeros(ComplexF64, dim, dim)
    for bra in 0:(dim - 1), ket in 0:(dim - 1)
        bm, km = _occ_modes(bra, m), _occ_modes(ket, m)
        length(bm) == length(km) || continue
        out[bra + 1, ket + 1] = isempty(bm) ? 1.0 + 0.0im :
                                det(S[bm .+ 1, km .+ 1])
    end
    return out
end

# ---------------------------------------------------------------------------
# native 腿集（每 flavor 一套；site 0-based）
# ---------------------------------------------------------------------------

struct Nat
    nl::Int          # 接口层数（= length(cells)）
    k::Int           # 节点数 = nl+2（Nt = nl+1）
    L::Int           # 4k（单 flavor）
end

"""node p（1..k）的腿 (br, stat) → 0-based site。腿序 per node: (bar⁺, eta⁺, bar⁻, eta⁻)。"""
function nsite(n::Nat, p::Int, br::Symbol, stat::Symbol)
    1 <= p <= n.k || throw(BoundsError(1:n.k, p))
    sub = (br === :+ ? 0 : 2) + (stat === :bar ? 1 : 2)
    return 4 * (p - 1) + sub - 1
end

"""接口 c（1..nl，= node c+1）的 IF 槽位（G 矩阵 1-based，接口局部索引）。"""
ifslot(n::Nat, c::Int, br::Symbol, stat::Symbol) =
    4 * (c - 1) + (br === :+ ? 1 : 3) + (stat === :bar ? 0 : 1)

"""接口 c 的 partner 槽位（测度 partner：bar↔eta 交叉）。"""
partner_slot(n::Nat, c::Int, br::Symbol, stat::Symbol) =
    ifslot(n, c, br, stat === :bar ? :eta : :bar)

"""bath pair（mode k, branch br, node p）的两个槽位（bar, eta），1-based。

bath pair 只存在于 node 1..nl+1（= k−1）：node 1 挂热边（t=0）、node nl+1
挂迹边（物理 T 端）；开放端 node k 无 cell 耦合、无 bath pair。"""
function bslot(n::Nat, m::Int, k::Int, br::Symbol, p::Int)
    nbp = n.k - 1
    base = 4 * n.nl + (k - 1) * 4 * nbp + (br === :+ ? 0 : 2 * nbp) + (p - 1) * 2
    return base + 1, base + 2
end

# ---------------------------------------------------------------------------
# 两扇区 G 装配（bath pair per (branch, node)，QuAPI exact-Δt 同构）
# ---------------------------------------------------------------------------

"""_add_mapped_q_matrix! 的移植（compact 版同配方）：Q kernel 把 (bar_imp,
eta_imp) IF 对与 bath 输出 bars/输入 etas 经 U 缝合，trace_bra 翻转输出行。"""
function q_add!(A::Matrix{ComplexF64}, bar_imp::Int, eta_imp::Int,
                bath_bars::Vector{Int}, bath_etas::Vector{Int}, U::AbstractMatrix;
                trace_bra::Bool=false)
    m = length(bath_bars)
    size(U) == (m + 1, m + 1) || throw(DimensionMismatch("Q 传播子维度不匹配"))
    put_G!(A, bar_imp, eta_imp, U[1, 1])
    for j in eachindex(bath_etas)
        put_G!(A, bar_imp, bath_etas[j], U[1, j + 1])
    end
    for i in eachindex(bath_bars)
        value = U[i + 1, 1]
        trace_bra && (value = -value)
        put_G!(A, bath_bars[i], eta_imp, value)
        for j in eachindex(bath_etas)
            bath_value = U[i + 1, j + 1]
            trace_bra && (bath_value = -bath_value)
            put_G!(A, bath_bars[i], bath_etas[j], bath_value)
        end
    end
    return A
end

"""耦合单粒子 h（mode 1 = 杂质，能级 0；mode k+1 = bath mode k）。"""
function coupled_h(V::Vector{ComplexF64}, Hbath::Matrix{ComplexF64})
    m = length(V)
    h = zeros(ComplexF64, m + 1, m + 1)
    for k in 1:m
        h[1, k+1] = h[k+1, 1] = V[k]
        h[k+1, k+1:k+1] .= Hbath[k, k]
    end
    return h
end

"""native 布局的 bath Gaussian G（接口 IF 槽位上，4·nl 维）。

cell c 连接 node c → c+1，锚在接口 c（= node c+1）的 (br) 对——与 compact
的 cell 端接口锚同构；bath pair 链 (+,1)→(+,2)→…→(+,nl+1)（forward）、
(−,nl+1)→(−,nl)→…→(−,1)（backward）；热边 (+,1)/(−,1)；迹边在 (·,nl+1)。"""
function native_G(Hbath::Matrix{ComplexF64}, V::Vector{ComplexF64}, β::Real,
                  cells::Vector{Float64}, n::Nat; trace_style::Symbol=:edge)
    m = length(V)
    nif = 4 * n.nl
    nbp = n.k - 1                      # bath pair 的 node 数（1..nl+1）
    nb = 4 * m * nbp
    A = zeros(ComplexF64, nif + nb, nif + nb)
    N = length(cells)
    (N == n.nl) || throw(DimensionMismatch("cells 数与层数不一致"))
    h = coupled_h(V, Hbath)
    Us = [exp(-1im * cell .* h) for cell in cells]
    bars(br, p) = [bslot(n, m, kk, br, p)[1] for kk in 1:m]
    etas(br, p) = [bslot(n, m, kk, br, p)[2] for kk in 1:m]
    # forward cells：输入 pair (+,c) etas → 输出 pair (+,c+1) bars，锚接口 c
    for c in 1:N
        flip = (trace_style === :flip && c == N)
        q_add!(A, ifslot(n, c, :+, :bar), ifslot(n, c, :+, :eta),
               bars(:+, c + 1), etas(:+, c), Us[c]; trace_bra = flip)
    end
    # backward cells：输入 pair (−,c+1) etas → 输出 pair (−,c) bars，锚接口 c
    for c in 1:N
        q_add!(A, ifslot(n, c, :-, :bar), ifslot(n, c, :-, :eta),
               bars(:-, c), etas(:-, c + 1), adjoint(Us[c]); trace_bra = false)
    end
    # trace 边（旋钮 :edge）：物理 T 端（node nl+1）的 − pair bar × + pair eta
    if trace_style === :edge
        for kk in 1:m
            put_G!(A, bslot(n, m, kk, :-, nbp)[1], bslot(n, m, kk, :+, nbp)[2], -1.0)
        end
    elseif trace_style !== :flip
        error("trace_style 必须是 :edge 或 :flip")
    end
    # thermal 边：(+,1) bars × (−,1) etas
    thermal = exp(-Float64(β) .* Hbath)
    for i in 1:m, j in 1:m
        put_G!(A, bslot(n, m, i, :+, 1)[1], bslot(n, m, j, :-, 1)[2], thermal[i, j])
    end
    # bath pair 测度 −γ̄γ（compact 版逐 pair 的 −1 对角同款）
    for kk in 1:m, br in (:+, :-), p in 1:nbp
        put_G!(A, bslot(n, m, kk, br, p)[1], bslot(n, m, kk, br, p)[2], -1.0)
    end
    E = A[1:nif, 1:nif]
    C = A[1:nif, nif+1:end]
    B = A[nif+1:end, nif+1:end]
    effective = E + C * (B \ transpose(C))
    G = 0.5 .* (effective .- transpose(effective))
    return G
end

"""归一化 ν = ZN/Z（bath 侧对象，与 compact 相同的公式）。"""
function bath_norm(V::Vector{ComplexF64}, Hbath::Matrix{ComplexF64}, β::Real,
                   cells::Vector{Float64})
    h = coupled_h(V, Hbath)
    Bprod = Matrix{ComplexF64}(I, length(V), length(V))
    for cell in cells
        Bprod = exp(-1im * cell .* h)[2:end, 2:end] * Bprod
    end
    R = exp(-Float64(β) .* Hbath)
    bw = Bprod * R * adjoint(Bprod)
    return tr(second_quantized(bw)) / tr(second_quantized(R))
end

"""消元核 D（接口腿上的 (bar,eta) 有向对），返回 pairs [(site_i, site_j, v)]。

接口 c 的轨迹腿 = node c+1；partner 槽位 = 接口 c 的 bar↔eta 交叉。"""
function native_D_pairs(n::Nat, G::AbstractMatrix; elim_sign::Int=-1)
    Ginv = inv(G)
    GA = (Ginv .- transpose(Ginv)) ./ 2
    legs = Tuple{Int,Int}[]   # (site, partner IF 槽位)
    for c in 1:n.nl, br in (:+, :-), stat in (:bar, :eta)
        push!(legs, (nsite(n, c + 1, br, stat), partner_slot(n, c, br, stat)))
    end
    sort!(legs, by = x -> x[1])
    pairs = Tuple{Int,Int,ComplexF64}[]
    for i in eachindex(legs), k in (i + 1):length(legs)
        v = elim_sign * GA[legs[i][2], legs[k][2]]
        iszero(v) && continue
        push!(pairs, (legs[i][1], legs[k][1], v))
    end
    return pairs
end

"""定向双线性因子：意图单项式 θ_row·θ_col（算符序），coeff；转为 site 序存储。"""
bilinear_pair(row::Int, col::Int, coeff::ComplexF64) =
    row < col ? (row, col, coeff) : (col, row, -coeff)

"""双线性对 → 单项式表项 ([i,j], c)。"""
bilinear_mono(pair::Tuple{Int,Int,ComplexF64}) = ([pair[1], pair[2]], pair[3])

"""native 骨架因子（单 flavor）：
  门核链（gap g=1..Nt，时长 links[g]）、ρ₀(node 1)、fold(node k，旋钮)、
  分辨率（:all 时逐 node×branch exp(−x̄x)；:none 时依赖全局积分定向）。"""
function native_skeleton_factors(n::Nat, links::Vector{Float64},
                                 Him::AbstractMatrix, rho0::AbstractMatrix;
                                 fold::Symbol=:neg, res::Symbol=:none)
    factors = Monomials[]
    # 门核链：forward ⟨x̄(g+1)⁺|A|x(g)⁺⟩，backward ⟨x̄(g)⁻|A†|x(g+1)⁻⟩
    for g in 1:(n.k - 1)
        push!(factors, kernel_monomials([nsite(n, g + 1, :+, :bar)],
                                        [nsite(n, g, :+, :eta)],
                                        exp((-1im * links[g]) .* Him)))
        push!(factors, kernel_monomials([nsite(n, g, :-, :bar)],
                                        [nsite(n, g + 1, :-, :eta)],
                                        exp((+1im * links[g]) .* Him)))
    end
    # bulkconnection 相邻重叠链（bare 核的测度补件，coeff +1，同 bulkconnection!）：
    # forward exp(x̄⁺(g+1)·x⁺(g))，backward exp(x̄⁻(g)·x⁻(g+1))
    for g in 1:(n.k - 1)
        push!(factors, Monomials([([], 1.0 + 0im),
                                  bilinear_mono(bilinear_pair(nsite(n, g + 1, :+, :bar),
                                                              nsite(n, g, :+, :eta),
                                                              1.0 + 0im))]))
        push!(factors, Monomials([([], 1.0 + 0im),
                                  bilinear_mono(bilinear_pair(nsite(n, g, :-, :bar),
                                                              nsite(n, g + 1, :-, :eta),
                                                              1.0 + 0im))]))
    end
    # ρ₀(node 1)：⟨x̄⁺(0)|ρ₀|x⁻(0)⟩
    push!(factors, kernel_monomials([nsite(n, 1, :+, :bar)],
                                    [nsite(n, 1, :-, :eta)], rho0))
    # fold（旋钮）：exp(∓x̄⁻(k)·x⁺(k)) / 无
    if fold !== :none
        c = (fold === :neg) ? -1.0 + 0im : 1.0 + 0im
        push!(factors, Monomials([([], 1.0 + 0im),
                                  bilinear_mono(bilinear_pair(nsite(n, n.k, :-, :bar),
                                                              nsite(n, n.k, :+, :eta),
                                                              c))]))
    end
    # 分辨率（旋钮）：
    #   :boundary —— 仅 node 1（接口分辨率已被交叉测度吸收，与 FormB 的
    #                 "接口无测度"约定一致；node 1 的 ρ₀=|0⟩⟨0| 核是常数，
    #                 bar⁺/eta⁻ 腿必须由分辨率消费，否则积分为零）；
    #   :all      —— 全部节点（过度计数接口分辨率，作对照）；
    #   :none     —— 无（诊断用）。
    if res === :boundary || res === :all
        nodes = (res === :boundary) ? [1] : collect(1:n.k)
        for p in nodes, br in (:+, :-)
            push!(factors, Monomials([([], 1.0 + 0im),
                                      bilinear_mono(bilinear_pair(nsite(n, p, br, :bar),
                                                                  nsite(n, p, br, :eta),
                                                                  -1.0 + 0im))]))
        end
    end
    return factors
end

"""node k（开放端）的 T 端分辨率：两对 (bar⁺,eta⁺)/(bar⁻,eta⁻) 的 exp(−x̄x)。"""
function nodek_res_factors(n::Nat)
    factors = Monomials[]
    for br in (:+, :-)
        push!(factors, Monomials([([], 1.0 + 0im),
                                  bilinear_mono(bilinear_pair(nsite(n, n.k, br, :bar),
                                                              nsite(n, n.k, br, :eta),
                                                              -1.0 + 0im))]))
    end
    return factors
end

# ---------------------------------------------------------------------------
# 收缩
# ---------------------------------------------------------------------------

"""开放腿 ρ 提取（open_rho 配方：bra 位 = bar⁺(k)，ket 位 = eta⁻(k)）。"""
function open_rho_native(integrated::Vector{ComplexF64}, n::Nat)
    barsites = [nsite(n, n.k, :+, :bar)]
    etasites = [nsite(n, n.k, :-, :eta)]
    rho = zeros(ComplexF64, 2, 2)
    for bra in 0:1, ket in 0:1
        loc = Int[]
        (bra == 1) && push!(loc, barsites[1])
        (ket == 1) && push!(loc, etasites[1])
        mask = 0
        for s in loc
            mask |= 1 << s
        end
        p = ket
        sgn = (isodd(mono_inversions(loc)) ? -1 : 1) *
              (isodd(p * (p - 1) ÷ 2) ? -1 : 1)
        rho[bra + 1, ket + 1] = sgn * integrated[mask + 1]
    end
    return rho
end

"""native 消元式稠密收缩。返回 (Z, ρ)。ρ 的开放腿 = (bar⁺(k), eta⁻(k))。"""
function native_case(; n::Nat, G::AbstractMatrix, pf::ComplexF64, nu::Number,
                     links::Vector{Float64}, Him::AbstractMatrix,
                     rho0::AbstractMatrix, elim_sign::Int=-1,
                     fold::Symbol=:neg, res::Symbol=:none, conv,
                     nodek_res::Bool=true)
    pairs = native_D_pairs(n, G; elim_sign)
    stack = Monomials[]
    append!(stack, bilinear_row_factors(pairs))
    append!(stack, native_skeleton_factors(n, links, Him, rho0; fold, res))
    L = n.L
    dense = embed(L, stack[1])
    dense .*= nu * pf
    buf = similar(dense)
    for k in 2:length(stack)
        mul_monomials!(buf, dense, L, stack[k])
        dense, buf = buf, dense
    end
    # ρ：保留开放腿
    keep = [nsite(n, n.k, :+, :bar), nsite(n, n.k, :-, :eta)]
    intsites = [s for s in 0:(L - 1) if all(s != kk for kk in keep)]
    poly = integrate_out(dense, L, intsites; step_sign = conv.step_sign,
                         order = conv.order)
    rho = open_rho_native(poly, n)
    # Z：全收缩。node-k 的 T 端分辨率（旋钮）+ 闭合（fold=:none 时补恒等闭合）
    Zpoly = dense
    if nodek_res
        for f in nodek_res_factors(n)
            mul_monomials!(buf, Zpoly, L, f)
            Zpoly, buf = buf, Zpoly
        end
    end
    if fold === :none
        cl = Monomials([([], 1.0 + 0im),
                        bilinear_mono(bilinear_pair(nsite(n, n.k, :+, :bar),
                                                    nsite(n, n.k, :-, :eta),
                                                    1.0 + 0im))])
        mul_monomials!(buf, Zpoly, L, cl)
        Zpoly, buf = buf, Zpoly
    end
    Z = integrate_out(Zpoly, L, collect(0:(L - 1)); step_sign = conv.step_sign,
                      order = conv.order)[1]
    return (Z = Z, rho = rho)
end

"""接口 IF 槽位集上的数值 Pfaffian（给定积分定向）。"""
function pf_mine_native(G::AbstractMatrix; step_sign::Int, order::Symbol)
    nif = size(G, 1)
    pairs = Tuple{Int,Int,ComplexF64}[]
    for i in 1:nif, j in (i + 1):nif
        iszero(G[i, j]) && continue
        push!(pairs, (i - 1, j - 1, G[i, j]))
    end
    factors = bilinear_row_factors(pairs)
    dense = embed(nif, factors[1])
    buf = similar(dense)
    for k in 2:length(factors)
        mul_monomials!(buf, dense, nif, factors[k])
        dense, buf = buf, dense
    end
    return integrate_out(dense, nif, collect(0:(nif - 1));
                         step_sign = step_sign, order = order)[1]
end

# ---------------------------------------------------------------------------
# 主流程：标定网格
# ---------------------------------------------------------------------------

const T_TOT_N = 0.4
const BETA_N = 1.3

function setup(; eps_d::Real)
    up = FiniteBath([0.4;;], [0.23 + 0.11im], BETA_N)
    imp = (eps = eps_d, rho0 = ComplexF64[1 0; 0 0])
    return up, imp
end

function main(; order::Int=1, N::Int=1)
    println("═"^76)
    println("Q2 标定：native 腿集两扇区消元 × 原生骨架（order $order, N=$N, 单自旋）")
    println("═"^76)
    conv0 = Conv(1, :inc, 1)
    sched = trotter_schedule(order, N, T_TOT_N / N)
    n = Nat(length(sched.cells), length(sched.cells) + 2, 4 * (length(sched.cells) + 2))
    up, _ = setup(eps_d = 0.0)
    Hbath = up.bath_hamiltonian
    V = up.hybridization
    cells = collect(Float64.(sched.cells))
    links = collect(Float64.(sched.links))
    nu = bath_norm(V, Hbath, BETA_N, cells)
    @printf("接口 IF 槽 %d 维；ν = %.6f%+.6fim；k = %d 节点\n",
            4 * n.nl, real(nu), imag(nu), n.k)

    # ---- Z 证书（ε_d = 0）：门核全恒等 ----
    _, imp0 = setup(eps_d = 0.0)
    Him0 = zeros(ComplexF64, 2, 2)
    best = Tuple{Float64,NamedTuple}[]
    for trace_style in (:edge, :flip), fold in (:neg, :pos, :none),
        res in (:boundary, :none, :all), nkr in (true, false)
        G = native_G(Hbath, V, BETA_N, cells, n; trace_style)
        pf = pf_mine_native(G; step_sign = conv0.step_sign, order = conv0.order)
        for norm in (:pf, :nupf)
            out = native_case(n = n, G = G, pf = pf,
                              nu = (norm === :pf ? 1.0 : nu), links = links,
                              Him = Him0, rho0 = imp0.rho0, fold = fold, res = res,
                              conv = conv0, nodek_res = nkr)
            errz = abs(out.Z - 1)
            @printf("  Z证书 trace=%-6s fold=%-5s res=%-9s nodek=%-5s norm=%-5s : Z = %.12f (|Z−1| = %.2e)\n",
                    trace_style, fold, res, nkr, norm, real(out.Z), errz)
            push!(best, (errz, (trace_style = trace_style, fold = fold, res = res,
                                nodek_res = nkr, norm = norm)))
        end
    end
    sort!(best, by = x -> x[1])
    winners = [b for b in best if b[1] < 1e-9]
    println("  → Z 通过组合数 = ", length(winners))
    foreach(b -> println("     ", b[2], "  |Z−1| = ", @sprintf("%.2e", b[1])), winners)

    # ---- ρ 判决（ε_d ≠ 0）： winners vs Fock 链 ----
    println()
    up, imp = setup(eps_d = 0.7)
    Him = ComplexF64[0 0; 0 imp.eps]
    ref = fock_chain_single(up, imp.eps, sched, imp.rho0)
    for (_, cfg) in winners
        G = native_G(up.bath_hamiltonian, V, BETA_N, cells, n;
                     trace_style = cfg.trace_style)
        pf = pf_mine_native(G; step_sign = conv0.step_sign, order = conv0.order)
        out = native_case(n = n, G = G, pf = pf,
                          nu = (cfg.norm === :pf ? 1.0 : nu), links = links,
                          Him = Him, rho0 = imp.rho0, fold = cfg.fold,
                          res = cfg.res, conv = conv0,
                          nodek_res = cfg.nodek_res)
        errr = maxabsdiff(out.rho, ref)
        @printf("  ρ判决 trace=%-6s fold=%-5s res=%-9s nodek=%-5s norm=%-5s : ρ-err = %.2e  %s\n",
                cfg.trace_style, cfg.fold, cfg.res, cfg.nodek_res, cfg.norm,
                errr, errr < 1e-9 ? "PASS" : "FAIL")
    end
    return winners
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
