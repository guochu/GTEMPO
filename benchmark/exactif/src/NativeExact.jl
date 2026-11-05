# =============================================================================
# src/NativeExact.jl —— 形态 B 消元式的 GTEMPO 原生管线（摆脱 IFGTEMPO）
#
# 把 FormB（消元式单扇区，与两扇区/Fock 逐位对齐）的构造移植到 GTEMPO 原生
# 实时管线：原生 lattice（只含杂质 Grassmann 变量）+ 原生 partialif 行扫 +
# 原生 sysdynamics/bulkconnection/boundarycondition/systhermalstate + 原生
# 观测量（integrate/gf/occupation）。
#
#   1. lattice  = GrassmannLattice(N = nl+1, δt, contour=:real, bands)
#                 node 1 = t=0（ρ₀ 层）；node 2..nl+1 = 接口 1..nl（D 源所在）；
#                 node k = nl+2 = 开放端（T：折返 + ρ(T) 观测）。
#   2. IF       = ExactDelta <: QuAPI.AbstractRealCorrelationFunction：
#                 D[λ,λ'] = −½(G⁻¹[i,k] − G⁻¹[k,i])（两扇区 bath Gaussian 的
#                 G⁻¹ 反对称化；λ,λ' = 接口轨迹腿，i,k = 其 bar↔eta 交叉
#                 partner 槽位）。getindex 按 ExactIFKernel 的读序约定
#                 （v = D[ā(i), a(j)]，跨支内嵌 Keldysh 负号）。
#                 由 GTEMPO hybriddynamics/partialif_hybrid 行扫消费。
#                 （稠密判决：消元核本身就是 (bar,eta) 有向对，行扫完整承载。）
#   3. sysdynamics_exact! = 调度版杂质门（sysdynamics_factors! 同款全窗口）+
#                 接口分辨率反测度（Δ 已吸收接口 resolution，乘幂零 (1+x̄x)
#                 抵消 GTEMPO 隐式 g_trace 的重复施加——FormB 同款技巧）+
#                 boundarycondition!（折返）+ systhermalstate!（ρ₀）。
#   4. 观测量   = 原生 integrate / gf / occupation，Z 含归一化 ν·Pf(G)。
#
# 两扇区 G 的 native 腿装配（QuAPI exact-Δt 同构的 pair 布局）：
#   bath pair per (mode, branch, node 1..nl+1)；cell c 连接 node c → c+1，
#   锚接口 c（= node c+1）；热边 (+,1)/(−,1)；迹边 (·,nl+1)。
#
# 参考：benchmark/exactif/form_b/（MATHEMATICS.md、formb_dense.jl、
# native_probe.jl 约定判决），README"形态 B"一节。
# =============================================================================
module NativeExact

using GTEMPO
using GTEMPO: TK
using LinearAlgebra
using QuAPI: AbstractRealCorrelationFunction
using Z2Tensors: TruncationScheme

export NativeSchedule, native_schedule, ExactDelta, exact_delta,
       sysdynamics_exact!, exactif_native_pipeline, open_rho_native

# ---------------------------------------------------------------------------
# Trotter 调度（IFGTEMPO trotter_schedule 的移植，GTEMPO 侧）
# ---------------------------------------------------------------------------

struct NativeSchedule
    cells::Vector{Float64}   # 二次型 cell 时长 ×nl
    links::Vector{Float64}   # 杂质门时长 ×(nl+1)
end

"""`order ∈ (1,2,4)` 的 Trotter 因子调度（与 IFGTEMPO trotter_schedule 逐位一致）。

- order 1: `A(δt)B(δt)` 同 pair 放置，cells = fill(δt, N)，links = [δt…δt, 0]
  （末位 0：开放端处于末时刻）；
- order 2: Strang 合并链 `A(δt/2)[B(δt)A(δt)]^{N−1}B(δt)A(δt/2)`，
  links = [δt/2, δt…, δt/2]；
- order 4: Suzuki `S4 = S2(aδt)²S2(bδt)`（2a+b=1，b<0），cells = repeat([a,a,b]δt,N)。
"""
function native_schedule(order::Integer, N::Integer, dt::Real)
    N >= 1 || throw(ArgumentError("number of macro steps must be at least one"))
    dtf = Float64(dt)
    if order == 1
        cells = fill(dtf, N)
        links = vcat(fill(dtf, N), 0.0)
    elseif order == 2
        cells = fill(dtf, N)
        links = vcat(dtf / 2, fill(dtf, N - 1), dtf / 2)
    elseif order == 4
        a = 1 / (2 - cbrt(2))
        b = 1 - 2a
        c = (a + b) / 2
        cells = repeat([a * dtf, a * dtf, b * dtf], N)
        links = vcat(a * dtf / 2, a * dtf, c * dtf)
        for _ in 2:N
            links = vcat(links, c * dtf, a * dtf, c * dtf)
        end
        links = vcat(links, b * dtf / 2)
    else
        throw(ArgumentError("unsupported Trotter order $order (use 1, 2 or 4)"))
    end
    @assert length(links) == length(cells) + 1
    return NativeSchedule(cells, links)
end

# ---------------------------------------------------------------------------
# 两扇区 bath Gaussian G（native 腿装配）+ 归一化 + Pfaffian
# ---------------------------------------------------------------------------

"""耦合单粒子 h（mode 1 = 杂质，能级 0；mode k+1 = bath mode k）。"""
function _coupled_h(V::Vector{ComplexF64}, Hbath::Matrix{ComplexF64})
    m = length(V)
    h = zeros(ComplexF64, m + 1, m + 1)
    for k in 1:m
        h[1, k+1] = h[k+1, 1] = V[k]
        h[k+1, k+1] = Hbath[k, k]
    end
    return h
end

"""反对称占位（1-based）：上三角 += v、下三角 −= v（i≠j）。"""
function _put_G!(A::Matrix{ComplexF64}, i::Int, j::Int, value)
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

"""接口 c（1..nl）的 IF 槽位（G 矩阵 1-based，接口局部索引）。

槽序同 native 腿序的接口投影：(bar⁺, eta⁺, bar⁻, eta⁻)。"""
_ifslot(c::Int, br::Symbol, stat::Symbol) =
    4 * (c - 1) + (br === :+ ? 1 : 3) + (stat === :bar ? 0 : 1)

"""接口 IF 槽的测度 partner（bar↔eta 交叉）。"""
_partner_slot(c::Int, br::Symbol, stat::Symbol) =
    _ifslot(c, br, stat === :bar ? :eta : :bar)

"""bath pair（mode k, branch br, node p）的两个槽位（bar, eta），1-based。

bath pair 只存在于 node 1..nl+1：node 1 挂热边（t=0），node nl+1 挂迹边
（物理 T 端，= 接口 nl）；开放端无 cell 耦合、无 bath pair。"""
function _bslot(nl::Int, m::Int, k::Int, br::Symbol, p::Int)
    nbp = nl + 1
    base = 4 * nl + (k - 1) * 4 * nbp + (br === :+ ? 0 : 2 * nbp) + (p - 1) * 2
    return base + 1, base + 2
end

"""Q kernel 缝合（IFGTEMPO _add_mapped_q_matrix! 的移植）：把 (bar_imp,
eta_imp) IF 对与 bath 输出 bars/输入 etas 经 U 缝合；trace_bra 翻转输出行
（compact 的 ⟨−ψ̄| 约定）。"""
function _q_add!(A::Matrix{ComplexF64}, bar_imp::Int, eta_imp::Int,
                 bath_bars::Vector{Int}, bath_etas::Vector{Int}, U::AbstractMatrix;
                 trace_bra::Bool=false)
    m = length(bath_bars)
    size(U) == (m + 1, m + 1) || throw(DimensionMismatch("Q 传播子维度不匹配"))
    _put_G!(A, bar_imp, eta_imp, U[1, 1])
    for j in eachindex(bath_etas)
        _put_G!(A, bar_imp, bath_etas[j], U[1, j + 1])
    end
    for i in eachindex(bath_bars)
        value = U[i + 1, 1]
        trace_bra && (value = -value)
        _put_G!(A, bath_bars[i], eta_imp, value)
        for j in eachindex(bath_etas)
            bath_value = U[i + 1, j + 1]
            trace_bra && (bath_value = -bath_value)
            _put_G!(A, bath_bars[i], bath_etas[j], bath_value)
        end
    end
    return A
end

"""
    contour_G(Hbath, V, β, cells; trace_style=:edge) -> (G, ν)

native 布局的两扇区 bath Gaussian `G`（接口 IF 槽位上，4·nl 维反对称）与
归一化 `ν = ZN/Z`。`Hbath`（bath 单粒子矩阵）、`V`（耦合幅）、`β`、`cells`
（二次型 cell 时长）。

装配（QuAPI exact-Δt 同构的 pair 布局 + FormB 的两扇区语义）：
cell c 连接 node c → c+1，锚接口 c 的 (br) 对；bath pair 链
(+,1)→(+,2)→…→(+,nl+1)（forward）、(−,nl+1)→…→(−,1)（backward）；
热边 (+,1) bars × (−,1) etas（exp(−β·Hbath)）；迹边（:edge = (−,nl+1) bar ×
(+,nl+1) eta 系数 −1；:flip = 末个 forward Q 输出行翻转）；每 bath pair
−γ̄γ 测度。Schur 消 bath 块后取反对称部分。
"""
function contour_G(Hbath::Matrix{ComplexF64}, V::Vector{ComplexF64}, β::Real,
                   cells::Vector{Float64}; trace_style::Symbol=:edge)
    m = length(V)
    nl = length(cells)
    nif = 4 * nl
    nbp = nl + 1
    nb = 4 * m * nbp
    A = zeros(ComplexF64, nif + nb, nif + nb)
    h = _coupled_h(V, Hbath)
    Us = [exp(-1im * cell .* h) for cell in cells]
    bars(br, p) = [_bslot(nl, m, kk, br, p)[1] for kk in 1:m]
    etas(br, p) = [_bslot(nl, m, kk, br, p)[2] for kk in 1:m]
    for c in 1:nl   # forward：输入 (+,c) etas → 输出 (+,c+1) bars，锚接口 c
        flip = (trace_style === :flip && c == nl)
        _q_add!(A, _ifslot(c, :+, :bar), _ifslot(c, :+, :eta),
                bars(:+, c + 1), etas(:+, c), Us[c]; trace_bra = flip)
    end
    for c in 1:nl   # backward：输入 (−,c+1) etas → 输出 (−,c) bars，锚接口 c
        _q_add!(A, _ifslot(c, :-, :bar), _ifslot(c, :-, :eta),
                bars(:-, c), etas(:-, c + 1), adjoint(Us[c]); trace_bra = false)
    end
    if trace_style === :edge
        for kk in 1:m
            _put_G!(A, _bslot(nl, m, kk, :-, nbp)[1],
                    _bslot(nl, m, kk, :+, nbp)[2], -1.0)
        end
    elseif trace_style !== :flip
        error("trace_style 必须是 :edge 或 :flip")
    end
    thermal = exp(-Float64(β) .* Hbath)
    for i in 1:m, j in 1:m
        _put_G!(A, _bslot(nl, m, i, :+, 1)[1], _bslot(nl, m, j, :-, 1)[2],
                thermal[i, j])
    end
    for kk in 1:m, br in (:+, :-), p in 1:nbp   # bath pair 测度 −γ̄γ
        _put_G!(A, _bslot(nl, m, kk, br, p)[1], _bslot(nl, m, kk, br, p)[2], -1.0)
    end
    E = A[1:nif, 1:nif]
    C = A[1:nif, nif+1:end]
    B = A[nif+1:end, nif+1:end]
    effective = E + C * (B \ transpose(C))
    G = 0.5 .* (effective .- transpose(effective))
    return G, bath_norm(V, Hbath, β, cells)
end

"""占据态模式列表。"""
_occ_modes(state::Int, m::Int) = [mode for mode in 0:(m - 1) if state & (1 << mode) != 0]

"""单粒子矩阵 S 的二次量子化（公开 Fock 基）。"""
function _second_quantized(S::AbstractMatrix)
    m = size(S, 1)
    dim = 1 << m
    out = zeros(ComplexF64, dim, dim)
    for bra in 0:(dim - 1), ket in 0:(dim - 1)
        bm, km = _occ_modes(bra, m), _occ_modes(ket, m)
        length(bm) == length(km) || continue
        out[bra + 1, ket + 1] = isempty(bm) ? 1.0 + 0.0im : det(S[bm .+ 1, km .+ 1])
    end
    return out
end

"""ν = ZN/Z：bath 链传播热权重的迹比（bath 配分归一；与 compact 公式一致）。"""
function bath_norm(V::Vector{ComplexF64}, Hbath::Matrix{ComplexF64}, β::Real,
                   cells::Vector{Float64})
    h = _coupled_h(V, Hbath)
    Bprod = Matrix{ComplexF64}(I, length(V), length(V))
    for cell in cells
        Bprod = exp(-1im * cell .* h)[2:end, 2:end] * Bprod
    end
    R = exp(-Float64(β) .* Hbath)
    bw = Bprod * R * adjoint(Bprod)
    return tr(_second_quantized(bw)) / tr(_second_quantized(R))
end

"""数值 Pfaffian（Parlett–Reid，O(n³)；FormB.jl 同款——勿用 2ⁿ 展开）。"""
function pfaffian(A0::Matrix{ComplexF64})
    n = size(A0, 1)
    (n % 2 == 0) || return zero(ComplexF64)
    A = copy(A0)
    pf = one(ComplexF64)
    for k in 1:2:n
        kp = k + 1
        kp > n && break
        if abs(A[kp, k]) < 1e-14
            r = 0
            for i in (kp + 1):n
                abs(A[i, k]) > 1e-14 && (r = i; break)
            end
            r == 0 && return zero(ComplexF64)
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

# ---------------------------------------------------------------------------
# ExactDelta：消元核 D 的 CorrelationFunction
# ---------------------------------------------------------------------------

"""
    ExactDelta <: QuAPI.AbstractRealCorrelationFunction

消元核 D（两扇区 bath Gaussian 的 G⁻¹ 反对称化，bar→eta 读矩阵）封装为
QuAPI 风格的实时关联函数，直接供 GTEMPO `hybriddynamics!` 的 partialif
行扫消费。字段：

- `pos`：读序腿位（GrassmannContourIndex → 列号；行/列都用 native node）；
- `M`：(ā, a) 读矩阵——`M[ā(node,branch), a(node,branch)] = GA[partner 槽]`，
  只在接口节点（2..nl+1）有源；
- `G`、`ν`、`pf`：两扇区 bath Gaussian、归一化、数值 Pfaffian
  （归一化 = ν·pf，由管线乘进 Z）。

getindex 按 ExactIFKernel 的读序约定：`v = M[ā(i 节点), a(j 节点)]`，
跨支块内嵌 Keldysh 负号（与 QuAPI 统一实时约定一致；partialif 行扫读数处
补偿回物理块）。
"""
mutable struct ExactDelta <: AbstractRealCorrelationFunction
    pos::Dict{GrassmannContourIndex,Int}
    M::Matrix{ComplexF64}
    G::Matrix{ComplexF64}
    ν::ComplexF64
    pf::ComplexF64
end
GTEMPO.TK.scalartype(::Type{<:ExactDelta}) = ComplexF64
branches(::ExactDelta) = (:+, :-)
normdelta(x::ExactDelta) = x.ν * x.pf

function Base.getindex(x::ExactDelta, i::GrassmannContourIndex, j::GrassmannContourIndex)
    (i.conj == j.conj) && return zero(ComplexF64)
    v = x.M[x.pos[GrassmannContourIndex(i.j, i.branch, i.band, true)],
            x.pos[GrassmannContourIndex(j.j, j.branch, j.band, false)]]
    return (i.branch != j.branch) ? -v : v
end

"""
    exact_delta(Hbath, V, β, sched, lattice; trace_style=:edge) -> ExactDelta

由两扇区 bath Gaussian 构造接口消元核。`sched`（NativeSchedule）给出 nl；
`lattice` 的 node p 承接 schedule 点 p−1（node 1 = t=0，接口 c = node c+1，
node k = 开放端）。D 只作用在接口腿上（node 2..nl+1）；读矩阵行 = ā 腿、
列 = a 腿（全节点），`M[ā(λ), a(λ')] = GA[partner(λ), partner(λ')]`。
"""
function exact_delta(Hbath::Matrix{ComplexF64}, V::Vector{ComplexF64}, β::Real,
                     sched::NativeSchedule, lattice::AbstractRealGrassmannLattice;
                     trace_style::Symbol=:edge)
    nl = length(sched.cells)
    (lattice.k == nl + 2) || throw(DimensionMismatch(
        "lattice 须含 nl+2 个节点（Nt = nl+1）；k = $(lattice.k)，nl = $nl"))
    G, ν = contour_G(Hbath, V, β, sched.cells; trace_style)
    GA = (inv(G) .- transpose(inv(G))) ./ 2
    k = lattice.k
    pos = Dict{GrassmannContourIndex,Int}()
    for band in 1:lattice.bands, br in (:+, :-), p in 1:k
        row = 4 * k * (band - 1) + 2 * k * (br === :+ ? 0 : 1) + 2 * (p - 1) + 1
        pos[GrassmannContourIndex(p, br, band, true)] = row    # 行（ā）
        pos[GrassmannContourIndex(p, br, band, false)] = row + 1   # 列（a）
    end
    M = zeros(ComplexF64, 4 * lattice.bands * k, 4 * lattice.bands * k)
    for band in 1:lattice.bands, c in 1:nl, br in (:+, :-),
        c2 in 1:nl, br2 in (:+, :-)
        # 行 = bar 腿（接口 c，node c+1），列 = eta 腿（接口 c2，node c2+1）
        v = GA[_partner_slot(c, br, :bar), _partner_slot(c2, br2, :eta)]
        iszero(v) && continue
        row = pos[GrassmannContourIndex(c + 1, br, band, true)]
        col = pos[GrassmannContourIndex(c2 + 1, br2, band, false)]
        M[row, col] = v
    end
    return ExactDelta(pos, M, G, ν, pfaffian(G))
end

# ---------------------------------------------------------------------------
# sysdynamics_exact!：调度门 + 接口反测度 + 折返 + ρ₀
# ---------------------------------------------------------------------------

"""
    sysdynamics_exact!(gmps, lattice, model, sched; bands, trunc, counter_c)
        -> GrassmannMPS

严格核管线的"系统侧"因子：把杂质动力学与积分骨架合并为一：

1. 调度版杂质门（`sysdynamics_factors!` 同款全窗口：gap g 上放置
   `fock_propagator(model, branch, sched.links[g])`；τ=0 恒等门不能跳过，
   其张量携带相干态分辨率重叠，是测度的一部分）；
2. 接口分辨率反测度：消元已把接口 resolution 吸收进 D，GTEMPO 隐式
   g_trace 会重复施加——逐接口（node 2..nl+1）逐支逐带乘幂零反测度
   `(1 + counter_c·x·x̄)`（(1−x̄x)(1+x̄x)=1 精确抵消；counter_c 默认 +1，
   与 FormB 的标定一致）；
3. `boundarycondition!`（Keldysh 折返闭合）；
4. `systhermalstate!`（ρ₀ 于 node 1）。
"""
function sysdynamics_exact!(gmps::GrassmannMPS, lattice::AbstractRealGrassmannLattice,
                            model::GTEMPO.AbstractImpurityHamiltonian, sched::NativeSchedule;
                            bands::Int=lattice.bands, trunc::TruncationScheme=GTEMPO.Defaults.alg_trunc_sys(),
                            counter_c::Int=1, β::Real=0.0)
    nl = length(sched.cells)
    (lattice.k == nl + 2) || throw(DimensionMismatch("lattice 须含 nl+2 个节点"))
    # 1. 调度门（全窗口，复用 common.jl 验证过的放置逻辑）
    gmps = _scheduled_gates!(gmps, lattice, model, sched; trunc)
    # 2. 接口反测度（node 2..nl+1，两支，各带）
    for band in 1:bands, node in 2:(nl + 1), br in (:+, :-)
        s_eta = lattice[GrassmannContourIndex(node, br, band, false)]
        s_bar = lattice[GrassmannContourIndex(node, br, band, true)]
        fac = GTEMPO.partialmpo(s_eta, [s_bar], [ComplexF64(counter_c)]) *
              GTEMPO.GrassmannMPS(ComplexF64, length(lattice))
        gmps, _ = GTEMPO.mult(gmps, fac; trunc = trunc)
    end
    # 3. 折返 + 4. ρ₀
    for band in 1:bands
        gmps = GTEMPO.boundarycondition!(gmps, lattice, band=band, trunc=trunc)
    end
    gmps = GTEMPO.systhermalstate!(gmps, lattice, model; trunc=trunc, β=β)
    return gmps
end

"""调度门放置（common.jl sysdynamics_factors! 的移植；τ=0 恒等门不跳过）。"""
function _scheduled_gates!(gmps::GrassmannMPS, lattice::AbstractRealGrassmannLattice,
                           model::GTEMPO.AbstractImpurityHamiltonian, sched::NativeSchedule;
                           trunc::TruncationScheme=GTEMPO.Defaults.alg_trunc_sys())
    (lattice.Nt == length(sched.cells) + 1) ||
        throw(DimensionMismatch("lattice.Nt 须为 nl+1"))
    cache = Dict{Tuple,Any}()
    for branch in (:+, :-), g in 1:(length(sched.links))
        fm = GTEMPO.fock_propagator(model, branch, sched.links[g])
        sparse = GTEMPO._propagator_sparsegmps(lattice, fm, g, branch, cache, nothing)
        gmps = GTEMPO.mult!(gmps, sparse, trunc=trunc)[1]
    end
    return gmps
end

# ---------------------------------------------------------------------------
# 完整管线
# ---------------------------------------------------------------------------

"""
    exactif_native_pipeline(bathp, model, N, δt; order=2, trunc, β, bands,
                            trace_style=:edge, counter_c=1)
        -> (lat, mpsK, mpsIs, deltas, norm)

GTEMPO 原生严格核实时管线（消元式）。`bathp`：每带的 (Hbath, V) 元组向量
（单带传 `(x,)`；两带传 `(up, dn)`，各为 NamedTuple(hb = Hbath, v = V)）。

返回 `lat`（原生 lattice）、`mpsK`（sysdynamics_exact 的系统因子）、
`mpsIs[band]`（ExactDelta 的 IF 因子，partialif 行扫生成）、`deltas[band]`
（ExactDelta 对象，含 G/ν/pf）、`norm = ∏_band ν·pf`（收缩 Z 须乘）。
观测量直接用 GTEMPO 的 integrate/gf/occupation（Z 乘 `norm`）。
"""
function exactif_native_pipeline(bathp::Union{Tuple,Vector}, model::GTEMPO.AbstractImpurityHamiltonian,
                                 N::Integer, δt::Real;
                                 order::Integer=2, trunc::TruncationScheme=GTEMPO.Defaults.alg_trunc_if(),
                                 β::Real=0.0, bands::Int=1,
                                 trace_style::Symbol=:edge, counter_c::Int=1)
    sched = native_schedule(order, N, δt)
    nl = length(sched.cells)
    lat = GTEMPO.GrassmannLattice(N=nl + 1, δt=δt, contour=:real, bands=bands)
    specs = normalize_bathspecs(bathp, bands)
    deltas = ExactDelta[]
    mpsIs = GTEMPO.GrassmannMPS[]
    norm = 1.0 + 0.0im
    for band in 1:bands
        d = exact_delta(specs[band].hb, specs[band].v, β, sched, lat;
                        trace_style)
        mpsI = GTEMPO.hybriddynamics(lat, d; band=band, trunc=trunc)
        push!(deltas, d)
        push!(mpsIs, mpsI)
        norm *= normdelta(d)
    end
    mpsK = sysdynamics_exact!(GTEMPO.vacuumstate(ComplexF64, lat), lat, model,
                              sched; bands, trunc, counter_c, β)
    return (lat = lat, mpsK = mpsK, mpsIs = mpsIs, deltas = deltas, norm = norm,
            sched = sched)
end

"""每带 bath 规格：(hb = Hbath 矩阵, v = 耦合幅向量)。单带可省元组层。"""
function normalize_bathspecs(bathp, bands::Int)
    specs = bathp isa Tuple ? collect(bathp) : collect(bathp)
    length(specs) == bands || throw(DimensionMismatch("bath 规格数须等于 bands"))
    out = NamedTuple[]
    for s in specs
        if s isa NamedTuple && haskey(s, :hb)
            push!(out, s)
        elseif s isa Tuple{Matrix{ComplexF64},Vector{ComplexF64}}
            push!(out, (hb = s[1], v = s[2]))
        else
            push!(out, (hb = s.bath_hamiltonian, v = s.hybridization))
        end
    end
    return out
end


# ---------------------------------------------------------------------------
# 开放腿 ρ(T) 读出（native lattice，node k 块；单带）
# ---------------------------------------------------------------------------
# 复合多项式（ExactDelta 行扫 IF × sysdynamics_exact 骨架）在"开放腿提取"
# 约定下精确（form_b/native_probe.jl 稠密判决；GTEMPO 的 GTerm 插入式读出
# 是另一种 T 端约定，不适用——见 README 约定记录）。单带实现：
#   1. 右扫 update_pair_right 闭合 pairs 3..L/2；
#   2. 边界块 sites 4..1 逐 site 吸收（两因子物理腿 g_fuse 合并为一）；
#      全部 4 腿保持开放（node-k 无分辨率，ρ 支路不需要——稠密判决）；
#   3. 折返腿 (site4 = bar⁻, site1 = eta⁺) 以 g_trace_phy 收缩（= 稠密
#      fold 因子 + 素积分的扫等价；FormB 标定）；
#   4. 从融合树 sector 直读 ρ（bra = site2 = bar⁺，ket = site3 = eta⁻）。
# 全局符号 ρ_sign 与 sc² 由 V=0 证书 + Fock 判决锁定。

using TensorOperations: @tensor
const GMPS = GrassmannMPS

"""右扫闭合边界块之外的 pairs（单带：pairs 3..L/2）。

env：{cod (d1, d2)（两因子 site-4 aux-right 双像）; dom (e,)}"""
function _right_sweep1(factors::Tuple{GMPS,GMPS})
    L = length(factors[1])
    tm = GrassmannTransferMatrix(factors...)
    env = r_RR(factors...)
    for j in div(L, 2):-1:3
        env = lmul!(scaling(tm), GTEMPO.update_pair_right(env, j, tm.states...))
    end
    return env
end

"""吸收 site s（两因子）：收缩 aux-right 双像，物理腿 g_fuse 合并放到
domain 尾部（e 之前）。

env 进：{cod (d1,d2); dom (opens…(ndom), e)}
env 出：{cod (l1,l2); dom (opens…, p, e)}"""
function _absorb_site2(env, T1::MPSTensor, T2::MPSTensor, ndom::Int)
    if ndom == 0
        @tensor t[l1, l2, p2, p1; e] := env[d2, d1, e] * T2[l2, p2, d2] *
                                        T1[l1, p1, d1]
    elseif ndom == 1
        @tensor t[l1, l2, p2, p1; o, e] := env[d2, d1, o, e] * T2[l2, p2, d2] *
                                           T1[l1, p1, d1]
    elseif ndom == 2
        @tensor t[l1, l2, p2, p1; o1, o2, e] := env[d2, d1, o1, o2, e] *
                                                T2[l2, p2, d2] * T1[l1, p1, d1]
    else
        @tensor t[l1, l2, p2, p1; o1, o2, o3, e] := env[d2, d1, o1, o2, o3, e] *
                                                    T2[l2, p2, d2] *
                                                    T1[l1, p1, d1]
    end
    t = GTEMPO.g_fuse(t, 3)                                   # fuse cod (p2, p1)
    # 重排：cod = (l1, l2)，dom = (opens…, p, e)
    if ndom == 0
        return TK.permute(t, (1, 2), (3, 4))
    elseif ndom == 1
        return TK.permute(t, (1, 2), (4, 3, 5))
    elseif ndom == 2
        return TK.permute(t, (1, 2), (4, 5, 3, 6))
    else
        return TK.permute(t, (1, 2), (4, 5, 6, 3, 7))
    end
end

"""单带开放腿 ρ(T)。factors = (mpsK, mpsIs[band])；返回 2×2 ρ。"""
function open_rho_native(factors::Tuple{GMPS,GMPS}; ρ_sign::Number=1.0 + 0im)
    env = _right_sweep1(factors)
    # 吸收 site 4 (bar⁻) → dom (p4, e)；site 3 (eta⁻) → (p4, p3, e)；
    # site 2 (bar⁺) → (p4, p3, p2, e)；site 1 (eta⁺) → (p4, p3, p2, p1, e)
    env = _absorb_site2(env, factors[1].data[4], factors[2].data[4], 0)
    env = _absorb_site2(env, factors[1].data[3], factors[2].data[3], 1)
    env = _absorb_site2(env, factors[1].data[2], factors[2].data[2], 2)
    env = _absorb_site2(env, factors[1].data[1], factors[2].data[1], 3)
    # 折返：重排 dom (p4, p3, p2, p1, e) → (p4, p1, p3, p2, e)，g_trace_phy
    # 总腿 = (l1, l2 | p4, p3, p2, p1, e)：p4 = 3, p1 = 6 → 相邻化后迹 (3,4)
    env = TK.permute(env, (1, 2), (3, 6, 4, 5, 7))
    env = GTEMPO.g_trace_phy(env, 3)
    env = lmul!(scaling(factors...)^2, env)
    # 末态：{cod (l1, l2); dom (p3 = ket, p2 = bra, e)}
    rho = zeros(ComplexF64, 2, 2)
    for (f1, f2) in TK.fusiontrees(env)
        unc = f2.uncoupled                              # dom = (p3, p2, e)
        (unc[3].n == 0) || continue                     # edge 扇区 0
        ket = unc[1].n                                  # eta⁻(k)
        bra = unc[2].n                                  # bar⁺(k)
        blk = env[f1, f2]
        rho[bra + 1, ket + 1] += ρ_sign * blk[1, 1]
    end
    return rho
end

end # module NativeExact
