# exactif.jl：严格核构造器（exactif_port_plan；E1 核心内核 + E4′/E6 契约）
#
# 定位（F6 用户裁决）：ExactIF 从"独立 IF 算法"转为"**严格核构造器**"（kernel
# 层）——`exactif_kernel` 产出 `ExactHybridizationKernel` 后，IF 算法层
# （PartialIF/PolarIF）零特例消费；ExactTTI 类型级拒绝（v2.3 算法支持集，
# 守卫在 influencefunctional.jl 的 hybriddynamics! 入口）。
#
# 移植源 = extensions/exactif/src/NativeExact.jl（1698 行；旧栈 Z2Tensors +
# QuAPI AbstractRealCorrelationFunction 行扫消费）。裁剪面（E5 用户裁决）：
# order ∈ {1, 2}（order 4 S4 门序列含负时间演化不移植）、等距离散
# （exact_Δt 非等距 cells 不移植）。跨支 Keldysh 负号在 blockification 期
# 折进 blocks（F2 符号约定收敛：构造期收编、读点直读）。
#
# 两条构造线（本文件 Phase B 落 QuAPI-Δ 线；compact 线随 Phase C 管线）：
#   QuAPI-Δ：IFD.exact_fermionic_Δt_cells 的补零 pair 网格（trotter_factors
#           调度）——Δ = 骨架相对 (ā,a) 点值矩阵、normdelta = prefactor
#           （Π n̄ₖ·det(BB)，幺正门相位下 Z ≡ 1 契约）；
#   compact：FormB 紧凑消元（两扇区 bath Gaussian 的 G⁻¹ 反对称化）——
#           normdelta = ν·Pf(G)。
#
# 腿-格点映射（新栈约定）：lattice 实时 Nt = npairs + 1 片（每 gap 一片 +
# **末位恒等窗**——旧栈 real 格 k = Nt+1 节点含末端开放点，恒等窗携带该窗的
# 相干态分辨率重叠、不可省），IF 腿 = Δ 点 1..Nt（点 p ↔ 片 p；旧栈点直接
# 映射节点号，新格片 = 旧节点）。门表 gates = [gate_durations; 0] 经折叠
# fock_sysdynamics 的 gates kwarg 装配（调度自管，q=1 反 Wick 测度、order=1
# 组合——Strang 拆分在表内首尾半步）。
#
# 归一化契约（E4′）：Z = integrate(K, Is) · kernel.normdelta——normdelta 经
# sysdynamics(kernel=…) 的链幅通道（logscale + 首站号）与分子同携带，比值
# 观测量严格相消；管线直用乘法则走本文件的显式 normdelta。

# ---------------------------------------------------------------------------
using LinearAlgebra: det, Diagonal, I

# Trotter 调度（IFGTEMPO trotter_schedule 的移植；order 4 裁剪——E5）
# ---------------------------------------------------------------------------

"""
    NativeSchedule

Trotter 调度表：`cells` = 二次型 cell 时长 ×nl、`links` = 杂质门时长 ×(nl+1)
（compact 线布局；QuAPI-Δ 线用 `trotter_factors` 的补零 pair 网格）。
"""
struct NativeSchedule
    cells::Vector{Float64}   # 二次型 cell 时长 ×nl
    links::Vector{Float64}   # 杂质门时长 ×(nl+1)
end

"""`order ∈ (1, 2)` 的 Trotter 因子调度（与 extension native_schedule 逐位
一致；order 4 不移植——S4 门序列含负时间演化，E5 用户裁决）。

- order 1: `A(δt)B(δt)` 同 pair 放置，cells = fill(δt, N)，links = [δt…δt, 0]
  （末位 0：开放端处于末时刻）；
- order 2: Strang 合并链 `A(δt/2)[B(δt)A(δt)]^{N−1}B(δt)A(δt/2)`，
  links = [δt/2, δt…, δt/2]。
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
    else
        throw(ArgumentError("unsupported Trotter order $order (exactif ports orders 1 and 2 only)"))
    end
    @assert length(links) == length(cells) + 1
    return NativeSchedule(cells, links)
end

"""order ∈ (1, 2) 的补零 pair 网格表：每 gap 一个 Trotter 因子（QuAPI-Δ 与
调度门管线共用）。

gate_durations[g] = gap g 的杂质门时长（0 = 恒等门，不跳过——恒等窗携带该
pair 的相干态分辨率重叠）；cell_durations[g] = gap g 的二次型 cell 时长
（0 = 无 cell；传入 `IFD.exact_fermionic_Δt_cells` 的逐 cell 表，点数 =
length+1）。order 2 = [δt/2, δt, …, δt, δt/2] 门 + N 个 δt cell 交错。
"""
function trotter_factors(order::Integer, N::Integer, δt::Real)
    N >= 1 || throw(ArgumentError("number of macro steps must be at least one"))
    dtf = Float64(δt)
    if order == 1
        npairs = N
        gates = fill(dtf, N)
        cells = fill(dtf, N)
        cells[end] = 0.0   # 末 gap 仅恒等门（开放端在末时刻）
    elseif order == 2
        npairs = 2N + 1
        gates = zeros(npairs); cells = zeros(npairs)
        links = vcat(dtf / 2, fill(dtf, N - 1), dtf / 2)
        gates[1:2:end] .= links
        cells[2:2:end] .= dtf
    else
        throw(ArgumentError("unsupported Trotter order $order (exactif ports orders 1 and 2 only)"))
    end
    return (npairs = npairs, gate_durations = gates, cell_durations = cells)
end

# ---------------------------------------------------------------------------
# bath Gaussian 机器（compact 线 + 归一化；extension 逐字移植）
# ---------------------------------------------------------------------------

"""耦合单粒子 h（mode 1 = 杂质，能级 0；mode k+1 = bath mode k）。"""
function _coupled_h(V::Vector{ComplexF64}, Hbath::Matrix{ComplexF64})
    m = length(V)
    # h[imp, bath] = conj(V)（d†c 幅）、h[bath, imp] = V（c†d 幅）——与
    # 参照链 fock_chain_single / IFGTEMPO bath_one_particle_hamiltonian 一致
    h = zeros(ComplexF64, m + 1, m + 1)
    for k in 1:m
        h[1, k+1] = conj(V[k])
        h[k+1, 1] = V[k]
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
# QuAPI-Δ 严格核数据源（extension quapi_delta 移植；IFD.exact_fermionic_Δt_cells）
# ---------------------------------------------------------------------------

"""
    StrictDelta

严格核中间载体（extension ExactDelta/ExactIFKernel 的数据面，去 QuAPI 容器
化）：`pos` = 读序腿位（GrassmannContourIndex → 矩阵行/列号）、`M` = (ā,a)
读矩阵（原始 Δ 条目；跨支负号在 `exactif_kernel` 的 blockification 折入——
F2 构造期符号收编）、`normdelta` = 带号 Z 常数（E4′ 单一存储）。
"""
struct StrictDelta
    pos::Dict{GrassmannContourIndex,Int}
    M::Matrix{ComplexF64}
    normdelta::Float64
end

"""腿-格点映射表（新栈约定）：lattice 实时片 p（1..Nt）↔ Δ 点 p；开放端点
Nt+1 无腿（丢弃）。键 = GrassmannContourIndex(片, 支, 带, conj)，行 = ā 腿、
列 = a 腿（(ā,a) 读矩阵约定）。"""
function _strict_pos(lattice::AbstractRealGrassmannLattice)
    k = _branch_len(lattice, :+)
    pos = Dict{GrassmannContourIndex,Int}()
    for band in 1:lattice.bands, br in (:+, :-), p in 1:k
        row = 4 * k * (band - 1) + 2 * k * (br === :+ ? 0 : 1) + 2 * (p - 1) + 1
        pos[GrassmannContourIndex(p, br, band, true)] = row
        pos[GrassmannContourIndex(p, br, band, false)] = row + 1
    end
    return pos
end

"""严格容器（IFD ExactRealCorrelationFunction 四分支块）→ legpos (ā,a) 读
矩阵（IFD 腿基序无损重构：块 = Δ[支 b1 ∧ conj 行, 支 b2 ∧ 非 conj 列]）"""
function _strict_legmatrix(corr, N::Int)
    basis = [GrassmannContourIndex(j0 + 1, br, 1, c)
             for j0 in 0:N for br in (:+, :-) for c in (false, true)]
    Δ = zeros(ComplexF64, length(basis), length(basis))
    for b1 in (:+, :-), b2 in (:+, :-)
        rows = [i for (i, l) in enumerate(basis) if l.branch == b1 && l.conj]
        cols = [j for (j, l) in enumerate(basis) if l.branch == b2 && !l.conj]
        Δ[rows, cols] .= IFD.branch(corr; b1=b1, b2=b2)
    end
    return basis, Δ
end

"""
    quapi_strictdelta(bath, cells, lattice) -> StrictDelta

QuAPI-Δ 严格核数据源（extension quapi_delta 的移植）：`IFD.exact_fermionic_Δt_cells`
的骨架相对 (ā,a) 读矩阵 + prefactor。**cell 补零约定**：Δ 点数 =
length(cells)+1；片 p（1..Nt）↔ 点 p，开放端点 Nt+1 无腿（丢弃）。Δ 原始块
= 物理指数系数（extension 实证标定：M 层无需附加跨支符号——跨支负号由
blockification 读取规则折入，见 `exactif_kernel`）。
"""
function quapi_strictdelta(bath, cells::Vector{Float64},
                           lattice::AbstractRealGrassmannLattice)
    # 纯浴对象契约（无 kicks——IF 组装不依赖杂质哈密顿量；杂质全项走 K 链门）
    corr, prefactor = IFD.exact_fermionic_Δt_cells(bath.spectrum; β=bath.β, cells=cells, μ=bath.μ)
    basis, Δ = _strict_legmatrix(corr, length(cells))
    pos = _strict_pos(lattice)
    k = _branch_len(lattice, :+)
    M = zeros(ComplexF64, 4 * lattice.bands * k, 4 * lattice.bands * k)
    qpos = Dict{GrassmannContourIndex,Int}(l => i for (i, l) in enumerate(basis))
    for (l1, i1) in qpos, (l2, i2) in qpos
        (l1.conj && !l2.conj) || continue
        p1, p2 = l1.j, l2.j
        (1 <= p1 <= k && 1 <= p2 <= k) || continue   # 开放端点无腿
        M[pos[GrassmannContourIndex(p1, l1.branch, l1.band, true)],
          pos[GrassmannContourIndex(p2, l2.branch, l2.band, false)]] = Δ[i1, i2]
    end
    return StrictDelta(pos, M, _real_normdelta(prefactor))
end

# ---------------------------------------------------------------------------
# exactif_kernel：StrictDelta → ExactHybridizationKernel（blockification）
# ---------------------------------------------------------------------------

"""
    exactif_kernel(bath, N, δt; order=2, bands=1) -> ExactHybridizationKernel
    exactif_kernel(delta::StrictDelta, lattice; bands=nothing) -> ExactHybridizationKernel

严格核构造器（E1/E4′/E6）：Trotter 补零 pair 网格（`trotter_factors`）调度 +
`IFD.exact_fermionic_Δt_cells` 严格消元 → (ā,a) 读矩阵按 (branch, band) 分块
→ **ExactHybridizationKernel**（dense Matrix blocks、带号 normdelta 单一
存储）。

**跨支负号折入**（F2）：StrictDelta.M 的读取规则 = 跨支块取负
（extension getindex 的 `(i.branch != j.branch) ? -v : v` 收进构造期）——
blocks[(b₁≠b₂)] = −M、blocks[(b₁==b₂)] = +M，核存储 = 物理终值（与
correlationfunction 的 QuAPI 存储约定同收敛），`fullfactor` 按键直读无补号。

`bands` = 带数（多带键共享同一块对象，零拷贝别名——correlationfunction 同
款）。`normdelta` = IFD prefactor（Π n̄ₖ·det(BB)；Z = ∫·normdelta 契约，
幺正门相位下 Z ≡ 1 机器精确）。IF 算法层零特例：`hybriddynamics(lat,
kernel, PartialIF()/PolarIF())` 直接消费；ExactTTI 类型级拒绝（v2.3）。
"""
function exactif_kernel(delta::StrictDelta, lattice::AbstractRealGrassmannLattice)
    k = _branch_len(lattice, :+)
    bands = lattice.bands
    T = ComplexF64
    blocks = Dict{Tuple{ContourBlockIndex,ContourBlockIndex},Union{Matrix{T},CorrelationMatrix{T}}}()
    for band in 1:bands, b1 in (:+, :-), b2 in (:+, :-)
        rows = [delta.pos[GrassmannContourIndex(p, b1, band, true)] for p in 1:k]
        cols = [delta.pos[GrassmannContourIndex(p, b2, band, false)] for p in 1:k]
        blk = delta.M[rows, cols]
        # 跨支负号折入（F2）：读取规则收进构造期
        blocks[(ContourBlockIndex(b1, band, true), ContourBlockIndex(b2, band, false))] =
            b1 === b2 ? blk : -blk
    end
    return ExactHybridizationKernel(blocks, delta.normdelta, collect(1:bands))
end

function exactif_kernel(bath, N::Integer, δt::Real; order::Integer=2, bands::Int=1)
    factors = trotter_factors(order, N, δt)
    lattice = GrassmannLattice(N=factors.npairs + 1, δt=δt, contour=:real, bands=bands)
    delta = quapi_strictdelta(bath, collect(Float64.(factors.cell_durations)), lattice)
    return exactif_kernel(delta, lattice)
end

# ---------------------------------------------------------------------------
# 裸 Berezin K 链机器（sysdynamics_exact——extension _scheduled_gates! +
# boundarycondition! + systhermalstate! 的 GM 因子机器换代移植）
#
# 设计定谳（plan §7 取证）：严格契约的 K 侧 = 旧栈裸机器——逐门**跨片二站点
# 传播子因子**（bra 变量在节点 idx+1、ket 在节点 idx；− 支 bra/ket 互换）+
# Keldysh 折返 exp(−ā⁻(N)a⁺(N)) + ρ₀ 因子（bra ā⁺@(+,1)、ket a⁻@(−,1)），
# 全部用 GM 因子原语（`lattice_monomial` 单项式 + `apply!` 楔积）表达。
# **不复用** `fock_sysdynamics` 的 S·W 权重胞——那是冻结核路线的
# anti-Wick/Weyl 测度象（O(δt)/O(δt²)，靠 Weyl 标量补偿、ED 锚 1e-4 容差
# 不可见），严格契约要机器精确。
#
# 节点-片映射（Weyl 胞布局定谳，plan §7.3）：片 j 的 a 腿 = 节点 j、ā 腿 =
# 节点 j+1 ⟹ 门 g（连接节点 g → g+1）的因子变量 = {ā@片 g+1, a@片 g}——
# 跨片二站点；IF Δ 点 p ↔ 片 p 的 (ā_p, a_p) 双腿（`_strict_pos`）。
# lattice 携带 Nt = npairs + 1 片（门表 + 末位恒等窗 = 末节点）。
# ---------------------------------------------------------------------------

"Fock 占据位 → 带（band 1 = MSB，FockMatrix 约定与 fock_charges 同源）"
_occ_bands(occ::Int, B::Int) = [b for b in 1:B if (occ >> (B - b)) & 1 == 1]

"""
    _fockop_factor(lat, M, bras, kets) -> GMPS

Fock 矩阵 M 的**精确**相干态单项式展开因子 Σ_{IJ} M[I,J]·m_I(bras)·m_J(kets)
（extension `_tosparsegmps` 展开的 GM 因子机器换代表达：逐单项式
`lattice_monomial`（GTerm 语义，给定序系数）→ `togmps` 物化 → GMPS 求和；
恒等项（I = J = 0 无变元）= 标量 × 真空链）。荷守恒矩阵契约：|I| ≠ |J| 的
矩阵元非零即抛（正常杂质传播子/密度矩阵全荷守恒——非中性单项式在 GM 因子
层不可表达，也物理为零）。单项式变元序 = c̄(bra, 带升序)·c(ket, 带降序)、
系数 = 矩阵元 × 序列反序奇偶（旧 `_fockpropagator_sites` 逐字；恒等门因子
1 − ā@g+1 a@g 自带该 pair 的相干态分辨率重叠——恒等门不可短路）。

用法 = `mult!(K, factor; trunc)` 一次压缩乘法上链（`truncrelerr(1e-14)` 与
旧机器 `Defaults.alg_trunc_sys` 同款——严格性在机器精度层保留）。
"""
function _fockop_factor(lat, M::AbstractMatrix,
                        bras::Vector{GrassmannContourIndex},
                        kets::Vector{GrassmannContourIndex})
    B = length(bras)
    (length(kets) == B) || throw(DimensionMismatch("bra/ket positions must match bands"))
    size(M) == (2^B, 2^B) || throw(DimensionMismatch(
        "expansion matrix must be $(2^B)×$(2^B), got $(size(M))"))
    dim = 2^B
    spaces = site_spaces(lat)
    fac = nothing
    for I in 0:dim-1, J in 0:dim-1
        v = ComplexF64(M[I+1, J+1])
        abs(v) < 1.0e-14 && continue
        occ_i = _occ_bands(I, B)
        occ_j = _occ_bands(J, B)
        length(occ_i) == length(occ_j) || throw(ArgumentError(
            "non-number-conserving Fock matrix (|M[$(I+1),$(J+1)]| = $(abs(v)) ≠ 0) — " *
            "the strict-route expansion only supports charge-conserving operators"))
        # 单项式序 = c̄(bra, 带升序) · c(ket, 带降序)，系数 = 矩阵元 × 序列反序
        # 奇偶（旧 _fockpropagator_sites 逐字；lattice_monomial 的给定序系数语义
        # 恰好承载同一对象：coeff·∏ 给定序 ⟹ 归一化反序号由其内部吸收）
        idx = vcat(bras[occ_i], reverse(kets[occ_j]))
        comp = isempty(idx) ? vacuumstate(ComplexF64, lat) :
               GM.togmps(ComplexF64, lattice_monomial(lat, Tuple(idx); coeff=v), spaces)
        fac = fac === nothing ? v * comp : fac + v * comp
    end
    return fac
end

"0..B−1 带的全部子集（升序带表；折返带子集展开用）"
_band_subsets(B::Int) = [[b for b in 1:B if (mask >> (b - 1)) & 1 == 1]
                          for mask in 0:(2^B - 1)]

"1..n 元素集的全部子集（升序元素表；逐点测度 2B-pair 子集展开用）"
_pair_subsets(n::Int) = [[k for k in 1:n if (mask >> (k - 1)) & 1 == 1]
                          for mask in 0:(2^n - 1)]

"""
    sysdynamics_exact(lattice, model, gates; initial_state=nothing) -> GMPS

严格核管线的**裸 Berezin** 系统因子（旧栈 `sysdynamics_exact!` 三件套的
GM 因子机器换代；`gates` = `trotter_factors` 的补零 pair 网格门时长表，
两支同表，Td 模型时间参数 = 门左端点累计、− 支 = T − Σᵢ≤g）：

1. 真空起始 + 逐门跨片传播子因子（+ 支 bra@节点 g+1/ket@节点 g、− 支
   bra@节点 g/ket@节点 g+1——`_propagator_sparsegmps` 的 ib/ik 约定逐字）；
2. Keldysh 折返 exp(Σ_b −ā⁻(N)a⁺(N))（`boundarycondition!` 逐字；多带 =
   带子集展开，系数 (−1)^{|S|}）；
3. ρ₀ 因子（`sysinitialstate!` 逐字：bra ā⁺@(+,1)、ket a⁻@(−,1)；缺省 =
   真空投影）。**无反测度**（QuAPI-Δ 相对测度 + 单位链骨架定义）。

乘法压缩 = `truncrelerr(1e-14)`（旧机器 `Defaults.alg_trunc_sys` 同款——
SVD 相对截断在机器精度层严格性保留，zidentity 证书仲裁）。⚠️ 恒等门不可
短路——其展开携带单位链 + 分辨率重叠（测度的一部分）。
"""
function sysdynamics_exact(lattice::AbstractRealGrassmannLattice, model::AbstractImpurityHamiltonian,
                           gates::Vector{Float64}; initial_state=nothing, counter_c::Int=0)
    npairs = length(gates)
    lattice.Nt == npairs + 1 || throw(DimensionMismatch(
        "strict-contour lattice carries Nt = npairs + 1 slices (gate table + " *
        "trailing identity window = final node), got Nt = $(lattice.Nt) vs gates = $npairs"))
    B = lattice.bands
    trunc = truncrelerr(1.0e-14)                              # 旧机器 alg_trunc_sys 同款
    spaces = site_spaces(lattice)
    K = vacuumstate(ComplexF64, lattice)
    T = sum(gates)
    vac = vacuumstate(ComplexF64, lattice)
    for g in 1:npairs
        for branch in (:+, :-)
            t = branch === :+ ? sum(gates[1:g-1]; init=0.0) :
                T - sum(gates[1:g]; init=0.0)
            M = Matrix(fock_propagator(model, branch, gates[g], t).data)
            # ⚠️ 恒等门不可跳过——恒等传播子展开 = 1 + ā@节点g+1 a@节点g（单位链
            # + 该 pair 的相干态分辨率重叠，extension「恒等门不能跳过」纪律）
            bras = branch === :+ ?
                [GrassmannContourIndex(g + 1, :+, b, true) for b in 1:B] :
                [GrassmannContourIndex(g, :-, b, true) for b in 1:B]
            kets = branch === :+ ?
                [GrassmannContourIndex(g, :+, b, false) for b in 1:B] :
                [GrassmannContourIndex(g + 1, :-, b, false) for b in 1:B]
            fac = _fockop_factor(lattice, M, bras, kets)
            K = mult!(K, fac; trunc=trunc)[1]
        end
    end
    # 1b. 接口分辨率反测度（counter_c ≠ 0：compact 核已吸收接口 resolution，
    #     逐接口片逐支逐带乘幂零 (1 + counter_c·x·x̄) 抵消隐式重复——extension
    #     sysdynamics_exact! counter_c 步骤；compact 路线标定 counter_c = −1）
    counter_c != 0 && for p in 2:(lattice.Nt - 1), br in (:+, :-), b in 1:B
        fac = one(ComplexF64) * vac +
              GM.togmps(ComplexF64, lattice_monomial(lattice,
                    (GrassmannContourIndex(p, br, b, true), GrassmannContourIndex(p, br, b, false));
                    coeff=ComplexF64(counter_c)), spaces)
        K = mult!(K, fac; trunc=trunc)[1]
    end
    # 2. Keldysh 折返 exp(Σ_b +ā⁻(N)a⁺(N))（带子集展开系数 (+1)^{|S|}；系数
    #    +1 由 zidentity 证书机器精度钉定——旧 boundarycondition! 的 GTerm −1
    #    在 lattice_monomial 给定序语义下换号为 +1）
    np = lattice.Nt
    fold = one(ComplexF64) * vac
    for mask in _band_subsets(B)
        isempty(mask) && continue
        idx = GrassmannContourIndex[]
        for b in mask
            push!(idx, GrassmannContourIndex(np, :-, b, true))
            push!(idx, GrassmannContourIndex(np, :+, b, false))
        end
        fold = fold + GM.togmps(ComplexF64, lattice_monomial(lattice, Tuple(idx); coeff=1.0), spaces)
    end
    K = mult!(K, fold; trunc=trunc)[1]
    # 3. ρ₀（bra ā⁺@(+,1)、ket a⁻@(−,1)）
    rho = initial_state === nothing ?
          Matrix{ComplexF64}(Diagonal(vcat([1.0], zeros(ComplexF64, 2^B - 1)))) :
          initial_state isa FockMatrix ? Matrix{ComplexF64}(initial_state.data) :
          Matrix{ComplexF64}(initial_state)
    size(rho) == (2^B, 2^B) || throw(DimensionMismatch(
        "initial state must be $(2^B)×$(2^B)"))
    rho ≈ Matrix{ComplexF64}(I, 2^B, 2^B) || begin
        bras = [GrassmannContourIndex(1, :+, b, true) for b in 1:B]
        kets = [GrassmannContourIndex(1, :-, b, false) for b in 1:B]
        K = mult!(K, _fockop_factor(lattice, rho, bras, kets); trunc=trunc)[1]
    end
    return K
end

# ---------------------------------------------------------------------------
# compact 消元核（extension compact_delta 逐字移植 + blockification）
# ---------------------------------------------------------------------------

"""
    compact_strictdelta(Hbath, V, β, sched, lattice; sgn=-1) -> StrictDelta

紧凑消元核（extension `compact_delta` 的逐字移植）：bath pair = trace_final +
forward(N) + backward(N) 共 2N+1 对；forward 末步 trace_bra 翻转（⟨−ψ̄| 约定）；
热边 first_forward.bars × first_backward.etas；每 pair 测度 −γ̄γ；Schur 消 bath
块取反对称；消元整体负号（FormB elim_sign = −1）。G 在接口 IF 槽位
（4·nl 维）；D = 消元核 = G⁻¹ 反对称化的 (bar, eta) partner 读矩阵；
normdelta = ν·Pf(G)。接口 c ∈ 1..nl ↔ lattice 片 c+1（whole-slice 节点映射）。
"""
function compact_strictdelta(Hbath::Matrix{ComplexF64}, V::Vector{ComplexF64}, β::Real,
                             sched::NativeSchedule, lattice::AbstractRealGrassmannLattice;
                             sgn::Int=-1)
    cells = sched.cells
    nl = length(cells)
    m = length(V)
    nif = 4 * nl
    total = nif + 2m * (2nl + 1)
    A = zeros(ComplexF64, total, total)
    h = _coupled_h(V, Hbath)
    Us = [exp(-1im * cell .* h) for cell in cells]
    cursor = nif
    pair_data = NamedTuple[]
    for label in (:trace_final, :forward, :backward)
        count = label === :trace_final ? 1 : nl
        for step in 0:(count - 1)
            bars = collect((cursor + 1):(cursor + m))
            etas = collect((cursor + m + 1):(cursor + 2m))
            push!(pair_data, (label = label, step = step, bars = bars, etas = etas))
            cursor += 2m
        end
    end
    trace_pair = pair_data[1]
    forward_pairs = pair_data[2:(nl + 1)]
    backward_pairs = pair_data[(nl + 2):end]
    for step in 0:(nl - 1)
        input_pair = forward_pairs[step + 1]
        output_pair = step == nl - 1 ? trace_pair : forward_pairs[step + 2]
        _q_add!(A, _ifslot(step + 1, :+, :bar), _ifslot(step + 1, :+, :eta),
                output_pair.bars, input_pair.etas, Us[step + 1];
                trace_bra = (step == nl - 1))
    end
    for step in 0:(nl - 1)
        output_pair = backward_pairs[step + 1]
        input_pair = step == nl - 1 ? trace_pair : backward_pairs[step + 2]
        _q_add!(A, _ifslot(step + 1, :-, :bar), _ifslot(step + 1, :-, :eta),
                output_pair.bars, input_pair.etas, adjoint(Us[step + 1]))
    end
    thermal = exp(-Float64(β) .* Hbath)
    for i in 1:m, j in 1:m
        _put_G!(A, forward_pairs[1].bars[i], backward_pairs[1].etas[j], thermal[i, j])
    end
    for pair in pair_data, mode in 1:m
        _put_G!(A, pair.bars[mode], pair.etas[mode], -1.0)
    end
    E = A[1:nif, 1:nif]
    C = A[1:nif, nif+1:end]
    B = A[nif+1:end, nif+1:end]
    effective = E + C * (B \ transpose(C))
    G = 0.5 .* (effective .- transpose(effective))
    ν = bath_norm(V, Hbath, β, cells)
    Ginv = inv(G)
    GA = (Ginv .- transpose(Ginv)) ./ 2
    k = lattice.Nt
    pos = Dict{GrassmannContourIndex,Int}()
    for band in 1:lattice.bands, br in (:+, :-), p in 1:k
        row = 4 * k * (band - 1) + 2 * k * (br === :+ ? 0 : 1) + 2 * (p - 1) + 1
        pos[GrassmannContourIndex(p, br, band, true)] = row
        pos[GrassmannContourIndex(p, br, band, false)] = row + 1
    end
    M = zeros(ComplexF64, 4 * lattice.bands * k, 4 * lattice.bands * k)
    for band in 1:lattice.bands, c in 1:nl, br in (:+, :-),
        c2 in 1:nl, br2 in (:+, :-)
        # 消元恒等式 exp(−½jᵀG⁻¹j) 的整体负号（FormB elim_sign = −1）
        v = sgn * GA[_partner_slot(c, br, :bar), _partner_slot(c2, br2, :eta)]
        row = pos[GrassmannContourIndex(c + 1, br, band, true)]
        col = pos[GrassmannContourIndex(c2 + 1, br2, band, false)]
        iszero(v) && continue
        M[row, col] = v
    end
    return StrictDelta(pos, M, _real_normdelta(ν * pfaffian(G)))
end

# ---------------------------------------------------------------------------
# 严格核管线（Phase B 主验收形态；extension exactif_native_quapi 换代移植）
# ---------------------------------------------------------------------------

"""
    exactif_pipeline(bath, model, N, δt; order=2, bands=1, initial_state=nothing)
        -> (lat, K, Is, kernel, factors)

严格核主流管线（extension exactif_native_quapi 的新栈换代；**无反测度**——
QuAPI-Δ 相对测度 + 单位链骨架定义，GM 收缩隐式测度正是其所要）：

1. `trotter_factors(order, N, δt)` 补零 pair 网格 → lattice（Nt = npairs 片）；
2. `exactif_kernel` → ExactHybridizationKernel → `hybriddynamics(PartialIF)`
   （IF 算法层零特例；多带 = 单链覆盖）；
3. K 链 = `sysdynamics_exact`（裸 Berezin 机器：逐门跨片传播子因子 +
   Keldysh 折返 + ρ₀/`initial_state`——见其 docstring；**不复用**
   `fock_sysdynamics` 的 S·W 权重胞）；
4. `Z = integrate(K, Is) · kernel.normdelta`（E4′：normdelta 与分子同携带；
   经 sysdynamics(kernel=…) 的链幅通道时自动携带，比值观测量相消）。

**纯浴对象契约**：核不携带任何杂质信息（IFD uniform API，无 kicks）——杂质全项（εd、U、时间依赖）走 K 链门（`sysdynamics_exact` 的全模型传播子），IF 组装对一般杂质模型无假设（extensions/exactif 目标定义）。

插入观测量（greater/lesser/occupation 级）在管线层重建：K′ =
`sysdynamics_exact(lat, model, pipe.gates; …)`——调度表与 Z 的 K 链同源
（公开测量入口的 sysdynamics 重建走均匀门，与 pair 网格调度不同源，
Exact 路线不用；插入算符重建面留 V2）。
"""
function exactif_pipeline(bath, model::AbstractImpurityHamiltonian, N::Integer, δt::Real;
                          order::Integer=2, bands::Int=1, initial_state=nothing)
    factors = trotter_factors(order, N, δt)
    gates = collect(Float64.(factors.gate_durations))              # 原始门表（npairs 条）
    lat = GrassmannLattice(N=factors.npairs + 1, δt=δt, contour=:real, bands=bands)
    delta = quapi_strictdelta(bath, collect(Float64.(factors.cell_durations)), lat)
    kernel = exactif_kernel(delta, lat)
    Is, _ = hybriddynamics(lat, kernel, PartialIF())
    K = sysdynamics_exact(lat, model, gates; initial_state=initial_state)
    return (lat = lat, K = K, Is = Is, kernel = kernel, factors = factors,
            gates = gates)
end

"严格核管线的配分函数（Z = integrate·normdelta；无插入分母）"
function exactif_Z(pipe)
    return integrate(pipe.K, pipe.Is) * pipe.kernel.normdelta
end

"严格核管线的插入观测量（管线层 K′ 重建——调度表与 Z 链同源）：
`exactif_expectation(pipe, model, insertions; initial_state)` = (num, Z)"
function exactif_expectation(pipe, model::AbstractImpurityHamiltonian,
                             insertions::Vector{OperatorInsertion}; initial_state=nothing)
    K′ = sysdynamics_exact(pipe.lat, model, pipe.gates; initial_state=initial_state)
    num = integrate(K′, pipe.Is)
    Z = integrate(pipe.K, pipe.Is) * pipe.kernel.normdelta
    return num / Z
end
