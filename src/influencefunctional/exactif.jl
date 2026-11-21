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
3. K 链 = 折叠 `fock_sysdynamics(gates=…)`（调度门自管：逐 gap
   gate_durations、两支同表、q=1 反 Wick、ρ₀/`initial_state` 折格 1、迹闭合
   格 Nt——旧 boundarycondition!/systhermalstate! 的折叠内联语义）；
4. `Z = integrate(K, Is) · kernel.normdelta`（E4′：normdelta 与分子同携带；
   经 sysdynamics(kernel=…) 的链幅通道时自动携带，比值观测量相消）。

插入观测量（greater/lesser/occupation 级）在管线层重建：K′ =
`fockmpo2gmps(fock_sysdynamics(lat, model; order=1, gates=…, insertions=…),
lat)`——调度表与 Z 的 K 链同源（公开测量入口的 sysdynamics 重建走均匀门，
与 pair 网格调度不同源，Exact 路线不用）。
"""
function exactif_pipeline(bath, model::AbstractImpurityHamiltonian, N::Integer, δt::Real;
                          order::Integer=2, bands::Int=1, initial_state=nothing)
    factors = trotter_factors(order, N, δt)
    gates = vcat(collect(Float64.(factors.gate_durations)), 0.0)   # 末位恒等窗（开放端）
    lat = GrassmannLattice(N=factors.npairs + 1, δt=δt, contour=:real, bands=bands)
    delta = quapi_strictdelta(bath, collect(Float64.(factors.cell_durations)), lat)
    kernel = exactif_kernel(delta, lat)
    Is, _ = hybriddynamics(lat, kernel, PartialIF())
    K = fockmpo2gmps(fock_sysdynamics(lat, model; order=1, gates=gates,
                                      initial_state=initial_state), lat)
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
    K′ = fockmpo2gmps(fock_sysdynamics(pipe.lat, model; order=1,
                                       gates=pipe.gates,
                                       insertions=insertions, initial_state=initial_state),
                      pipe.lat)
    num = integrate(K′, pipe.Is)
    Z = integrate(pipe.K, pipe.Is) * pipe.kernel.normdelta
    return num / Z
end
