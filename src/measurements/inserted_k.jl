# inserted_k.jl：inserted-K 边界语义观测量（D5 两轨之一；secondorderif
# measurement.jl 同约定）。
#
# 全部经「带插入重建 K → GM 引擎 integrate → /Z」：
#   occupation(branch j; band)      n̂（逐带或总）在 (branch, j) 边界 :left
#   nn(i, j; b1, b2, band)          n̂ 双插入密度-密度（ED 验证 ✓）
#   greater/lesser/matsubara        ✓ 对拍 U1GTEMPO2 机器精度（2026-11-15）。
#                                   五处口径差全部定案（historical；插入构造已
#                                   符号化——plan/fockmpo_dress_selfread_plan.md
#                                   2026-11-20）：荷现在从符号 FockOp 按位派生
#                                   （V3——d=−q、ddag=+q 内建，无手工标签）；
#                                   ①builder prefix 不再消费插入元数据（converter
#                                   从穿戴链 δ 读数自读）；②push 序 = cutleft
#                                   foldr 首 push 最左（同切位 0⁺：greater =
#                                   d·ddag、lesser/电流 = ddag·d）；③matsubara
#                                   时序因子 (i≥j ? −1 : +1) 且混合轮廓 τ 支放行；
#                                   ④order/coupling_norm 显式 kwargs 必须转发；
#                                   ⑤total-n̂ = Σ_b n_b 不可表为单 FockOp（V6
#                                   推论）——逐带求和（期望线性、数值恒等）。
#   electriccurrent/heatcurrent     GF 求和精确版；raw_lead_correlator 同锚定
#                                   （平衡态净电流离散伪影 O(δt)，旧栈同款量级）。
# 边界 j ∈ 0..N（0 = ρ₀ 切位，⚠️ :- 支 j = 0 = 链末切位）；rt 顶点 branch=:+、
# τ 顶点 branch=:τ。τ 链 G(β,0) = +⟨d†d⟩（热环等时同号；j==0 平移路由属
# secondorderif M≥2 的 Z2 栈绕尾技巧，GM 栈 F≡0 前提下不需要——D5 验收时按
# golden 裁决）。
#
# 扫描型测量：K 重建走 fock_sysdynamics（穿戴内联）+ fockmpo2gmps（荷账自读，
# 2026-11-20 起唯一路径）；制备对象摊销已随 V2 删除——扫描成本回到逐点重建
# （FockMPO 级缓存/摊销归 V2 议题）。

"插入 K 重建的公共尾程：K′ = sysdynamics(...; insertions) → integrate(K′, Is...) / Z
（K 构造 = fock_sysdynamics（穿戴）+ fockmpo2gmps（自读）唯一路径）。
`kernel` 转发 sysdynamics——Z 常数按核类型分派（Frozen = Weyl 标量、
Exact = normdelta 链幅携带，exactif_port_plan E6）。"
function _inserted_expectation(lattice, model, Is, Z, insertions; order=1, coupling_norm=0.0,
                               kernel=nothing, initial_state=nothing, spec=nothing)
    K = sysdynamics(lattice, model; order, initial_state, coupling_norm, insertions,
                    kernel=kernel,
                    spec=spec === nothing ? TK.sectortype(lattice) : spec)
    num = integrate(K, Is...)
    return num / Z
end

"""
    occupation(lattice, model, kernel, Is...; j, branch, band=nothing,
               order=1, coupling_norm=0.0, initial_state=nothing, Z) -> Complex

⟨n̂(t_j)⟩（`band` 缺省 = 总占据 Σ_b n̂_b——逐带求和，期望线性数值恒等；j ≥ 1
在 (branch, j) 边界，j = 0 折 ρ₀/热迹）。`Is` = 各带 IF 链；`Z` = 无插入分母。
"""
function occupation(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian,
                    kernel::Union{Nothing,AbstractHybridizationKernel}, Is...; j::Integer, branch::Symbol=default_branch(lattice),
                    band::Union{Nothing,Int}=nothing, order::Integer=1,
                    coupling_norm::Real=0.0, initial_state=nothing, Z)
    if band === nothing
        acc = 0.0im
        for b in 1:lattice.bands
            ins = [OperatorInsertion(branch, Int(j), _single_band_op(:n, b, lattice.bands))]
            acc += _inserted_expectation(lattice, model, Is, Z, ins; order=Int(order),
                                         coupling_norm=coupling_norm, kernel=kernel, initial_state=initial_state)
        end
        return acc
    end
    ins = [OperatorInsertion(branch, Int(j), _single_band_op(:n, Int(band), lattice.bands))]
    return _inserted_expectation(lattice, model, Is, Z, ins; order=Int(order),
                                 coupling_norm=coupling_norm, kernel=kernel, initial_state=initial_state)
end

"缺省分支：虚时 :τ、实时/混合 rt 顶点 :+"
default_branch(lattice::AbstractGrassmannLattice) =
    lattice isa AbstractImagGrassmannLattice ? :τ : :+

"""
    raw_correlation(lattice, model, kernel, Is...; i, j=0, kind=:greater, ...) -> Complex

两点函数 ⟨d(tᵢ) ddag(tⱼ)⟩（greater）/ ⟨ddag(tⱼ) d(tᵢ)⟩（lesser）；实时/混合
专用（虚时用 matsubara）。0⁺ push 序约定见文件头。奇算符插入的荷从符号派生
（V3）——奇插入链总荷非零时由 converter 的 Σδ 守卫拒绝（ArgumentError）。
"""
function raw_correlation(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian,
                         kernel::Union{Nothing,AbstractHybridizationKernel}, Is...; i::Integer, j::Integer=0,
                         kind::Symbol=:greater, order::Integer=1,
                         coupling_norm::Real=0.0, initial_state=nothing)
    (lattice isa AbstractRealGrassmannLattice || lattice isa AbstractMixedGrassmannLattice) ||
        throw(ArgumentError("实时/混合轮廓专用；纯虚时用 matsubara"))
    bands = lattice.bands
    a = _single_band_op(:a, 1, bands)
    adag = _single_band_op(:adag, 1, bands)
    ins = OperatorInsertion[]
    if kind === :greater
        # ⟨d(tᵢ) ddag(tⱼ)⟩：荷 = 算符自荷（符号派生——:a = −q、:adag = +q）。
        # push 序 = 算符左乘序（cutleft foldr 首 push 最左、最后作用）：同切位
        # 0⁺ greater = d·ddag ⟹ d 先 push
        push!(ins, OperatorInsertion(:+, Int(i), a))
        push!(ins, OperatorInsertion(:+, Int(j), adag))
    elseif kind === :lesser
        # ⟨ddag(tⱼ) d(tᵢ)⟩：同切位 0⁺ = ddag·d（a†a）⟹ ddag 先 push
        push!(ins, OperatorInsertion(:-, Int(j), adag))
        push!(ins, OperatorInsertion(:+, Int(i), a))
    else
        throw(ArgumentError("kind must be :greater or :lesser"))
    end
    return _inserted_expectation(lattice, model, Is, 1.0, ins; order=Int(order),
                                 coupling_norm=coupling_norm, kernel=kernel, initial_state=initial_state)
end

"greater(tᵢ, tⱼ=0) = −i⟨d(tᵢ)ddag(tⱼ)⟩（归一化走 occupation 同款 Z 分母）"
function greater(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian,
                 kernel::Union{Nothing,AbstractHybridizationKernel}, Is...; i::Integer, j::Integer=0, Z,
                 kwargs...)
    (lattice isa AbstractRealGrassmannLattice || lattice isa AbstractMixedGrassmannLattice) ||
        throw(ArgumentError("虚时轮廓用 matsubara"))
    return -1.0im * _normalized_two_point(lattice, model, kernel, Is, Z, i, j, :greater; kwargs...)
end

"lesser(tᵢ, tⱼ=0) = +i⟨ddag(tⱼ)d(tᵢ)⟩"
function lesser(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian,
                kernel::Union{Nothing,AbstractHybridizationKernel}, Is...; i::Integer, j::Integer=0, Z,
                kwargs...)
    (lattice isa AbstractRealGrassmannLattice || lattice isa AbstractMixedGrassmannLattice) ||
        throw(ArgumentError("虚时轮廓用 matsubara"))
    return 1.0im * _normalized_two_point(lattice, model, kernel, Is, Z, i, j, :lesser; kwargs...)
end

"两点归一化公共体（greater/lesser 宿主）：order/coupling_norm 显式转发
（缺省随 kernel 元数据；hybridization 的 diagnostics 未登记 order 时返回
1——**golden 教训：忽略调用方显式值会以 order=1 链对 order=2 的 Z 归一**）"
function _normalized_two_point(lattice, model, kernel, Is, Z, i, j, kind;
                               order::Integer=1,
                               coupling_norm::Real=0.0, kwargs...)
    bands = lattice.bands
    a = _single_band_op(:a, 1, bands)
    adag = _single_band_op(:adag, 1, bands)
    ins = OperatorInsertion[]
    if kind === :greater
        # push 序 = 左乘序（同切位 0⁺ greater = d·ddag）
        push!(ins, OperatorInsertion(:+, Int(i), a))
        push!(ins, OperatorInsertion(:+, Int(j), adag))
    else
        # lesser = ddag·d（a†a）
        push!(ins, OperatorInsertion(:-, Int(j), adag))
        push!(ins, OperatorInsertion(:+, Int(i), a))
    end
    return _inserted_expectation(lattice, model, Is, Z, ins;
                                 order=Int(order), coupling_norm=coupling_norm, kernel=kernel,
                                 initial_state=get(kwargs, :initial_state, nothing))
end

"""
    matsubara(lattice, model, kernel, Is...; i, j=0, ...) -> Complex

虚时格林函数 G(τᵢ − τⱼ)：d@(τ,i)、ddag@(τ,j) 均 :left。τ 链等时
G(β,0) = +⟨d†d⟩（热环等时同号）。
"""
function matsubara(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian,
                   kernel::Union{Nothing,AbstractHybridizationKernel}, Is...; i::Integer, j::Integer=0, Z,
                   order::Integer=1, coupling_norm::Real=0.0,
                   initial_state=nothing)
    # τ↔τ 插入在混合（Kadanoff）轮廓同样合法（虚时支 = 热态构建腿）
    lattice isa AbstractImagGrassmannLattice || lattice isa AbstractMixedGrassmannLattice ||
        throw(ArgumentError("matsubara 需虚时或混合轮廓"))
    bands = lattice.bands
    ins = [OperatorInsertion(:τ, Int(i), _single_band_op(:a, 1, bands)),
           OperatorInsertion(:τ, Int(j), _single_band_op(:adag, 1, bands))]
    # Matsubara 时序约定（U1GTEMPO2 matsubara 同款）：G(τᵢ,τⱼ) =
    # (i ≥ j ? −1 : +1)·⟨d(τᵢ) ddag(τⱼ)⟩（τ 有序排列 ⟨T_τ dd†⟩ 的分支号）
    return (i >= j ? -1 : 1) * _inserted_expectation(lattice, model, Is, Z, ins; order=Int(order),
                                                     coupling_norm=coupling_norm, kernel=kernel, initial_state=initial_state)
end

"""
    nn(lattice, model, kernel, Is...; i, j, b1, b2, band, order, coupling_norm,
       initial_state, Z) -> Complex

密度-密度关联 ⟨n̂(b1,i; band1) n̂(b2,j; band2)⟩（非连通），inserted-K 双插入
精确版；同切位同带（n̂² = n̂）退化为单插入。`band` = Int 或 (band1, band2)。
"""
function nn(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian,
            kernel::Union{Nothing,AbstractHybridizationKernel}, Is...; i::Integer, j::Integer,
            b1::Symbol=default_branch(lattice), b2::Symbol=b1,
            band::Union{Int,Tuple{Int,Int}}=(1, 1), order::Integer=1,
            coupling_norm::Real=0.0, initial_state=nothing, Z)
    band1, band2 = band isa Int ? (band, band) : band
    if i == j && b1 === b2 && band1 == band2
        ins = [OperatorInsertion(b1, Int(i), _single_band_op(:n, band1, lattice.bands))]
    else
        ins = [OperatorInsertion(b1, Int(i), _single_band_op(:n, band1, lattice.bands)),
               OperatorInsertion(b2, Int(j), _single_band_op(:n, band2, lattice.bands))]
    end
    return _inserted_expectation(lattice, model, Is, Z, ins; order=Int(order),
                                 coupling_norm=coupling_norm, kernel=kernel, initial_state=initial_state)
end

# kernel 常数取数点（E4′ 单一取数：kernel 的 Z 常数按类型取——Frozen
# `.coupling_norm` = Σα、Exact `.normdelta` = 带号 Z 常数；nothing 路径 = 0）。
# 注意：sysdynamics 的常数消费不再经此默认值（测量层转发 `kernel` 本体、按型
# 分派 `_constant_lscale`——本函数供管线 Z 常数取数用）
kernel_coupling_norm(k::FrozenHybridizationKernel) = k.coupling_norm
kernel_coupling_norm(k::ExactHybridizationKernel) = k.normdelta
kernel_coupling_norm(::Nothing) = 0.0

"电流 η 系数的块读取（band 键化；η⁺⁻ = −块值——电流算符的 η = 无号物理块，
与指数系数差 Keldysh 基一个相对号，公式层常数，见 correlationfunction.jl 注）"
function _kernel_eta(kernel::AbstractHybridizationKernel, b::Symbol, band::Int)
    band in kernel.bands || throw(BoundsError(kernel.bands, band))
    return kernel.blocks[(ContourBlockIndex(:+, band, true), ContourBlockIndex(b, band, false))]
end

# ---------------------------------------------------------------------------
# Currents（Phase 5c 精确版：GF 求和公式 + 带荷双插入；旧 currents.jl 同公式）
# ---------------------------------------------------------------------------

"""
    raw_lead_correlator(lattice, model, kernel, Is...; k, j, b2=:+, band=1, ...) -> Complex

⟨ā(+,k; band) a(b2,j; band)⟩——current 公式的原始 GF（带荷双插入，荷从符号
派生；约定与 golden 验证过的 greater/lesser 同款：:left 乘）。
"""
function raw_lead_correlator(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian,
                             kernel::Union{Nothing,AbstractHybridizationKernel}, Is...; k::Integer, j::Integer,
                             b2::Symbol=:+, band::Int=1, order::Integer=1,
                             coupling_norm::Real=0.0, initial_state=nothing, Z)
    bands = lattice.bands
    # push 序 = 左乘序：同切位 (k == j) 电流取 a†a 0⁺ 极限 ⟹ ddag 先 push
    ins = [OperatorInsertion(:+, Int(k), _single_band_op(:adag, band, bands)),
           OperatorInsertion(b2, Int(j), _single_band_op(:a, band, bands))]
    return _inserted_expectation(lattice, model, Is, Z, ins; order=Int(order),
                                 coupling_norm=coupling_norm, kernel=kernel, initial_state=initial_state)
end

"""
    electriccurrent(lattice, model, kernel, Is...; k, band=1, Z, ...) -> Complex

边界 t_k 的电流（GF 求和精确版）：`(2/δt)·Σ_{j<k}[η⁺⁺[k,j]·⟨ā(+,k)a(+,j)⟩ +
η⁺⁻[k,j]·⟨ā(+,k)a(−,j)⟩]`，η⁺⁺ = kernel 块直读（raw）、η⁺⁻ =
−blocks[(:+,:−)]。⚠️ η 块禁走顶点移位 accessor（旧栈 ED 电流对比 4/8 vs 8/8
教训在案，旧 currents.jl 注释原文）。k 固定、j 变化的求和形态正是制备缓存的
摊销区。
"""
function electriccurrent(lattice::AbstractRealGrassmannLattice, model::AbstractImpurityHamiltonian,
                         kernel::Union{Nothing,AbstractHybridizationKernel}, Is...; k::Integer, band::Int=1, Z,
                         order::Integer=1, coupling_norm::Real=0.0,
                         initial_state=nothing)
    ηpp = Matrix(_kernel_eta(kernel, :+, band))   # 循环底板一次性物化（惰性视图 → 稠密）
    ηpm = -Matrix(_kernel_eta(kernel, :-, band))
    curr = 0.0im
    for j in 1:k-1
        curr += ηpp[k, j] * raw_lead_correlator(lattice, model, kernel, Is...; k=k, j=j, b2=:+,
                                                band=band, Z=Z, order=Int(order),
                                                coupling_norm=coupling_norm, initial_state=initial_state)
        curr += ηpm[k, j] * raw_lead_correlator(lattice, model, kernel, Is...; k=k, j=j, b2=:-,
                                                band=band, Z=Z, order=Int(order),
                                                coupling_norm=coupling_norm, initial_state=initial_state)
    end
    return 2 * curr / lattice.δt
end

"""
    electriccurrent_fast(lattice, model, kernel, Is...; k, band=1, Z, ...) -> Complex

electriccurrent 的流算符 GMPS 快速版（备份 electriccurrent_fast 同型）：流算符
本身构造成 GMPS（current_operator，胞坐标 0⁺ 单项式直和），`J_k ∧ K′` 一次
分子收缩 + 单发 integrate 替代逐 (k,j) 的 K′ 重建 + 积分。⚠️ 算符 = 胞坐标
0⁺ 语义（备份 build_current_mpo 同款），与边界插入的 electriccurrent 差
O(胞内接触项)——边界语义锚定路径仍是 electriccurrent；两轨快速版
（cached_electriccurrent_fast，同一算符走快照缓存）值等价。
"""
function electriccurrent_fast(lattice::AbstractRealGrassmannLattice, model::AbstractImpurityHamiltonian,
                              kernel::Union{Nothing,AbstractHybridizationKernel}, Is...; k::Integer, band::Int=1, Z,
                              order::Integer=1, coupling_norm::Real=0.0,
                              initial_state=nothing)
    Int(k) == 1 && return 0.0im
    op = current_operator(lattice, kernel, Int(k); band)
    K = sysdynamics(lattice, model; order=Int(order), initial_state, coupling_norm, kernel=kernel)
    W = op * K
    num = integrate(W, Is...)
    return 2 * num / (Z * lattice.δt)
end

"浴谱 J(ω) → ωJ(ω)（热流权重；DiracDelta = α→ωα，一般有界谱经 similar 加权）"
function _weighted_bath(bath)
    spec = bath.spectrum
    spec isa DiracDelta && return similar(bath, f=similar(spec, α=spec.ω * spec.α))
    # DiscreteSpectrum 无函数形式 similar（离散点线性插值容器）：权重逐点乘 ω
    spec isa DiscreteSpectrum && return similar(bath, f=similar(spec, spec.ws .* spec.ws, spec.ws .* spec.fs))
    return similar(bath, f=similar(spec, w -> w * spec.f(w)))
end

"""
    heatcurrent(lattice, model, bath, Is...; k, band=1, Z, order, coupling_norm,
                initial_state) -> Complex

热流 = ω 加权谱的电流（heatcorrelationfunction 同语义：J(ω) → ωJ(ω)）。
η 系数取自加权核；关联函数按调用方原 K 的 order/coupling_norm 重建。
"""
function heatcurrent(lattice::AbstractRealGrassmannLattice, model::AbstractImpurityHamiltonian,
                     bath::AbstractNormalBath, Is...; k::Integer, band::Int=1, Z,
                     order::Integer=1, coupling_norm::Real=0.0, initial_state=nothing)
    kernel_η = correlationfunction(_weighted_bath(bath), lattice)
    ηpp = Matrix(_kernel_eta(kernel_η, :+, band))   # 循环底板一次性物化（惰性视图 → 稠密）
    ηpm = -Matrix(_kernel_eta(kernel_η, :-, band))
    curr = 0.0im
    for j in 1:k-1
        curr += ηpp[k, j] * raw_lead_correlator(lattice, model, nothing, Is...; k=k, j=j, b2=:+,
                                                band=band, Z=Z, order=Int(order),
                                                coupling_norm=coupling_norm, initial_state=initial_state)
        curr += ηpm[k, j] * raw_lead_correlator(lattice, model, nothing, Is...; k=k, j=j, b2=:-,
                                                band=band, Z=Z, order=Int(order),
                                                coupling_norm=coupling_norm, initial_state=initial_state)
    end
    return 2 * curr / lattice.δt
end
