# inserted_k.jl：inserted-K 边界语义观测量（D5 两轨之一；secondorderif
# measurement.jl 同约定）。
#
# 全部经「带插入重建 K → GM 引擎 integrate → /Z」：
#   occupation(branch j; band)      n̂（逐带或总）在 (branch, j) 边界 :left
#   nn(i, j; b1, b2, band)          n̂ 双插入密度-密度（ED 验证 ✓）
#   greater/lesser/matsubara        ✓ 已修复并对拍 U1GTEMPO2 机器精度（2026-11-15）。
#                                   五处口径差全部定案：①builder prefix 必须以
#                                   cut-0 插入荷为初值并逐胞累计（sysdynamics.jl；奇
#                                   插入链迹归零的真凶——偶算符 charge=nothing
#                                   不受影响）；②荷方向 = 算符自荷（d=−q、
#                                   ddag=+q，correlation_insertions 同款，此前全
#                                   反向）；③push 序 = cutleft foldr 首 push 最左
#                                   （同切位 0⁺：greater = d·ddag、lesser/电流 =
#                                   ddag·d）；④matsubara 时序因子 (i≥j ? −1 : +1)
#                                   且混合轮廓 τ 支放行；⑤order/coupling_norm 显
#                                   式 kwargs 必须转发（_normalized_two_point 曾
#                                   忽略 → 以 order=1 链对 order=2 的 Z 归一）。
#                                   双栈对拍：实时 greater/lesser 3 点、虚时
#                                   matsubara 4 点、occupation 4 点、lead 关联
#                                   6 点（含同切位）全部 Δ ≤ 5e-15。
#   electriccurrent/heatcurrent     GF 求和精确版；raw_lead_correlator 同上锚定
#                                   （平衡态净电流离散伪影 O(δt)，旧栈同款量级）。
# 边界 j ∈ 0..N（0 = ρ₀ 切位）；rt 顶点 branch=:+、τ 顶点 branch=:τ。
# τ 链 G(β,0) = +⟨d†d⟩（热环等时同号；j==0 平移路由属 secondorderif M≥2 的
# Z2 栈绕尾技巧，GM 栈 F≡0 前提下不需要——D5 验收时按 golden 裁决）。

"""
    _fock_ops(bands) -> Vector{Matrix{ComplexF64}}

各带 JW 湮灭算符（band 1 = MSB；FockMatrix 约定）。
"""
_fock_ops(bands::Int) = _fock_annihilators(bands)

"n̂ 算符（band=nothing → 全带总占据）"
function _number_op(bands::Int, band::Union{Nothing,Int})
    d = _fock_ops(bands)
    if band === nothing
        acc = d[1]' * d[1]
        for b in 2:bands; acc += d[b]' * d[b]; end
        return acc
    end
    (1 <= band <= bands) || throw(BoundsError(1:bands, band))
    return d[band]' * d[band]
end

"插入 K 重建的公共尾程：K′ = sysdynamics(...; insertions) → integrate(K′, Is...) / Z"
function _inserted_expectation(lattice, model, Is, Z, insertions; order=1, coupling_norm=0.0,
                               initial_state=nothing, spec=nothing)
    K = sysdynamics(lattice, model; order, initial_state, coupling_norm, insertions,
                    spec=spec === nothing ? TK.sectortype(lattice) : spec)
    num = integrate(K, Is...)
    return num / Z
end

"""
    occupation(lattice, model, kernel, Is...; j, branch, band=nothing,
               order=1, coupling_norm=0.0, initial_state=nothing, Z) -> Complex

⟨n̂(t_j)⟩（`band` 缺省 = 总占据；j ≥ 1 在 (branch, j) 边界，j = 0 折 ρ₀/
热迹）。`Is` = 各带 IF 链；`Z` = 无插入分母。
"""
function occupation(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian,
                    kernel::Union{Nothing,HybridizationKernel}, Is...; j::Integer, branch::Symbol=default_branch(lattice),
                    band::Union{Nothing,Int}=nothing, order::Integer=kernel_order(kernel),
                    coupling_norm::Real=kernel_coupling_norm(kernel), initial_state=nothing, Z)
    ins = [OperatorInsertion(branch, Int(j), _number_op(lattice.bands, band))]
    return _inserted_expectation(lattice, model, Is, Z, ins; order=Int(order),
                                 coupling_norm=coupling_norm, initial_state=initial_state)
end

"缺省分支：虚时 :τ、实时/混合 rt 顶点 :+"
default_branch(lattice::AbstractGrassmannLattice) =
    lattice isa AbstractImagGrassmannLattice ? :τ : :+

"""
    raw_correlation(lattice, model, kernel, Is...; i, j=0, kind=:greater, ...) -> Complex

两点函数 ⟨d(tᵢ) ddag(tⱼ)⟩（greater）/ ⟨ddag(tⱼ) d(tᵢ)⟩（lesser）；实时/混合
专用（虚时用 matsubara）。0⁺ push 序约定见文件头。⚠️ 奇算符插入必须显式带荷
（charge=var_charge/qneg）：charge=nothing 时 builder 的 prefix 不移位，奇单项
式的电荷流不可表示 → K' 切分产出零站 → 链归零（偶算符 n̂ 不受影响）。
"""
function raw_correlation(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian,
                         kernel::Union{Nothing,HybridizationKernel}, Is...; i::Integer, j::Integer=0,
                         kind::Symbol=:greater, order::Integer=kernel_order(kernel),
                         coupling_norm::Real=kernel_coupling_norm(kernel), initial_state=nothing)
    (lattice isa AbstractRealGrassmannLattice || lattice isa AbstractMixedGrassmannLattice) ||
        throw(ArgumentError("实时/混合轮廓专用；纯虚时用 matsubara"))
    spec = TK.sectortype(lattice)
    qc = var_charge(spec, 1, lattice.bands)
    d = _fock_ops(lattice.bands)[1]; ddag = d'
    ins = OperatorInsertion[]
    if kind === :greater
        # ⟨d(tᵢ) ddag(tⱼ)⟩：荷 = 算符自荷（U1GTEMPO2 correlation_insertions 同款
        # ——d 湮灭 = −q、ddag 产生 = +q）。push 序 = 算符左乘序（cutleft foldr
        # 首 push 最左、最后作用）：同切位 0⁺ greater = d·ddag ⟹ d 先 push
        push!(ins, OperatorInsertion(:+, Int(i), d; charge=qneg(qc)))
        push!(ins, OperatorInsertion(:+, Int(j), ddag; charge=qc))
    elseif kind === :lesser
        # ⟨ddag(tⱼ) d(tᵢ)⟩：同切位 0⁺ = ddag·d（a†a）⟹ ddag 先 push
        push!(ins, OperatorInsertion(:-, Int(j), ddag; charge=qc))
        push!(ins, OperatorInsertion(:+, Int(i), d; charge=qneg(qc)))
    else
        throw(ArgumentError("kind must be :greater or :lesser"))
    end
    return _inserted_expectation(lattice, model, Is, 1.0, ins; order=Int(order),
                                 coupling_norm=coupling_norm, initial_state=initial_state)
end

"greater(tᵢ, tⱼ=0) = −i⟨d(tᵢ)ddag(tⱼ)⟩（归一化走 occupation 同款 Z 分母）"
function greater(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian,
                 kernel::Union{Nothing,HybridizationKernel}, Is...; i::Integer, j::Integer=0, Z,
                 kwargs...)
    (lattice isa AbstractRealGrassmannLattice || lattice isa AbstractMixedGrassmannLattice) ||
        throw(ArgumentError("虚时轮廓用 matsubara"))
    return -1.0im * _normalized_two_point(lattice, model, kernel, Is, Z, i, j, :greater; kwargs...)
end

"lesser(tᵢ, tⱼ=0) = +i⟨ddag(tⱼ)d(tᵢ)⟩"
function lesser(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian,
                kernel::Union{Nothing,HybridizationKernel}, Is...; i::Integer, j::Integer=0, Z,
                kwargs...)
    (lattice isa AbstractRealGrassmannLattice || lattice isa AbstractMixedGrassmannLattice) ||
        throw(ArgumentError("虚时轮廓用 matsubara"))
    return 1.0im * _normalized_two_point(lattice, model, kernel, Is, Z, i, j, :lesser; kwargs...)
end

"两点归一化公共体（greater/lesser 宿主）：order/coupling_norm 显式转发
（缺省随 kernel 元数据；hybridization 的 diagnostics 未登记 order 时返回
1——**golden 教训：忽略调用方显式值会以 order=1 链对 order=2 的 Z 归一**）"
function _normalized_two_point(lattice, model, kernel, Is, Z, i, j, kind;
                               order::Integer=kernel_order(kernel),
                               coupling_norm::Real=kernel_coupling_norm(kernel), kwargs...)
    spec = TK.sectortype(lattice)
    qc = var_charge(spec, 1, lattice.bands)
    d = _fock_ops(lattice.bands)[1]; ddag = d'
    ins = OperatorInsertion[]
    if kind === :greater
        # push 序 = 左乘序（同切位 0⁺ greater = d·ddag）
        push!(ins, OperatorInsertion(:+, Int(i), d; charge=qneg(qc)))
        push!(ins, OperatorInsertion(:+, Int(j), ddag; charge=qc))
    else
        # lesser = ddag·d（a†a）
        push!(ins, OperatorInsertion(:-, Int(j), ddag; charge=qc))
        push!(ins, OperatorInsertion(:+, Int(i), d; charge=qneg(qc)))
    end
    return _inserted_expectation(lattice, model, Is, Z, ins;
                                 order=Int(order), coupling_norm=coupling_norm,
                                 initial_state=get(kwargs, :initial_state, nothing))
end

"""
    matsubara(lattice, model, kernel, Is...; i, j=0, ...) -> Complex

虚时格林函数 G(τᵢ − τⱼ)：d@(τ,i)、ddag@(τ,j) 均 :left。τ 链等时
G(β,0) = +⟨d†d⟩（热环等时同号）。
"""
function matsubara(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian,
                   kernel::Union{Nothing,HybridizationKernel}, Is...; i::Integer, j::Integer=0, Z,
                   order::Integer=kernel_order(kernel), coupling_norm::Real=kernel_coupling_norm(kernel),
                   initial_state=nothing)
    # τ↔τ 插入在混合（Kadanoff）轮廓同样合法（虚时支 = 热态构建腿）
    lattice isa AbstractImagGrassmannLattice || lattice isa AbstractMixedGrassmannLattice ||
        throw(ArgumentError("matsubara 需虚时或混合轮廓"))
    d = _fock_ops(lattice.bands)[1]; ddag = d'
    qc = var_charge(TK.sectortype(lattice), 1, lattice.bands)
    ins = [OperatorInsertion(:τ, Int(i), d; charge=qneg(qc)),
           OperatorInsertion(:τ, Int(j), ddag; charge=qc)]
    # Matsubara 时序约定（U1GTEMPO2 matsubara 同款）：G(τᵢ,τⱼ) =
    # (i ≥ j ? −1 : +1)·⟨d(τᵢ) ddag(τⱼ)⟩（τ 有序排列 ⟨T_τ dd†⟩ 的分支号）
    return (i >= j ? -1 : 1) * _inserted_expectation(lattice, model, Is, Z, ins; order=Int(order),
                                                     coupling_norm=coupling_norm, initial_state=initial_state)
end

"""
    nn(lattice, model, kernel, Is...; i, j, b1, b2, band, order, coupling_norm,
       initial_state, Z) -> Complex

密度-密度关联 ⟨n̂(b1,i; band1) n̂(b2,j; band2)⟩（非连通），inserted-K 双插入
精确版；同切位同带（n̂² = n̂）退化为单插入。`band` = Int 或 (band1, band2)。
"""
function nn(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian,
            kernel::Union{Nothing,HybridizationKernel}, Is...; i::Integer, j::Integer,
            b1::Symbol=default_branch(lattice), b2::Symbol=b1,
            band::Union{Int,Tuple{Int,Int}}=(1, 1), order::Integer=kernel_order(kernel),
            coupling_norm::Real=kernel_coupling_norm(kernel), initial_state=nothing, Z)
    band1, band2 = band isa Int ? (band, band) : band
    if i == j && b1 === b2 && band1 == band2
        ins = [OperatorInsertion(b1, Int(i), _number_op(lattice.bands, band1))]
    else
        ins = [OperatorInsertion(b1, Int(i), _number_op(lattice.bands, band1)),
               OperatorInsertion(b2, Int(j), _number_op(lattice.bands, band2))]
    end
    return _inserted_expectation(lattice, model, Is, Z, ins; order=Int(order),
                                 coupling_norm=coupling_norm, initial_state=initial_state)
end

# kernel 元数据便捷（耦合范数/order 随 kernel 携带）；nothing（显式 order/
# coupling_norm 路径）安全——kwarg 默认值在 kernel=nothing 时也求值
kernel_coupling_norm(k::HybridizationKernel) = k.coupling_norm === nothing ? 0.0 : k.coupling_norm
kernel_coupling_norm(::Nothing) = 0.0
kernel_order(k::HybridizationKernel) = get(k.diagnostics, :order, 1)
kernel_order(::Nothing) = 1

# ---------------------------------------------------------------------------
# Currents（Phase 5c 精确版：GF 求和公式 + 带荷双插入；旧 currents.jl 同公式）
# ---------------------------------------------------------------------------

"""
    raw_lead_correlator(lattice, model, kernel, Is...; k, j, b2=:+, band=1, ...) -> Complex

⟨ā(+,k; band) a(b2,j; band)⟩——current 公式的原始 GF（带荷双插入，约定与
golden 验证过的 greater/lesser 同款：charge=nothing、:left 乘）。
"""
function raw_lead_correlator(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian,
                             kernel::Union{Nothing,HybridizationKernel}, Is...; k::Integer, j::Integer,
                             b2::Symbol=:+, band::Int=1, order::Integer=kernel_order(kernel),
                             coupling_norm::Real=kernel_coupling_norm(kernel), initial_state=nothing, Z)
    qc = var_charge(TK.sectortype(lattice), band, lattice.bands)
    d = _fock_ops(lattice.bands)[band]
    ddag = d'
    # push 序 = 左乘序：同切位 (k == j) 电流取 a†a 0⁺ 极限 ⟹ ddag 先 push
    ins = [OperatorInsertion(:+, Int(k), ddag; charge=qc),
           OperatorInsertion(b2, Int(j), d; charge=qneg(qc))]
    return _inserted_expectation(lattice, model, Is, Z, ins; order=Int(order),
                                 coupling_norm=coupling_norm, initial_state=initial_state)
end

"""
    electriccurrent(lattice, model, kernel, Is...; k, band=1, Z, ...) -> Complex

边界 t_k 的电流（GF 求和精确版）：`(2/δt)·Σ_{j<k}[η⁺⁺[k,j]·⟨ā(+,k)a(+,j)⟩ +
η⁺⁻[k,j]·⟨ā(+,k)a(−,j)⟩]`，η⁺⁺ = kernel.blocks[(:+,:+)]（raw）、η⁺⁻ =
−blocks[(:+,:−)]。⚠️ η 块禁走顶点移位 accessor（旧栈 ED 电流对比 4/8 vs 8/8
教训在案，旧 currents.jl 注释原文）。
"""
function electriccurrent(lattice::AbstractRealGrassmannLattice, model::AbstractImpurityHamiltonian,
                         kernel::Union{Nothing,HybridizationKernel}, Is...; k::Integer, band::Int=1, Z,
                         order::Integer=kernel_order(kernel), coupling_norm::Real=kernel_coupling_norm(kernel),
                         initial_state=nothing)
    ηpp = Matrix(kernel.blocks[(:+, :+)])       # 循环底板一次性物化（惰性视图 → 稠密）
    ηpm = -Matrix(kernel.blocks[(:+, :-)])
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
    kernel_η = hybridization(lattice, _weighted_bath(bath))
    ηpp = Matrix(kernel_η.blocks[(:+, :+)])     # 循环底板一次性物化（惰性视图 → 稠密）
    ηpm = -Matrix(kernel_η.blocks[(:+, :-)])
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
