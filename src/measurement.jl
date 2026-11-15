# measurement.jl：观测量（inserted-K 重建路线，D5 两轨之一；secondorderif
# measurement.jl 同约定）。
#
# 全部经「带插入重建 K → GM 引擎 integrate → /Z」：
#   occupation(branch j; band)      n̂（逐带或总）在 (branch, j) 边界 :left
#   nn(i, j; b1, b2, band)          n̂ 双插入密度-密度（ED 验证 ✓）
#   greater/lesser/matsubara        ⚠️ 未决（诊断已收窄）：奇算符（d/ddag）插
#                                   入的 K' 链迹恒零/垃圾。已排除：①builder——
#                                   忠实复刻（含 cut-0 折叠 ρ₀ + prefix 机制）
#                                   逐胞结构与 U1GTEMPO2 参照完全一致；②荷方向
#                                   （与 U1GTEMPO2 correlation_insertions 一致）；
#                                   ③宇并特异（:u1 整数荷同样失败）；④引擎 vs
#                                   多项式——独立 dense oracle 与引擎双双确认
#                                   多项式迹精确为零，而参照栈同构胞数据给出
#                                   正确值 0.66−0.07i。⇒ **决定性实验已完
#                                   成（2026-11-15）**：U1GTEMPO2 直接依赖本
#                                   GM（同 uuid），参照栈 pairize 出的奇插入
#                                   K' 链喂同一 GM.Zvalue 机器精度正确
#                                   （0.66081−0.06945i）——**GM 引擎无罪，
#                                   port 侧 split_cell/from_core/pairize 之一
#                                   为凶**（逐胞 nnz 结构与参照一致但张量内容
#                                   层有分歧；:u1 整数荷忠实链同样垃圾，非宇
#                                   并特异）。**下会议首步**：双栈各 dump cell-1
#                                   切分出的两个变量张量全块结构（扇区键+值）
#                                   逐块 diff（:u1 荷标签可直接对应他们的
#                                   ChargeCharge 首槽）。
#   electriccurrent/heatcurrent     公式正确（旧 currents.jl 同式），依赖上
#                                   述两点 GF 路线——路线修复前勿用。
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
                    spec=spec === nothing ? symmetry_spec(lattice.symmetry) : spec)
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
                         kernel::HybridizationKernel, Is...; i::Integer, j::Integer=0,
                         kind::Symbol=:greater, order::Integer=kernel_order(kernel),
                         coupling_norm::Real=kernel_coupling_norm(kernel), initial_state=nothing)
    (lattice isa AbstractRealGrassmannLattice || lattice isa AbstractMixedGrassmannLattice) ||
        throw(ArgumentError("实时/混合轮廓专用；纯虚时用 matsubara"))
    spec = symmetry_spec(lattice.symmetry)
    qc = var_charge(spec, 1, lattice.bands)
    d = _fock_ops(lattice.bands)[1]; ddag = d'
    ins = OperatorInsertion[]
    if kind === :greater
        # ⟨d(tᵢ) ddag(tⱼ)⟩：同切位 0⁺ = ⟨dd†⟩ ⟹ ddag 先 push、d 后 push
        push!(ins, OperatorInsertion(:+, Int(j), ddag; charge=qneg(qc)))
        push!(ins, OperatorInsertion(:+, Int(i), d; charge=qc))
    elseif kind === :lesser
        # ⟨ddag(tⱼ) d(tᵢ)⟩：d@(+,i) left；ddag@(−,j) right
        push!(ins, OperatorInsertion(:+, Int(i), d; charge=qc))
        push!(ins, OperatorInsertion(:-, Int(j), ddag; charge=qneg(qc), side=:right))
    else
        throw(ArgumentError("kind must be :greater or :lesser"))
    end
    return _inserted_expectation(lattice, model, Is, 1.0, ins; order=Int(order),
                                 coupling_norm=coupling_norm, initial_state=initial_state)
end

"greater(tᵢ, tⱼ=0) = −i⟨d(tᵢ)ddag(tⱼ)⟩（归一化走 occupation 同款 Z 分母）"
function greater(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian,
                 kernel::HybridizationKernel, Is...; i::Integer, j::Integer=0, Z,
                 kwargs...)
    (lattice isa AbstractRealGrassmannLattice || lattice isa AbstractMixedGrassmannLattice) ||
        throw(ArgumentError("虚时轮廓用 matsubara"))
    return -1.0im * _normalized_two_point(lattice, model, kernel, Is, Z, i, j, :greater; kwargs...)
end

"lesser(tᵢ, tⱼ=0) = +i⟨ddag(tⱼ)d(tᵢ)⟩"
function lesser(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian,
                kernel::HybridizationKernel, Is...; i::Integer, j::Integer=0, Z,
                kwargs...)
    (lattice isa AbstractRealGrassmannLattice || lattice isa AbstractMixedGrassmannLattice) ||
        throw(ArgumentError("虚时轮廓用 matsubara"))
    return 1.0im * _normalized_two_point(lattice, model, kernel, Is, Z, i, j, :lesser; kwargs...)
end

function _normalized_two_point(lattice, model, kernel, Is, Z, i, j, kind; kwargs...)
    spec = symmetry_spec(lattice.symmetry)
    qc = var_charge(spec, 1, lattice.bands)
    d = _fock_ops(lattice.bands)[1]; ddag = d'
    ins = OperatorInsertion[]
    if kind === :greater
        push!(ins, OperatorInsertion(:+, Int(j), ddag; charge=qneg(qc)))
        push!(ins, OperatorInsertion(:+, Int(i), d; charge=qc))
    else
        push!(ins, OperatorInsertion(:+, Int(i), d; charge=qc))
        push!(ins, OperatorInsertion(:-, Int(j), ddag; charge=qneg(qc), side=:right))
    end
    return _inserted_expectation(lattice, model, Is, Z, ins;
                                 order=kernel_order(kernel), coupling_norm=kernel_coupling_norm(kernel),
                                 initial_state=get(kwargs, :initial_state, nothing))
end

"""
    matsubara(lattice, model, kernel, Is...; i, j=0, ...) -> Complex

虚时格林函数 G(τᵢ − τⱼ)：d@(τ,i)、ddag@(τ,j) 均 :left。τ 链等时
G(β,0) = +⟨d†d⟩（热环等时同号）。
"""
function matsubara(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian,
                   kernel::HybridizationKernel, Is...; i::Integer, j::Integer=0, Z,
                   order::Integer=kernel_order(kernel), coupling_norm::Real=kernel_coupling_norm(kernel),
                   initial_state=nothing)
    lattice isa AbstractImagGrassmannLattice || throw(ArgumentError("虚时轮廓专用"))
    d = _fock_ops(lattice.bands)[1]; ddag = d'
    qc = var_charge(symmetry_spec(lattice.symmetry), 1, lattice.bands)
    ins = [OperatorInsertion(:τ, Int(i), d; charge=qc),
           OperatorInsertion(:τ, Int(j), ddag; charge=qneg(qc))]
    return _inserted_expectation(lattice, model, Is, Z, ins; order=Int(order),
                                 coupling_norm=coupling_norm, initial_state=initial_state)
end

# kernel 元数据便捷（耦合范数/order 随 kernel 携带）
kernel_coupling_norm(k::HybridizationKernel) = k.coupling_norm === nothing ? 0.0 : k.coupling_norm
kernel_order(k::HybridizationKernel) = get(k.diagnostics, :order, 1)

# ---------------------------------------------------------------------------
# Phase 5b：算符路线（GM ExpectationCache 快照 + SparseGMPS 单项式）——
# 旧栈 cached_gf/TwosideExpectationCache 同构：缓存一次、任意多观测量只扫
# 算符支撑；与 inserted-K 重建路线（上方）在 golden 对拍下等价（D5 两轨）。
# ---------------------------------------------------------------------------

"""
    lattice_monomial(lattice, indices::NTuple{N,GrassmannContourIndex}; coeff=1) -> GM.SparseGMPS

轮廓格点处的 Grassmann 单项式算符（旧 GTerm 同语义）：变量位 = lattice 格点
映射，表系数 = coeff·(−1)^{给定序逆序数}（升序变量基）；站张量键恒平凡、
每站一个 digit 槽置位（vac/pair/η/η̄），支撑 = 触及站的闭包。
"""
function lattice_monomial(lattice::AbstractGrassmannLattice, indices::NTuple{N,GrassmannContourIndex};
                          coeff::Number=1) where {N}
    spec = symmetry_spec(lattice.symmetry)
    G = sector_type(spec)
    positions = [lattice[ix] for ix in indices]   # 变量位（给定序）
    inv = 0
    for i in eachindex(positions), j in (i + 1):length(positions)
        positions[i] > positions[j] && (inv += 1)
    end
    c0 = coeff * (-1)^inv
    touched = sort!(unique((p + 1) ÷ 2 for p in positions))
    first_s, last_s = first(touched), last(touched)
    spaces = site_spaces(lattice)
    charges = _site_charges(spec, lattice)
    z = zero_charge(spec)
    data = Any[]; q = z   # q = 已吞并变量的累积荷（键穿线；末站须回到平凡）
    for s in first_s:last_s
        o1 = (2s - 1) in positions; o2 = (2s) in positions
        E = spaces[s]
        T = TK.zeros(ComplexF64, qspace(spec, [q]) ⊗ E ← qspace(spec, [q]))
        T[(sector(spec, q), sector(spec, z), TK.dual(sector(spec, q)))][:, 1, :] .= 1
        if o1 || o2
            basis = GM.grassmann_basis(E)
            tag = o1 ? (o2 ? basis.pair : basis.eta) : basis.etabar
            qout = qadd(qadd(q, o1 ? charges[2s - 1] : z), o2 ? charges[2s] : z)
            T = TK.zeros(ComplexF64, qspace(spec, [q]) ⊗ E ← qspace(spec, [qout]))
            T[(sector(spec, q), tag.sector, TK.dual(sector(spec, qout)))][:, tag.slot, :] .=
                s == first_s ? c0 : 1
            q = qout
        end
        push!(data, T)
    end
    q == z || throw(ArgumentError("operator monomial is not charge-neutral (open legs)"))
    return GM.SparseGMPS(convert(Vector{typeof(first(data))}, data), collect(first_s:last_s))
end

"""
    environments_cache(K, Is...) -> GM.ExpectationCache

多链联合收缩缓存（K 为第一条链 = 算符目标链，GM A/Bs 约定）；Z 与全部
cached_* 观测量共用，构建一次 O(L)，快照不可变、无限复用。
"""
environments_cache(K::GM.GMPS, Is::GM.GMPS...) = GM.environments(K, Is...)

"GM ScaledScalar → 包皮肤 ScaledComplex"
_scaled(v::GM.ScaledScalar) = ScaledComplex(v.logabs, ComplexF64(v.phase))

"缓存 Z（无插入分母；与 integrate(K, Is...) 同值、读左端 O(1)）"
cached_Z(cache::GM.ExpectationCache) = unscale(_scaled(GM.Zvalue(cache)))

"""
    cached_gf(lattice, indices, K, Is...) -> Complex

轮廓格点单项式的原始 Grassmann 关联函数 ⟨Π g⟩（已归一化）。

实现 = `op ∧ K` 稀疏楔积（只触碰算符支撑站）+ 引擎单发全迹
（`GM.Zvalue`，O(L) 左扫、无右扫补偿通道）——相对 `raw_correlation`
省掉逐点 K 重建（M=2 实测单点 ~88% 成本）。

⚠️ 语义（D5 裁决）：Chrono 格的胞变量位于 Weyl 胞两端（a = 胞左端 t_{c−1}、
ā = 胞右端 t_c），单项式关联是**胞坐标 0⁺ 点**的量，与 `occupation`/
`greater`/`lesser` 的边界矩阵插入点差 O(胞内接触项)（δt=0.15 时 ~1e-3 档，
随 δt 收敛）——边界观测量请用 inserted-K 路线（上方）。

⚠️ 未用 `GM.expectation(op, env)` 快照扫掠：其右扫对数幅值补偿存在缺口
（实测 ⟨op⟩ 差一个链常数标度 ~1e-5 档，见 GrassmannMPS berezin.jl 右扫
rightlogscales 记账），修复前本包装走精确路径。
"""
function cached_gf(lattice::AbstractGrassmannLattice, indices::NTuple{N,GrassmannContourIndex},
                   K::GM.GMPS, Is::GM.GMPS...) where {N}
    op = lattice_monomial(lattice, indices)
    W = op * K
    return unscale(_scaled(GM.Zvalue([W, Is...])))
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

"""
    cached_correlation(lattice, K, Is...; i, j=0, kind=:greater, band=1, cache) -> Complex

原始两点单项式关联（胞坐标 0⁺ 点，语义与容差见 cached_gf 的 D5 说明）：
  greater：⟨a(+,i) ā(+,j)⟩  lesser：⟨ā(−,i) a(+,j)⟩  matsubara：⟨a(τ,i) ā(τ,j)⟩
"""
function cached_correlation(lattice::AbstractGrassmannLattice, K::GM.GMPS, Is::GM.GMPS...;
                            i::Integer, j::Integer=0, kind::Symbol=:greater, band::Int=1,
                            cache::GM.ExpectationCache=GM.environments(K, Is...))
    if kind === :matsubara
        lattice isa AbstractImagGrassmannLattice || throw(ArgumentError("matsubara 需虚时轮廓"))
        idx = (GrassmannContourIndex(Int(i), :τ, band, false),
               GrassmannContourIndex(Int(j), :τ, band, true))
    elseif kind === :greater
        idx = (GrassmannContourIndex(Int(i), :+, band, false),
               GrassmannContourIndex(Int(j), :+, band, true))
    elseif kind === :lesser
        idx = (GrassmannContourIndex(Int(i), :-, band, true),
               GrassmannContourIndex(Int(j), :+, band, false))
    else
        throw(ArgumentError("kind must be :greater, :lesser or :matsubara"))
    end
    return cached_gf(lattice, idx, K, Is...; cache=cache)
end

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
    qc = var_charge(symmetry_spec(lattice.symmetry), band, lattice.bands)
    d = _fock_ops(lattice.bands)[band]
    ddag = d'
    ins = [OperatorInsertion(b2, Int(j), d; charge=qc),
           OperatorInsertion(:+, Int(k), ddag; charge=qneg(qc))]
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
                         kernel::HybridizationKernel, Is...; k::Integer, band::Int=1, Z,
                         order::Integer=kernel_order(kernel), coupling_norm::Real=kernel_coupling_norm(kernel),
                         initial_state=nothing)
    ηpp = kernel.blocks[(:+, :+)]
    ηpm = -kernel.blocks[(:+, :-)]
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
    ηpp = kernel_η.blocks[(:+, :+)]
    ηpm = -kernel_η.blocks[(:+, :-)]
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

# ---------------------------------------------------------------------------
# 公开皮肤
# ---------------------------------------------------------------------------

export occupation, greater, lesser, matsubara, default_branch, nn,
       raw_lead_correlator, electriccurrent, heatcurrent
