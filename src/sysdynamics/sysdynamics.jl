# sysdynamics.jl：杂质系统因子 K 构造（槽算符/插入/切位/存储序 + 公开入口）
#
# 公开 sysdynamics 全权委托 fock_sysdynamics +
# fockmpo2gmps（本目录 fockmpo.jl；plan/sysdynamics_fock_mpo_plan.md）：
#   Fock 直积 MPO 开链装配 → 切位 ≥ 1 插入稀疏乘法 → 转换（ρ₀/切位 0
#   三明治 + prefix + oddlater + 迹闭合 + split_cell/pairize 分站 + reorder）。
# 旧径（U1GTEMPO2 移植的 dense-core ring automaton + SystemFactorPreparation
# 两层缓存 + backend 开关）已于 2026-11-18 随 V2 删除（对拍闸门 203 项绿后）。
#
# 本文件保留的机器（新径与测量层共用）：
#   - local_weights：槽算符积展开（q=1 anti-Wick / q=1/2 Weyl；W₄ gaussian）；
#   - OperatorInsertion：插入算符（结构化 ==/hash）；
#   - fock_charges / _cell_step2 / _model_propagator / insertion_cut_index /
#     _storage_cells / _cell_tensors（split_cell 分站装置）。
#
# 语义（INTEGRATION.md §4 对齐，R2：以 U1GTEMPO2 为准）：
#   一阶（q=1）：mats = S·M（+支）/ M·S（−支）——Lie 顺序；
#   二阶（q=1/2，Weyl）：mats = S·M·S（对称半步），Weyl 标量
#   lscale = coupling_norm/4·Σ step²（coupling_norm = 全带总耦合 Σ_b α_b）。
#   ρ₀ 折进时序单元 1（(+,1) / τ₁），迹闭合折进末单元（(−,1) / τ_N）——
#   boundarycondition! 对 K 链不再需要（D6）。
#
# 费米符号：构造张量扇区全 bosonic，机器符号由 GM.@grassmann/GrassmannBackend
# 注入；oddlater Koszul 符号显式携带（同 U1GTEMPO2）。

# ---------------------------------------------------------------------------
# Fock 空间算符（band 1 = 最高有效位，FockMatrix 约定；JW 穿越更低索引带）
# ---------------------------------------------------------------------------

function _fock_annihilators(M::Int)
    dim = 2^M
    d = [zeros(ComplexF64, dim, dim) for _ in 1:M]
    for b in 1:M
        bit = 1 << (M - b)
        mask_cross = (dim - 1) & ~(bit - 1) & ~bit   # 更低索引带的位
        for st in 0:dim-1
            (st & bit) != 0 || continue
            d[b][(st ⊻ bit) + 1, st + 1] = iseven(count_ones(st & mask_cross)) ? 1.0 : -1.0
        end
    end
    return d
end

"""
    local_weights(spec, bands, order) -> Vector{Matrix{ComplexF64}}

槽算符积展开（4^bands 个 2^bands×2^bands 矩阵）。每带 4 槽
[q·I−n, −d†, d, W₄]：q = 1（一阶 anti-Wick）/ 1/2（Weyl）；
W₄ = n（一阶）/ q·I+n（二阶，Gaussian 测度适配）。带间 graded⊗ 符号
(−1)^{Σ_{b<b'} popcount·popcount}（U1GTEMPO2 管线约定，R2）。
槽编码 m 的每带 2bit：0=(无) 1=(ā) 2=(a) 3=(aā)，band 1 = 低位组（与
split_cell 的 charges 序 LSB 剥离一致；Fock 态矩阵基仍为 band 1 = MSB，
两者独立——m 只枚举槽积、不寻址 Fock 基）。
M≥2 的带间 graded⊗ 符号因子经 `graded_sign_factor` 注入（R2 标定点）。
"""
graded_sign_factor(bands::Int) = 1   # 默认 OFF（golden 通道裁决：对 :z2 结果不敏感）
measure_mode() = :gaussian           # W4 = qI+n（gaussian）/ I（bare）——R2 标定点

function local_weights(::Type{<:Sector}, bands::Int, order::Int)
    order in (1, 2) || throw(ArgumentError("order must be 1 or 2"))
    d = _fock_annihilators(bands)
    dim = 2^bands
    q = order == 2 ? 0.5 : 1.0
    locals = [let nb = d[b]' * d[b]
                  W4 = order == 2 ? (measure_mode() === :gaussian ? q * Matrix{ComplexF64}(I, dim, dim) + nb :
                                      Matrix{ComplexF64}(I, dim, dim)) :
                       Matrix{ComplexF64}(nb)
                  [Matrix{ComplexF64}(q * I - nb), -d[b]', d[b], W4]
              end for b in 1:bands]
    weights = Vector{Matrix{ComplexF64}}(undef, 4^bands)
    for m in 0:(4^bands)-1
        # m 的每带 2bit 组：band 1 = 低位组（U1GTEMPO2 m = Σ slot_b·4^{b−1} 同款
        # 布局）——split_cell 按 charges 序 LSB→MSB 剥离物理位，charges 以 band 1
        # 在前，两组约定必须一致，否则 bands≥2 时带内容互换（M=1 平凡不敏感）
        slots = [(m >> (2 * (b - 1))) & 3 for b in 1:bands]
        O = Matrix{ComplexF64}(I, dim, dim)
        sign = 1
        for b1 in 1:bands, b2 in (b1+1):bands
            sign *= (-1)^(count_ones(slots[b1]) * count_ones(slots[b2]))
        end
        # R2/T1 标定：graded⊗ 符号按带数启停（M=1 不敏感恒 1；M≥2 的标定
        # 由 golden 通道裁决——见 local_weights docstring）
        sign *= graded_sign_factor(bands)
        for b in 1:bands
            O = O * locals[b][slots[b] + 1]
        end
        weights[m+1] = sign * O
    end
    return weights
end

"Fock 基态（0..dim−1）的荷：占据带的 var_charge 之和"
function fock_charges(::Type{G}, bands::Int) where {G<:Sector}
    dim = 2^bands
    z = one(G)
    out = Vector{typeof(z)}(undef, dim)
    for st in 0:dim-1
        q = z
        for b in 1:bands
            ((st >> (bands - b)) & 1) == 1 && (q = qadd(q, var_charge(G, b, bands)))
        end
        out[st + 1] = q
    end
    return out
end

# ---------------------------------------------------------------------------
# 插入算符（测量与 ρ₀；Phase 5 消费）
# ---------------------------------------------------------------------------

"""
    OperatorInsertion(branch, boundary, matrix; charge=nothing, side=:left)

物理算符插入：在 `branch` 分支的边界 `boundary`（(+,0)/τ 0 = ρ₀ 切位）乘进链。
`matrix` 为 2^bands×2^bands Fock 矩阵；`charge` = 算符的荷标签（缺省 = 零荷）。
`side`（secondorderif WeylInsertion 同语义）：`:left` 乘在单元左（O·M，
轮廓更晚读法）、`:right` 乘在单元右（M·O）。同切位多插入的 0⁺ 语义由 push
序承载（先 push = 轮廓更早）。
"""
struct OperatorInsertion
    branch::Symbol
    boundary::Int
    matrix::Matrix{ComplexF64}
    charge::Any
    side::Symbol
    function OperatorInsertion(branch::Symbol, boundary::Int, matrix::AbstractMatrix;
                               charge=nothing, side::Symbol=:left)
        dim = size(matrix, 1)
        (dim == size(matrix, 2) && dim <= 16 && ispow2(dim)) ||
            throw(ArgumentError("operator must be a square 2^bands×2^bands matrix"))
        all(isfinite, matrix) || throw(ArgumentError("finite insertion matrix required"))
        side in (:left, :right) || throw(ArgumentError("side must be :left or :right"))
        new(branch, boundary, Matrix{ComplexF64}(matrix), charge, side)
    end
end

# 结构化相等/hash（制备对象格张量缓存键的先决；用户每次测量调用新构造插入
# 对象，=== 不成立）。字段逐项 + 矩阵元素相等（dim ≤ 16，成本可忽略）；
# 矩阵 hash 按列主序线性迭代，与 == 的元素序一致。charge 用 isequal 保持
# hash/== 一致（nothing / Sector 值类型）。
function Base.:(==)(a::OperatorInsertion, b::OperatorInsertion)
    a.branch === b.branch && a.boundary == b.boundary && a.side === b.side &&
        isequal(a.charge, b.charge) && a.matrix == b.matrix
end

function Base.hash(x::OperatorInsertion, h::UInt)
    h = hash(x.branch, h)
    h = hash(x.boundary, h)
    h = hash(x.side, h)
    h = hash(x.charge, h)
    for v in x.matrix
        h = hash(v, h)
    end
    return h
end

# ---------------------------------------------------------------------------
# 共享装配工具（split_cell 分站装置 + 传播子/切位/存储序助手）
# ---------------------------------------------------------------------------

"单元数据 → 变量粒度裸张量列（split_cell 逐格 SVD 切分；无数值截断）"
function _cell_tensors(spec, core, leftlabels, rightlabels, charges)
    return split_cell(spec, core, leftlabels, rightlabels, charges)
end

"实时格 Δt²、τ 格 δτ²（Weyl 标量的 Σ step²）"
_cell_step2(lattice, b) = (b === :τ ? lattice.δτ^2 : -lattice.δt^2)

"模型传播子（一阶全步 / 二阶半步；Td 模型带区间左端点 t）"
function _model_propagator(model, b::Symbol, dt::Real, t::Union{Nothing,Real})
    if model isa ConstImpurityHamiltonian
        return Matrix(fock_propagator(model, b, dt).data)
    elseif model isa AbstractTdImpurityHamiltonian
        return Matrix(fock_propagator(model, b, dt, t === nothing ? 0.0 : t).data)
    else
        throw(ArgumentError("unsupported impurity Hamiltonian $(typeof(model))"))
    end
end

"""
    insertion_cut_index(lattice, ins) -> Int

插入切位 → 时序格序索引（U1GTEMPO2 insertion_cut 同式；cut k 的插入左乘进
时序格 k 的 mats，cut 0 = ρ₀ 边界）：τ j → j；+ j → nimag + j；
− j → nimag + 2nreal − j。boundary ∈ 0..格数（该分支）。
"""
function insertion_cut_index(lattice::AbstractGrassmannLattice, ins::OperatorInsertion)
    ni = lattice isa AbstractRealGrassmannLattice ? 0 :
         lattice isa AbstractImagGrassmannLattice ? lattice.Nτ : lattice.Nτ
    nt = lattice isa AbstractRealGrassmannLattice ? lattice.Nt : 
         lattice isa AbstractMixedGrassmannLattice ? lattice.Nt : 0
    ins.branch === :τ && return ins.boundary
    ins.branch === :+ && return ni + ins.boundary
    ins.branch === :- && return ni + 2nt - ins.boundary
    throw(ArgumentError("branch must be :τ, :+ or :-"))
end

"存储序格表（lattice _index 的块序反查）：[τ Nτ…1][+ 1..Nt][− Nt…1] 等"
function _storage_cells(lattice::AbstractGrassmannLattice)
    cells = Tuple{Symbol,Int}[]
    lattice isa AbstractImagGrassmannLattice && append!(cells, [(:τ, j) for j in lattice.kτ:-1:1])
    if lattice isa AbstractRealGrassmannLattice
        append!(cells, [(:+, j) for j in 1:lattice.k])
        append!(cells, [(:-, j) for j in lattice.k:-1:1])
    end
    if lattice isa AbstractMixedGrassmannLattice
        append!(cells, [(:τ, j) for j in 1:lattice.Nτ])
        append!(cells, [(:+, j) for j in 1:lattice.Nt])
        append!(cells, [(:-, j) for j in lattice.Nt:-1:1])
    end
    return cells
end


# ---------------------------------------------------------------------------
# 公开皮肤（旧 API 名保持；plan §3.3；2026-11-18 起全权委托 fock_sysdynamics
# + fockmpo2gmps——旧径 SystemFactorPreparation/_system_factor_reference 与
# backend 开关随 V2 删除）
# ---------------------------------------------------------------------------

"""
    sysdynamics(lattice, model; order=1, initial_state=nothing,
                coupling_norm=0.0, insertions=OperatorInsertion[], trunc) -> GMPS

杂质系统因子 K（order=1 相干 anti-Wick 语义，q=1 局部权重）。实现 =
`fock_sysdynamics`（Fock 直积 MPO 开链装配）→ 切位 ≥ 1 插入稀疏乘法 →
`fockmpo2gmps`（ρ₀/切位 0 三明治 + prefix + oddlater + 迹闭合 +
split_cell/pairize 分站 + reorder）。ρ₀/迹闭合已折进单元
（boundarycondition! 不再需要）；`initial_state`（FockMatrix/矩阵）仅实时
轮廓生效（缺省 = 真空投影）。
"""
function sysdynamics(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian;
                     order::Integer=1, kwargs...)
    order in (1, 2) || throw(ArgumentError("order must be 1 or 2"))
    return _system_factor_fock(lattice, model; order=Int(order),
                               spec=TK.sectortype(lattice), kwargs...)
end
