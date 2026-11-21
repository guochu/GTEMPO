# sysdynamics.jl：杂质系统因子 K 构造（槽算符/插入/切位/存储序 + 公开入口）
#
# 公开 sysdynamics 全权委托 fock_sysdynamics +
# fockmpo2gmps（本目录 fockmpo.jl；plan/sysdynamics_fock_mpo_plan.md +
# plan/fockmpo_dress_selfread_plan.md）：
#   fock_sysdynamics 一次装配产出**穿戴完毕**的 Fock 链（开链装配 → 切位 ≥ 1
#   插入稀疏乘法 → ρ₀/切位 0 三明治折首站 → oddlater 逐站号 → 符号荷 vs 链 δ
#   交叉自检）；fockmpo2gmps 纯翻译 + 荷账自读（槽名义荷表 + 站读数 δ——
#   不消费任何插入/初态元数据）+ 迹闭合 + split_cell/pairize 分站 + reorder。
# 旧径（U1GTEMPO2 移植的 dense-core ring automaton + SystemFactorPreparation
# 两层缓存 + backend 开关）已于 2026-11-18 随 V2 删除（对拍闸门 203 项绿后）。
#
# 本文件保留的机器（新径与测量层共用）：
#   - local_weights：槽算符积展开（q=1 anti-Wick / q=1/2 Weyl；W₄ gaussian）；
#   - OperatorInsertion：插入算符（符号 FockOp；荷按位派生；结构化 ==/hash）；
#   - fock_charges / _cell_step2 / _model_propagator / insertion_cut_index /
#     _storage_cells / _cell_tensors（split_cell 分站装置）。
# 符号算符层（LocalOp/FockOp）在 localops.jl。
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
    OperatorInsertion(branch, j, op::FockOp; side=:left)

物理算符插入（全程符号态——V2/V7）：在 `branch` 分支的第 `j` 时间步乘进链。
**j = 物理时间步数**（(+,j) 与 (−,j) 同一物理时间点；虚时 j 同口径）；j = 0 =
ρ₀ 切位（三明治），⚠️ `:-` 支 j = 0 = 链末切位（合法物理插入）——过滤按
`insertion_cut_index` 的切位映射而非 j 字段。`op` = 符号 FockOp，荷从符号按位
派生（V3——无手工荷标签）；`side`（secondorderif WeylInsertion 同语义）：
`:left` 乘在单元左（O·M，轮廓更晚读法）、`:right` 乘在单元右（M·O）。同切位
多插入的 0⁺ 语义由 push 序承载（先 push = 轮廓更早）。构造守卫：j ≥ 0；
分支格数上限（越界）在 `insertion_cut_index` 提前抛。
"""
struct OperatorInsertion
    branch::Symbol
    j::Int                       # 原 boundary：时间步数（V1）
    op::FockOp                   # 原 matrix + charge：符号态，荷按位派生（V2/V3）
    side::Symbol
    function OperatorInsertion(branch::Symbol, j::Int, op::FockOp; side::Symbol=:left)
        j >= 0 || throw(ArgumentError("insertion time step must be ≥ 0, got j = $j"))
        side in (:left, :right) || throw(ArgumentError("side must be :left or :right"))
        new(branch, j, op, side)
    end
end

# 结构化相等/hash（制备对象格张量缓存键的先决；用户每次测量调用新构造插入
# 对象，=== 不成立）。字段逐项 + FockOp 符号逐位（Vector == 逐元素）；物化同
# 矩阵不同符号形视为不同对象——缓存多算一份、数值不受影响。
function Base.:(==)(a::OperatorInsertion, b::OperatorInsertion)
    a.branch === b.branch && a.j == b.j && a.side === b.side && a.op == b.op
end

function Base.hash(x::OperatorInsertion, h::UInt)
    h = hash(x.branch, h)
    h = hash(x.j, h)
    h = hash(x.side, h)
    return hash(x.op, h)
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
− j → nimag + 2nreal − j（反向）。j ∈ 0..该分支格数（越界/轮廓缺支提前抛）。
"""
function insertion_cut_index(lattice::AbstractGrassmannLattice, ins::OperatorInsertion)
    ni = lattice isa AbstractRealGrassmannLattice ? 0 :
         lattice isa AbstractImagGrassmannLattice ? lattice.Nτ : lattice.Nτ
    nt = lattice isa AbstractRealGrassmannLattice ? lattice.Nt :
         lattice isa AbstractMixedGrassmannLattice ? lattice.Nt : 0
    if ins.branch === :τ
        ni > 0 || throw(ArgumentError("branch :τ requires an imag or mixed contour lattice"))
        0 <= ins.j <= ni || throw(ArgumentError("τ insertion j = $(ins.j) out of range 0..$ni"))
        return ins.j
    elseif ins.branch === :+
        nt > 0 || throw(ArgumentError("branch :+ requires a real or mixed contour lattice"))
        0 <= ins.j <= nt || throw(ArgumentError(":+ insertion j = $(ins.j) out of range 0..$nt"))
        return ni + ins.j
    elseif ins.branch === :-
        nt > 0 || throw(ArgumentError("branch :- requires a real or mixed contour lattice"))
        0 <= ins.j <= nt || throw(ArgumentError(":- insertion j = $(ins.j) out of range 0..$nt"))
        return ni + 2nt - ins.j
    end
    throw(ArgumentError("branch must be :τ, :+ or :-"))
end

"存储序格表（lattice _index 的块序反查，2026-11-20 统一链存储序 + v7 虚轴翻转）：
imag [τ 1…Nτ]（存储 = 时序）；real [(+,1),(−,1),(+,2),(−,2),…（± 成对、时间升序）]；
mixed [τ Nτ…1][(+,1),(−,1),…（± 成对）]"
function _storage_cells(lattice::AbstractGrassmannLattice)
    cells = Tuple{Symbol,Int}[]
    lattice isa AbstractImagGrassmannLattice && append!(cells, [(:τ, j) for j in 1:lattice.kτ])
    if lattice isa AbstractRealGrassmannLattice
        for j in 1:lattice.k
            push!(cells, (:+, j))
            push!(cells, (:-, j))
        end
    end
    if lattice isa AbstractMixedGrassmannLattice
        append!(cells, [(:τ, j) for j in lattice.Nτ:-1:1])
        for j in 1:lattice.Nt
            push!(cells, (:+, j))
            push!(cells, (:-, j))
        end
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
                coupling_norm=0.0, kernel=nothing, insertions=OperatorInsertion[],
                trunc) -> GMPS

杂质系统因子 K（order=1 相干 anti-Wick 语义，q=1 局部权重）。实现 =
`fock_sysdynamics`（穿戴内联：Fock 直积 MPO 装配 + 切位 ≥ 1 插入稀疏乘法 +
ρ₀/切位 0 三明治折首站 + oddlater + 符号荷交叉自检）→ `fockmpo2gmps`（荷账
自读 + 迹闭合 + split_cell/pairize 分站 + reorder）。ρ₀/迹闭合已折进单元
（boundarycondition! 不再需要）；`initial_state`（FockMatrix/矩阵）仅实时
轮廓生效（缺省 = 真空投影）。

**Z 常数消费**（exactif_port_plan E6）：`kernel::Union{Nothing,
AbstractHybridizationKernel}` 按类型分派——`kernel === nothing` = 现行
`coupling_norm` → Weyl 公式（Frozen 缺省语义）；`kernel !== nothing` =
`_constant_lscale(kernel, lattice, order)`（Frozen = Weyl 公式，显式
`coupling_norm` 优先、0 = 未给哨兵取 kernel.coupling_norm；Exact =
log|normdelta| + 号折首站，Weyl 通道不适用）。`logscale` 累加语义（E2）不变。
"""
function sysdynamics(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian;
                     order::Integer=1, kwargs...)
    order in (1, 2) || throw(ArgumentError("order must be 1 or 2"))
    return _system_factor_fock(lattice, model; order=Int(order),
                               spec=TK.sectortype(lattice), kwargs...)
end
