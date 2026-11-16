# sysdynamics.jl：杂质系统因子 K 的 cell automaton 构造（plan §3.3/D2）
#
# 移植自 U1GTEMPO2 system.jl/prepared_system.jl（charge 通道机器精度验证），
# 泛化点（D1 对称性手选）：
#   - 荷标签/扇区映射走 bridge 的 spec 钩子（Int / NTuple{2,Int}）；
#   - S = fock_propagator(model, branch, dt[, t])（GTEMPO 模型接口直供传播子，
#     一阶全步 / 二阶半步 dt/2）；
#   - local_weight 按 M 带积展开（U1GTEMPO2 的 graded⊗ 符号约定泛化为
#     Σ_{b<b'} popcount(slot_b)·popcount(slot_b')）。
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

# ---------------------------------------------------------------------------
# 参考装配（dense core + split_cell；golden 锚定路径）
# ---------------------------------------------------------------------------

"单元数据 → 变量粒度裸张量列（split_cell 逐格 SVD 切分；无数值截断）"
function _cell_tensors(spec, core, leftlabels, rightlabels, charges)
    return split_cell(spec, core, leftlabels, rightlabels, charges)
end

"""
    _system_factor_reference(lattice, model; order=1, initial_state=nothing,
                             coupling_norm=0.0, insertions=OperatorInsertion[],
                             trunc=nothing, truncerrs=nothing, spec=TK.Z2Irrep) -> GMPS

时序 ring automaton（U1GTEMPO2 system_factor_reference 同型）：每格
(b,j) 的 core 经 split_cell 切成变量粒度张量列，全部格拼接后 pairize，
再 reorder! 到 lattice 存储序。ρ₀ = `initial_state`（FockMatrix 或
2^bands×2^bands 矩阵；实时缺省 = 真空投影 δ₀₀、τ/混合 = 热迹 I）折进单元 1；
迹闭合折进末单元。
"""
function _system_factor_reference(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian;
                                  order::Int=1, initial_state=nothing, coupling_norm::Real=0.0,
                                  insertions=OperatorInsertion[], trunc=nothing,
                                  truncerrs=nothing, spec::Type{<:Sector}=TK.Z2Irrep)
    order in (1, 2) || throw(ArgumentError("order must be 1 or 2"))
    order == 2 && !(coupling_norm >= 0) && throw(ArgumentError("order=2 requires coupling_norm ≥ 0 (total Σ_b α_b)"))
    bands = lattice.bands
    bands_guard(spec, bands)
    num_bands(model) == bands || throw(ArgumentError("model bands $(num_bands(model)) ≠ lattice bands $bands"))
    dim = 2^bands
    z = one(spec)
    cells = chronological_cells(lattice)
    L = length(cells)
    L > 0 || throw(ArgumentError("empty lattice"))
    # ρ₀：实时 = initial_state（缺省 = 真空投影 δ₀₀，旧 GTEMPO 真空初态语义）；
    # τ/混合 = 热迹 I（τ 段自身建热态，initial_state 忽略）
    rho = Matrix{ComplexF64}(I, dim, dim)
    if lattice isa AbstractRealGrassmannLattice
        rho = Matrix{ComplexF64}(Diagonal(vcat([1.0], zeros(dim - 1))))
        initial_state === nothing || (rho = initial_state isa FockMatrix ?
            Matrix{ComplexF64}(initial_state.data) : Matrix{ComplexF64}(initial_state))
    end
    # 插入切位（时序格序索引；cut 0 = ρ₀ 边界）
    cuts = Dict{Int,Vector{OperatorInsertion}}()
    for ins in insertions
        k = insertion_cut_index(lattice, ins)
        push!(get!(cuts, k, OperatorInsertion[]), ins)
    end
    total_q = z
    for ins in insertions
        total_q = qadd(total_q, ins.charge === nothing ? z : ins.charge)
    end
    total_q == z || throw(ArgumentError("total insertion charge must be zero"))
    # 切位算符（side 语义：left = O·M 折左（foldr——后 push 更晚 = 更左）；
    # right = M·O 折右（push 序））
    cutleft(k) = foldr((ins, A) -> ins.matrix * A, filter(i -> i.side === :left, get(cuts, k, OperatorInsertion[]));
                       init=Matrix{ComplexF64}(I, dim, dim))
    cutright(k) = foldl((A, ins) -> A * ins.matrix, filter(i -> i.side === :right, get(cuts, k, OperatorInsertion[]));
                        init=Matrix{ComplexF64}(I, dim, dim))
    # ρ₀ 切位（cut 0）：left = O·ρ、right = ρ·O
    rho = cutleft(0) * rho * cutright(0)
    weights = local_weights(spec, bands, order)
    fqs = fock_charges(spec, bands)
    pairs = [(s, t) for t in 1:dim for s in 1:dim]
    data = nothing
    # prefix = 已折入状态的累积插入荷（U1GTEMPO2 同款）：cut 0 的插入已折 ρ₀
    # （上行 rho = cutleft(0)·rho·cutright(0)），其荷必须计入初值——否则奇插入
    # （荷≠0）首张量右键荷不随 ρ₀ 的荷移位 → 块路由错位、链迹归零；逐胞累计
    # 而非重置（leftlabels[k] 相对全部已过切位）。偶算符 charge=nothing 不进
    # prefix，occupation/nn 路径不受影响。
    prefix = foldl((p, ins) -> ins.charge === nothing ? p : qadd(p, ins.charge),
                   get(cuts, 0, OperatorInsertion[]); init=z)
    labels = Any[]
    charges_per_cell = cell_charges(spec, bands)
    for (k, (b, j)) in enumerate(cells)
        # 传播子：一阶全步 / 二阶半步；Td 模型带区间左端点 t
        dt = b === :τ ? lattice.δτ : lattice.δt
        t = b === :+ ? (j - 1) * dt : (b === :- ? (lattice.Nt - j) * dt : nothing)
        Shalf = order == 2
        S = _model_propagator(model, b, Shalf ? dt / 2 : dt, t)
        mats = order == 2 ? [S * M * S for M in weights] :
               b === :- ? [M * S for M in weights] : [S * M for M in weights]
        Ol = cutleft(k); Or = cutright(k)
        mats = [Ol * M * Or for M in mats]
        # 奇插入穿越全部更早 Grassmann 系数的 Koszul 符号
        oddlater = 0
        for (cut, items) in cuts, ins in items
            cut >= k || continue
            q = ins.charge === nothing ? z : ins.charge
            oddlater += GM.fermionparity(q) ? 1 : 0
        end
        for m in 0:dim^2-1
            mats[m+1] .*= (-1)^(oddlater * count_ones(m))
        end
        leftlabels = k == 1 ? [z] : [qsub(qsub(fqs[s], fqs[t]), prefix) for (s, t) in pairs]
        for (cut, items) in cuts, ins in items
            cut == k || continue
            ins.charge === nothing || (prefix = qadd(prefix, ins.charge))
        end
        rightlabels = k == L ? [z] : [qsub(qsub(fqs[s], fqs[t]), prefix) for (s, t) in pairs]
        core = zeros(ComplexF64, length(leftlabels), dim^2, length(rightlabels))
        for m in 0:dim^2-1
            M = mats[m+1]
            if L == 1
                core[1, m+1, 1] = tr(M * rho)
            elseif k == 1
                for (r, (s, t)) in enumerate(pairs)
                    core[1, m+1, r] = (M * rho)[s, t]
                end
            elseif k == L
                for (l, (s, t)) in enumerate(pairs)
                    core[l, m+1, 1] = M[t, s]
                end
            else
                for (l, (s, t)) in enumerate(pairs), (r, (v, w)) in enumerate(pairs)
                    t == w && (core[l, m+1, r] = M[v, s])
                end
            end
        end
        all(iszero, core) && return monomial(spec, repeat(charges_per_cell, L); coeff=0)
        ts = _cell_tensors(spec, core, leftlabels, rightlabels, charges_per_cell)
        data === nothing ? (data = ts) : append!(data, ts)
        append!(labels, [(b, j, v) for v in 1:2bands])
    end
    lscale = order == 2 ? real(coupling_norm / 4 * sum(_cell_step2(lattice, b) for (b, _) in cells)) : 0.0
    out = pairize(spec, data)
    GM.setlogscale!(out, lscale)
    # 时序构造 → 存储序：每格 bands 站成块移动（站粒度置换——GM.reorder! 的
    # perm 语义 = 新链第 k 站 = 原链第 perm[k] 站；上游 GM.swap! 的费米换序号
    # 已修复，多占据组态保真）
    chrono = Dict(cell => k for (k, cell) in enumerate(cells))
    perm = Int[]
    for cell in _storage_cells(lattice)
        k = chrono[cell]
        append!(perm, [(k - 1) * bands .+ (1:bands)...])
    end
    GM.reorder!(out, perm; trunc=GM.truncrelerr(1e-14))
    if trunc !== nothing
        _, err = GM.truncate!(out; trunc)
        truncerrs === nothing || push!(truncerrs, err)
    end
    return out
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

# spec symbol 反查（守卫报错信息用）

# ---------------------------------------------------------------------------
# 公开皮肤（旧 API 名保持；plan §3.3）
# ---------------------------------------------------------------------------

"""
    sysdynamics(lattice, model; order=1, initial_state=nothing,
                coupling_norm=0.0, insertions=OperatorInsertion[], trunc) -> GMPS

杂质系统因子 K（cell automaton；order=1 相干 anti-Wick 语义，q=1 局部权重）。
ρ₀/迹闭合已折进单元（boundarycondition! 不再需要）；`initial_state`
（FockMatrix/矩阵）仅实时轮廓生效（缺省 = 真空投影）。
"""
function sysdynamics(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian;
                     order::Integer=1, kwargs...)
    order in (1, 2) || throw(ArgumentError("order must be 1 or 2"))
    return _system_factor_reference(lattice, model; order=Int(order),
                                    spec=TK.sectortype(lattice), kwargs...)
end

"""
    sysdynamics_weyl2order(lattice, model; coupling_norm, kwargs...) -> GMPS

Weyl 二阶 K（order=2：q=1/2 对称半步 S·M·S + Weyl 标量
exp(coupling_norm/4·Σ step²)；`coupling_norm` = 全带总耦合 Σ_b α_b，必填 ≥ 0）。
IF 与收缩零改动（INTEGRATION.md §2 F1）。
"""
function sysdynamics_weyl2order(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian;
                                coupling_norm::Real, kwargs...)
    return _system_factor_reference(lattice, model; order=2, coupling_norm=coupling_norm,
                                    spec=TK.sectortype(lattice), kwargs...)
end
