# fockmpo.jl：Fock 表象直积 MPO + fockmpo2gmps 转换（荷账自读）
# （plan/sysdynamics_fock_mpo_plan.md；2026-11-20 起穿戴内联 fock_sysdynamics、
#  converter 纯翻译——plan/fockmpo_dress_selfread_plan.md：无独立 dressing 层）
#
# 架构（两段闭环）：
#   fock_sysdynamics → FockMPO（**穿戴完毕**：开链装配 → 切位 ≥ 1 稀疏乘法 →
#     ρ₀/切位 0 三明治折首站（在稀疏乘法之后！）→ oddlater 逐站 slot 号 →
#     符号荷 vs 链 δ 读数交叉自检（P5））；仍开链——迹闭合归 converter
#   fockmpo2gmps → 荷账自读（槽名义荷表 + 站读数 δ，不消费任何插入/初态元
#     数据）+ 迹闭合 + 键重封块 → split_cell + pairize（基线分站装置，符号
#     规范逐位继承）→ reorder!（v7 虚轴翻转后 imag perm 恒等 = 无操作）→ GMPS
#
# V1 定谳（原型实证，见 debug/fockmpo_proto/）：
#   1. **分站符号规范不可局部化**：[band1-first] GMPS 的系数带 GM 楔积规范相位
#      σ(occ)（ED 验证、与 IF 链配对），对胞级/带级局部因子不可分解（GF(2) 拟合
#      10 候选全矛盾）。直构无 SVD 路线无法继承该规范（bands=2 N≥3 符号翻转
#      24/4096，σ 谱与拟合脚本留 debug/fockmpo_proto/）。
#   2. 因此 fockmpo2gmps 的分站步骤复用基线自己的 split_cell+pairize（逐位一致的
#      符号规范）；SVD = 逐胞满秩精确的纯符号规范装置，代价与基线相同、无回归。
#      「无 SVD 直构」为 V2 议题（plan §6 的效率本就推迟）。
#   3. **FockMPO 键腿 = 非分格（平凡扇区）**：对角单项式腿携带荷时，MPO 张量的
#      融合约束（q_bL ⊗ q_m ≅ q_m ⊗ q_bR ⟺ q_bL = q_bR）与值结构
#      （q(a_out) = q(a_in) + q_m）矛盾——带荷对角 MPO 腿不可分格表达。
#      扇区结构由 fockmpo2gmps 重建（split_cell 的荷差键 + pairize 的分格站）。
#
# 语义（闸门 = 手工闭环对拍 + observables golden U1GTEMPO2 机器精度锚）：
#   一阶（q=1）：τ/+ 支 mats = S·M、− 支 M·S；二阶（Weyl）：S·M·S + Weyl 标量
#   lscale = coupling_norm/4·Σ step²（_system_factor_fock/sysdynamics 层
#   setlogscale!——coupling_norm 不进 fock_sysdynamics/fockmpo2gmps）。
#   ρ₀/切位 0 插入三明治（Ol₀·ρ·Or₀）由 fock_sysdynamics 折进首站（在切位 ≥ 1
#   稀疏乘法之后）；切位 ≥ 1 插入走 sparsefockmpo `*`；oddlater Koszul 号
#   (−1)^{oddlater·popcount(m)}（oddlater = 切位 ≥ k 的奇插入数，谓词 cut ≥ k）
#   由 fock_sysdynamics 逐站折进链。插入荷账由 converter 从穿戴链自读：
#   δ_k = f_q(out) − f_q(in) − q_m（未 dressing 链 δ ≡ 0）、prefix(键 k) =
#   Σ_{j≤k} δ_j（与旧 cut 累计逐位恒等——plan §2.2 验证例）。
#   单项式腿 4^bands 态的槽编码与 local_weights 一致（band b = 2bit 低位组，
#   0=无 1=ā(−d†) 2=a(d) 3=aā(W₄)）。

# ---------------------------------------------------------------------------
# 类型（纯容器；腿位型 = MPOTensor：(键L ⊗ 单项式出) ← (单项式入 ⊗ 键R)）
# ---------------------------------------------------------------------------

"Fock 键 / 单项式腿空间（V1 非分格：全平凡扇区，理由见文件头 3）"
_fock_bondspace(G::Type{<:Sector}, dim::Int) = Vect[G](one(G) => dim)   # one(G)：ProductSector 系不接受 Int 键

"""
    FockMPO{T<:MPOTensor}

Fock 表象直积 MPO（纯容器，与 GMPS 同构的类型参数化）：每物理格一 site，
bands 合并进 site 内部（Fock 维 2^bands、单项式腿 4^bands，site 内块对角——
slot 指标过站不变，每单项式一个 (dim,dim) 算符块）。站位解释由 lattice 携带：
site k ↔ chronological_cells(lattice)[k]。**链 = 穿戴态**（`fock_sysdynamics`
产出：ρ₀/切位 0 三明治已折、插入已乘、oddlater 号已折）、仍开链（迹闭合归
`fockmpo2gmps`，可复用于其它收缩口径）。logscale 不挂 FockMPO——Weyl 标量在
fockmpo2gmps 产出 GMPS 后由 _system_factor_fock/sysdynamics 层 setlogscale!。
"""
struct FockMPO{T<:GM.MPOTensor}
    data::Vector{T}
    function FockMPO(data::Vector{T}) where {T<:GM.MPOTensor}
        isempty(data) && throw(ArgumentError("empty FockMPO"))
        new{T}(data)
    end
end

Base.length(f::FockMPO) = length(f.data)
Base.getindex(f::FockMPO, i::Int) = f.data[i]
Base.copy(f::FockMPO) = FockMPO([copy(t) for t in f.data])   # 深拷贝：张量为引用语义，就地穿戴/乘法路线的调用方需要独立链

"""
    SparseFockMPO{T<:MPOTensor}

稀疏算符链（镜像 GM SparseGMPS）：positions = 支撑链位（chrono 格序 1..L），
data = 每支撑位的总算符（同位多算符的列乘积，push 序 = 0⁺ 序：先 push 最左）。
一个 SparseFockMPO 只装一个 side（:left 走 `A*K`、:right 走 `K*A`；greater/lesser
双位插入 = 同 side 两位，左右各建一个）。乘法 = 纯张量收缩，无 graded 符号。
"""
struct SparseFockMPO{T<:GM.MPOTensor}
    positions::Vector{Int}
    data::Vector{T}
    function SparseFockMPO(positions::Vector{Int}, data::Vector{T}) where {T<:GM.MPOTensor}
        length(positions) == length(data) ||
            throw(DimensionMismatch("positions and data must have equal length"))
        new{T}(positions, data)
    end
end

# ---------------------------------------------------------------------------
# fock_sysdynamics：装配（纯 Fock 无 Grassmann 内容；开链）
# ---------------------------------------------------------------------------

"模型传播子（一阶全步 / 二阶半步；Td 模型带区间左端点 t）——sysdynamics.jl 同款"
function _fock_propagator_at(lattice::AbstractGrassmannLattice, model, b::Symbol, j::Int, order::Int)
    dt = b === :τ ? lattice.δτ : lattice.δt
    t = b === :+ ? (j - 1) * dt : (b === :- ? (lattice.Nt - j) * dt : nothing)
    return Matrix(_model_propagator(model, b, order == 2 ? dt / 2 : dt, t))
end

"""
    _fock_sysdynamics_chronological(lattice, model; order=1, insertions=OperatorInsertion[],
                     initial_state=nothing) -> FockMPO

（deprecated 时序装配体——供 fockmpo2gmps_deprecated 与闸门对照使用；
现行装配见 foldedmpo.jl 的折叠 fock_sysdynamics。）时序 Fock 直积 MPO，
**一次装配产出穿戴完毕的链**（P1——无独立 dressing 层）：①开链装配（逐物理格一 site，site 内块对角，单元算符与
sysdynamics 同口径：一阶 τ/+ 支 S·M、− 支 M·S，二阶 S·M·S）；②切位 ≥ 1 插入
稀疏乘法（`sparsefockmpo *`，过滤按 `insertion_cut_index` 切位）；③ρ₀/切位 0
三明治折首站（**在稀疏乘法之后**——切位 1 插入作用在站 1 时语义 =
(Ol₁·M₁·Or₁)·(Ol₀·ρ·Or₀)；边界态内联：imag/mixed = 热迹 I、real = 真空投影 P，
`initial_state`（FockMatrix/矩阵）仅实时轮廓覆盖）；④oddlater 逐站 slot 对角号
（谓词 cut ≥ k，奇偶源 = 符号荷）；⑤P5 交叉自检：符号荷累计 vs 链 δ 读数逐站
核对（不一致 = 物化 bug——站内一致性检查抓不住的"一致但错"缺口，如 JW 符号
错折 → ArgumentError）。

**不含** 迹闭合/Weyl 标量——前者归 `fockmpo2gmps`，后者由调用方
（_system_factor_fock/sysdynamics 层）setlogscale! 补上。
"""
function _fock_sysdynamics_chronological(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian;
                          order::Integer=1, insertions=OperatorInsertion[], initial_state=nothing)
    order in (1, 2) || throw(ArgumentError("order must be 1 or 2"))
    bands = lattice.bands
    G = TK.sectortype(lattice)
    bands_guard(G, bands)
    num_bands(model) == bands ||
        throw(ArgumentError("model bands $(num_bands(model)) ≠ lattice bands $bands"))
    cells = chronological_cells(lattice)
    isempty(cells) && throw(ArgumentError("empty lattice"))
    weights = local_weights(G, bands, Int(order))
    dim = 2^bands
    focksp = _fock_bondspace(G, dim)
    monosp = _fock_bondspace(G, dim^2)
    z = one(G)
    data = [begin
        S = _fock_propagator_at(lattice, model, b, j, Int(order))
        mats = order == 2 ? [S * W * S for W in weights] :
               b === :- ? [W * S for W in weights] : [S * W for W in weights]
        t = TK.zeros(ComplexF64, focksp ⊗ monosp ← monosp ⊗ focksp)
        blk = t[(z, z, TK.dual(z), TK.dual(z))]
        for m in 0:dim^2-1
            blk[:, m+1, m+1, :] .= mats[m+1]   # 单项式腿对角（slot 过站不变）
        end
        t
    end for (b, j) in cells]
    f = FockMPO(data)

    # ② 切位 ≥ 1 稀疏乘法（左右各一 SparseFockMPO；切位 0 归三明治——过滤按
    # 切位映射而非 j 字段：:- 支 j = 0 = 链末切位）
    ge1 = filter(i -> insertion_cut_index(lattice, i) >= 1, insertions)
    Wl = sparsefockmpo(lattice, filter(i -> i.side === :left, ge1))
    Wr = sparsefockmpo(lattice, filter(i -> i.side === :right, ge1))
    f = Wr * (Wl * f)

    # ③ ρ₀/切位 0 三明治折首站（在稀疏乘法之后——红线 1；push 序 = 0⁺ 序）。
    # 边界态内联（P3 恒折）：imag/mixed = 热迹 I、real = 真空投影 P；
    # `initial_state`（FockMatrix/矩阵）仅实时轮廓覆盖
    rho = if lattice isa AbstractRealGrassmannLattice
        initial_state === nothing ? Matrix{ComplexF64}(Diagonal(vcat([1.0], zeros(dim - 1)))) :
        (initial_state isa FockMatrix ? Matrix{ComplexF64}(initial_state.data) :
         Matrix{ComplexF64}(initial_state))
    else
        Matrix{ComplexF64}(I, dim, dim)
    end
    size(rho) == (dim, dim) || throw(DimensionMismatch("initial state must be $dim×$dim"))
    rho = foldr((i, A) -> _fockop_matrix(i.op, bands) * A,
                filter(i -> insertion_cut_index(lattice, i) == 0 && i.side === :left, insertions);
                init=rho)
    rho = foldl((A, i) -> A * _fockop_matrix(i.op, bands),
                filter(i -> insertion_cut_index(lattice, i) == 0 && i.side === :right, insertions);
                init=rho)
    blk1 = f.data[1][(z, z, TK.dual(z), TK.dual(z))]
    for m in 0:dim^2-1
        blk1[:, m+1, m+1, :] .= blk1[:, m+1, m+1, :] * rho   # 右乘进全部 slot 块
    end

    # ④ oddlater 逐站 slot 对角号（谓词 cut ≥ k 逐字继承；标量因子对线性核
    # 构造可提出，折进链与折进核逐位恒等）
    bycut = Dict{Int,Vector{OperatorInsertion}}()
    for ins in insertions
        push!(get!(bycut, insertion_cut_index(lattice, ins), OperatorInsertion[]), ins)
    end
    for (k, t) in enumerate(f.data)
        oddlater = 0
        for (cut, items) in bycut, ins in items
            cut >= k || continue
            oddlater += GM.fermionparity(_fockop_charge(ins.op, G, bands)) ? 1 : 0
        end
        oddlater == 0 && continue
        blk = t[(z, z, TK.dual(z), TK.dual(z))]
        for m in 0:dim^2-1
            blk[:, m+1, m+1, :] .*= (-1)^(oddlater * count_ones(m))
        end
    end

    # ⑤ P5 交叉自检：符号荷累计 vs 链 δ 读数逐站核对（converter 侧独立再读，
    # 双路径互证）。切位 0 插入随三明治折进站 1（prefix 恒等式：δ₁ 含切位 0/1）。
    fqs1 = fock_charges(G, bands)
    qms = _slot_charges(G, bands)
    for (k, t) in enumerate(f.data)
        expected = z
        for (cut, items) in bycut, ins in items
            (cut == k || (cut == 0 && k == 1)) || continue
            expected = qadd(expected, _fockop_charge(ins.op, G, bands))
        end
        actual = _site_charge_shift(t[(z, z, TK.dual(z), TK.dual(z))], fqs1, qms, bands)
        actual == expected || throw(ArgumentError(
            "fock_sysdynamics: dressing charge mismatch at site $k — symbol charge " *
            "$expected vs chain δ reading $actual (materialization bug?)"))
    end
    return f
end

# ---------------------------------------------------------------------------
# 荷账自读装置（P2/P5 共用：槽名义荷表 + 站读数 δ）
# ---------------------------------------------------------------------------

"""
槽名义荷表（先验，纯函数）：槽编码 m 每带 2bit 组 `slots[b] = (m >> 2(b−1)) & 3`
（band 1 = 低位组，与 local_weights 逐位同源）：0 = q·I−n（0）、1 = −d†（+q_b）、
2 = d（−q_b）、3 = W₄（0）——`q_m = Σ_b ([1]−[2])·var_charge(G, b, bands)`。
带间 graded⊗ 符号与 JW ±1 均纯标量，不影响稀疏性（荷读数）。
"""
function _slot_charges(::Type{G}, bands::Int) where {G<:Sector}
    out = Vector{typeof(one(G))}(undef, 4^bands)
    for m in 0:4^bands-1
        q = one(G)
        for b in 1:bands
            slot = (m >> (2 * (b - 1))) & 3
            slot == 1 && (q = qadd(q, var_charge(G, b, bands)))
            slot == 2 && (q = qadd(q, qneg(var_charge(G, b, bands))))
        end
        out[m+1] = q
    end
    return out
end

"""
站荷移位读数（站内自验证三件套之 1/3 在此）：站张量腿序 =
(键L ⊗ 单项式出) ← (单项式入 ⊗ 键R)，槽对角；条目 `blk[s, m+1, m+1, u] ≠ 0`
（容差 = 1e-12 × 站内 max|blk|）⟹ `f_q(s) − f_q(u) = q_m + δ_k`，即
**δ_k = f_q(out) − f_q(in) − q_m**（任一非零条目；未 dressing 链 δ ≡ 0——
机器可检不变量，test/api/fockmpo.jl 锚）。全条目低于容差 → δ = 0（零链，
split_cell 语义）。

三件套：①站内一致性——同站不同条目/槽的 δ 读数全等，不一致 = 装配电荷不相干
（电荷非对角初态，A1′）→ ArgumentError；③slot 对角断言——离对角 slot 条目
≈ 0（读数前置检查）。②（Σδ = 0 环闭合）在 fockmpo2gmps。
"""
function _site_charge_shift(blk, fqs1, qms::Vector{G}, bands::Int) where {G<:Sector}
    dim = size(blk, 1)
    tol = 1e-12 * maximum(abs, blk)
    # ③ slot 对角断言（读数前置）：离对角 slot 条目 ≈ 0
    for m1 in 1:dim^2, m2 in 1:dim^2, u in 1:dim, s in 1:dim
        m1 == m2 && continue
        abs(blk[s, m1, m2, u]) > tol && throw(ArgumentError(
            "monomial leg not slot-diagonal (off-diagonal |entry| > $tol)"))
    end
    # ① 站内一致性读数：非零条目的 δ 全等
    δ = nothing
    for m in 1:dim^2, u in 1:dim, s in 1:dim
        abs(blk[s, m, m, u]) > tol || continue
        d = qsub(qsub(fqs1[s], fqs1[u]), qms[m])
        δ === nothing ? (δ = d) :
        δ == d || throw(ArgumentError(
            "site charge bookkeeping inconsistent (charge-non-diagonal dressing?): " *
            "read δ = $δ and $d"))
    end
    return δ === nothing ? one(G) : δ
end

# ---------------------------------------------------------------------------
# 插入稀疏算符链（切位 ≥ 1；切位 0 随边界态三明治内联在 fock_sysdynamics）
# ---------------------------------------------------------------------------

"""
    sparsefockmpo(lattice, insertions) -> SparseFockMPO

    插入算符 → 稀疏算符链（V7 装配点物化：FockOp 在此物化为 Fock 矩阵）：
切位 ≥ 1 的插入按 `insertion_cut_index` 映射到链位，同位算符按 push 序折成
列积（foldr/foldl，先 push 最左，0⁺ 序）。**切位 0 插入在此过滤掉**（归
`fock_sysdynamics` 的 ρ₀/切位 0 三明治）；全部插入须同 side（左右各建一个
SparseFockMPO）。⚠️ 过滤按 cut index 而非 `j` 字段——`:-` 支 j = 0 = 链末
切位（cut = 2Nt），不是 ρ₀ 切位（lesser(j=0) 路径）。"""
function sparsefockmpo(lattice::AbstractGrassmannLattice, insertions::Vector{OperatorInsertion})
    ins1 = filter(i -> insertion_cut_index(lattice, i) >= 1, insertions)
    side = isempty(ins1) ? :left : first(ins1).side
    all(i -> i.side === side, ins1) ||
        throw(ArgumentError("mixed sides: build one SparseFockMPO per side (:left / :right)"))
    G = TK.sectortype(lattice)
    bands = lattice.bands
    dim = 2^bands
    focksp = _fock_bondspace(G, dim)
    onesp = _fock_bondspace(G, 1)
    z = one(G)
    I0 = Matrix{ComplexF64}(I, dim, dim)
    bycut = Dict{Int,Vector{OperatorInsertion}}()
    for ins in ins1
        push!(get!(bycut, insertion_cut_index(lattice, ins), OperatorInsertion[]), ins)
    end
    positions = sort!(collect(keys(bycut)))
    data = [begin
        # 左：foldr（先 push 最左）；右：foldl（先 push 最左）——基线 cutleft/cutright 同款
        OL = foldr((i, A) -> _fockop_matrix(i.op, bands) * A,
                   filter(i -> i.side === :left, bycut[p]); init=I0)
        OR = foldl((A, i) -> A * _fockop_matrix(i.op, bands),
                   filter(i -> i.side === :right, bycut[p]); init=I0)
        t = TK.zeros(ComplexF64, focksp ⊗ onesp ← onesp ⊗ focksp)
        t[(z, z, TK.dual(z), TK.dual(z))][:, 1, 1, :] .= OL * OR
        t
    end for p in positions]
    return SparseFockMPO(positions, data)
end

# 左作用：A_k·K_k = 插在 cell k 前（现行 :left）
function Base.:*(A::SparseFockMPO, K::FockMPO)
    L = length(K)
    z = one(TK.sectortype(first(K.data)))
    newdata = Vector{eltype(K.data)}(K.data)
    for (p, tA) in zip(A.positions, A.data)
        1 <= p <= L || throw(BoundsError(K, p))
        tK = K.data[p]
        t = TK.zeros(TK.scalartype(tK), TK.space(tK, 1) ⊗ TK.space(tK, 2) ← TK.space(tK, 3) ⊗ TK.domain(tK)[2])
        bA = tA[(z, z, TK.dual(z), TK.dual(z))]   # (dim, 1, 1, dim)
        bK = tK[(z, z, TK.dual(z), TK.dual(z))]   # (dim, 4^B, 4^B, dim)
        bO = t[(z, z, TK.dual(z), TK.dual(z))]
        d = size(bK, 1)
        for m in 1:size(bK, 2)
            bO[:, m, m, :] .= reshape(bA, d, d) * bK[:, m, m, :]
        end
        newdata[p] = t
    end
    return FockMPO(newdata)
end

# 右作用：K_k·A_k = 插在 cell k 后（现行 :right）
function Base.:*(K::FockMPO, A::SparseFockMPO)
    L = length(K)
    z = one(TK.sectortype(first(K.data)))
    newdata = Vector{eltype(K.data)}(K.data)
    for (p, tA) in zip(A.positions, A.data)
        1 <= p <= L || throw(BoundsError(K, p))
        tK = K.data[p]
        t = TK.zeros(TK.scalartype(tK), TK.space(tK, 1) ⊗ TK.space(tK, 2) ← TK.space(tK, 3) ⊗ TK.domain(tK)[2])
        bA = tA[(z, z, TK.dual(z), TK.dual(z))]
        bK = tK[(z, z, TK.dual(z), TK.dual(z))]
        bO = t[(z, z, TK.dual(z), TK.dual(z))]
        d = size(bK, 1)
        for m in 1:size(bK, 2)
            bO[:, m, m, :] .= bK[:, m, m, :] * reshape(bA, d, d)
        end
        newdata[p] = t
    end
    return FockMPO(newdata)
end

# ---------------------------------------------------------------------------
# fockmpo2gmps：迹闭合 + ρ₀/切位 0 三明治 + 键重封块 + oddlater
#               → split_cell + pairize（基线分站装置，符号逐位继承）→ GMPS
# ---------------------------------------------------------------------------

"""
    fockmpo2gmps_deprecated(f::FockMPO, lattice; trunc=nothing, truncerrs=nothing) -> GMPS

（deprecated——时序装配 + 末端 reorder! 重排；现行装配见 foldedmpo.jl 的折叠
fockmpo2gmps。语义等价性由 todense 闸门锚定。v7 虚轴翻转后 imag 的
_storage_cells = 时序 ⟹ 本函数 perm 恒等、reorder! 退化为无操作，仅作
跨布局兼容保留。）FockMPO → GMPS（迹闭合 + 分站；
**P2 纯翻译 + 荷账自读**）：`f` 必须是
`fock_sysdynamics` 产出的**穿戴链**（ρ₀/切位 0 三明治、切位 ≥ 1 插入、oddlater
号已折进链）——本函数不消费任何插入/初态元数据（`rho0`/`insertions` kwargs
已随 2026-11-20 换代删除）。逐格自读站荷移位 δ_k（槽名义荷表 + 站读数：
δ_k = f_q(out) − f_q(in) − q_m；未 dressing 链 δ ≡ 0）→ 键标签
prefix(键 k) = Σ_{j≤k} δ_j（与旧 cut 累计逐位恒等——plan §2.2 验证例）→
(s,t) 对差键核（线程注入：首格直读、中间 = M[v,s]δ_{t,w}、末格 =
M[t,s] 右边界真空、L==1 = 迹闭合）→ split_cell（基线分站装置：荷差键 + 逐格
SVD 切分，符号规范逐位继承）→ pairize → reorder!（Chrono→存储序；imag 下
v7 翻转后 perm 恒等）→ GMPS。

站内自验证三件套：①站内一致性（δ 读数全等——不一致 = 电荷不相干 dressing，
如电荷非对角初态 → ArgumentError）；②Σ_k δ_k = 0 环闭合守卫（≡ 旧总插入荷
守卫，源 = 链矩阵）；③slot 对角断言——①③在 `_site_charge_shift`。

kwargs：`trunc`/`truncerrs` = reorder 后的 GM.truncate! 透传（缺省不截断）。
**Weyl 标量不在本函数**——调用方（_system_factor_fock/sysdynamics 层）持
coupling_norm，对产出的 GMPS 做 `logscale += lscale`（⚠️ 累加非覆盖：
truncate!/swap! 的压缩归一化因子已吸收进 logscale，E2 语义）。
"""
function fockmpo2gmps_deprecated(f::FockMPO, lattice::AbstractGrassmannLattice;
                      trunc=nothing, truncerrs=nothing)
    G = TK.sectortype(lattice)
    bands = lattice.bands
    dim = 2^bands
    z = one(G)
    cells = chronological_cells(lattice)
    L = length(cells)
    length(f.data) == L ||
        throw(DimensionMismatch("FockMPO has $(length(f.data)) sites, lattice has $L cells"))

    # 荷账自读（P2）：站读数 δ（站内一致性 + slot 对角校验在 _site_charge_shift）
    fqs1 = fock_charges(G, bands)
    qms = _slot_charges(G, bands)
    deltas = [_site_charge_shift(f.data[k][(z, z, TK.dual(z), TK.dual(z))], fqs1, qms, bands)
              for k in 1:L]
    total_q = foldl(qadd, deltas; init=z)
    total_q == z ||
        throw(ArgumentError("total chain charge must be zero (Σδ = $total_q ≠ 0)"))

    # prefix 恒等式：prefix(键 k) = Σ_{j≤k} δ_j（≡ 旧切位 0..k 插入荷累计）
    prefixes = Vector{typeof(z)}(undef, L)
    p = z
    for k in 1:L
        p = qadd(p, deltas[k])
        prefixes[k] = p
    end

    pairs_st = [(s, t) for t in 1:dim for s in 1:dim]   # 1-based，基线同款枚举序
    chpc = convert(Vector{eltype(fqs1)}, cell_charges(G, bands))

    data = nothing
    for (k, (b, j)) in enumerate(cells)
        first_cell = k == 1
        last_cell = k == L
        leftlabels = first_cell ? [z] :
                     [qsub(qsub(fqs1[s], fqs1[t]), prefixes[k-1]) for (s, t) in pairs_st]
        rightlabels = last_cell ? [z] :
                      [qsub(qsub(fqs1[s], fqs1[t]), prefixes[k]) for (s, t) in pairs_st]
        core = zeros(ComplexF64, length(leftlabels), dim^2, length(rightlabels))
        blk = f.data[k][(z, z, TK.dual(z), TK.dual(z))]   # (dim, 4^B, 4^B, dim)
        if L == 1
            for m in 0:dim^2-1
                core[1, m+1, 1] = tr(blk[:, m+1, m+1, :])   # 迹闭合（ρ₀ 已折进链）
            end
        elseif first_cell
            for m in 0:dim^2-1, (r, (s, t)) in enumerate(pairs_st)
                core[1, m+1, r] = blk[s, m+1, m+1, t]   # 首格直读（无 ρ 收缩、无 sgn）
            end
        elseif last_cell
            for m in 0:dim^2-1, (l, (s, t)) in enumerate(pairs_st)
                core[l, m+1, 1] = blk[t, m+1, m+1, s]
            end
        else
            for m in 0:dim^2-1, (l, (s, t)) in enumerate(pairs_st), (r, (v, w)) in enumerate(pairs_st)
                t == w && (core[l, m+1, r] = blk[v, m+1, m+1, s])
            end
        end
        column = _cell_tensors(G, core, leftlabels, rightlabels, chpc)
        isempty(column) && return monomial(G, repeat(chpc, L); coeff=0)
        data === nothing ? (data = copy(column)) : append!(data, column)
    end
    out = pairize(G, data)
    # 时序构造 → 存储序：每格 bands 站成块移动（GM.reorder! perm 语义 = 新链第 k 站
    # = 原链第 perm[k] 站；swap! 费米换序号已修复）。v7 虚轴翻转后 imag 的
    # _storage_cells = 时序 ⟹ perm 恒等、reorder! 退化为无操作（保留以兼容
    # real/mixed 布局与未来布局变更；graded 换序号语义见 plan §7）。
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

# ---------------------------------------------------------------------------
# gmps2fockmpo：每胞合并（Fock 层再操作 / 跨口径调试对拍；非闭形数学逆）
# ---------------------------------------------------------------------------

"""
    gmps2fockmpo(g::GMPS, lattice) -> FockMPO

新形 GMPS → 每胞合并 FockMPO（**V2**）：胞内 bands 个带站收缩为单 MPOTensor，
单项式腿 = 各带站物理空间的楔积积空间。**当前未实施**——合并的两个前提待
定谳：①楔积积空间的 (sector, slot) 基底索引映射（TensorKit ⊗ 融合约定）；
②合并秩回归 MPOTensor (2,2) 的腿融合装置。迹闭合烙在边界，本函数亦非
fockmpo2gmps 的数学逆（用途 = Fock 层再操作/跨口径调试对拍）。
"""
function gmps2fockmpo(g::GMPS, lattice::AbstractGrassmannLattice)
    throw(ArgumentError("gmps2fockmpo：V2 未实施（plan/sysdynamics_fock_mpo_plan.md §3.2；楔积基底映射与合并秩待定谳）"))
end

# ---------------------------------------------------------------------------
# 闭环装配：sysdynamics 公开入口的实现体（穿戴内联 fock_sysdynamics——2026-11-20
# 起切位 ≥ 1 乘法/ρ₀ 三明治/oddlater 全在 fock_sysdynamics 内完成，无独立
# dressing 层；plan/fockmpo_dress_selfread_plan.md P1）
# ---------------------------------------------------------------------------

"Weyl 标量 lscale（order=2 = coupling_norm/4·Σ step²，一阶 0）——旧径同式同序"
_fock_weyl_lscale(lattice::AbstractGrassmannLattice, order::Int, coupling_norm::Real) =
    order == 2 ? real(coupling_norm / 4 *
                      sum(_cell_step2(lattice, b) for (b, _) in chronological_cells(lattice))) :
    0.0

# Z 常数类型分派（exactif_port_plan E6）：
# - Frozen = Weyl 公式随 order。`coupling_norm` kwarg = 显式值优先（0 = 未给
#   哨兵 → kernel.coupling_norm）；**不做与 kernel 值的一致性抛错**——多带装配
#   的既有惯例是 kernel 全带 cn ≠ K 侧显式标量（如 two-band ED 锚：单 bath 核
#   全带 cn = 2α、K 侧传 α；比值观测量下常数分子分母同因子严格相消，两者各自
#   一致即可，绝对 Z 调用方自管一致性）。
# - Exact = log|normdelta| 与 order 无关——严格核无冻结顶点 O(δt²) 测度宽度
#   可补偿；号不进 log，由 _system_factor_fock 的 lmul! 折进首站。Weyl 通道
#   对 Exact 无意义：显式非零 coupling_norm 响抛（恒错禁止：normdelta 当
#   耦合范数用）。
_constant_lscale(k::FrozenHybridizationKernel, lattice::AbstractGrassmannLattice, order::Int;
                 coupling_norm::Real=k.coupling_norm) =
    _fock_weyl_lscale(lattice, order, coupling_norm == 0 ? k.coupling_norm : coupling_norm)
_constant_lscale(k::ExactHybridizationKernel, ::AbstractGrassmannLattice, ::Int) =
    log(abs(k.normdelta))

"""
    _system_factor_fock(lattice, model; order=1, coupling_norm=0.0, kernel=nothing,
                        initial_state=nothing, insertions=OperatorInsertion[],
                        trunc=nothing, truncerrs=nothing, spec=nothing) -> GMPS

sysdynamics 的实现体（三行闭环）：`fock_sysdynamics`（穿戴内联：切位 ≥ 1 稀疏
乘法 + ρ₀/切位 0 三明治 + oddlater + P5 自检，`initial_state` 仅实时轮廓生效）
→ `fockmpo2gmps`（荷账自读 + 迹闭合 + 分站）→ Z 常数 `logscale +=`。
**常数按类型分派**（exactif_port_plan E6）：`kernel === nothing` → 现行 Weyl
公式（`coupling_norm` 输入，order=2 要求 ≥ 0）；`kernel !== nothing` →
`_constant_lscale(kernel, lattice, order)`（Frozen：显式 coupling_norm 优先、
0 = 未给哨兵取 kernel.coupling_norm，无一致性抛错——多带惯例见
`_constant_lscale` 注；Exact：log|normdelta| + 号折首站（`lmul!`），显式非零
coupling_norm 响抛）。⚠️ 累加非覆盖：logscale 是链幅的 log 场，E2 语义下
truncate!/swap! 会把压缩归一化因子吸收进去，覆盖会抹掉这些因子；乘
exp(lscale) 语义 = 累加。`spec` 兼容旧径 kwarg（须 = sectortype(lattice)）。
"""
function _system_factor_fock(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian;
                             order::Int=1, coupling_norm::Real=0.0,
                             kernel::Union{Nothing,AbstractHybridizationKernel}=nothing,
                             initial_state=nothing,
                             insertions=OperatorInsertion[], trunc=nothing, truncerrs=nothing,
                             spec=nothing)
    order in (1, 2) || throw(ArgumentError("order must be 1 or 2"))
    kernel isa ExactHybridizationKernel && coupling_norm != 0 && throw(ArgumentError(
        "Exact kernel chains carry their Z constant via normdelta (chain-amplitude " *
        "dispatch); explicit coupling_norm (Weyl scalar) does not apply"))
    spec === nothing || spec === TK.sectortype(lattice) ||
        throw(ArgumentError("spec $spec ≠ sectortype(lattice) $(TK.sectortype(lattice))"))
    Kf = _fock_sysdynamics_chronological(lattice, model; order=order, insertions=insertions,
                                         initial_state=initial_state)
    out = fockmpo2gmps_deprecated(Kf, lattice; trunc=trunc, truncerrs=truncerrs)
    if kernel === nothing
        order == 2 && !(coupling_norm >= 0) &&
            throw(ArgumentError("order=2 requires coupling_norm ≥ 0 (total Σ_b α_b)"))
        GM.setlogscale!(out, GM.logscale(out) + _fock_weyl_lscale(lattice, order, coupling_norm))
    else
        GM.setlogscale!(out, GM.logscale(out) + _constant_lscale(kernel, lattice, order;
                                                                 coupling_norm=coupling_norm))
        kernel isa ExactHybridizationKernel && GM.lmul!(sign(kernel.normdelta), out)
    end
    return out
end


