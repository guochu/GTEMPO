# fockmpo.jl：Fock 表象直积 MPO + fockmpo2gmps 转换
# （plan/sysdynamics_fock_mpo_plan.md；2026-11-18 起为 sysdynamics 唯一路径——
#  旧径 _system_factor_reference/SystemFactorPreparation 随 V2 删除）
#
# 架构（三段闭环）：
#   fock_sysdynamics → FockMPO（Fock 键直积 MPO，开链：无 ρ₀/迹闭合/插入）
#   sparsefockmpo `*` → 插入算符乘进链（左/右作用；纯张量收缩，无 graded 符号）
#   fockmpo2gmps → 迹闭合 + ρ₀/切位 0 三明治 + 键重封块 + oddlater → split_cell
#   + pairize（基线分站装置，符号规范逐位继承）→ reorder! → GMPS
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
# 语义（迁移期经 203 项对拍闸门对齐旧径逐位；旧径删除后由闸门接线对拍
# + observables golden U1GTEMPO2 机器精度锚承接）：
#   一阶（q=1）：τ/+ 支 mats = S·M、− 支 M·S；二阶（Weyl）：S·M·S + Weyl 标量
#   lscale = coupling_norm/4·Σ step²（_system_factor_fock/sysdynamics 层
#   setlogscale!——2026-11-18 起 coupling_norm 不进 fock_sysdynamics/fockmpo2gmps）。
#   ρ₀/切位 0 插入三明治（Ol₀·ρ·Or₀）由 fockmpo2gmps 的 rho0/insertions kwargs
#   折进首格核；切位 ≥ 1 插入走 sparsefockmpo `*`；oddlater Koszul 号
#   (−1)^{oddlater·popcount(m)} 由 fockmpo2gmps 按插入元数据统一施加（两次乘法
#   可组合：指数相加）。插入总荷必须为零（奇插入链的 Grassmann 荷流闭合）。
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
site k ↔ chronological_cells(lattice)[k]。logscale 不挂 FockMPO——Weyl 标量在
fockmpo2gmps 产出 GMPS 后由 _system_factor_fock/sysdynamics 层 setlogscale!。
ρ₀/插入在链上以
块乘法/转换 kwargs 表达，链本体保持开链（无 ρ₀/迹闭合，可复用于其它收缩口径）。
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
Base.copy(f::FockMPO) = FockMPO([copy(t) for t in f.data])   # 深拷贝：fock_initstate! 原地折块，浅拷贝会穿透

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
    fock_sysdynamics(lattice, model; order=1) -> FockMPO

杂质系统因子 K 的 Fock 直积 MPO（开链，纯 Fock 无 Grassmann 内容）：逐物理格
一 site，site 内块对角（每单项式 m 一个 (dim,dim) 算符块）。单元算符与
sysdynamics 同口径：一阶 τ/+ 支 S·M、− 支 M·S；二阶 S·M·S。**不含** ρ₀/
迹闭合/插入/oddlater/Weyl 标量——ρ₀ 与切位 0 插入走 fockmpo2gmps 的
rho0/insertions kwargs，切位 ≥ 1 插入走 sparsefockmpo `*`，oddlater 由
fockmpo2gmps 施加，Weyl 标量由调用方（_system_factor_fock/sysdynamics 层）
setlogscale! 补上。
"""
function fock_sysdynamics(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian;
                          order::Integer=1)
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
    return FockMPO(data)
end

# ---------------------------------------------------------------------------
# ρ₀ / 插入
# ---------------------------------------------------------------------------

"""
    fock_vacuumstate(lattice; initial_state=nothing) -> Matrix{ComplexF64}

Fock 侧初态/ρ₀（sysdynamics 的 ρ₀ 逻辑镜像）：imag/mixed = 热迹 I；real =
真空投影 δ₀₀，`initial_state`（FockMatrix/矩阵）给定时为该态。
"""
function fock_vacuumstate(lattice::AbstractGrassmannLattice; initial_state=nothing)
    dim = 2^lattice.bands
    if lattice isa AbstractRealGrassmannLattice
        rho = Matrix{ComplexF64}(Diagonal(vcat([1.0], zeros(dim - 1))))
        initial_state === nothing || (rho = initial_state isa FockMatrix ?
            Matrix{ComplexF64}(initial_state.data) : Matrix{ComplexF64}(initial_state))
        return rho
    end
    return Matrix{ComplexF64}(I, dim, dim)
end

"""
    fock_initstate!(lattice, f::FockMPO, rho0=nothing) -> FockMPO

    ρ₀ 折叠（原地）：ρ₀ 右乘进首 site 全部单项式块（=(M·ρ₀) 同语义）。
`rho0` 缺省 = fock_vacuumstate(lattice)。**切位 0 插入不经本函数**——
Ol₀·ρ·Or₀ 三明治序由 fockmpo2gmps 的 rho0/insertions kwargs 保证
（先 fold 后乘会得到 Ol₀·O₁·ρ 的错误序）；走 kwargs 路线时不要再调本函数。
"""
function fock_initstate!(lattice::AbstractGrassmannLattice, f::FockMPO, rho0=nothing)
    dim = 2^lattice.bands
    rho = rho0 === nothing ? fock_vacuumstate(lattice) :
          (rho0 isa FockMatrix ? Matrix{ComplexF64}(rho0.data) : Matrix{ComplexF64}(rho0))
    size(rho) == (dim, dim) || throw(DimensionMismatch("rho0 must be $dim×$dim"))
    z = one(TK.sectortype(lattice))
    blk = f.data[1][(z, z, TK.dual(z), TK.dual(z))]
    for m in 0:dim^2-1
        blk[:, m+1, m+1, :] .= blk[:, m+1, m+1, :] * rho
    end
    return f
end

"""
    sparsefockmpo(lattice, insertions) -> SparseFockMPO

    插入算符 → 稀疏算符链：切位 ≥ 1 的插入按 `insertion_cut_index` 映射到链位，
同位算符按 push 序折成列积（foldr/foldl，先 push 最左，0⁺ 序）。**切位 0 插入
在此过滤掉**（归 fockmpo2gmps 的 rho0/insertions 三明治）；全部插入须同 side
（左右各建一个 SparseFockMPO）。⚠️ 过滤按 cut index 而非 `boundary` 字段——
`:-` 支 boundary 0 = 链末切位（cut = 2Nt），不是 ρ₀ 切位（lesser(j=0) 路径）。"""
function sparsefockmpo(lattice::AbstractGrassmannLattice, insertions::Vector{OperatorInsertion})
    ins1 = filter(i -> insertion_cut_index(lattice, i) >= 1, insertions)
    side = isempty(ins1) ? :left : first(ins1).side
    all(i -> i.side === side, ins1) ||
        throw(ArgumentError("mixed sides: build one SparseFockMPO per side (:left / :right)"))
    G = TK.sectortype(lattice)
    dim = 2^lattice.bands
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
        OL = foldr((i, A) -> i.matrix * A, filter(i -> i.side === :left, bycut[p]);
                   init=I0)
        OR = foldl((A, i) -> A * i.matrix, filter(i -> i.side === :right, bycut[p]);
                   init=I0)
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
    fockmpo2gmps(f::FockMPO, lattice; rho0=nothing, insertions=OperatorInsertion[],
                 trunc=nothing, truncerrs=nothing) -> GMPS

FockMPO → GMPS（迹闭合 + 分站）：逐格把 FockMPO site 的 (dim,dim) 单项式块重排
进 (s,t) 对差键核（线程注入：cell1 = (M₁ρ₀)[s,t]、中间 = M[v,s]δ_{t,w}、末格 =
M[t,s] 右边界真空）→ split_cell（基线分站装置：荷差键 + 逐格 SVD 切分，符号规范
逐位继承）→ pairize → reorder!（Chrono→存储序）→ GMPS。

kwargs：
- `rho0`：ρ₀（缺省 = fock_vacuumstate(lattice)；imag/mixed = I、real = 真空投影）。
  切位 0 插入在此与 ρ₀ 三明治（Ol₀·ρ·Or₀，push 序 = 0⁺ 序）后折进首格核——
  与基线 (M₁·Ol₀·ρ·Or₀) 逐位同序。
- `insertions`：**全部**插入（含切位 0）——仅作标签元数据：prefix 荷累计（左键
  标签移位）+ oddlater Koszul 号 (−1)^{oddlater·popcount(m)}（oddlater = 切位 ≥ k
  的奇插入数）+ 总荷守卫（Σ charge = 0）。插入的矩阵本体应已乘进链（sparsefockmpo
  `*`），本函数不重复乘。
- `trunc`/`truncerrs`：reorder 后的 GM.truncate! 透传（缺省不截断）。

**Weyl 标量不在本函数**——调用方（_system_factor_fock/sysdynamics 层）持
coupling_norm，对产出的 GMPS 做 `logscale += lscale`（⚠️ 累加非覆盖：
truncate!/swap! 的压缩归一化因子已吸收进 logscale，E2 语义）。
"""
function fockmpo2gmps(f::FockMPO, lattice::AbstractGrassmannLattice;
                      rho0=nothing, insertions=OperatorInsertion[],
                      trunc=nothing, truncerrs=nothing)
    G = TK.sectortype(lattice)
    bands = lattice.bands
    dim = 2^bands
    z = one(G)
    cells = chronological_cells(lattice)
    L = length(cells)
    length(f.data) == L ||
        throw(DimensionMismatch("FockMPO has $(length(f.data)) sites, lattice has $L cells"))

    # 插入元数据：切位映射、总荷守卫
    cuts = Dict{Int,Vector{OperatorInsertion}}()
    total_q = z
    for ins in insertions
        push!(get!(cuts, insertion_cut_index(lattice, ins), OperatorInsertion[]), ins)
        total_q = qadd(total_q, ins.charge === nothing ? z : ins.charge)
    end
    total_q == z || throw(ArgumentError("total insertion charge must be zero"))

    # ρ₀ + 切位 0 三明治（Ol₀·ρ·Or₀；push 序 = 0⁺ 序）——过滤按 cut index：
    # `:-` 支 boundary 0 映射链末切位（cut = 2Nt），不是 ρ₀ 切位（与旧径 cuts
    # 键同口径）
    rho = rho0 === nothing ? fock_vacuumstate(lattice) :
          (rho0 isa FockMatrix ? Matrix{ComplexF64}(rho0.data) : Matrix{ComplexF64}(rho0))
    size(rho) == (dim, dim) || throw(DimensionMismatch("rho0 must be $dim×$dim"))
    rho = foldr((i, A) -> i.matrix * A,
                filter(i -> insertion_cut_index(lattice, i) == 0 && i.side === :left, insertions);
                init=rho)
    rho = foldl((A, i) -> A * i.matrix,
                filter(i -> insertion_cut_index(lattice, i) == 0 && i.side === :right, insertions);
                init=rho)

    # prefix 荷（标签口径 = leftlabels；rightlabels 由 (prefix, 切位荷) 唯一决定）
    prefix = foldl((p, ins) -> ins.charge === nothing ? p : qadd(p, ins.charge),
                   get(cuts, 0, OperatorInsertion[]); init=z)
    pairs_st = [(s, t) for t in 1:dim for s in 1:dim]   # 1-based，基线同款枚举序
    fqs1 = fock_charges(G, bands)
    chpc = convert(Vector{eltype(fqs1)}, cell_charges(G, bands))

    data = nothing
    for (k, (b, j)) in enumerate(cells)
        first_cell = k == 1
        last_cell = k == L
        leftlabels = first_cell ? [z] : [qsub(qsub(fqs1[s], fqs1[t]), prefix) for (s, t) in pairs_st]
        prefix_left = prefix
        for (cut, items) in cuts, ins in items
            cut == k || continue
            ins.charge === nothing || (prefix = qadd(prefix, ins.charge))
        end
        rightlabels = last_cell ? [z] : [qsub(qsub(fqs1[s], fqs1[t]), prefix) for (s, t) in pairs_st]
        # oddlater：切位 ≥ k 的奇插入数（切位 0 计入——基线同款）
        oddlater = 0
        for (cut, items) in cuts, ins in items
            cut >= k || continue
            q = ins.charge === nothing ? z : ins.charge
            oddlater += GM.fermionparity(q) ? 1 : 0
        end
        core = zeros(ComplexF64, length(leftlabels), dim^2, length(rightlabels))
        blk = f.data[k][(z, z, TK.dual(z), TK.dual(z))]   # (dim, 4^B, 4^B, dim)
        if L == 1
            for m in 0:dim^2-1
                core[1, m+1, 1] = tr(blk[:, m+1, m+1, :] * rho) *
                                  (-1)^(oddlater * count_ones(m))
            end
        elseif first_cell
            for m in 0:dim^2-1
                sgn = (-1)^(oddlater * count_ones(m))
                for (r, (s, t)) in enumerate(pairs_st)
                    # (M₁·ρ)[s, t]：ρ 已含切位 0 三明治
                    core[1, m+1, r] = sum(blk[s, m+1, m+1, u] * rho[u, t] for u in 1:dim) * sgn
                end
            end
        elseif last_cell
            for m in 0:dim^2-1
                sgn = (-1)^(oddlater * count_ones(m))
                for (l, (s, t)) in enumerate(pairs_st)
                    core[l, m+1, 1] = blk[t, m+1, m+1, s] * sgn
                end
            end
        else
            for m in 0:dim^2-1
                sgn = (-1)^(oddlater * count_ones(m))
                for (l, (s, t)) in enumerate(pairs_st), (r, (v, w)) in enumerate(pairs_st)
                    t == w && (core[l, m+1, r] = blk[v, m+1, m+1, s] * sgn)
                end
            end
        end
        column = _cell_tensors(G, core, leftlabels, rightlabels, chpc)
        isempty(column) && return monomial(G, repeat(chpc, L); coeff=0)
        data === nothing ? (data = copy(column)) : append!(data, column)
    end
    out = pairize(G, data)
    # 时序构造 → 存储序：每格 bands 站成块移动（GM.reorder! perm 语义 = 新链第 k 站
    # = 原链第 perm[k] 站；swap! 费米换序号已修复）
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
# 闭环装配：sysdynamics 公开入口的实现体（2026-11-18 起唯一路径——旧径随 V2
# 删除；plan/sysdynamics_fock_mpo_plan.md）
# ---------------------------------------------------------------------------

"切位 ≥ 1 插入乘进 Fock 链（左右各一 SparseFockMPO，切位 0 归 fockmpo2gmps
的 rho0/insertions 三明治）——test/api/fockmpo.jl 闸门 _wf_build 同配方。
乘法非原地：产物在支撑位新建张量、其余透传别名，输入链（含缓存链）不被污染。"
function _fockmpo_with_insertions(lattice::AbstractGrassmannLattice, Kf::FockMPO,
                                  insertions)
    ge1 = filter(i -> insertion_cut_index(lattice, i) >= 1, insertions)
    Wl = sparsefockmpo(lattice, filter(i -> i.side === :left, ge1))
    Wr = sparsefockmpo(lattice, filter(i -> i.side === :right, ge1))
    return Wr * (Wl * Kf)
end

"Weyl 标量 lscale（order=2 = coupling_norm/4·Σ step²，一阶 0）——旧径同式同序"
_fock_weyl_lscale(lattice::AbstractGrassmannLattice, order::Int, coupling_norm::Real) =
    order == 2 ? real(coupling_norm / 4 *
                      sum(_cell_step2(lattice, b) for (b, _) in chronological_cells(lattice))) :
    0.0

"""
    _system_factor_fock(lattice, model; order=1, coupling_norm=0.0, initial_state=nothing,
                        insertions=OperatorInsertion[], trunc=nothing, truncerrs=nothing,
                        spec=nothing) -> GMPS

sysdynamics 的实现体（`fock_sysdynamics` → 切位 ≥ 1 插入稀疏乘法 →
`fockmpo2gmps`：ρ₀/切位 0 三明治 + 迹闭合 + oddlater + prefix + reorder）。
**coupling_norm（Weyl 标量）只在此层输入**：order=2 要求 ≥ 0，经
`_fock_weyl_lscale` 算 lscale 后**累加**到产出 GMPS 的 logscale
（`logscale += lscale`——logscale 是链幅的 log 场，E2 语义下 truncate!/swap!
会把压缩归一化因子吸收进去，覆盖会抹掉这些因子；乘 exp(lscale) 语义 =
累加）。语义经 test/api/fockmpo.jl 对拍闸门锚定：
`initial_state` 仅实时轮廓生效（虚时/混合 = 热迹 I——此处按轮廓门控，
`fockmpo2gmps` 的 rho0 直通）；`spec` 兼容旧径 kwarg（须 = sectortype(lattice)）。
"""
function _system_factor_fock(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian;
                             order::Int=1, coupling_norm::Real=0.0, initial_state=nothing,
                             insertions=OperatorInsertion[], trunc=nothing, truncerrs=nothing,
                             spec=nothing)
    order in (1, 2) || throw(ArgumentError("order must be 1 or 2"))
    order == 2 && !(coupling_norm >= 0) &&
        throw(ArgumentError("order=2 requires coupling_norm ≥ 0 (total Σ_b α_b)"))
    spec === nothing || spec === TK.sectortype(lattice) ||
        throw(ArgumentError("spec $spec ≠ sectortype(lattice) $(TK.sectortype(lattice))"))
    Kf = _fockmpo_with_insertions(lattice,
                                  fock_sysdynamics(lattice, model; order=order),
                                  insertions)
    rho0 = lattice isa AbstractRealGrassmannLattice ? initial_state : nothing
    out = fockmpo2gmps(Kf, lattice; rho0=rho0, insertions=insertions,
                       trunc=trunc, truncerrs=truncerrs)
    GM.setlogscale!(out, GM.logscale(out) + _fock_weyl_lscale(lattice, order, coupling_norm))
    return out
end


