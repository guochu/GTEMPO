# foldedmpo.jl：折叠装配（plan/unified_chain_ordering_plan.md §7；2026-11-20）
#
# 装配内序（无 reorder!）：折叠恒等式
#   Tr[A⁺₁⋯A⁺_N A⁻_N⋯A⁻₁] = Tr[A⁻_N⋯A⁻₁ · ρ₀ · A⁺₁⋯A⁺_N] = Tr[ρ_N],
#   ρ_j = A⁻_j ρ_{j−1} A⁺_j
# 链站序 = 时间升序 = 实时段存储序（± 成对）。
#
# 折叠 FockMPO 结构：
#   real  : [𝓛₁ … 𝓛_Nt]——合并超算子格，键 = pair 态 (s,t)（维 dim²，
#           索引 ps = (t−1)·dim+s，与 pairs_st 枚举同序），单项式腿 =
#           双支槽积组合（维 4^{2B}，m = m⁻ + 4^B·m⁺），
#           blk[(s,t), m, m, (v,w)] = A⁻_{m⁻}[s,v]·A⁺_{m⁺}[w, t]；
#           ρ₀（initial_state/cut-0 三明治）折入格 1（rank-3：blk[(s,t),m,1] =
#           (A⁻ρ₀A⁺)[s,t]），Tr[ρ_Nt] 闭合在格 Nt（Nt=1 时 Tr[A⁻ρ₀A⁺]）。
#   mixed : [τ_N…τ₂（bra-side 格：blk[(s,t),m,m,(v,w)] = δ_{s,v}·A_τ[m][w,t]），
#            τ₁（rank-3 vec 头：blk[(s,t),m,1] = (A_τ₁ρ_s)[s,t]，ρ_s = cut-0
#            三明治）、𝓛₁ … 𝓛_Nt]——恒等式 vec(τ₁⋯τ_N) = (I⊗τ_Nᵀ)⋯(I⊗τ₂ᵀ)
#            (τ₁⊗I)|I⟩⟩：τ₁ ket-side（vec 头）、τ₂..τ_N bra-side（X·τ_j 追加），
#            τ 段输出 bond = vec(τP) = ρ₀ 槽（seam）。
#   imag  : 时序 per-支格（deprecated 同构）——v7 虚轴翻转后存储 = 时序，
#           converter 时序核直出（无 reorder!）。
#
# 键荷账：pair 差标签 q(s)−q(t)−prefix（prefix = 插入荷沿线程累计）；
# 合并格荷移 δ = q(m⁻)+q(m⁺)+Σq_ins、τ 格 δ = q(m)+Σq_ins（q_ins 两支同号
# 进入 pair 差）。Σδ = 0 环闭合守卫。

# ---------------------------------------------------------------------------
# 折叠荷账装置
# ---------------------------------------------------------------------------

"合并格 pair 槽荷表：m = m⁻ + 4^B·m⁺（低 bands 组 = − 支、高 = + 支），
q_m = q(m⁻) + q(m⁺)"
function _pair_slot_charges(::Type{G}, bands::Int) where {G<:Sector}
    n = 4^(2bands)
    out = Vector{typeof(one(G))}(undef, n)
    for m in 0:n-1
        q = one(G)
        for b in 1:2bands
            bb = b <= bands ? b : b - bands
            slot = (m >> (2 * (b - 1))) & 3
            slot == 1 && (q = qadd(q, var_charge(G, bb, bands)))
            slot == 2 && (q = qadd(q, qneg(var_charge(G, bb, bands))))
        end
        out[m+1] = q
    end
    return out
end

"pair 差荷：fp[(s,t)] = q(s) − q(t)（ps = (t−1)dim+s 枚举序）"
_pair_diffs(fqs1, dim) = [qsub(fqs1[s], fqs1[t]) for t in 1:dim for s in 1:dim]

"""
    _pair_site_charge_shift(blk, fp, qms, rank4::Bool) -> Sector

pair 版站荷移位读数：rank4 = true 时 `blk[ps, m, m, pw] ≠ 0` ⟹
fp[ps] − fp[pw] − q_m = δ；rank3（头格）时 `blk[ps, m, m, 1] ≠ 0` ⟹
fp[ps] − q_m = δ。站内一致性 + slot 对角断言同单线版。
"""
function _pair_site_charge_shift(blk, fp::Vector{G}, qms::Vector{G}, rank4::Bool) where {G<:Sector}
    D = length(fp)
    tol = 1e-12 * maximum(abs, blk)
    M = length(qms)
    for m1 in 1:M, m2 in 1:M, pw in (rank4 ? (1:D) : (1:1)), ps in 1:D
        m1 == m2 && continue
        abs(blk[ps, m1, m2, pw]) > tol && throw(ArgumentError(
            "monomial leg not slot-diagonal (off-diagonal |entry| > $tol)"))
    end
    δ = nothing
    for m in 1:M, pw in (rank4 ? (1:D) : (1:1)), ps in 1:D
        abs(blk[ps, m, m, pw]) > tol || continue
        d = rank4 ? qsub(qsub(fp[ps], fp[pw]), qms[m]) : qsub(fp[ps], qms[m])
        δ === nothing ? (δ = d) :
        δ == d || throw(ArgumentError(
            "pair charge bookkeeping inconsistent (charge-non-diagonal dressing?): " *
            "read δ = $δ and $d"))
    end
    return δ === nothing ? one(G) : δ
end

"折叠插入乘子：:left → O·A（foldr 首 push 最左）、:right → A·O（foldl）；
返回 Dict j => (OL, OR)"
function _folded_insertion_ops(ins_list, bands, dim)
    ol = Dict{Int,Matrix{ComplexF64}}()
    or = Dict{Int,Matrix{ComplexF64}}()
    I0 = Matrix{ComplexF64}(I, dim, dim)
    for ins in ins_list
        if ins.side === :left
            ol[ins.j] = haskey(ol, ins.j) ? _fockop_matrix(ins.op, bands) * ol[ins.j] :
                        _fockop_matrix(ins.op, bands) * I0
        else
            or[ins.j] = haskey(or, ins.j) ? or[ins.j] * _fockop_matrix(ins.op, bands) :
                        I0 * _fockop_matrix(ins.op, bands)
        end
    end
    return ol, or
end

# ---------------------------------------------------------------------------
# 折叠 fock_sysdynamics
# ---------------------------------------------------------------------------

"""
    fock_sysdynamics(lattice, model; order=1, insertions=OperatorInsertion[],
                     initial_state=nothing, gates=nothing) -> FockMPO

杂质系统因子的**折叠** Fock 直积 MPO（装配内序 = 存储序，converter 无
reorder!）：real = 合并超算子格链 [𝓛₁…𝓛_Nt]（ρ₀ 折入格 1、Tr[ρ_Nt] 闭合
格 Nt）；mixed = [τ_N…τ₂(bra-side) τ₁(vec 头) 𝓛₁…𝓛_Nt]（τ 段输出 bond =
ρ₀ 槽）；imag = 时序 per-支格（v7 翻转后时序 = 存储序，converter 直出）。
切位语义不变：(±,j) → 合并格 j、(τ,j) → τ 格 j、cut-0 → ρ₀ 三明治（real）/
τ₁ 格（mixed, imag）。

**调度门**（`gates::Vector{Float64}`，exactif_port_plan E5——严格核管线的
Trotter 调度自管，exactif.jl 传 `trotter_factors` 的 gate_durations 表）：
real/mixed 格逐合并格 j 的门时长 = `gates[j]`（两支同表；Td 模型时间参数 =
门左端点累计 Σᵢ<ⱼ gates[i]，− 支 = T − Σᵢ≤ⱼ gates[i]）。调度表自带 Strang
拆分（首尾半步）⟹ `order` 必须 = 1（q=1 反 Wick 测度、S·W / W·S 组合——
逐格内再劈半步会与表双重劈分）；imag 格不支持（热环无调度概念）。
"""
function fock_sysdynamics(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian;
                          order::Integer=1, insertions=OperatorInsertion[], initial_state=nothing,
                          gates::Union{Nothing,Vector{Float64}}=nothing)
    order in (1, 2) || throw(ArgumentError("order must be 1 or 2"))
    if gates !== nothing
        lattice isa AbstractImagGrassmannLattice && throw(ArgumentError(
            "scheduled gates (gates kwarg) apply to real/mixed contours only"))
        order == 1 || throw(ArgumentError(
            "scheduled gates carry their own Strang splitting (half-step table) — " *
            "order must be 1 (the gate table is the Trotter order)"))
        length(gates) == lattice.Nt || throw(DimensionMismatch(
            "gates must have Nt = $(lattice.Nt) entries, got $(length(gates))"))
    end
    bands = lattice.bands
    G = TK.sectortype(lattice)
    bands_guard(G, bands)
    num_bands(model) == bands ||
        throw(ArgumentError("model bands $(num_bands(model)) ≠ lattice bands $bands"))
    cells = chronological_cells(lattice)
    isempty(cells) && throw(ArgumentError("empty lattice"))
    weights = local_weights(G, bands, Int(order))
    dim = 2^bands
    z = one(G)
    nplus = 4^bands
    ins_real = filter(i -> i.branch !== :τ, insertions)
    ins_tau = filter(i -> i.branch === :τ, insertions)
    # 调度门传播子（real/mixed 合并格；时长 gates[j]、Td 时间参数 = 门左端点；
    # − 支时间参数 = T − Σᵢ≤ⱼ gates[i]（区间后端点，与均匀路由 (Nt−j)δt 同语义））
    scheduled_S = function (b::Symbol, j::Int)
        gates === nothing && return _fock_propagator_at(lattice, model, b, j, Int(order))
        T = sum(gates)
        t = b === :+ ? sum(gates[1:j-1]; init=0.0) :
            b === :- ? T - sum(gates[1:j]; init=0.0) : nothing
        return Matrix(_model_propagator(model, b, gates[j], t))
    end

    if lattice isa AbstractImagGrassmannLattice
        # ---- imag：时序 per-支格（deprecated 同构；ρ₀=I + cut-0 三明治折 (τ,1)）----
        focksp = _fock_bondspace(G, dim)
        monosp = _fock_bondspace(G, dim^2)
        rho_s = foldr((i, A) -> _fockop_matrix(i.op, bands) * A,
            filter(i -> insertion_cut_index(lattice, i) == 0 && i.side === :left, ins_tau);
            init=Matrix{ComplexF64}(I, dim, dim))
        rho_s = foldl((A, i) -> A * _fockop_matrix(i.op, bands),
            filter(i -> insertion_cut_index(lattice, i) == 0 && i.side === :right, ins_tau);
            init=rho_s)
        olt, ort = _folded_insertion_ops(filter(i -> insertion_cut_index(lattice, i) >= 1, ins_tau), bands, dim)
        data = [begin
            S = _fock_propagator_at(lattice, model, b, j, Int(order))
            mats = order == 2 ? [S * W * S for W in weights] : [S * W for W in weights]
            mats = haskey(olt, j) || haskey(ort, j) ?
                   [get(olt, j, Matrix{ComplexF64}(I, dim, dim)) * M *
                    get(ort, j, Matrix{ComplexF64}(I, dim, dim)) for M in mats] : mats
            mats = j == 1 ? [M * rho_s for M in mats] : mats   # ρ₀ 三明治后于 (τ,1) 插入（红线 1）
            oddlater = count(ins -> insertion_cut_index(lattice, ins) >= kk &&
                                 GM.fermionparity(_fockop_charge(ins.op, G, bands)), insertions)
            t = TK.zeros(ComplexF64, focksp ⊗ monosp ← monosp ⊗ focksp)
            blk = t[(z, z, TK.dual(z), TK.dual(z))]
            for m in 0:dim^2-1
                blk[:, m+1, m+1, :] .= oddlater == 0 ? mats[m+1] :
                                       mats[m+1] * (-1)^(oddlater * count_ones(m))
            end
            t
        end for (kk, (b, j)) in enumerate(cells)]
        return FockMPO(data)
    end

    # ---- real/mixed：折叠装配 ----
    dim2 = dim^2
    fockp2 = _fock_bondspace(G, dim2)
    monop2 = _fock_bondspace(G, nplus^2)
    monosp = _fock_bondspace(G, nplus)
    onesp = _fock_bondspace(G, 1)
    I0 = Matrix{ComplexF64}(I, dim, dim)

    cells_by = Dict{Tuple{Symbol,Int},Vector{Matrix{ComplexF64}}}()
    for (b, j) in cells
        S = scheduled_S(b, j)
        cells_by[(b, j)] = order == 2 ? [S * W * S for W in weights] :
                           b === :- ? [W * S for W in weights] : [S * W for W in weights]
    end
    olp, orp = _folded_insertion_ops(filter(i -> i.j >= 1, ins_real), bands, dim)
    olm, orm = _folded_insertion_ops(filter(i -> i.j >= 1, ins_real), bands, dim)
    getp = (d, j) -> get(d, j, I0)

    "合并格 𝓛_j：blk[(s,t),m,m,(v,w)] = A⁻[s,v]·A⁺[w,t]；head=true 时 rank-3
    ρ₀ 收缩（仅 real 格 1）"
    function merged_site(j::Int, head::Bool, rho0::Matrix{ComplexF64})
        Am = [get(olm, j, I0) * M * get(orm, j, I0) for M in cells_by[(:-, j)]]
        Ap = [get(olp, j, I0) * M * get(orp, j, I0) for M in cells_by[(:+, j)]]
        if head
            t = TK.zeros(ComplexF64, fockp2 ⊗ monop2 ← monop2 ⊗ onesp)
            blk = t[(z, z, TK.dual(z), TK.dual(one(G)))]
            for mn in 0:nplus-1, mp in 0:nplus-1
                m = mn + nplus * mp + 1
                B = Am[mn+1] * rho0 * Ap[mp+1]     # (A⁻ρ₀A⁺)[s,t]
                for t2 in 1:dim, s in 1:dim
                    blk[(t2-1)*dim+s, m, m, 1] = B[s, t2]
                end
            end
        else
            t = TK.zeros(ComplexF64, fockp2 ⊗ monop2 ← monop2 ⊗ fockp2)
            blk = t[(z, z, TK.dual(z), TK.dual(z))]
            for mn in 0:nplus-1, mp in 0:nplus-1
                m = mn + nplus * mp + 1
                An = Am[mn+1]; Apn = Ap[mp+1]
                for v in 1:dim, w in 1:dim, t2 in 1:dim, s in 1:dim
                    blk[(t2-1)*dim+s, m, m, (w-1)*dim+v] = An[s, v] * Apn[w, t2]
                end
            end
        end
        return t
    end

    if lattice isa AbstractRealGrassmannLattice
        # ρ₀（initial_state 仅实时覆盖）+ cut-0 三明治（在半支乘法之后——红线 1）
        rho0 = initial_state === nothing ?
               Matrix{ComplexF64}(Diagonal(vcat([1.0], zeros(dim - 1)))) :
               (initial_state isa FockMatrix ? Matrix{ComplexF64}(initial_state.data) :
                Matrix{ComplexF64}(initial_state))
        size(rho0) == (dim, dim) || throw(DimensionMismatch("initial state must be $dim×$dim"))
        rho0 = foldr((i, A) -> _fockop_matrix(i.op, bands) * A,
            filter(i -> i.j == 0 && i.side === :left, ins_real); init=rho0)
        rho0 = foldl((A, i) -> A * _fockop_matrix(i.op, bands),
            filter(i -> i.j == 0 && i.side === :right, ins_real); init=rho0)
        Nt = lattice.Nt
        f = FockMPO([merged_site(j, j == 1, rho0) for j in 1:Nt])
    else
        # ---- mixed ----
        Nτ = lattice.Nτ; Nt = lattice.Nt
        data = Any[]
        # τ_N…τ₂：bra-side 格（X·τ_j 追加；blk = δ_{s,v}·A_τ[m][w,t]）
        for jj in Nτ:-1:2
            At = [get(olt, jj, I0) * M * get(ort, jj, I0) for M in cells_by[(:τ, jj)]]
            t = TK.zeros(ComplexF64, fockp2 ⊗ monosp ← monosp ⊗ fockp2)
            blk = t[(z, z, TK.dual(z), TK.dual(z))]
            for mn in 0:nplus-1
                At_m = At[mn+1]     # (τ_jᵀ 语义：A[w,t] 提供 bra 线追加)
                for v in 1:dim, w in 1:dim, t2 in 1:dim, s in 1:dim
                    s == v || continue
                    blk[(t2-1)*dim+s, mn+1, mn+1, (w-1)*dim+v] = At_m[w, t2]
                end
            end
            push!(data, t)
        end
        # τ₁：rank-3 vec 头（cut-0 三明治 ρ_s = O_l·I·O_r 折入：vec((A₁ρ_s)[s,t])）
        rho_s = foldr((i, A) -> _fockop_matrix(i.op, bands) * A,
            filter(i -> i.j == 0 && i.side === :left, ins_tau); init=I0)
        rho_s = foldl((A, i) -> A * _fockop_matrix(i.op, bands),
            filter(i -> i.j == 0 && i.side === :right, ins_tau); init=rho_s)
        A1 = [get(olt, 1, I0) * M * get(ort, 1, I0) for M in cells_by[(:τ, 1)]]
        t = TK.zeros(ComplexF64, fockp2 ⊗ monosp ← monosp ⊗ onesp)
        blk = t[(z, z, TK.dual(z), TK.dual(one(G)))]
        for mn in 0:nplus-1
            B = A1[mn+1] * rho_s      # (A₁ρ_s)[s,t]
            for t2 in 1:dim, s in 1:dim
                blk[(t2-1)*dim+s, mn+1, mn+1, 1] = B[s, t2]
            end
        end
        push!(data, t)
        # 合并格 𝓛₁…𝓛_Nt（seam bond = τ 段输出 = ρ₀ 槽；head=false 全 rank-4）
        append!(data, [merged_site(j, false, I0) for j in 1:Nt])
        f = FockMPO([t for t in data])
    end

    # ---- P5 交叉自检：符号荷 vs pair δ 读数 ----
    fqs1 = fock_charges(G, bands)
    fp = _pair_diffs(fqs1, dim)
    qm_pair = _pair_slot_charges(G, bands)
    qm_single = _slot_charges(G, bands)
    is_mixed = lattice isa AbstractMixedGrassmannLattice
    Nτ = is_mixed ? lattice.Nτ : 0
    "插入 → 线程位置（切位 0 → 头格）"
    function pos_of(ins)
        if lattice isa AbstractRealGrassmannLattice
            ins.j == 0 ? 1 : ins.j
        else
            ins.branch === :τ ? (ins.j == 0 ? Nτ : Nτ - ins.j + 1) : Nτ + ins.j
        end
    end
    for (k, t) in enumerate(f.data)
        expected = z
        for ins in insertions
            pos_of(ins) == k || continue
            expected = qadd(expected, _fockop_charge(ins.op, G, bands))
        end
        blk = t[(z, z, TK.dual(z), TK.dual(z))]
        if lattice isa AbstractImagGrassmannLattice
            actual = _site_charge_shift(blk, fqs1, qm_single, bands)
        elseif lattice isa AbstractRealGrassmannLattice
            actual = _pair_site_charge_shift(blk, fp, qm_pair, k != 1)
        else
            actual = k == 1 ? _pair_site_charge_shift(blk, fp, qm_single, false) :
                     _pair_site_charge_shift(blk, fp, k <= Nτ ? qm_single : qm_pair, true)
        end
        actual == expected || throw(ArgumentError(
            "fock_sysdynamics: dressing charge mismatch at site $k — symbol charge " *
            "$expected vs chain δ reading $actual (materialization bug?)"))
    end
    return f
end

# ---------------------------------------------------------------------------
# 折叠 fockmpo2gmps：装配内序直出（无 reorder!、无 trunc——纯翻译；v7 起
# 三轮廓全直出，2026-11-21 用户裁决删 trunc 输入）
# ---------------------------------------------------------------------------

"""
    fockmpo2gmps(f::FockMPO, lattice) -> GMPS

折叠 FockMPO → GMPS（**装配内序直出，无 reorder!、无截断后处理**）：`f` 必须是
折叠 `fock_sysdynamics` 的产物（real/mixed 为 pair 态键折叠链、imag 为时序
per-支链——v7 虚轴翻转后 imag 时序 = 存储序）。逐格直出核（real/mixed 按
pair 态键线程；imag 时序核直出）→ split_cell + pairize → GMPS。converter 是
纯翻译：无 trunc 输入（截断属 IF 算法层/SVDCompression 与调用方，2026-11-21
用户裁决删除——套件从不经此管道传 trunc）。

站内自验证三件套（pair 版）：①站内 δ 读数全等（电荷不相干 dressing →
ArgumentError）；②Σδ = 0 环闭合；③slot 对角断言。
**Weyl 标量不在本函数**——调用方 setlogscale!（累加，E2 语义）。
"""
function fockmpo2gmps(f::FockMPO, lattice::AbstractGrassmannLattice)
    G = TK.sectortype(lattice)
    bands = lattice.bands
    dim = 2^bands
    z = one(G)
    fqs1 = fock_charges(G, bands)
    fp = _pair_diffs(fqs1, dim)
    qm_pair = _pair_slot_charges(G, bands)
    qm_single = _slot_charges(G, bands)
    chpc = convert(Vector{eltype(fqs1)}, cell_charges(G, bands))
    chpc2 = vcat(chpc, chpc)                      # 合并格变量序：+ 支全带 → − 支全带

    if lattice isa AbstractImagGrassmannLattice
        # ---- imag：时序 = 存储序（v7 虚轴翻转），时序核直出无 reorder! ----
        L = length(f.data)
        cells = chronological_cells(lattice)
        deltas = [_site_charge_shift(f.data[k][(z, z, TK.dual(z), TK.dual(z))],
                                     fqs1, qm_single, bands) for k in 1:L]
        total_q = foldl(qadd, deltas; init=z)
        total_q == z || throw(ArgumentError("total chain charge must be zero (Σδ ≠ 0)"))
        prefixes = Vector{typeof(z)}(undef, L)
        pp = z
        for k in 1:L
            pp = qadd(pp, deltas[k])
            prefixes[k] = pp
        end
        pairs = [(s, t) for t in 1:dim for s in 1:dim]
        data = nothing
        for (k, (b, j)) in enumerate(cells)
            first_cell = k == 1
            last_cell = k == L
            leftlabels = first_cell ? [z] :
                         [qsub(qsub(fqs1[s], fqs1[t]), prefixes[k-1]) for (s, t) in pairs]
            rightlabels = last_cell ? [z] :
                          [qsub(qsub(fqs1[s], fqs1[t]), prefixes[k]) for (s, t) in pairs]
            core = zeros(ComplexF64, length(leftlabels), dim^2, length(rightlabels))
            blk = f.data[k][(z, z, TK.dual(z), TK.dual(z))]
            if L == 1
                for m in 0:dim^2-1
                    core[1, m+1, 1] = tr(blk[:, m+1, m+1, :])
                end
            elseif first_cell
                for m in 0:dim^2-1, (r, (s, t)) in enumerate(pairs)
                    core[1, m+1, r] = blk[s, m+1, m+1, t]
                end
            elseif last_cell
                for m in 0:dim^2-1, (l, (s, t)) in enumerate(pairs)
                    core[l, m+1, 1] = blk[t, m+1, m+1, s]
                end
            else
                for m in 0:dim^2-1, (l, (s, t)) in enumerate(pairs), (r, (v, w)) in enumerate(pairs)
                    t == w && (core[l, m+1, r] = blk[v, m+1, m+1, s])
                end
            end
            column = _cell_tensors(G, core, leftlabels, rightlabels, chpc)
            isempty(column) && return monomial(G, repeat(chpc, L); coeff=0)
            data === nothing ? (data = copy(column)) : append!(data, column)
        end
        out = pairize(G, data)
        return out
    end
    # ---- real/mixed：折叠 pair 线程直出 ----
    is_mixed = lattice isa AbstractMixedGrassmannLattice
    Nτ = is_mixed ? lattice.Nτ : 0
    L = length(f.data)
    deltas = let
        out = Vector{typeof(z)}(undef, L)
        for k in 1:L
            blk = f.data[k][(z, z, TK.dual(z), TK.dual(z))]
            rank4 = k != 1                               # rank-3 头（mixed τ_N vec / real ρ₀）
            qms = is_mixed ? (k == 1 ? qm_single : (k <= Nτ ? qm_single : qm_pair)) : qm_pair
            out[k] = _pair_site_charge_shift(blk, fp, qms, rank4)
        end
        out
    end
    total_q = foldl(qadd, deltas; init=z)
    total_q == z || throw(ArgumentError("total chain charge must be zero (Σδ ≠ 0)"))
    prefixes = Vector{typeof(z)}(undef, L)
    p = z
    for k in 1:L
        p = qadd(p, deltas[k])
        prefixes[k] = p
    end
    pair_labels(j) = [qsub(qsub(fqs1[s], fqs1[t]), prefixes[j]) for (s, t) in pairs_st(dim)]
    cores = Any[]
    for k in 1:L
        blk = f.data[k][(z, z, TK.dual(z), TK.dual(z))]
        lo = k == 1 ? [z] : pair_labels(k - 1)
        ro = k == L ? [z] : pair_labels(k)
        τ_site = is_mixed && k <= Nτ
        ch = is_mixed ? (τ_site ? chpc : chpc2) : chpc2
        head = k == 1
        tail = k == L
        rank3 = head   # rank-3 头（mixed τ_N vec / real ρ₀ 收缩）
        S = size(blk, 2)                                     # 站 slot 维（合并 4^{2B} / τ 4^B）
        if rank3 && tail                                         # L==1：Tr[A⁻ρ₀A⁺]（real）或 Tr[ρ₁] 闭合
            core = zeros(ComplexF64, 1, S, 1)
            for m in 1:S, s in 1:dim
                core[1, m, 1] = blk[(s-1)*dim+s, m, m, 1]
            end
        elseif rank3                                             # rank-3 头直读
            core = zeros(ComplexF64, 1, S, length(ro))
            for m in 1:S, (r, (s, t)) in enumerate(pairs_st(dim))
                core[1, m, r] = blk[(t-1)*dim+s, m, m, 1]
            end
        elseif tail                                              # Tr[ρ_L] 闭合
            core = zeros(ComplexF64, length(lo), S, 1)
            for m in 1:S, (l, (v, w)) in enumerate(pairs_st(dim)), s in 1:dim
                core[l, m, 1] += blk[(s-1)*dim+s, m, m, (w-1)*dim+v]
            end
        else                                                     # rank-4 直读
            core = zeros(ComplexF64, length(lo), S, length(ro))
            for m in 1:S, (l, (v, w)) in enumerate(pairs_st(dim)), (r, (s, t)) in enumerate(pairs_st(dim))
                core[l, m, r] = blk[(t-1)*dim+s, m, m, (w-1)*dim+v]
            end
        end
        push!(cores, (core, lo, ro, ch))
        end

    data = nothing
    for (core, lo, ro, ch) in cores
        column = _cell_tensors(G, core, lo, ro, ch)
        isempty(column) && return monomial(G, repeat(chpc, length(lattice) ÷ (2bands)); coeff=0)
        data === nothing ? (data = copy(column)) : append!(data, column)
    end
    out = pairize(G, data)
    return out
end

"pair 态枚举（ps = (t−1)·dim+s，与核行/列序一致）"
pairs_st(dim) = [(s, t) for t in 1:dim for s in 1:dim]
