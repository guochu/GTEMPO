# ============================================================================
# system.jl — Weyl（二阶）/ 一阶系统因子 K
#
# 算法（benchmark/weyl-secondorder.html；对齐 U1GTEMPO2 system_factor_reference）：
#
#   二阶：U₂(λ) = exp(λ²c/4) · e^{−λH_S/2} · Q_W[e^{−λH_Bᴳ}] · e^{−λH_S/2}
#   一阶：U₁(λ) = e^{−λH_S} · R₁(λ)（Lie 顺序，反支为 M·S）
#
# 局部权重用推导的裸 Berezin 测度版本（GTEMPO 的 integrate 是裸配对：一阶管线
# 中相干态高斯权重放在传播子 overlap 因子里，而 Weyl K 是无 overlap 的
# automaton 链），因此无需 U1GTEMPO2 的 Gaussian 测度适配（aā → I/2+n）：
#
#   M_q = q·I − n − a·d† + ā·d + aā·I,   q = 1/2（Weyl）、1（一阶）
#
# 折叠链（实时每顶点折叠正反支，"X → A₊XA₋"；虚时无反支）：
#   K = Tr_Fock[ A_N⋯A_1 ρ₀ B₁⋯B_N ]
# 第一成员沿 site 顺序（顶点 N → 1）穿线 A_j = s·S₊M_jS₊；第二成员反向穿线
# B_j = s·S₋M_jS₋（与 U1GTEMPO2 相同的无转置形式；实 Fock 数据下
# Mᵀ = M†、S₋ᵀ = S₋ 保证与 𝒰† 严格等价）。单元 1 折叠 ρ₀，单元 N 闭合迹。
# ============================================================================

# ---- Fock 空间算符（band 1 = 最高有效位，FockMatrix 约定） ----------------
# JW 约定（2026-11-05 与 U1GTEMPO2 fock_operators 对齐）：带 b 的算符穿越
# 所有更低索引带（b' < b）——band1 无符号、band2 带 (−1)^{n₁}。旧实现穿越
# 更高索引带（恰好反转），单带不受影响、M≥2 与 U1 槽矩阵冲突。
function _fock_annihilators(M::Int)
    dim = 2^M
    d = [zeros(ComplexF64, dim, dim) for _ in 1:M]
    for b in 1:M
        bit = 1 << (M - b)
        mask_cross = (dim - 1) & ~(bit - 1) & ~bit   # 更低索引带的位
        for st in 0:dim-1
            (st & bit) != 0 || continue
            d[b][(st ⊻ bit)+1, st+1] = iseven(count_ones(st & mask_cross)) ? 1.0 : -1.0
        end
    end
    return d
end

"""
    local_weight(bands, order; mode=:bare) -> Vector{Matrix{ComplexF64}}

槽算符：每带 4 个候选 [q·I−n, −d†, d, W₄]（q = 1/2 Weyl / 1 一阶），按带序
乘积展开成 4^bands 个 2^bands×2^bands 矩阵。槽编码 m 的每带 2bit：
0=(无), 1=(a), 2=(ā), 3=(aā)，高位带在前。

`mode` / `variant` 选择测度适配（ED 校准，见 validation/calibrate.jl）：
- `mode=:gaussian`（默认）：W₄ = I/2+n（Weyl）/ n（一阶）—— 与 U1GTEMPO2
  相同的 Gaussian 测度适配。
- `variant=:dag`（默认）：耦合槽 a ↦ +d†、ā ↦ +d。HTML 推导的 M_W 记号
  （a 系数 = −d†）对应 U1GTEMPO2 的 Gaussian Berezin 配对；GTEMPO 的
  `g_trace_phy` 分级 trace 配对对奇宇称扇区多出一个 −1，等价于把 a-单
  monomial 的系数翻号。`:std`/`:swap` 等其余变体仅作对照。
"""
function local_weight(bands::Int, order::Int; mode::Symbol = :gaussian,
                      variant::Symbol = :dag)
    order in (1, 2) || throw(ArgumentError("order 必须为 1 或 2"))
    mode in (:bare, :gaussian) || throw(ArgumentError("mode 必须为 :bare 或 :gaussian"))
    d = _fock_annihilators(bands)
    dim = 2^bands
    Id = Matrix{ComplexF64}(I, dim, dim)
    q = order == 2 ? 0.5 : 1.0
    weight = Vector{Matrix{ComplexF64}}(undef, 4^bands)
    for m in 0:(4^bands)-1
        O = Id
        for b in 1:bands
            s = (m >> (2 * (bands - b))) & 3
            nb = d[b]' * d[b]
            W4 = mode === :bare ? Id : (order == 2 ? (0.5 * Id + nb) : nb)
            # a ↦ +d†（GTEMPO 分级 trace 配对的适配，见上方 docstring）；
            # variant kwarg 仅供校准对照，默认 :dag 即正确约定
            Mb = if s == 0
                q * Id - nb
            elseif s == 1
                variant === :std ? (-d[b]') : d[b]'
            elseif s == 2
                variant === :std ? d[b] : d[b]
            else
                W4
            end
            O = O * Mb
        end
        weight[m+1] = O
    end
    return weight
end

# ---- 插入算符 -------------------------------------------------------------
"""
    WeylInsertion(branch, boundary, matrix; side=:left)

物理算符插入：在 `branch` 分支的整数边界 `boundary`（0..N）处乘进链。
`side=:left`（默认）乘在 ρ₀/A 侧（`O·X`），`side=:right` 乘在 B 侧
（`X·O`，即同一物理边界的反向支读法）。对齐 U1GTEMPO2 的 OperatorInsertion
（基准用 2^bands 维 Fock 矩阵）。
"""
struct WeylInsertion
    branch::Symbol
    boundary::Int
    matrix::Matrix{ComplexF64}
    side::Symbol
    function WeylInsertion(branch::Symbol, boundary::Int, matrix::AbstractMatrix;
                           side::Symbol = :left)
        side in (:left, :right) || throw(ArgumentError("side 必须为 :left 或 :right"))
        new(branch, boundary, Matrix{ComplexF64}(matrix), side)
    end
end

# ---- 单元核装配 -----------------------------------------------------------
# 变量 site 顺序（与格点内顺序一致）：
#   实时：(a⁺₁,ā⁺₁,…,a⁺_M,ā⁺_M, a⁻₁,ā⁻₁,…,a⁻_M,ā⁻_M)
#   虚时：(a₁,ā₁,…,a_M,ā_M)
# 带 b 的槽值 = a_bit + 2·ā_bit，槽编码 = Σ_b slot_b · 4^(M−b)。
function _code_from_ys(ys, off::Int, M::Int)
    c = 0
    for b in 1:M
        c += (ys[off + 2b - 1] + 2 * ys[off + 2b]) << (2 * (M - b))
    end
    return c
end

function _slot_operators(M::Int, order::Int; mode::Symbol = :gaussian)
    # 每带 4 槽基矩阵（slot ∈ {0: qI−n, 1: a†(:dag), 2: a, 3: aā}）。
    # 注意：不能用 2-bit 槽码去索引 local_weight 的 4^M 组合码列表——单带时
    # 二者恰好重合，M≥2 会取错矩阵并重复计入 qI−n 因子（2026-11-05 ③ 定位）。
    # 结构 = Π_b locals[b][slot_b]（按带序普通张量积）。
    # 2026-11-05 裁决：U1GTEMPO2 的带间 graded 符号 (−1)^{|a||b|} 属于其
    # Gaussian 测度配对 + 取向反转管线（其 IF = e^{−Δᵀ}），与本仓库裸 Berezin
    # 管线不匹配——照搬即符号重复计入。U=0 因子化检验（k_factor_probe：
    # 联合链 vs K₁∧K₁ 逐支撑 0 失配）实证本管线正确组合 =
    # JW 穿越（_fock_annihilators，保留）+ 无 graded⊗ 符号（此处移除）。
    # M=1 时该符号本就不生效，单带基线不受影响。
    dim = 2^M
    mode in (:bare, :gaussian) || throw(ArgumentError("mode 必须为 :bare 或 :gaussian"))
    d = _fock_annihilators(M)
    q = order == 2 ? 0.5 : 1.0
    locals = [let nb = d[b]' * d[b]
                  W4 = mode === :bare ? Matrix{ComplexF64}(I, dim, dim) :
                       Matrix{ComplexF64}(order == 2 ? 0.5 * I + nb : nb)
                  [Matrix{ComplexF64}(q * I - nb), d[b]', d[b], W4]
              end for b in 1:M]
    ops = Matrix{ComplexF64}[]
    for code in 0:(4^M)-1
        O = Matrix{ComplexF64}(I, dim, dim)
        for b in 1:M
            sb = (code >> (2 * (M - b))) & 3
            O = O * locals[b][sb + 1]
        end
        push!(ops, O)
    end
    return ops
end

"""
    system_factor(c::AbstractWeylContour, model; order=2, coupling_norm=0,
                  rho0=nothing, insertions=WeylInsertion[]) -> GrassmannMPS

在 Weyl 轮廓 `c` 上构造系统因子 K（GTEMPO GrassmannMPS，与
`influence_functional` 的 IF 一起交给 GTEMPO 的 `integrate` 收缩）。

- `model`：GTEMPO `ConstImpurityHamiltonian`（ToulouseIM / AndersonIM 等），
  通过 `fock_propagator` 取精确 Fock 传播子。
- `order=2`：Weyl 权重 + 对称半步 + 标量补偿 exp(λ²c/4)；`order=1`：
  anti-Wick 权重 + Lie 全步。一、二阶共用同一 IF。
- `rho0`：实时初态 Fock 矩阵（默认真空投影子）；虚时恒为 I（热迹）。
- `insertions`：观测量分子的物理算符插入。
"""
function system_factor(c::AbstractWeylContour, model; order::Int = 2,
                       coupling_norm::Real = 0.0,
                       rho0::Union{Nothing,AbstractMatrix} = nothing,
                       insertions::AbstractVector{WeylInsertion} = WeylInsertion[],
                       mode::Symbol = :gaussian, variant::Symbol = :dag,
                       oddlater_mode::Symbol = :none,
                       reorder_sign::Bool = true)
    order in (1, 2) || throw(ArgumentError("order 必须为 1 或 2"))
    if c isa WeylRealTime
        return _system_factor_rt(c, model; order, coupling_norm,
                                 rho0, insertions, mode, variant,
                                 oddlater_mode, reorder_sign)
    elseif c isa WeylMixedTime
        return _system_factor_mixed(c, model; order, coupling_norm,
                                    rho0, insertions, mode, variant,
                                    oddlater_mode, reorder_sign)
    else
        return _system_factor_tau(c, model; order, coupling_norm,
                                  rho0, insertions, mode, variant)
    end
end

"""
    _system_factor_rt — 实时 K：U1GTEMPO2 `system_factor_reference` 的时序链
    逐行翻译（Z2 宇称投影）。

时序单元：(+,1), …, (+,N), (−,N), …, (−,1)；切位 (+,j) → j、(−,j) → 2N−j+1。
每单元（VERBATIM U1 语义）：
    mats = order==2 ? [S·M·S] : (−支 ? [M·S] : [S·M])，S = e^{−lam·H/(2 或 1)}
    cutop(k)：该切位插入矩阵左乘；ρ₀ ← cutop(0)·ρ₀
    oddlater = 切位 ≥ k 的奇插入数；mats[m+1] ×= (−1)^{oddlater·popcount(m)}
    核：k==1：[1, m+1, (s,t)] = (M·ρ₀)[s,t]（左真空）
        k==L：[(s,t), m+1, 1] = M[t,s]（transpose-form 闭合，右真空）
        内部：t == w && core[(s,t), m+1, (v,w)] = M[v,s]（第二成员恒等穿线）
键宇称 = U1 电荷标签 (s−t−prefix) 的宇称投影 = pairpar(s,t) ⊕ 累积奇插入
宇称 F（prefix₀ = 切位 0 插入）。
单元张量直接铺入 vacuum GMPS 的格点位置（tile_cells!；与虚时路径同构）。
"""
function _system_factor_rt(c::WeylRealTime, model; order::Int = 2,
                           coupling_norm::Real = 0.0,
                           rho0::Union{Nothing,AbstractMatrix} = nothing,
                           insertions::AbstractVector{WeylInsertion} = WeylInsertion[],
                           mode::Symbol = :gaussian, variant::Symbol = :dag,
                           oddlater_mode::Symbol = :none,
                           reorder_sign::Bool = true)
    N = c.Nt
    M = c.bands
    dim = 2^M
    dim2 = dim * dim
    L = 2N                                   # 时序单元数
    weight = local_weight(M, order; mode = mode)
    slotops = _slot_operators(M, order; mode = mode)
    lattice = weyl_lattice(c)
    h = c.dt
    sc = order == 2 ? exp(-real(h^2) * coupling_norm * M / 4) : 1.0   # exp(λ²c/4), λ=±i h；c = 总耦合 Σ_b α_b

    # ---- 插入按切位分桶（U1 insertion_cut：+j → j；−j → 2N−j+1）----
    cutbucket = Dict{Int,Vector{Matrix{ComplexF64}}}()
    for O in insertions
        cut = O.branch === :+ ? O.boundary : 2N - O.boundary   # U1: ni+2Nt−boundary
        push!(get!(cutbucket, cut, Matrix{ComplexF64}[]), Matrix{ComplexF64}(O.matrix))
    end
    isoddop(Mx) = any(!iszero(Mx[i, j]) &&
                      paritybits(i - 1) != paritybits(j - 1)
                      for i in 1:size(Mx, 1), j in 1:size(Mx, 2))
    oddcount(cc::Int) = count(isoddop, get(cutbucket, cc, Matrix{ComplexF64}[]))
    function cutop(cc::Int)
        O = Matrix{ComplexF64}(I, dim, dim)
        for Mx in get(cutbucket, cc, Matrix{ComplexF64}[])
            O = Mx * O
        end
        return O
    end

    # ---- ρ₀（实时：初始概率对角；切位 0 插入左乘；prefix₀ = 切位 0 奇宇称）----
    rho = isnothing(rho0) ? diagm(0 => ComplexF64[1.0, fill(0.0, dim - 1)...]) :
                            Matrix{ComplexF64}(rho0)
    size(rho) == (dim, dim) || throw(DimensionMismatch("rho0 必须为 $dim×$dim"))
    rho = cutop(0) * rho
    F₀ = oddcount(0) % 2                    # prefix₀ 的宇称

    # ---- 逐时序格准备（mats / 宇称标签，时序角色不变）----
    cell_branch(k::Int) = k <= N ? :+ : :-
    cell_boundary(k::Int) = k <= N ? k : 2N - k + 1
    pairs = [(s, t) for t in 1:dim for s in 1:dim]   # t 大股，s 最快（U1 约定）
    _pairpar(i) = paritybits((i - 1) % dim) ⊕ paritybits((i - 1) ÷ dim)

    mats_k = Vector{Vector{Matrix{ComplexF64}}}(undef, L)
    inpar_k = Vector{Vector{Int}}(undef, L)
    outpar_k = Vector{Vector{Int}}(undef, L)
    F = F₀                                   # prefix 宇称（累积到当前切位之前）
    for k in 1:L
        b = cell_branch(k)
        hfrac = order == 2 ? 0.5 : 1.0
        S = Matrix{ComplexF64}(GTEMPO.fock_propagator(model, b, h * hfrac).data)
        mats = [order == 2 ? S * slotops[c2+1] * S :
                (b === :- ? slotops[c2+1] * S : S * slotops[c2+1]) for c2 in 0:(4^M)-1]
        order == 2 && (mats .*= sc)
        # cutop(k)：该切位插入矩阵左乘（U1 VERBATIM）
        Mk0 = cutop(k)
        mats = [Mk0 * m for m in mats]
        # oddlater 手工符号：默认 :none——退耦极限标定（α=1e-10，gt/lt/occ
        # 全部与解析值 1e-14/15 一致）：分级收缩已全权供给，手工再加即错。
        if oddlater_mode !== :none
            oddlater = 0
            lo = oddlater_mode === :incl ? k : k + 1
            for k2 in lo:L
                oddlater += oddcount(k2)
            end
            isodd(oddlater) && (mats = [paritybits(c2) == 1 ? (-m) : m
                                        for (c2, m) in enumerate(mats)])
        end
        f_k = oddcount(k) % 2                # 本格奇插入（含 cutop 左乘的）
        mats_k[k] = mats
        inpar_k[k] = (k == 1) ? [0] : [_pairpar(i) ⊕ F for i in 1:dim2]
        outpar_k[k] = (k == L) ? [0] : [_pairpar(i) ⊕ (F ⊕ f_k) for i in 1:dim2]
        F = F ⊕ f_k
    end

    # ---- 时序单元核装配（VERBATIM U1 语义；ChronoKeldysh 下时序序 = 存储序）----
    tensors = TensorMap[]
    for k in 1:L
        mats = mats_k[k]
        inpar = inpar_k[k]
        outpar = outpar_k[k]
        insize = (k == 1) ? 1 : dim2
        outsize = (k == L) ? 1 : dim2
        V = 2M
        C = zeros(ComplexF64, insize, ntuple(Returns(2), V)..., outsize)
        for m in 0:(4^M)-1
            Mk = mats[m+1]
            ysv = ntuple(v -> (m >> (v - 1)) & 1, V)
            # 格间重排簿记补偿（2026-11-04 K 系数交叉核对定位）：本构造相对
            # U1 多出 (−1)^{C(m,2)}，m = 奇选取格数；按 Σ_q [n_q odd]·[前置
            # 变量宇称奇] 逐格施加。reorder_sign=false 恢复旧行为（对照用）。
            cellodd = reorder_sign ? Int(paritybits(m)) : 0
            if k == 1
                for (r, (s, t)) in enumerate(pairs)
                    C[1, ntuple(x -> ysv[x] + 1, V)..., r] = (Mk * rho)[s, t]
                end
            elseif k == L
                for (l, (s, t)) in enumerate(pairs)
                    C[l, ntuple(x -> ysv[x] + 1, V)..., 1] =
                        (cellodd * inpar[l] == 1 ? -1 : 1) * Mk[t, s]
                end
            else
                for (l, (s, t)) in enumerate(pairs), (r, (v, w)) in enumerate(pairs)
                    t == w && (C[l, ntuple(x -> ysv[x] + 1, V)..., r] =
                        (cellodd * inpar[l] == 1 ? -1 : 1) * Mk[v, s])
                end
            end
        end
        all(iszero, C) && error("实时 K：单元 $k 的核全零")
        append!(tensors, split_core(CellCore(C, inpar, outpar)))
    end

    # ---- 时序直铺（零 permute!：存储序 = 链的本征序，MPS 即时序链本身）----
    # 每格 2M 个位置、(a,ā) 奇偶对齐 ⟹ AC 积分器的位置对 = 共轭对且不跨格。
    Tc = typeof(first(tensors))
    return GrassmannMPS(convert(Vector{Tc}, tensors))
end

"""
    _system_factor_mixed — 混合（Kadanoff–Baym）轮廓 K：τ 段 + 实时 +/− 段

时序单元序：τ₁…τ_Nτ → (+,1)…(+,Nt) → (−,Nt)…(−,1)，单 (s,t) 线穿线、真空
开边界（τ₁ 时序-in = 全 Fock 热迹折叠 ρ₀=I；(−,Nt)…… 时序末格（−,1）迹闭合）。
τ 格用虚时权重（order 2：sc_τ·S·M·S，S = e^{−hτH/2}）；实时格用实时权重
（order 2：sc_rt·S·M·S，S = e^{∓ihtH/2}；order 1：τ/+(S·M)、−(M·S)）。
Weyl 标量按段：τ 格 exp(hτ²·Σα/4)、rt 格 exp(−ht²·Σα/4)（order 2）。
ChronoMixed 时序直铺（零 permute!）。
"""
function _system_factor_mixed(c::WeylMixedTime, model; order::Int = 2,
                              coupling_norm::Real = 0.0,
                              rho0::Union{Nothing,AbstractMatrix} = nothing,
                              insertions::AbstractVector{WeylInsertion} = WeylInsertion[],
                              mode::Symbol = :gaussian, variant::Symbol = :dag,
                              oddlater_mode::Symbol = :none,
                              reorder_sign::Bool = true)
    Nτ = c.Nτ
    Nt = c.Nt
    M = c.bands
    dim = 2^M
    dim2 = dim * dim
    hτ = c.dτ
    ht = c.dt
    scτ = order == 2 ? exp(hτ^2 * coupling_norm * M / 4) : 1.0
    scrt = order == 2 ? exp(-real(ht^2) * coupling_norm * M / 4) : 1.0
    slotops = _slot_operators(M, order; mode = mode)
    lattice = weyl_lattice(c)

    # ---- 插入按切位分桶（切位 = 时序序号：τ j → j；+ j → Nτ+j；− j → Nτ+2Nt−j+1）----
    cutbucket = Dict{Int,Vector{Matrix{ComplexF64}}}()
    for O in insertions
        cut = O.branch === :τ ? O.boundary :
              O.branch === :+ ? Nτ + O.boundary : Nτ + 2Nt - O.boundary + 1
        push!(get!(cutbucket, cut, Matrix{ComplexF64}[]), Matrix{ComplexF64}(O.matrix))
    end
    isoddop(Mx) = any(!iszero(Mx[i, j]) &&
                      paritybits(i - 1) != paritybits(j - 1)
                      for i in 1:size(Mx, 1), j in 1:size(Mx, 2))
    oddcount(cc::Int) = count(isoddop, get(cutbucket, cc, Matrix{ComplexF64}[]))
    function cutop(cc::Int)
        O = Matrix{ComplexF64}(I, dim, dim)
        for Mx in get(cutbucket, cc, Matrix{ComplexF64}[])
            O = Mx * O
        end
        return O
    end

    # ρ₀ = I（τ 链起点 = 全 Fock 热迹；混合轮廓的热平衡由 τ 段本身建立）；
    # 切位 0（τ 边界 0）插入折叠进 R；prefix₀ = 其奇宇称。
    R = Matrix{ComplexF64}(I, dim, dim)
    for O in get(cutbucket, 0, Matrix{ComplexF64}[])
        R = O * R
    end
    F₀ = oddcount(0) % 2

    hfrac = order == 2 ? 0.5 : 1.0
    pairs = [(s, t) for t in 1:dim for s in 1:dim]   # t 大股，s 最快（U1 约定）
    _pairpar(i) = paritybits((i - 1) % dim) ⊕ paritybits((i - 1) ÷ dim)

    # ---- 时序单元定义（分支 / 边界 / 传播子 / 权重标量）----
    L = Nτ + 2Nt
    cell_branch(k::Int) = k <= Nτ ? :τ : (k <= Nτ + Nt ? :+ : :-)
    cell_boundary(k::Int) = k <= Nτ ? k : (k <= Nτ + Nt ? k - Nτ : L - k + 1)

    gmps = vacuumstate(ComplexF64, lattice)
    F = F₀
    for k in 1:L
        b = cell_branch(k)
        j = cell_boundary(k)
        hfrac_local = order == 2 ? 0.5 : 1.0
        if b === :τ
            S = Matrix{ComplexF64}(GTEMPO.fock_propagator(model, :τ, hτ * hfrac_local).data)
            sc = scτ
        else
            S = Matrix{ComplexF64}(GTEMPO.fock_propagator(model, b, ht * hfrac_local).data)
            sc = scrt
        end
        mats = [order == 2 ? S * slotops[c2+1] * S :
                (b === :− || b === :- ? slotops[c2+1] * S : S * slotops[c2+1]) for c2 in 0:(4^M)-1]
        order == 2 && (mats .*= sc)
        Mk0 = cutop(k)
        mats = [Mk0 * m for m in mats]
        if oddlater_mode !== :none
            oddlater = 0
            lo = oddlater_mode === :incl ? k : k + 1
            for k2 in lo:L
                oddlater += oddcount(k2)
            end
            isodd(oddlater) && (mats = [paritybits(c2) == 1 ? (-m) : m
                                        for (c2, m) in enumerate(mats)])
        end
        f_k = oddcount(k) % 2
        # τ 首格：ρ₀ 折叠（真空 in）；混合末格（−,1）：迹闭合（真空 out）；
        # 其余：dim2 穿线。
        is_fold = (k == 1)
        is_close = (k == L)
        inpar = is_fold ? [0] : [_pairpar(i) ⊕ F for i in 1:dim2]
        outpar = is_close ? [0] : [_pairpar(i) ⊕ (F ⊕ f_k) for i in 1:dim2]
        insize = is_fold ? 1 : dim2
        outsize = is_close ? 1 : dim2
        V = 2M
        C = zeros(ComplexF64, insize, ntuple(Returns(2), V)..., outsize)
        for m in 0:(4^M)-1
            Mk = mats[m+1]
            ysv = ntuple(v -> (m >> (v - 1)) & 1, V)
            # 格间重排簿记补偿（2026-11-04 rt 定位 + 2026-11-05 mixed 扫描实证）：
            # (−1)^{C(m,2)} 属于分级装配序本身（不分支）：混合链的全部单元
            # （τ 格与 rt 格）均施加，reorder_sign=false 恢复旧行为（对照用）。
            cellodd = reorder_sign ? Int(paritybits(m)) : 0
            if is_fold
                for (r, (s, t)) in enumerate(pairs)
                    C[1, ntuple(x -> ysv[x] + 1, V)..., r] = (Mk * R)[s, t]
                end
            elseif is_close
                for (l, (s, t)) in enumerate(pairs)
                    C[l, ntuple(x -> ysv[x] + 1, V)..., 1] =
                        (cellodd * inpar[l] == 1 ? -1 : 1) * Mk[t, s]
                end
            else
                for (l, (s, t)) in enumerate(pairs), (r, (v, w)) in enumerate(pairs)
                    t == w && (C[l, ntuple(x -> ysv[x] + 1, V)..., r] =
                        (cellodd * inpar[l] == 1 ? -1 : 1) * Mk[v, s])
                end
            end
        end
        all(iszero, C) && error("混合 K：单元 $k ($b,$j) 核全零")
        ts = split_core(CellCore(C, inpar, outpar))
        block = b === :τ ? j + 3 : (b === :+ ? Nτ + 3 + j : Nτ + 3 + Nt + (Nt - j + 1))
        for v in 1:2M
            gmps.data[(block - 1) * 2M + v] = ts[v]
        end
        F = F ⊕ f_k
    end
    GTEMPO.unset_svectors!(gmps)
    return gmps
end


"""
    _system_factor_tau — 虚时 K（U1 语义 + F 宇称路由；已验证）。

单元 1..N 沿 τ 顶点 2..N+1。折叠（j=1，第二轴 = 右真空 [0]，ρ₀ 在第二
成员接合）；闭合（j=N，第一轴 = 左真空 [0]，A[tR,sR] transpose-闭合）；
内部：第一轴 = F_j、第二轴 = F_{j−1}（F = 切位奇插入累积宇称），
|ys| = f_j 自动筛选。
"""

function _system_factor_tau(c::WeylImagTime, model; order::Int = 2,
                            coupling_norm::Real = 0.0,
                            rho0::Union{Nothing,AbstractMatrix} = nothing,
                            insertions::AbstractVector{WeylInsertion} = WeylInsertion[],
                            mode::Symbol = :gaussian, variant::Symbol = :dag,
                           oddlater_mode::Symbol = :none)
    N = c.Nτ
    M = c.bands
    dim = 2^M
    pairdim = dim * dim
    slotops = _slot_operators(M, order; mode = mode)
    lattice = weyl_lattice(c)
    h = c.dτ
    sc = order == 2 ? exp(h^2 * coupling_norm * M / 4) : 1.0   # Weyl 标量 = 总耦合 Σ_b α_b（2026-11-05 Z 因子化实证）

    isoddop(Mx) = any(!iszero(Mx[i, j]) &&
                      paritybits(i - 1) != paritybits(j - 1)
                      for i in 1:size(Mx, 1), j in 1:size(Mx, 2))

    insmap = Dict{Tuple{Symbol,Int},Vector{Matrix{ComplexF64}}}()
    for O in insertions
        push!(get!(insmap, (O.branch, O.boundary), Matrix{ComplexF64}[]),
              Matrix{ComplexF64}(O.matrix))
    end

    # ρ₀ = I；boundary-0 插入折叠（U1 cutop(0)·ρ₀）；F₀ = 其奇宇称
    R = Matrix{ComplexF64}(I, dim, dim)
    for O in get(insmap, (:τ, 0), Matrix{ComplexF64}[])
        R = O * R
    end
    F0 = count(isoddop, get(insmap, (:τ, 0), Matrix{ComplexF64}[])) % 2

    hfrac = order == 2 ? 0.5 : 1.0
    positions = Int[]
    tensors = TensorMap[]
    _pairpar(i) = paritybits((i - 1) % dim) ⊕ paritybits((i - 1) ÷ dim)
    F_before = F0                          # prefix₀ 计入首格的 f₁
    for j in 1:N
        S = Matrix{ComplexF64}(GTEMPO.fock_propagator(model, :τ, h * hfrac).data)
        Amats = Vector{Matrix{ComplexF64}}(undef, 4^M)
        for code in 0:(4^M)-1
            Mop = slotops[code+1]
            Amats[code+1] = order == 2 ? sc * (S * Mop * S) : S * Mop
        end
        # cutop(j)：该单元的插入矩阵左乘
        for O in get(insmap, (:τ, j), Matrix{ComplexF64}[])
            Amats = [O * m for m in Amats]
        end
        # 注意：不施 U1 的 oddlater 显式符号——分级收缩已供给全部奇算符
        # 交叉符号；显式再加会重复计入（退耦浴标定：:none 与 ED 精确一致，
        # :u1 反而带 δt 无关的 τ 线性偏移）。
        # f_j 保持纯净（仅切位 j 的奇插入宇称）；boundary-0 的 F₀ 已在
        # F_before 初值中，切位 j 的键扇区 = pairpar ⊕ (F_before ⊕ f_j) = pairpar ⊕ F_j。
        f_j = count(isoddop, get(insmap, (:τ, j), Matrix{ComplexF64}[])) % 2

        insize = (j == N) ? 1 : pairdim
        outsize = (j == 1) ? 1 : pairdim
        # F 路由（宇称轴，v15 已验证方向）：第一轴（朝晚端）= pairpar ⊕ F_j；
        # 第二轴（朝早端）= pairpar ⊕ F_{j−1}。边界：折叠右端 / 闭合左端 = 真空 [0]。
        inpar = (j == N) ? [0] : [_pairpar(i) ⊕ (F_before ⊕ f_j) for i in 1:insize]
        outpar = (j == 1) ? [0] : [_pairpar(i) ⊕ F_before for i in 1:outsize]

        V = 2M
        C = zeros(ComplexF64, insize, ntuple(Returns(2), V)..., outsize)
        for ci in 0:(2^V)-1
            ys = ntuple(k2 -> (ci >> (k2 - 1)) & 1, V)
            A = Amats[_code_from_ys(ys, 0, M)+1]
            # 2026-11-05 裁决：旧 τ 格间楔补偿 (−1)^{Σ_{j≥2}Nⱼ}（M≥2）已移除——
            # 与 IF 预翻同批引入、同源（未带置换宇称的 CSV 对比假象）。
            # U=0 因子化检验（z_factor_probe.jl）：移除后 K₂ = K₁∧K₁ 逐支撑
            # 精确成立；保留补偿反而破坏因子化（Z₂/Z₁² 0.992→0.854）。
            for so in 1:outsize, si in 1:insize
                sL = (si - 1) % dim + 1; tL = (si - 1) ÷ dim + 1
                sR = (so - 1) % dim + 1; tR = (so - 1) ÷ dim + 1
                v = zero(ComplexF64)
                if j == 1 && j == N
                    # 退化 N=1：两侧皆边界 —— 全 Fock 迹闭合 Tr[A·R]
                    #（真空闭合会丢失热迹；2026-11-04 calibrate ① N=1 暴露）
                    for s in 1:dim, s0 in 1:dim
                        v += A[s, s0] * R[s0, s]
                    end
                elseif j == 1
                    for s0 in 1:dim
                        v += A[sL, s0] * R[s0, tL]
                    end
                elseif j == N
                    v = A[tR, sR]
                else
                    v = (tL == tR) ? A[sL, sR] : zero(ComplexF64)
                end
                iszero(v) || (C[si, ntuple(x -> ys[x] + 1, V)..., so] = v)
            end
        end
        all(iszero, C) && error("虚时 K：单元 $j 核全零")
        append!(tensors, split_core(CellCore(C, inpar, outpar)))
        v = cell_vertex(c, j)
        for b in 1:M, cj in (false, true)
            push!(positions, lattice[GrassmannContourIndex(v, :τ, b, cj)])
        end
        F_before = F_before ⊕ f_j
    end

    gmps = vacuumstate(ComplexF64, lattice)
    return tile_cells!(gmps, positions, tensors)
end
