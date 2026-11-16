# bridge.jl：GrassmannMPS 后端适配层（对称性泛型化）
#
# 布局映射（plan/gtempo_grassmannmps_backend_plan.md §2，D1 裁决）：
#   构造路径（sysdynamics/influencefunctional）保持变量粒度（每变量一站、2 维
#   物理腿、bosonic 扇区），经 `pairize` 相邻同味站对楔积融合成 pair 布局 GMPS
#   （每胞 2 站：up 形/down 形 grassmannspace，4 维站）后进入 GrassmannMPS。
#
# 对称性用户手选（2026-11-13 D1 裁决，§3.4；2026-11-16 A1/A2 裁决：spec =
# TK 扇区类型本身、荷标签 = TK 扇区——Julia 类型即标签惯用法）。荷标签
# （构造层纯数据）类型随 spec：
#   Z2Irrep      :z2   各带宇并合并——旧 GTEMPO 语义；sector = mod 2
#   Z2Z2         :z2z2 每自旋宇并：奇带 = 各轨道 up、偶带 = down，偶数带
#   U1Irrep      :u1   各带总粒子数合并；sector = 原值
#   ChargeCharge :u1u1 每自旋粒子数：奇带 = up、偶带 = down，偶数带；U1GTEMPO2 同款
# 费米符号不在标签里——由 GM.GrassmannBackend 显式注入（@grassmann 宏）。
# 参考实现 = U1GTEMPO2 bridge.jl（ChargeCharge 通道机器精度验证过）；
# 其余三通道靠 GM abelian 电池 + 本包跨对称性一致电池 + ED 三层兜底（R9）。

# ---------------------------------------------------------------------------
# 对称性 spec（扇区类型注册面：本文件是唯一允许按对称性分支的位置之一，
# 与 GM/grassmannspaces.jl 的注册表对表）
# ---------------------------------------------------------------------------

const SYMMETRY_SPECS = IdDict{Symbol,Type}(
    :z2 => TK.Z2Irrep, :z2z2 => Z2Z2, :u1 => TK.U1Irrep, :u1u1 => ChargeCharge)

"Symbol → spec（扇区类型；未知 Symbol 显式拒绝，禁静默降级——D1 裁决）"
symmetry_spec(s::Symbol)::Type{<:Sector} = get(SYMMETRY_SPECS, s, nothing) === nothing ?
    throw(ArgumentError("unknown symmetry :$s (allowed: :z2, :z2z2, :u1, :u1u1)")) :
    SYMMETRY_SPECS[s]

# 荷标签代数（TK 扇区直用）：融合 = first(⊗)（GM 自身惯用法）、负 = dual、
# 减 = 加负。宇并形 mod 2 语义由扇区融合自然承载（+1 与 −1 同扇区不压重）。

@inline qadd(a::Sector, b::Sector) = first(⊗(a, b))
@inline qneg(a::Sector) = TK.dual(a)
@inline qsub(a::Sector, b::Sector) = qadd(a, TK.dual(b))

"荷排序键（A1 裁决：复现旧标签序——Int 升序 / 元组字典序。TK 扇区 isless
是 |q| 序，直接 sort 会重排键扇区布局破 golden，故所有排序走 sortkey）"
sortkey(s::TK.Z2Irrep) = Int(s.n)
sortkey(s::TK.U1Irrep) = Int(s.charge)
sortkey(s::Z2Z2) = (Int(s[1].n), Int(s[2].n))
sortkey(s::ChargeCharge) = (Int(s[1].charge), Int(s[2].charge))

"荷分辨形 bands 守卫（:z2z2/:u1u1）：多轨道完全支持——轨道交错约定
（奇带 = 各轨道 up、偶带 = down；bands = 2·轨道数），唯一结构要求 = bands
为偶数；语义要求（杂质 Hamiltonian 与 IF 确实携带该对称性，如 KanamoriIM）
由模型/调用方保证。单荷形 :z2/:u1 无要求。"
function bands_guard(::Type{G}, bands::Int) where {G<:Sector}
    if (G === Z2Z2 || G === ChargeCharge) && isodd(bands)
        throw(ArgumentError(
            "symmetry $G resolves per-spin flavor and requires an even number of bands " *
            "(orbital interleaving: odd bands = up, even bands = down), got bands = $bands"))
    end
    return bands
end

"变量荷：band b 上 η 生成元的荷（η̄ = qneg）
——宇并形算术 mod 2 后 η̄ ≡ η（同一扇区），带内 (a,ā) 区分在站形槽位，
不在荷上；荷分辨形 band→槽位映射照旧。"
var_charge(::Type{TK.Z2Irrep}, band::Int, bands::Int) = TK.Z2Irrep(1)
var_charge(::Type{TK.U1Irrep}, band::Int, bands::Int) = TK.U1Irrep(1)
# 轨道交错：奇带 = 各轨道 up（第一味）、偶带 = down（第二味）——多轨道
# （bands = 2·轨道数）与双带同式；旧「band == 1 ? up : down」写法在 bands ≥ 3
# 时把轨道 2 的 up 错标成 down（2026-11-16 修正，随 bands_guard 放宽）
var_charge(::Type{Z2Z2}, band::Int, bands::Int) =
    isodd(band) ? Z2Z2((TK.Z2Irrep(1), TK.Z2Irrep(0))) : Z2Z2((TK.Z2Irrep(0), TK.Z2Irrep(1)))
var_charge(::Type{ChargeCharge}, band::Int, bands::Int) =
    isodd(band) ? ChargeCharge(TK.U1Irrep(1), TK.U1Irrep(0)) : ChargeCharge(TK.U1Irrep(0), TK.U1Irrep(1))

"站形：band b 的物理站空间（单荷形全站同形；双荷形按带奇偶交替朝向——
:u1u1 = U1GTEMPO2 同款 up/down 交替站形；朝向序 = GM 规范形注册表
canonical_grassmann_spaces 的 [up, down] 约定）"
function station_space(::Type{G}, band::Int) where {G<:Sector}
    spaces = GM.canonical_grassmann_spaces(G)
    return length(spaces) == 1 ? only(spaces) : spaces[isodd(band) ? 1 : 2]
end

"变量荷 → 站朝向（_convert_site 用；单荷形无朝向 → 恒 :up 合法）"
station_flavor(::Type{TK.Z2Irrep}, q) = :up
station_flavor(::Type{TK.U1Irrep}, q) = :up
station_flavor(::Type{Z2Z2}, q) = iszero(q[2].n) ? :up : :down
station_flavor(::Type{ChargeCharge}, q) = iszero(q[2].charge) ? :up : :down

"键荷序列 → 分级空间（荷计数聚合；排序 = sortkey 复现旧标签升序保证扇区
插入序确定——A1）"
function qspace(::Type{G}, labels::AbstractVector) where {G<:Sector}
    TK.Vect[G](q => count(==(q), labels) for q in sort(unique(labels); by = sortkey))
end
"2 维物理腿空间：vacuum ⊕ generator q"
pspace(::Type{G}, q::Sector) where {G<:Sector} = qspace(G, [one(G), q])
"物理腿的非平凡扇区（2 维腿只有一个生成元扇区）"
occupied(p) = only(filter(!=(one(TK.sectortype(p))), collect(TK.sectors(p))))

"每胞变量荷序列（Chrono 排布每格 = [(band b: a, ā)] 逐带）"
function cell_charges(::Type{G}, bands::Int) where {G<:Sector}
    out = Sector[]
    for b in 1:bands
        qc = var_charge(G, b, bands)
        push!(out, qc, qneg(qc))
    end
    return out
end

"BCS 守卫（D7 裁决）：BCS 反常关联破坏粒子数守恒、宇称恒守恒 → 仅 :z2"
bcs_symmetry_guard(sym::Symbol) =
    sym === :z2 || throw(ArgumentError(
        "BCS bath supports only symmetry=:z2 (anomalous pairing breaks particle-number " *
        "conservation; only parity is conserved), got :$sym"))

# ---------------------------------------------------------------------------
# 链几何：变量 ↔ 站映射（Chrono 排布每格 2·bands 变量 = bands 站）
# ---------------------------------------------------------------------------

"链站空间序列：站 s 的 band = mod1(s, bands)"
function site_spaces(::Type{G}, nstations::Int, bands::Int) where {G<:Sector}
    [station_space(G, mod1(s, bands)) for s in 1:nstations]
end

"变量位（1..2·nstations）→ 站位"
@inline station_of(varpos::Int) = (varpos + 1) ÷ 2

# ---------------------------------------------------------------------------
# 截断（A3 裁决 2026-11-16：Truncation = MAK.TruncationStrategy，无状态
# findtruncated_svd 重载——top-D 全局预算 ∧ 相对尾容差；误差只走 GM
# truncate!/mult 返回的 truncerr 标准通道。荷排序键 = sortkey，跨荷全局预算）
# ---------------------------------------------------------------------------

struct Truncation <: MAK.TruncationStrategy
    D::Int                                # maxdim → D（对齐 GM TruncateDimCutoff.D / MAK truncrank(D) 惯例）
    rtol::Float64
    function Truncation(; D::Integer=128, rtol::Real=1e-12)
        D > 0 && isfinite(rtol) && 0 <= rtol < 1 ||
            throw(ArgumentError("invalid truncation"))
        new(D, rtol)
    end
end

"Truncation 的 svd_trunc 扩展（= 旧 LegacyTruncation 去 cutoff 去记账）：排序
键 (-|v|, sortkey, j) ≡ 旧 (-v, 标签, j)；保留判据 top-D ∧ 相对尾能量
tails[k+1] > rtol²·total（rtol 是丢弃尾巴的相对 L2 能量判据，与 GM
TruncateDimCutoff.ϵ 的逐值噪声地板无换算关系；rtol=0 与全零谱每键保 1 条）"
function MAK.findtruncated_svd(values::TK.SectorVector, s::Truncation)
    entries = sort!(collect((abs(v), q, j) for (q, vs) in pairs(values)
                            for (j, v) in enumerate(vs));
                    by = e -> (-e[1], sortkey(e[2]), e[3]))
    tails = zeros(length(entries) + 1)
    for i in length(entries):-1:1
        tails[i] = tails[i+1] + abs2(entries[i][1])
    end
    k = 0
    for (_, _, _) in entries
        k >= s.D && break
        k > 0 && tails[k+1] <= s.rtol^2 * tails[1] && break
        k += 1
    end
    keep = Dict{keytype(values),Int}()
    for (_, q, j) in entries[1:k]
        keep[q] = max(get(keep, q, 0), j)
    end
    return TK.SectorDict(q => (1:get(keep, q, 0)) for q in keys(values))
end

"结构完整性检查（旧 check_charge 语义）：非空、logscale 有限、左边界真空、
键匹配、物理站为规范同站形"
function _check_structure(x::GM.GMPS)
    isempty(x) && throw(ArgumentError("empty GMPS"))
    isfinite(GM.logscale(x)) || throw(ArgumentError("nonfinite GMPS logscale"))
    TK.space(x[1], 1) == oneunit(TK.space(x[1], 1)) || throw(ArgumentError("left boundary must be vacuum"))
    for i in 1:length(x)-1
        TK.domain(x[i])[1] == GM.space(x[i+1], 1) || throw(ArgumentError("bond mismatch between sites $i and $(i+1)"))
    end
    foreach(i -> GM.grassmann_basis(GM.physpace(x, i)), 1:length(x))
    true
end

# 截断/乘法退役（A3）：正则化压缩 = GM.truncate!、压缩乘法 = GM.mult——
# 误差走两者返回的 truncerr 标量（逐键相对截断误差之和），本层不再自造
# compress!/multiply 机器。

# ---------------------------------------------------------------------------
# pairize：变量粒度裸张量链 → pair 布局 GMPS（对称性泛型版）
# ---------------------------------------------------------------------------

"站形（flavor）：物理生成元荷 → 站朝向（单荷形恒 :up）"
_station_space(::Type{G}, q) where {G<:Sector} =
    station_space(G, station_flavor(G, q) === :up ? 1 : 2)

"变量粒度站张量（cod=(键出,物理 2 维), dom=(键入)）→ pair 布局 4 维站张量：
**键腿恒等**（荷扇区原样保留——荷分辨对称性的效率即在于此，键两端恒等
映射故必一致）；物理腿按味嵌入各自站形的槽——vac→unit 槽 1、η 变量 →
grassmann_basis.eta 标签、η̄ 变量 → grassmann_basis.etabar 标签（异扇区形 =
扇区分辨、自对偶形 = 同扇区分槽——U1GTEMPO2 原版只支持异扇区，自对偶槽位
是本泛化版的关键补丁），pair 槽留空（逐对 _fuse_adjacent 楔积时填充）。
`conj` = 变量的共轭旗（pairize 内奇位 = a/η、偶位 = ā/η̄）。"
function _convert_site(::Type{<:Sector}, t, conj::Bool)
    q = occupied(GM.space(t, 2))
    E = _station_space(TK.sectortype(t), q)
    basis = GM.grassmann_basis(E)
    tag = conj ? basis.etabar : basis.eta
    G = TK.sectortype(t)
    out = TK.zeros(TK.scalartype(t), TK.space(t, 1) ⊗ E ← TK.domain(t)[1])
    for (f1, f2) in TK.fusiontrees(t)
        if f1.uncoupled[2] == one(G)
            # 物理 (0,0) 块（vac）写进站 (0,0) 槽 1
            out[(f1.uncoupled[1], one(G), TK.dual(only(f2.uncoupled)))][:, 1, :] .= t[f1, f2][:, 1, :]
        else
            # 生成元块写进站形的目标味标签（sector + slot）
            out[(f1.uncoupled[1], tag.sector, TK.dual(only(f2.uncoupled)))][:, tag.slot, :] .= t[f1, f2][:, 1, :]
        end
    end
    out
end

"相邻两站（同味变量对）楔积融合：**收缩中间键** + wedge_vertex 站级楔积顶点；
机器符号由 GM.@grassmann 承担"
function _fuse_adjacent(a, b)
    GM.domain(a)[1] == GM.space(b, 1) ||
        throw(ArgumentError("adjacent-site fusion requires matching middle bond"))
    mv = GM.wedge_vertex(TK.scalartype(a), GM.space(a, 2))
    GM.@grassmann out[l p; r] := a[l sa; m] * b[m sb; r] * mv[p; sa sb]
    out
end

"""
    pairize(spec, tensors) -> GMPS

变量粒度裸张量链（偶数长）→ pair 布局 GMPS：逐对 (θ站, θ̄站) 嵌入后相邻融合
（中间键收缩 + 楔积）。链长减半、物理维 2×4。
"""
function pairize(::Type{G}, tensors::AbstractVector) where {G<:Sector}
    iseven(length(tensors)) || error("pairize requires an even variable count (θ/θ̄ paired): got $(length(tensors))")
    interp = [_convert_site(G, t, iseven(i)) for (i, t) in enumerate(tensors)]
    data = [_fuse_adjacent(interp[2j-1], interp[2j]) for j in 1:(length(tensors) ÷ 2)]
    GM.GMPS(data)
end

# ---------------------------------------------------------------------------
# 链构造 API（from_core/split_cell/monomial；全部停留在变量粒度层、输出裸
# 张量——pairize 后才进 GrassmannMPS。U1GTEMPO2 同名函数的对称性泛型版）
# ---------------------------------------------------------------------------

function from_core(::Type{G}, data, leftlabels, rightlabels, c) where {G<:Sector}
    t = TK.zeros(ComplexF64, qspace(G, leftlabels) ⊗ pspace(G, c) ← qspace(G, rightlabels))
    z = one(G)
    for q in unique(leftlabels), (si, s) in ((1, z), (2, c))
        r = qadd(q, s)
        r in rightlabels || continue
        li = findall(==(q), leftlabels); ri = findall(==(r), rightlabels)
        t[(q, s, TK.dual(r))] .= data[li, si:si, ri]
    end
    t
end

# Split only fixed-size unit cells, never the full temporal polynomial.
function split_cell(::Type{G}, data::AbstractArray{ComplexF64,3}, leftlabels, rightlabels, charges) where {G<:Sector}
    work = Array(data); labels = leftlabels
    z = one(G)
    T0 = typeof(TK.zeros(ComplexF64, qspace(G, [z]) ⊗ pspace(G, charges[1]) ← qspace(G, [z])))
    tensors = T0[]
    for c in charges[1:end-1]
        rest = size(work, 2) ÷ 2; mat = reshape(work, 2length(labels), :)
        rowqs = [qadd(q, s == z ? z : c) for s in (z, c) for q in labels]
        newqs = G[]; pieces = []
        for q in sort(unique(rowqs); by = sortkey)
            rows = findall(==(q), rowqs)
            Uc, S, V = MAK.svd_compact(mat[rows, :])
            # No numerical singular-value cutoff on factor construction.
            keep = count(>(0), diag(S)); keep == 0 && continue
            push!(pieces, (rows, Uc[:, 1:keep], S[1:keep, 1:keep] * V[1:keep, :]))
            append!(newqs, fill(q, keep))
        end
        isempty(newqs) && throw(ArgumentError("zero cell"))
        core = zeros(ComplexF64, length(labels), 2, length(newqs)); next = zeros(ComplexF64, length(newqs), rest, length(rightlabels)); offset = 0
        for (rows, Uc, R) in pieces
            k = size(Uc, 2); dst = offset+1:offset+k
            reshape(core, 2length(labels), :)[rows, dst] = Uc
            next[dst, :, :] = reshape(R, k, rest, length(rightlabels)); offset += k
        end
        push!(tensors, from_core(G, core, labels, newqs, c)); work = next; labels = newqs
    end
    push!(tensors, from_core(G, work, labels, rightlabels, charges[end])); tensors
end

"""
    monomial(spec, charges, positions=Int[]; coeff=1) -> GMPS

单项式链：`positions` 上的变量带各自生成元荷、其余真空；逆序号吸收进 coeff。
空 charges 抛错（零时刻轮廓走 solve 的零时路径）。
"""
function monomial(::Type{G}, charges::AbstractVector, positions::AbstractVector{<:Integer}=Int[]; coeff::Number=1) where {G<:Sector}
    isempty(charges) && throw(ArgumentError("use the zero-time path for empty contours"))
    all(i -> 1 <= i <= length(charges), positions) || throw(BoundsError(charges, positions))
    length(unique(positions)) == length(positions) || throw(ArgumentError("repeated Grassmann variable"))
    all(isfinite, (coeff,)) || throw(ArgumentError("nonfinite coefficient"))
    inv = sum(positions[i] > positions[j] for i in eachindex(positions) for j in i+1:length(positions); init=0)
    z = one(G)
    q = z
    data = map(eachindex(charges)) do i
        s = i in positions ? charges[i] : z; r = qadd(q, s)
        t = TK.zeros(ComplexF64, qspace(G, [q]) ⊗ pspace(G, charges[i]) ← qspace(G, [r]))
        t[(q, s, TK.dual(r))] .= i == 1 ? coeff * (-1)^inv : 1
        q = r; t
    end
    pairize(G, data)
end

# ---------------------------------------------------------------------------
# 系数提取（站基矢指标 0..3 经 grassmann_basis 定位，不硬编码索引；
# 占据映射：变量序 (2k−1,2k) → 站 k 的 η/η̄ 槽）
# ---------------------------------------------------------------------------

function _occupation_index(p, o1::Int, o2::Int)
    b = GM.grassmann_basis(p)
    tag = o1 == 0 ? (o2 == 0 ? b.vac : b.etabar) : (o2 == 0 ? b.eta : b.pair)
    basis = [(s, j) for s in collect(TK.sectors(p)) for j in 1:TK.dim(p, s)]
    findfirst(isequal((tag.sector, tag.slot)), basis) - 1
end

# 系数抽查直接消费 GM.coefficient_scaled/coefficient（站级 digit 约定：每站一个
# 0..3 基指标，vac0/pair1/eta2/etabar3）；变量级 (o1,o2)→digit 的映射助手留作
# 测量层（Phase 5）私用。
function _var_coefficient_scaled(x::GM.GMPS, occupations)
    L = length(x)
    length(occupations) == 2L || throw(DimensionMismatch("one binary occupation per variable (2 per station)"))
    all(s -> s in (0, 1), occupations) || throw(ArgumentError("binary occupations required"))
    occ = [_occupation_index(GM.physpace(x, k), occupations[2k-1], occupations[2k]) for k in 1:L]
    z = GM.coefficient_scaled(x, occ)
    z
end

# ---------------------------------------------------------------------------
# 结构查询（旧导出 API 保名）
# ---------------------------------------------------------------------------

charge_dimensions(x::GM.GMPS) = [Dict(q => TK.dim(s, q) for q in TK.sectors(s)) for s in GM.bondspaces(x)]
check_charge(x::GM.GMPS) = _check_structure(x)

# ---------------------------------------------------------------------------
# 收缩引擎接线（U1GTEMPO2 §2.5：integrate_scaled→Zvalue、EnvironmentCache→
# environments、中切口分母条件数在本层）
# ---------------------------------------------------------------------------

function integrate_scaled(states::GM.GMPS...)
    # 上限 = GM 引擎的均匀循环链数（K ≤ 5：K + 4 IF ⟹ norb=2 荷分辨形满员；
    # 更高轨道数走合并宇并/总荷形 :z2/:u1——单 IF 链不限轨道数）
    1 <= length(states) <= 5 || throw(ArgumentError("integrate supports one to five independent factors (GM engine chain bound K ≤ 5)"))
    z = GM.Zvalue(collect(states))
    z
end
integrate(states::GM.GMPS...) = value(integrate_scaled(states...))

"""
快照 + GrassmannMPS ExpectationCache 包装 + 中切口分母条件数（旧 middle-cut
cross 语义：|boundary_scalar(L_mid·R_mid)|/norms，站粒度切割取链半长）。
"""
struct EnvironmentCache
    cache::GM.ExpectationCache
    denominator_condition::Float64
end
function EnvironmentCache(states::Tuple)
    1 <= length(states) <= 5 || throw(ArgumentError("one to five factors required (GM engine chain bound)"))
    cache = GM.environments(states...)
    mid = cache.L ÷ 2                      # 站粒度半链切割（旧对粒度 mid 同语义）
    lmid = cache.leftenv[mid+1]; rmid = cache.rightenv[mid+1]
    cond = abs(GM._boundary_pair(lmid, rmid)) / max(norm(lmid) * norm(rmid), floatmin(Float64))
    EnvironmentCache(cache, min(cond, cache.denominator_condition))
end
integrate_scaled(c::EnvironmentCache) = integrate_scaled(c.cache.states...)
integrate(c::EnvironmentCache) = value(integrate_scaled(c))
