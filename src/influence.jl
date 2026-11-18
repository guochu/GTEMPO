# influence.jl：浴影响泛函 IF（plan §Phase 3，D3 首批 = PartialIF + ExactTTI）
#
# 移植自 U1GTEMPO2 influence.jl（charge 通道机器精度、与 secondorderif 端到端
# ≤1.4e-11 验证），泛化点（D1）：
#   - spin::Symbol(:up/:down) → band::Int（band b 的生成元荷 = var_charge(spec,b)；
#     双荷形奇带=up/偶带=down（多轨道轨道交错），单荷形全带同荷、带=站位）；
#   - kernel blocks 为 canonical Δ[(a,b)][i,j]（IF = exp(Σ ā_i Δ_ij a_j)），
#     由 hybridization() 从 QuAPI corr 提取（含实时跨支/混合 τ↔rt 符号定标）。
#
# TDVPIF/XTRGIF/TTIIF 暂停（D3 裁决：第二批 TTIIF/XTRGIF 走 GM.SchurGMPS/
# timeevo/expmult! 另行移植）。

# ---------------------------------------------------------------------------
# kernel 定义与构造已迁 correlationfunction.jl（plan/if_expmult W2：
# HybridizationKernel{T} + ContourBlockIndex + correlationfunction 双构造）；
# 本文件只保留 IF 装配消费面。
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# ⚠️ deprecated（W3 语义换代冻结档案，交叉验证期后退役）：逐行秩一因子积
# PartialIF——现行装配走 influencefunctional/partialif.jl 的 fullfactor +
# expmult!×PartialExp（同一算符的整形式表示）；本段保留供每因子 dense 对拍。
# ---------------------------------------------------------------------------

"""
    partial_row(spec, charges, anchor, columns, values) -> GMPS

第 anchor 行的偏 IF 因子 exp(ā_anchor Σ_j Δ_j a_j) 的链表示
（U1GTEMPO2 partial_row 的对称性泛型版）。
"""
function partial_row(::Type{G}, charges, anchor, columns, values) where {G<:Sector}
    qbar = charges[anchor]; q = qneg(qbar); z = one(G)
    ws = Dict(zip(columns, values))
    data = map(eachindex(charges)) do i
        ls = i == 1 ? [z] : [z, i <= anchor ? q : qbar]
        rs = i == length(charges) ? [z] : [z, i < anchor ? q : qbar]
        t = TK.zeros(ComplexF64, qspace(G, ls) ⊗ pspace(G, charges[i]) ← qspace(G, rs))
        t[(z, z, z)] .= 1
        if i == anchor
            q in ls && (t[(q, qbar, z)] .= 1)
            qbar in rs && (t[(z, qbar, TK.dual(qbar))] .= 1)
        else
            pending = i < anchor ? q : qbar
            pending in ls && pending in rs &&
                (t[(pending, z, TK.dual(pending))] .= 1)
            c = get(ws, i, 0im)
            if !iszero(c)
                if i < anchor
                    t[(z, q, TK.dual(q))] .= -c
                else
                    t[(qbar, q, z)] .= c
                end
            end
        end
        t
    end
    pairize(G, data)
end

"""
    partial_if(lattice, kernel; band=1, trunc=Truncation(), truncerrs) -> GMPS

PartialIF：IF = Π_row exp(ā_row Σ_j Δ_row,j a_j) 的逐行因子积（精确）。
`band` 选择 IF 所在带（生成元荷/站位随带）。行/列序 = 链存储格序
（_storage_cells），values 与 cols 逐位对齐。`truncerrs` = 逐行乘法压缩的
GM truncerr 记账（标量通道）。
"""
function partial_if(lattice::AbstractGrassmannLattice, kernel::HybridizationKernel;
                    band::Int=1, trunc::Truncation=Truncation(), truncerrs::Vector{Float64}=Float64[])
    spec = TK.sectortype(lattice)
    out = vacuumstate(lattice)
    charges = _site_charges(spec, lattice)
    cells = _storage_cells(lattice)              # 存储格序
    cols = [_cell_var_index(lattice, cell, band, false) for cell in cells]   # a 变量位
    ncells = length(cells)
    for (m, (b, j)) in enumerate(cells)
        # 行 = 格 (b,j)（branch 局部下标 j）；列值按存储格序逐位取核
        values = [kernel.blocks[(ContourBlockIndex(b, band), ContourBlockIndex(cell[1], band))][j, cell[2]] for cell in cells]
        all(iszero, values) && continue
        anchor = _cell_var_index(lattice, (b, j), band, true)   # ā 变量位
        out, err = GM.mult(out, partial_row(spec, charges, anchor, cols, values); trunc)
        push!(truncerrs, err)
    end
    out
end

"IF 变量荷序列（全链 = nblocks 格 × 每格 2·bands 变量，存储序）"
function _site_charges(::Type{G}, lattice) where {G<:Sector}
    per_cell = cell_charges(G, lattice.bands)
    return repeat(per_cell, nblocks(lattice))
end

"单元 (b,j) 的带 band 变量位：conj=false → a（奇位）、true → ā（偶位）"
function _cell_var_index(lattice, cell::Tuple{Symbol,Int}, band::Int, conj::Bool)
    return block_offset(lattice, cell_block(lattice, cell)) + 2 * (band - 1) + Int(conj)
end


# 矩阵铅笔拟合 fit_exponentials 已迁 influencefunctional/exacttti.jl（W4 拟合层）。

# ---------------------------------------------------------------------------
# 外代数记忆元胞（U1GTEMPO2 memory_cell/embed_cells/contact_factor 原样）
# ---------------------------------------------------------------------------

function wedge_sign(a, b)
    a & b != 0 && return 0
    isodd(sum(count_ones(b & ((1 << i) - 1)) for i in 0:63 if a & (UInt64(1) << i) != 0; init=0)) ? -1 : 1
end
function poly_mul(a, b)
    out = Dict{Int,ComplexF64}()
    for (i, x) in a, (j, y) in b
        s = wedge_sign(i, j); s == 0 && continue
        out[i|j] = get(out, i|j, 0im) + s * x * y
    end
    filter!(p -> !iszero(last(p)), out); out
end
linear_poly(v, offset=0) = Dict(1 << (i - 1 + offset) => ComplexF64(x) for (i, x) in enumerate(v) if !iszero(x))

"Exterior-algebra memory transfer. Openings are anchored at lag one, so q=0 is valid."
function memory_cell(roots, opens, closes)
    n = length(first(opens)); m = length(roots)
    opening = Dict(0 => 1.0 + 0im); subs = Dict{Int,ComplexF64}[]
    for a in 1:m
        xi = Dict(1 << (n + a - 1) => 1.0 + 0im)
        opening = poly_mul(opening, merge(Dict(0 => 1.0 + 0im), poly_mul(linear_poly(opens[a]), xi)))
        push!(subs, merge(linear_poly(closes[a]), Dict(1 << (n + a - 1) => ComplexF64(roots[a]))))
    end
    data = zeros(ComplexF64, 1 << m, 1 << n, 1 << m)
    for input in 0:(1 << m)-1
        p = Dict(0 => 1.0 + 0im)
        for a in 1:m; input & (1 << (a - 1)) == 0 || (p = poly_mul(p, subs[a])); end
        p = poly_mul(p, opening)
        for (mask, value) in p; data[input+1, (mask & ((1 << n) - 1)) + 1, (mask >> n) + 1] += value; end
    end
    data
end

function embed_cells(spec, lattice, localdata, memorycharges, cellvars; band=1)
    charges = _site_charges(spec, lattice)
    L = 2 * nblocks(lattice) * lattice.bands
    data = Vector{Any}(undef, L)
    labels = [reduce(qadd, (memorycharges[a] for a in eachindex(memorycharges) if mask & (1 << (a - 1)) != 0);
                      init=one(spec)) for mask in 0:(1 << length(memorycharges))-1]
    for (k, vars) in enumerate(cellvars)
        left = k == 1 ? [1] : collect(eachindex(labels)); right = k == length(cellvars) ? [1] : collect(eachindex(labels))
        ts = split_cell(spec, localdata[left, :, right], labels[left], labels[right], charges[vars])
        for (v, t) in zip(vars, ts); data[v] = t; end
    end
    # Spectator sites propagate the current virtual space with vacuum coefficient.
    active = Set(v for vars in cellvars for v in vars)
    for i in 1:L
        i in active && continue
        v = i == 1 ? qspace(spec, [one(spec)]) : TK.domain(data[i-1])[1]
        t = TK.zeros(ComplexF64, v ⊗ pspace(spec, charges[i]) ← v)
        G = TK.sectortype(t)
        for q in TK.sectors(v)
            b = t[(q, one(G), TK.dual(q))]
            reshape(b, size(b, 1), size(b, 3)) .= I(size(b, 1))
        end
        data[i] = t
    end
    out = pairize(spec, data); check_charge(out); out
end

# validate_tti 已删（2026-11-16 kernel.blocks 物化退役）：blocks 是 QuAPI corr
# 的惰性 branch 视图，Toeplitz 平移不变是 CorrelationMatrix 存储的结构性保证
# （统一约定/Keldysh 恒等式由 QuAPI 构造器校验）——旧校验退化为恒真。

"A whole segment's equal-cell exponential is a product of fixed-size local cells."
function contact_factor(spec, lattice, band, bs, D)
    n = bs[1] === :τ ? lattice.kτ : (lattice isa AbstractRealGrassmannLattice ? lattice.k : lattice.Nt)
    nv = 2length(bs)
    poly = Dict(0 => 1.0 + 0im)
    for i in 1:length(bs), j in 1:length(bs)
        bar = 2i - 1; ann = 2j - 2
        mask = (1 << bar) | (1 << ann)
        value = (bar < ann ? 1 : -1) * D[i, j]
        iszero(value) || (poly = poly_mul(poly, Dict(0 => 1.0 + 0im, mask => ComplexF64(value))))
    end
    core = zeros(ComplexF64, 1, 1 << nv, 1)
    for (mask, value) in poly; core[1, mask + 1, 1] = value; end
    if bs == (:+, :-)
        # 实时 ± 合并接触：虚拟布局构建 + reorder 到 Chrono 站序（同 real_factors）
        vars = [_virtual_merged_vars(lattice, band, m) for m in 1:n]
        fac = embed_cells(spec, lattice, core, Any[], vars; band)
        GM.reorder!(fac, _virtual_to_chrono_perm(lattice, n); trunc=GM.truncrelerr(1e-14))
        return fac
    end
    js = bs[1] === :τ ? (1:n) : (n:-1:1)
    vars = [[_cell_var_index(lattice, (b, j), band, bar) for b in bs for bar in (false, true)] for j in js]
    embed_cells(spec, lattice, core, Any[], vars; band)
end

# ---------------------------------------------------------------------------
# ⚠️ deprecated（W4 换代冻结档案）：ExactTTI 旧三通道装配（memory_cell/
# contact_factor/embed_cells 发射；tti_* 新构造见 influencefunctional/
# exacttti.jl）。混合 cross 通道的 rank-one 系数向量落在 a 型腿（cross_factors
# rows/cols 同 conj=false）= 潜在味错位（无测试覆盖、新路由由互换对拍裁决）。
#
# ExactTTI：实/虚/交叉三族指数因子（U1GTEMPO2 real/imag/cross_factors 同型）
#
# 布局适配：memory_cell/contact 的合并胞转移按 U1GTEMPO2 存储 [± 相邻、j 降序]
# 设计。Chrono 实时存储 = [(+,1)…(+,n)][(−,n)…(−,1)] ⟹ real 因子在**虚拟
# ± 相邻布局**上构建（embed_cells 原样），再 `GM.reorder!` 置换到 Chrono 站序
# （站粒度置换、费米换序号上游已修）。Chrono 虚时存储 [格 Nτ…1] 与 compact
# 热通道链（降序物理时间、initial 在 τ=β）**天然同序** ⟹ 不 reorder 直嵌。
# ---------------------------------------------------------------------------

"虚拟合并胞 m（块 2m−1 = (+,n−m+1)、块 2m = (−,n−m+1)）的 band 带变量位
（顺序 (a+,ā+,a−,ā−) = memory_cell opens/closes 腿序）"
function _virtual_merged_vars(lattice, band, m)
    B = lattice.bands
    plus0 = (2m - 2) * 2B + 2(band - 1) + 1
    minus0 = (2m - 1) * 2B + 2(band - 1) + 1
    return [plus0, plus0 + 1, minus0, minus0 + 1]
end

"虚拟 ± 相邻布局 → Chrono 存储（τ 头原位、实时部分 [+ 升序][− 降序]）的
站置换：GM.reorder! 语义 = 新链第 k 站 = 原链第 perm[k] 站"
function _virtual_to_chrono_perm(lattice, n)
    B = lattice.bands
    off = lattice isa AbstractMixedGrassmannLattice ? lattice.Nτ : 0
    blockperm(k) = k <= off ? k :
                   k <= off + n ? 2(n - (k - off)) + 1 :
                   2(k - off - n)
    return [(blockperm(k) - 1) * B + s for k in 1:(off + 2n) for s in 1:B]
end

function real_factors(emit, spec, lattice, k::HybridizationKernel, band, alg, diagnostics)
    lattice isa Union{AbstractRealGrassmannLattice,AbstractMixedGrassmannLattice} || return nothing
    n = lattice.k
    blocks(x, y) = k.blocks[(ContourBlockIndex(x, band), ContourBlockIndex(y, band))]
    residuals = Float64[]; nfactors = 0
    g = ComplexF64[blocks(:+, :+)[i + 1, 1] for i in 1:n-1]
    l = ComplexF64[blocks(:+, :-)[i + 1, 1] for i in 1:n-1]
    # Keldysh tails audit（structured Lag 块近似：直读 kernel 的 Toeplitz 性）
    locations = [(max(d, 0) + 1, max(-d, 0) + 1) for d in -(n - 1):n-1]
    error2 = 0.0; total2 = 0.0
    for (i, j) in locations
        i == j && continue
        v = i > j ? (g[i-j], l[i-j], -g[i-j], -l[i-j]) : (-conj(l[j-i]), conj(l[j-i]), -conj(g[j-i]), conj(g[j-i]))
        weight = n - abs(i - j)
        for (pair, x) in zip(((:+, :+), (:+, :-), (:-, :+), (:-, :-)), v)
            error2 += weight * abs2(x - blocks(pair...)[i, j]); total2 += weight * abs2(blocks(pair...)[i, j])
        end
    end
    err = sqrt(error2) / max(sqrt(total2), eps())
    err <= max(alg.rtol, 1e-12) || throw(ArgumentError("ExactTTI real Keldysh relations fail: $err"))
    diagnostics[:keldysh_residual] = err
    q = var_charge(spec, band, lattice.bands)
    vars = [_virtual_merged_vars(lattice, band, m) for m in 1:n]   # 虚拟 ± 相邻布局
    perm = _virtual_to_chrono_perm(lattice, n)
    for (component, y) in ((:greater, g), (:lesser, l))
        # rightanchors 默认 true：超临界极点精确折叠为 :right，由下方守卫拒绝；
        # 传 false 会把增长极点静默投影到单位圆，架空该守卫。
        f = fit_exponentials(y; rtol=alg.rtol, maxterms=alg.maxterms)
        push!(residuals, f.relative_error)
        expansions = get!(diagnostics, :real_expansions, Dict{Symbol,Any}())
        expansions[component] = (; f..., fit_interval=(1, n - 1), convention=:integrated_lag_one)
        for (w, r, a) in zip(f.weights, f.nodes, f.anchors)
            a == :left || throw(ArgumentError("growing real pole unsupported; use PartialIF"))
            z = sqrt(w); zc = sqrt(conj(w))
            opens = [ComplexF64[0, z, 0, -z], ComplexF64[zc, 0, -zc, 0]]
            closes = component === :greater ? [ComplexF64[z, 0, 0, 0], ComplexF64[0, 0, 0, zc]] :
                     [ComplexF64[0, 0, z, 0], ComplexF64[0, zc, 0, 0]]
            data = memory_cell([r, conj(r)], opens, closes)
            fac = embed_cells(spec, lattice, data, [qneg(q), q], vars; band)
            GM.reorder!(fac, perm; trunc=GM.truncrelerr(1e-14))
            emit(fac); nfactors += 1
        end
    end
    diagnostics[:real_fit_residuals] = residuals; diagnostics[:real_pole_factors] = nfactors
    D = [blocks(a, b)[1, 1] for a in (:+, :-), b in (:+, :-)]
    all(iszero, D) || emit(contact_factor(spec, lattice, band, (:+, :-), D))
    nothing
end

function imag_factors(emit, spec, lattice, k::HybridizationKernel, band, alg, diagnostics)
    lattice isa Union{AbstractImagGrassmannLattice,AbstractMixedGrassmannLattice} || return nothing
    z = one(spec)
    n = lattice.kτ
    M = k.blocks[(ContourBlockIndex(:τ, band), ContourBlockIndex(:τ, band))]   # Toeplitz 平移不变由 corr 存储结构性保证（QuAPI 构造器校验）
    lower = ComplexF64[M[j+1, 1] for j in 1:n-1]; upper = ComplexF64[M[1, j+1] for j in 1:n-1]
    thermal = norm(upper + reverse(lower)) / max(norm(lower), eps())
    thermal <= max(alg.rtol, 1e-12) || throw(ArgumentError("imaginary kernel violates thermal pairing"))
    f = fit_exponentials(lower; rtol=alg.rtol, maxterms=alg.maxterms)
    diagnostics[:imag_expansion] = (; f..., fit_interval=(1, n - 1), convention=:integrated_lag_one)
    shift = 0im
    qcharge = var_charge(spec, band, lattice.bands)
    vars = [[_cell_var_index(lattice, (:τ, j), band, bar) for bar in (false, true)] for j in n:-1:1]
    for (w, q, a) in zip(f.weights, f.nodes, f.anchors)
        Z = 1 + q^n; abs(Z) >= 1e-8 || throw(ArgumentError("nearly singular thermal pole"))
        amp = sqrt(w * Z); d = ComplexF64[1 0; 0 q]; F = ComplexF64[0 1; 0 0]
        bar = (d, amp * transpose(F)); ann = (Matrix{ComplexF64}(I, 2, 2), amp * F)
        core = zeros(ComplexF64, 4, 4, 4)
        for x in 0:1, y in 0:1
            mat = a === :left ? bar[x+1] * ann[y+1] : (-1)^(x * y) * bar[y+1] * ann[x+1]
            core[:, 1+x+2y, :] .= kron(Matrix{ComplexF64}(I, 2, 2), mat)
        end
        labels = a === :left ? [z, qneg(qcharge), qcharge, z] :
                 [z, qcharge, qneg(qcharge), z]
        ts = nothing; initial = ComplexF64[1, 0, 0, 1]; terminal = initial / Z
        for j in 1:n
            ci = core; ls = labels; rs = labels
            if j == 1; ci = reshape(transpose(initial) * reshape(ci, 4, :), 1, 4, 4); ls = [z]; end
            if j == n; ci = reshape(reshape(ci, :, 4) * terminal, size(ci, 1), 4, 1); rs = [z]; end
            part = split_cell(spec, ci, ls, rs, [qcharge, qneg(qcharge)])
            ts === nothing ? (ts = part) : append!(ts, part)
        end
        # Chrono 虚时存储 [格 Nτ…1]（块 b = 格 n−b+1）与 compact 链（降序物理
        # 时间、initial 在 τ=β = 链首）**天然同序**——不再 reorder，直嵌：
        # single 站 m（= 格 n−m+1）→ 块 m → 全链站 (m−1)·bands + band
        single = pairize(spec, ts)
        emit(embed_single(spec, lattice, single, [(m - 1) * lattice.bands + band for m in 1:n]))
        # Channel contact is w*q^(n-1) for left and its negative for right.
        shift += (a === :left ? 1 : -1) * w * q^(n - 1)
    end
    contact = M[1, 1] + shift
    iszero(contact) || emit(contact_factor(spec, lattice, band, (:τ,), reshape([contact], 1, 1)))
    diagnostics[:imag_fit_residual] = f.relative_error; diagnostics[:thermal_residual] = thermal
    nothing
end

"把紧致单带链（pair 布局）嵌回全链：station_positions[m] = single 第 m 站在全链
的站位；缺席站位填键一致的真空穿透站；logscale 原样携带"
function embed_single(spec, lattice, single, station_positions)
    issorted(station_positions) || throw(ArgumentError("embedding must preserve order"))
    spaces = site_spaces(lattice)
    L = length(lattice) ÷ 2
    data = Vector{typeof(single[1])}(undef, L); filled = Set{Int}()
    for (m, p) in enumerate(station_positions)
        TK.space(single[m], 2) == spaces[p] || throw(TK.SpaceMismatch("embedded station $m does not match site $p"))
        data[p] = single[m]; push!(filled, p)
    end
    for i in 1:L
        i in filled && continue
        v = i == 1 ? TK.oneunit(spaces[1]) : TK.domain(data[i-1])[1]
        t = TK.zeros(ComplexF64, v ⊗ spaces[i] ← v)
        G = TK.sectortype(t)
        for q in TK.sectors(v)
            b = t[(q, one(G), TK.dual(q))]
            for kk in 1:size(b, 1); b[kk, 1, kk] = 1; end   # 仅 vac 槽（行 1）
        end
        data[i] = t
    end
    out = GMPS(data); GM.setlogscale!(out, GM.logscale(single)); out
end

"General interleaved rank-one exponential 1 + A*B, four automaton states."
function rankone_factor(spec, charges, a, b)
    z = one(spec)
    qa = only(unique(charges[findall(!iszero, a)]))
    qb = qneg(qa)
    all(i -> charges[i] == qb, findall(!iszero, b)) || throw(ArgumentError("rank-one charges differ"))
    labels = [z, z, qa, qb]
    data = map(eachindex(charges)) do i
        C0 = Matrix{ComplexF64}(I, 4, 4); C1 = zeros(ComplexF64, 4, 4)
        C1[1, 3] = a[i]; C1[1, 4] = b[i]; C1[3, 2] = b[i]; C1[4, 2] = -a[i]
        left = i == 1 ? reshape(ComplexF64[1, 0, 0, 0], 1, 4) : Matrix{ComplexF64}(I, 4, 4)
        right = i == length(charges) ? reshape(ComplexF64[1, 1, 0, 0], 4, 1) : Matrix{ComplexF64}(I, 4, 4)
        core = cat(reshape(left * C0 * right, size(left, 1), 1, size(right, 2)),
                   reshape(left * C1 * right, size(left, 1), 1, size(right, 2)); dims=2)
        from_core(spec, core, i == 1 ? [z] : labels, i == length(charges) ? [z] : labels, charges[i])
    end
    pairize(spec, data)
end

function cross_factors(emit, spec, lattice, k::HybridizationKernel, band, alg, diagnostics)
    lattice isa AbstractMixedGrassmannLattice && lattice.Nt > 0 || return nothing
    ranks = Int[]; errors = Float64[]
    for (rowbs, colbs) in (((:τ,), (:+, :-)), ((:+, :-), (:τ,)))
        rows = [_cell_var_index(lattice, (b, j), band, false) for b in rowbs for j in 1:lattice.Nτ]
        cols = [_cell_var_index(lattice, (b, j), band, false) for b in colbs for j in 1:lattice.Nt]
        M = vcat([hcat([k.blocks[(ContourBlockIndex(a, band), ContourBlockIndex(b, band))] for b in colbs]...) for a in rowbs]...)
        U, S, V = MAK.svd_compact(M); sv = diag(S); tot = sum(abs2, sv); r = length(sv)
        while r > 0 && sum(abs2, sv[r:end]) <= alg.cross_rtol^2 * tot; r -= 1; end
        push!(ranks, r); push!(errors, iszero(tot) ? 0.0 : sqrt(sum(abs2, sv[r+1:end]) / tot))
        for j in 1:r
            a = zeros(ComplexF64, length(_site_charges(spec, lattice))); b = zero(a)
            a[rows] .= sqrt(sv[j]) * U[:, j]; b[cols] .= sqrt(sv[j]) * V[j, :]
            emit(rankone_factor(spec, _site_charges(spec, lattice), a, b))
        end
    end
    diagnostics[:cross_ranks] = ranks; diagnostics[:cross_residuals] = errors; nothing
end

# IF 算法类型/InfluenceFunctional/influence_functional 入口已迁
# influencefunctional/influencefunctional.jl（W1 拆分 + W3/W4 后端统一）。

# ---------------------------------------------------------------------------
# 公开皮肤：hybriddynamics（corr → kernel → IF 链；W2 裁决：band/orbital
# 信息在 kernel 中——签名收窄为 (lattice, kernel, alg)，多带 kernel 返回
# 逐带 IF 序（kernel.bands 序），单带 kernel 返回单 IF）
# ---------------------------------------------------------------------------

"""
    hybriddynamics(lattice, kernel, alg=PartialIF(); trunc) -> GMPS | Tuple{GMPS...}
    hybriddynamics(lattice, kernel; alg=PartialIF(), trunc) = hybriddynamics(lattice, kernel, alg; trunc)

逐带 IF 链（influence_functional 便捷入口）：band 覆盖 = `kernel.bands`。
单带 kernel 返回单 GMPS（与 `Is...` 变参消费兼容）；多带 kernel 返回
`Tuple`（逐带序）。多带 IF 请配合 `correlationfunction(...)` 全带核一次装配
（免逐带重建 corr）。
"""
function hybriddynamics(lattice::AbstractGrassmannLattice, kernel::HybridizationKernel,
                        alg::InfluenceFunctionalAlgorithm;
                        trunc::Truncation=Truncation())
    all(b -> 1 <= b <= lattice.bands, kernel.bands) || throw(BoundsError(1:lattice.bands, kernel.bands))
    if length(kernel.bands) == 1
        return influence_functional(lattice, kernel; alg, band=kernel.bands[1], trunc).state
    end
    return ntuple(b -> influence_functional(lattice, kernel; alg, band=b, trunc).state,
                  length(kernel.bands))
end

function hybriddynamics(lattice::AbstractGrassmannLattice, kernel::HybridizationKernel;
                        alg::InfluenceFunctionalAlgorithm=PartialIF(),
                        trunc::Truncation=Truncation())
    return hybriddynamics(lattice, kernel, alg; trunc)
end

"integrate(lattice, x0, x1...) —— Berezin 迹（链 = (K, Is...)，GM 引擎）"
integrate(lattice::AbstractGrassmannLattice, x0::GM.GMPS, x1::GM.GMPS...) = integrate(x0, x1...)
integrate_scaled(lattice::AbstractGrassmannLattice, x0::GM.GMPS, x1::GM.GMPS...) = integrate_scaled(x0, x1...)

# ---------------------------------------------------------------------------
# BCS bath：反常通道 partial IF（D7 裁决：仅 :z2——反常关联破坏粒子数守恒、
# 宇称恒守恒；旧栈 bcsinfluencefunctional/partialif.jl 的 GM 重写）
# ---------------------------------------------------------------------------

"BCS 指数的行/列腿（旧 _bcscolumns 同约定）：行 = (band1,ā)/(band2,ā)/(band2,a)，
列 = (band1,a)/(band2,a)/(band2,ā)，各分支全节点；行同节点列由装配处剔除
（同变量自配对 ξ² = 0 且系数位不可重复）"
function _bcs_legs(lattice::AbstractGrassmannLattice, band1::Int, band2::Int)
    rows = GrassmannContourIndex[]
    cols = GrassmannContourIndex[]
    for b in branches(lattice)
        for j in 1:(b === :τ ? lattice.kτ : lattice.kt)
            push!(cols, GrassmannContourIndex(Int(j), b, band1, false))
            push!(cols, GrassmannContourIndex(Int(j), b, band2, false))
            push!(cols, GrassmannContourIndex(Int(j), b, band2, true))
        end
    end
    for b in branches(lattice)
        for j in 1:(b === :τ ? lattice.kτ : lattice.kt)
            push!(rows, GrassmannContourIndex(Int(j), b, band1, true))
            push!(rows, GrassmannContourIndex(Int(j), b, band2, true))
            push!(rows, GrassmannContourIndex(Int(j), b, band2, false))
        end
    end
    return rows, cols
end

"""
    bcs_partial_if(lattice, kernel; trunc=Truncation(), truncerrs) -> GMPS

BCS 反常 partial IF = Π_row exp(行腿 × Σ 列腿系数) 的逐行因子积（精确）。
系数直读 kernel 块（有向指数系数语义、无 conj 翻转；Keldysh 相对号已在
`correlationfunction` 构造期折进块值——实时跨支块）。行/列腿覆盖指数四类
单项式 (ā↑,a↑)/(ā↓,a↓)/(ā↑,ā↓)/(a↓,a↑)——腿的 (band, conj) 组合按
`_bcs_leg_conjs` 味表选择对应块键。**仅 :z2**（对称性守卫 D7）。
"""
function bcs_partial_if(lattice::AbstractGrassmannLattice, kernel::HybridizationKernel;
                        trunc::Truncation=Truncation(),
                        truncerrs::Vector{Float64}=Float64[])
    spec = TK.sectortype(lattice)
    spec === TK.Z2Irrep || throw(ArgumentError(
        "BCS IF supports only symmetry=$Z2Symmetry (anomalous pairing breaks particle-number " *
        "conservation; only parity is conserved), got $(symmetry(lattice))"))
    length(kernel.bands) == 2 || throw(ArgumentError(
        "BCS kernel requires an orbital (band pair), got bands $(kernel.bands)"))
    band1, band2 = kernel.bands
    charges = _site_charges(spec, lattice)
    rows, cols = _bcs_legs(lattice, band1, band2)
    col_pos = [lattice[c] for c in cols]
    out = vacuumstate(lattice)
    for r in rows
        anchor = lattice[r]
        ws = Dict{Int,ComplexF64}()
        for (c, p) in zip(cols, col_pos)
            c == r && continue
            key = (ContourBlockIndex(branch(r), r.band), ContourBlockIndex(branch(c), c.band))
            # 键存在且腿 conj 与味表一致（(ā↑,a↓) 型误配对无块、跳过）
            haskey(kernel.blocks, key) &&
                (r.conj, c.conj) == _bcs_leg_conjs(r.band, c.band) || continue
            v = kernel.blocks[key][r.j, c.j]
            # 块值已带构造期约定（实时跨支 Keldysh 号已折入），直读即物理终值；
            # 同变量自配对 ξ² = 0 且系数位不可重复（上方 c == r 剔除）
            iszero(v) || (ws[p] = v)
        end
        isempty(ws) && continue
        ps = sort!(collect(keys(ws)))
        out, err = GM.mult(out, partial_row(spec, charges, anchor, ps, [ws[p] for p in ps]); trunc)
        push!(truncerrs, err)
    end
    out
end

"""
    bcs_hybriddynamics(lattice, bath::AbstractBCSBath, alg=PartialIF(); orbital=1, trunc) -> GMPS

BCS bath 的 IF 链便捷入口：`orbital` 选 Nambu 对（band1 = 2·orbital−1、
band2 = 2·orbital），对称性守卫 :z2（D7）。kernel 构造 =
`correlationfunction(bath, lattice; orbital)`。
"""
function bcs_hybriddynamics(lattice::AbstractGrassmannLattice, bath::AbstractBCSBath,
                            alg::InfluenceFunctionalAlgorithm=PartialIF();
                            orbital::Int=1, trunc::Truncation=Truncation())
    alg isa PartialIF || throw(ArgumentError(
        "BCS IF supports only PartialIF on the GM stack (D7), got $(typeof(alg))"))
    kernel = correlationfunction(bath, lattice; orbital)
    return influence_functional(lattice, kernel; alg, band=first(kernel.bands), trunc).state
end
