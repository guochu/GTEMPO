# influence.jl：浴影响泛函 IF（plan §Phase 3，D3 首批 = PartialIF + ExactTTI）
#
# 移植自 U1GTEMPO2 influence.jl（charge 通道机器精度、与 secondorderif 端到端
# ≤1.4e-11 验证），泛化点（D1）：
#   - spin::Symbol(:up/:down) → band::Int（band b 的生成元荷 = var_charge(spec,b)；
#     双荷形 bands≤2 时 band↔槽位，单荷形全带同荷、带=站位）；
#   - kernel blocks 为 canonical Δ[(a,b)][i,j]（IF = exp(Σ ā_i Δ_ij a_j)），
#     由 hybridization() 从 QuAPI corr 提取（含实时跨支/混合 τ↔rt 符号定标）。
#
# TDVPIF/XTRGIF/TTIIF 暂停（D3 裁决：第二批 TTIIF/XTRGIF 走 GM.SchurGMPS/
# timeevo/expmult! 另行移植）。

# ---------------------------------------------------------------------------
# kernel（canonical 块 + QuAPI 提取）
# ---------------------------------------------------------------------------

"""
    HybridizationKernel

单元级杂化核（canonical 约定）：`blocks[(b1,b2)][i,j]` = IF 指数系数
ā_i^{b1} Δ a_j^{b2}（i, j = 单元下标 1..N）。`coupling_norm` = 全带总耦合
Σ_b α_b（Weyl 标量补偿用；order=2 必需）。
"""
struct HybridizationKernel{L}
    lattice::L
    blocks::Dict{Tuple{Symbol,Symbol},Matrix{ComplexF64}}
    coupling_norm::Union{Nothing,Float64}
    diagnostics::Dict{Symbol,Any}
end

function HybridizationKernel(lattice, blocks; coupling_norm=nothing, diagnostics=Dict{Symbol,Any}())
    out = Dict{Tuple{Symbol,Symbol},Matrix{ComplexF64}}()
    for a in branches(lattice), b in branches(lattice)
        haskey(blocks, (a, b)) || throw(ArgumentError("missing kernel block $a/$b"))
        m = blocks[(a, b)]
        na = a === :τ ? lattice.kτ : (lattice isa AbstractRealGrassmannLattice ? lattice.k : lattice.Nt)
        nb = b === :τ ? lattice.kτ : (lattice isa AbstractRealGrassmannLattice ? lattice.k : lattice.Nt)
        size(m) == (na, nb) || throw(DimensionMismatch("kernel shape $a/$b: $(size(m)) ≠ ($na,$nb)"))
        all(isfinite, m) || throw(ArgumentError("nonfinite kernel block $a/$b"))
        out[(a, b)] = ComplexF64.(m)
    end
    coupling_norm === nothing || (isfinite(coupling_norm) && coupling_norm >= 0) ||
        throw(ArgumentError("invalid coupling norm"))
    HybridizationKernel(lattice, out, coupling_norm === nothing ? nothing : Float64(coupling_norm),
                        Dict{Symbol,Any}(diagnostics))
end

"""
    hybridization(lattice, corr; coupling_norm=nothing) -> HybridizationKernel

从 QuAPI 离散化函数直接提取 canonical 单元级核（IF = exp(Σ ā_i Δ_ij a_j)）——
secondorderif kernels.jl 同款三路径（ED 验证过）：
- 虚时：Δτ 的单元矩阵 [1:Nτ, 1:Nτ]；
- 实时：Δt 谱路径（离散谱精确极点装配），四块 [1:Nt, 1:Nt] 直读
  （块值 = 单元级全宽积分；QuAPI 内嵌的 Keldysh 跨支号与旧消费端 gate
  在净语义中相消 ⟹ canonical = raw）；
- 混合：Δm 九块，rt [1:Nt]/τ [1:Nτ]；Δm 的 τ↔rt 交叉块相对物理核差一个
  整体负号（KernelCorrView 预翻 = 净读数）⟹ canonical 取负。
`coupling_norm` 由调用方传入（全带总耦合 Σ_b α_b）。
"""
function hybridization(lattice::AbstractGrassmannLattice, bath::AbstractBath; coupling_norm=nothing)
    blocks = Dict{Tuple{Symbol,Symbol},Matrix{ComplexF64}}()
    if lattice isa AbstractImagGrassmannLattice
        corr = QuAPI.Δτ(bath; N=lattice.Nτ, δτ=lattice.δτ)
        blocks[(:τ, :τ)] = Matrix{ComplexF64}(corr.data[1:lattice.Nτ, 1:lattice.Nτ])
    elseif lattice isa AbstractRealGrassmannLattice
        qt = QuAPI.Δt(bath; N=lattice.Nt, t=lattice.Nt * lattice.δt)
        for (b1, b2) in ((:+, :+), (:+, :-), (:-, :+), (:-, :-))
            Mq = QuAPI.branch(qt; b1=b1, b2=b2)
            blocks[(b1, b2)] = Matrix{ComplexF64}(Mq[1:lattice.Nt, 1:lattice.Nt])
        end
    else
        qm = QuAPI.Δm(bath; Nτ=lattice.Nτ, t=lattice.Nt * lattice.δt, Nt=lattice.Nt, δτ=lattice.δτ)
        ri(c2) = c2 === :τ ? (1:lattice.Nτ) : (1:lattice.Nt)
        for b1 in (:+, :-, :τ), b2 in (:+, :-, :τ)
            Mq = QuAPI.branch(qm; b1=b1, b2=b2)
            M = Matrix{ComplexF64}(Mq[ri(b1), ri(b2)])
            (b1 === :τ) !== (b2 === :τ) && (M = -M)   # Δm τ↔rt 交叉块内嵌负号
            blocks[(b1, b2)] = M
        end
    end
    return HybridizationKernel(lattice, blocks; coupling_norm)
end

function hybridization(lattice::AbstractGrassmannLattice, corr; coupling_norm=nothing)
    throw(ArgumentError("hybridization 需要浴对象（QuAPI 直接提取）；经 correlationfunction 的 accessor 路线已随 T7 标定退役"))
end

# ---------------------------------------------------------------------------
# PartialIF：逐行秩一因子积（精确行指数；U1GTEMPO2 partial_if 同型）
# ---------------------------------------------------------------------------

"""
    partial_row(spec, charges, anchor, columns, values) -> GMPS

第 anchor 行的偏 IF 因子 exp(ā_anchor Σ_j Δ_j a_j) 的链表示
（U1GTEMPO2 partial_row 的对称性泛型版）。
"""
function partial_row(spec::AbstractSymmetrySpec, charges, anchor, columns, values)
    qbar = charges[anchor]; q = qneg(qbar); z = zero_charge(spec)
    ws = Dict(zip(columns, values))
    data = map(eachindex(charges)) do i
        ls = i == 1 ? [z] : [z, i <= anchor ? q : qbar]
        rs = i == length(charges) ? [z] : [z, i < anchor ? q : qbar]
        t = TK.zeros(ComplexF64, qspace(spec, ls) ⊗ pspace(spec, charges[i]) ← qspace(spec, rs))
        t[(sector(spec, z), sector(spec, z), sector(spec, z))] .= 1
        if i == anchor
            q in ls && (t[(sector(spec, q), sector(spec, qbar), sector(spec, z))] .= 1)
            qbar in rs && (t[(sector(spec, z), sector(spec, qbar), TK.dual(sector(spec, qbar)))] .= 1)
        else
            pending = i < anchor ? q : qbar
            pending in ls && pending in rs &&
                (t[(sector(spec, pending), sector(spec, z), TK.dual(sector(spec, pending)))] .= 1)
            c = get(ws, i, 0im)
            if !iszero(c)
                if i < anchor
                    t[(sector(spec, z), sector(spec, q), TK.dual(sector(spec, q)))] .= -c
                else
                    t[(sector(spec, qbar), sector(spec, q), sector(spec, z))] .= c
                end
            end
        end
        t
    end
    pairize(spec, data)
end

"""
    partial_if(lattice, kernel; band=1, trunc=Truncation(), reports) -> GMPS

PartialIF：IF = Π_row exp(ā_row Σ_j Δ_row,j a_j) 的逐行因子积（精确）。
`band` 选择 IF 所在带（生成元荷/站位随带）。行/列序 = 链存储格序
（_storage_cells），values 与 cols 逐位对齐。
"""
function partial_if(lattice::AbstractGrassmannLattice, kernel::HybridizationKernel;
                    band::Int=1, trunc::Truncation=Truncation(), reports=CompressionReport[])
    spec = symmetry_spec(lattice.symmetry)
    out = vacuumstate(lattice)
    charges = _site_charges(spec, lattice)
    cells = _storage_cells(lattice)              # 存储格序
    cols = [_cell_var_index(lattice, cell, band, false) for cell in cells]   # a 变量位
    ncells = length(cells)
    for (m, (b, j)) in enumerate(cells)
        # 行 = 格 (b,j)（branch 局部下标 j）；列值按存储格序逐位取核
        values = [kernel.blocks[(b, cell[1])][j, cell[2]] for cell in cells]
        all(iszero, values) && continue
        anchor = _cell_var_index(lattice, (b, j), band, true)   # ā 变量位
        out = multiply(out, partial_row(spec, charges, anchor, cols, values); trunc, reports)
    end
    out
end

"IF 变量荷序列（全链 = nblocks 格 × 每格 2·bands 变量，存储序）"
function _site_charges(spec::AbstractSymmetrySpec, lattice)
    per_cell = cell_charges(spec, lattice.bands)
    return repeat(per_cell, nblocks(lattice))
end

"单元 (b,j) 的带 band 变量位：conj=false → a（奇位）、true → ā（偶位）"
function _cell_var_index(lattice, cell::Tuple{Symbol,Int}, band::Int, conj::Bool)
    return block_offset(lattice, cell_block(lattice, cell)) + 2 * (band - 1) + Int(conj)
end

# ---------------------------------------------------------------------------
# 矩阵铅笔拟合（U1GTEMPO2 fit_exponentials 原样）
# ---------------------------------------------------------------------------

"Two-end anchored matrix pencil: all returned propagation nodes have modulus ≤ 1."
function fit_exponentials(y; rtol=1e-10, maxterms=40, pencil_rows=128, rightanchors::Bool=true)
    0 < rtol < 1 && maxterms >= 1 && pencil_rows >= 1 || throw(ArgumentError("invalid fit settings"))
    y = ComplexF64.(y); L = length(y)
    all(isfinite, y) || throw(ArgumentError("nonfinite fit samples"))
    (L == 0 || all(iszero, y)) && return (; weights=ComplexF64[], nodes=ComplexF64[], anchors=Symbol[], relative_error=0.0)
    L == 1 && return (; weights=[y[1]], nodes=ComplexF64[0], anchors=[:left], relative_error=0.0)
    h = L ÷ 2
    rows = h <= pencil_rows ? collect(0:h-1) : unique(round.(Int, range(0, h - 1; length=pencil_rows)))
    H0 = [y[i+j] for i in rows, j in 1:L-h]; H1 = [y[i+j+1] for i in rows, j in 1:L-h]
    U, S, V = MAK.svd_compact(H0); best = Inf
    for r in 1:min(maxterms, length(diag(S)))
        S[r, r] > 0 || break
        A = (U[:, 1:r]' * H1 * V[1:r, :]') * Diagonal(1 ./ diag(S)[1:r])
        nodes = ComplexF64.(MAK.eig_vals(A)); all(isfinite, nodes) || continue
        nodes = [!rightanchors && abs(q) > 1 ? q / abs(q) : q for q in nodes]
        anchors = [rightanchors && abs(q) > 1 ? :right : :left for q in nodes]
        nodes = [a == :right ? inv(q) : q for (a, q) in zip(anchors, nodes)]
        B = ComplexF64[q^(a == :left ? d - 1 : L - d) for d in 1:L, (q, a) in zip(nodes, anchors)]
        u, s, v = MAK.svd_compact(B); sv = diag(s); keep = sv .> maximum(sv) * 1e-14
        w = v[keep, :]' * ((u[:, keep]' * y) ./ sv[keep]); err = norm(B * w - y) / norm(y); best = min(best, err)
        err <= rtol && return (; weights=w, nodes, anchors, relative_error=err)
    end
    # Short finite windows can be underdetermined for a pencil. An exact bounded
    # exponential interpolation is a valid finite-window expansion, not a fit
    # fallback to another IF algorithm.
    if L <= maxterms
        nodes = ComplexF64[0.8 * cis(2pi * j / L) for j in 0:L-1]
        B = [q^(d - 1) for d in 1:L, q in nodes]; w = B \ y
        err = norm(B * w - y) / norm(y)
        err <= rtol && return (; weights=w, nodes, anchors=fill(:left, L), relative_error=err)
    end
    throw(ArgumentError("matrix-pencil fit residual $best exceeds $rtol with maxterms=$maxterms"))
end

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
                      init=zero_charge(spec)) for mask in 0:(1 << length(memorycharges))-1]
    for (k, vars) in enumerate(cellvars)
        left = k == 1 ? [1] : collect(eachindex(labels)); right = k == length(cellvars) ? [1] : collect(eachindex(labels))
        ts = split_cell(spec, localdata[left, :, right], labels[left], labels[right], charges[vars])
        for (v, t) in zip(vars, ts); data[v] = t; end
    end
    # Spectator sites propagate the current virtual space with vacuum coefficient.
    active = Set(v for vars in cellvars for v in vars)
    for i in 1:L
        i in active && continue
        v = i == 1 ? qspace(spec, [zero_charge(spec)]) : TK.domain(data[i-1])[1]
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

validate_tti(M::AbstractMatrix; tol) = _validate_tti_matrix(M, tol)
function _validate_tti_matrix(M, tol)
    n = size(M, 1); size(M, 2) == n || throw(DimensionMismatch())
    err = maximum((abs(M[i, j] - M[i >= j ? i - j + 1 : 1, i >= j ? 1 : j - i + 1]) for i in 1:n, j in 1:n); init=0.0)
    err <= tol * max(norm(M), eps()) || throw(ArgumentError("kernel is not TTI: residual $err")); err
end

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
    blocks = k.blocks
    residuals = Float64[]; nfactors = 0
    for a in (:+, :-), b in (:+, :-); validate_tti(blocks[(a, b)]; tol=alg.rtol); end
    g = ComplexF64[blocks[(:+, :+)][i + 1, 1] for i in 1:n-1]
    l = ComplexF64[blocks[(:+, :-)][i + 1, 1] for i in 1:n-1]
    # Keldysh tails audit（structured Lag 块近似：直读 kernel 的 Toeplitz 性）
    locations = [(max(d, 0) + 1, max(-d, 0) + 1) for d in -(n - 1):n-1]
    error2 = 0.0; total2 = 0.0
    for (i, j) in locations
        i == j && continue
        v = i > j ? (g[i-j], l[i-j], -g[i-j], -l[i-j]) : (-conj(l[j-i]), conj(l[j-i]), -conj(g[j-i]), conj(g[j-i]))
        weight = n - abs(i - j)
        for (pair, x) in zip(((:+, :+), (:+, :-), (:-, :+), (:-, :-)), v)
            error2 += weight * abs2(x - blocks[pair][i, j]); total2 += weight * abs2(blocks[pair][i, j])
        end
    end
    err = sqrt(error2) / max(sqrt(total2), eps())
    err <= max(alg.rtol, 1e-12) || throw(ArgumentError("ExactTTI real Keldysh relations fail: $err"))
    diagnostics[:keldysh_residual] = err
    q = var_charge(spec, band, lattice.bands)
    vars = [_virtual_merged_vars(lattice, band, m) for m in 1:n]   # 虚拟 ± 相邻布局
    perm = _virtual_to_chrono_perm(lattice, n)
    for (component, y) in ((:greater, g), (:lesser, l))
        f = fit_exponentials(y; rtol=alg.rtol, maxterms=alg.maxterms, rightanchors=false)
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
    D = [blocks[(a, b)][1, 1] for a in (:+, :-), b in (:+, :-)]
    all(iszero, D) || emit(contact_factor(spec, lattice, band, (:+, :-), D))
    nothing
end

function imag_factors(emit, spec, lattice, k::HybridizationKernel, band, alg, diagnostics)
    lattice isa Union{AbstractImagGrassmannLattice,AbstractMixedGrassmannLattice} || return nothing
    z = zero_charge(spec)
    n = lattice.kτ
    M = k.blocks[(:τ, :τ)]; validate_tti(M; tol=alg.rtol)
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
    z = zero_charge(spec)
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
        M = vcat([hcat([k.blocks[(a, b)] for b in colbs]...) for a in rowbs]...)
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

# ---------------------------------------------------------------------------
# IF 算法类型与入口（公开皮肤保名）
# ---------------------------------------------------------------------------

abstract type InfluenceFunctionalAlgorithm end

"""
    PartialIF(; maxdim=128, rtol=1e-12, cutoff=0)

逐行秩一因子积 IF（精确行指数 + 乘间压缩；SciPost Phys. Core 7, 063 (2024)）。
"""
Base.@kwdef struct PartialIF <: InfluenceFunctionalAlgorithm
    maxdim::Int = 128
    rtol::Float64 = 1e-12
    cutoff::Float64 = 0.0
end
_trunc(alg::PartialIF) = Truncation(maxdim=alg.maxdim, rtol=alg.rtol, cutoff=alg.cutoff)

"""
    ExactTTI(; rtol=1e-10, maxterms=40, cross_rtol=rtol, factor_rtol=nothing)

矩阵铅笔拟合 IF（U1GTEMPO2 ExactTTI 同型）：实轴 greater/lesser 双通道 +
虚时热通道 + 混合交叉通道，各展成指数因子精确楔积。
"""
struct ExactTTI <: InfluenceFunctionalAlgorithm
    rtol::Float64
    maxterms::Int
    cross_rtol::Float64
    factor_rtol::Union{Nothing,Float64}
    function ExactTTI(; rtol::Real=1e-10, maxterms::Integer=40, cross_rtol::Real=rtol,
                      factor_rtol::Union{Nothing,Real}=nothing)
        Float64(rtol); Int(maxterms)
        0 < rtol < 1 && 0 < cross_rtol < 1 && maxterms >= 1 ||
            throw(ArgumentError("invalid IF fit settings"))
        factor_rtol === nothing || (isfinite(factor_rtol) && 0 <= factor_rtol < 1) ||
            throw(ArgumentError("invalid factor cleanup tolerance"))
        new(Float64(rtol), Int(maxterms), Float64(cross_rtol),
            factor_rtol === nothing ? nothing : Float64(factor_rtol))
    end
end

struct InfluenceFunctional{K,S}
    kernel::K
    state::S
    band::Int
    diagnostics::Dict{Symbol,Any}
    reports::Vector{CompressionReport}
end

"""
    influence_functional(lattice, kernel; alg=PartialIF(), band=1, trunc=Truncation()) -> InfluenceFunctional

IF 构造入口（Phase 3 首批：PartialIF / ExactTTI）。
"""
function influence_functional(lattice::AbstractGrassmannLattice, kernel::HybridizationKernel;
                              alg::InfluenceFunctionalAlgorithm=PartialIF(), band::Int=1,
                              trunc::Truncation=Truncation())
    1 <= band <= lattice.bands || throw(BoundsError(1:lattice.bands, band))
    reports = CompressionReport[]
    diag = copy(kernel.diagnostics)
    diag[:algorithm] = alg; diag[:truncation] = trunc
    spec = symmetry_spec(lattice.symmetry)
    if alg isa PartialIF
        state = partial_if(lattice, kernel; band, trunc=_trunc(alg), reports)
    elseif alg isa ExactTTI
        state = vacuumstate(lattice)
        cleanup_reports = CompressionReport[]
        for builder in (real_factors, imag_factors, cross_factors)
            builder(function (f)
                        if alg.factor_rtol !== nothing
                            push!(cleanup_reports, compress!(f; trunc=Truncation(maxdim=typemax(Int), rtol=alg.factor_rtol)))
                        end
                        state = multiply(state, f; trunc, reports)
                    end, spec, lattice, kernel, band, alg, diag)
        end
        diag[:factor_cleanup_reports] = cleanup_reports
    else
        throw(ArgumentError("unsupported IF algorithm $(typeof(alg))"))
    end
    InfluenceFunctional(kernel, state, band, diag, reports)
end

# ---------------------------------------------------------------------------
# 公开皮肤：hybriddynamics（旧 API 名；corr → kernel → IF 链）
# ---------------------------------------------------------------------------

"""
    hybriddynamics(lattice, corr, alg=PartialIF(); band=1, coupling_norm=nothing) -> GMPS

从 QuAPI 关联函数构建 `band` 带的影响泛函链（influence_functional 便捷入口）。
`coupling_norm`（全带总耦合 Σ_b α_b）仅在配合 Weyl 二阶 K 时需要，随 kernel
携带给 sysdynamics_weyl2order。
"""
function hybriddynamics(lattice::AbstractGrassmannLattice, bath::AbstractBath,
                        alg::InfluenceFunctionalAlgorithm=PartialIF();
                        band::Int=1, coupling_norm::Union{Nothing,Real}=nothing,
                        trunc::Truncation=Truncation())
    kernel = hybridization(lattice, bath; coupling_norm=coupling_norm)
    return influence_functional(lattice, kernel; alg, band, trunc).state
end

"integrate(lattice, x0, x1...) —— Berezin 迹（链 = (K, Is...)，GM 引擎）"
integrate(lattice::AbstractGrassmannLattice, x0::GM.GMPS, x1::GM.GMPS...) = integrate(x0, x1...)
integrate_scaled(lattice::AbstractGrassmannLattice, x0::GM.GMPS, x1::GM.GMPS...) = integrate_scaled(x0, x1...)
