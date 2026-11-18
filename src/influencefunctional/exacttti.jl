# exacttti.jl：ExactTTI——tti_* 三通道因子构造 + GM expmult!×PolarExp 装配
# （plan/if_expmult W4；拟合层不动 = ExpExp anchored pencil 供给）
#
# 通道 → GM 因子容器映射（按块结构的可表达性裁定，dense 互换对拍裁决）：
#   - 虚时热通道：ExpDecaySumSum（锚点列拟合池）——纯虚时格的 τ/τ 环核 =
#     上下三角两条几何级联 + 等时接触修正（ExpDecaySum 的对角 = Σαₖ 语义，
#     (T,F)/(F,T) 两取向的对角相消，修正子 λ=0 只携接触）——极点压缩保留；
#     混合格 τ 段 → Pairing（τ 站显式表 + 全矩阵 C，ExpDecay cell 平铺无法
#     限制到分支子格）。
#   - 实时双通道（greater/lesser）：Pairing × 4（±± 同支、±∓ 跨支）——
#     ± 链拼接（− 支存储反序）使跨支 Toeplitz 结构不可用统一 cell 级联表达；
#     C = 拟合重构（Keldysh 关系审计通过的 (g,l) 拟合列 ⟹ 与现行 pole 因子
#     同一算符）+ 原始接触；混合格 → 原始块 C（无重构）。
#   - 混合交叉通道：Pairing × 2（τ→rt / rt→τ；块值直读；现行 rankone 发射
#     的 GM 原生对应物，SVD 秩一压缩由 PolarExp 承接 = cross_rtol 口径）。
#
# 审计保留：热配对（虚时）/ Keldysh 尾（实时）残差随 audit 返回，并入
# InfluenceFunctional.diagnostics（kernel diagnostics 字段已随 W2 删除）。

"""
    fit_exponentials(y; rtol=1e-10, maxterms=40, pencil_rows=128)
        -> (; weights, nodes, anchors, relative_error)

双端锚定矩阵铅笔（Hua–Sarkar：Hankel 行子采样 + 一次 SVD 最小阶扫描）——核心 =
`ExpExp.anchored_matrix_pencil_sweep`（`anchored_matrix_pencil` 的单 SVD 秩扫
变体）。`|q| > 1` 极点一律经恒等式 `q^(d−1) = q^(L−1)·(1/q)^(L−d)` 精确折叠为
`:right` 锚（权重吸收常数），返回节点全部模长 ≤ 1。

IF 域策略层（留本函数，ExpExp 不收——F6）：零序列 / `L == 1` 早退；铅笔扫完
不达 `rtol` 且 `L ≤ maxterms` 时以单位根 `q_j = 0.8·cis(2πj/L)` 做精确有界
插值（合法有限窗展开，非回退到其它 IF 算法）；仍不达则透传 sweep 的
ArgumentError（消息带最优残差）。
"""
function fit_exponentials(y; rtol=1e-10, maxterms=40, pencil_rows=128)
    0 < rtol < 1 && maxterms >= 1 && pencil_rows >= 1 || throw(ArgumentError("invalid fit settings"))
    y = ComplexF64.(y); L = length(y)
    all(isfinite, y) || throw(ArgumentError("nonfinite fit samples"))
    (L == 0 || all(iszero, y)) && return (; weights=ComplexF64[], nodes=ComplexF64[], anchors=Symbol[], relative_error=0.0)
    L == 1 && return (; weights=[y[1]], nodes=ComplexF64[0], anchors=[:left], relative_error=0.0)
    try
        w, nodes, anchors, err =
            anchored_matrix_pencil_sweep(y; tol=rtol, n=maxterms, pencil_rows=pencil_rows)
        return (; weights=w, nodes=nodes, anchors=anchors, relative_error=err)
    catch e
        e isa ArgumentError || rethrow()
        # Short finite windows can be underdetermined for a pencil. An exact bounded
        # exponential interpolation is a valid finite-window expansion, not a fit
        # fallback to another IF algorithm.
        L <= maxterms || rethrow()
        nodes = ComplexF64[0.8 * cis(2pi * j / L) for j in 0:L-1]
        B = [q^(d - 1) for d in 1:L, q in nodes]; w = B \ y
        err = norm(B * w - y) / norm(y)
        err <= rtol && return (; weights=w, nodes=nodes, anchors=fill(:left, L), relative_error=err)
        rethrow()
    end
end

# ---------------------------------------------------------------------------
# 站表与列读取助手
# ---------------------------------------------------------------------------

"链站号：格 (b,j) 的 band 带站 = (链块号 − 1)·bands + band"
_station(lattice::AbstractGrassmannLattice, cell::Tuple{Symbol,Int}, band::Int) =
    (cell_block(lattice, cell) - 1) * lattice.bands + band

_kernel_block(k::HybridizationKernel, b1::Symbol, b2::Symbol, band::Int) =
    k.blocks[(ContourBlockIndex(b1, band), ContourBlockIndex(b2, band))]

"拟合列的 lag 重建（ExactTTI 因子锚 = :left，:right 由增长极点守卫拒绝）"
_fitcolumn(f, m::Int) = sum(w * r^(m - 1) for (w, r) in zip(f.weights, f.nodes); init=0.0im)

# ---------------------------------------------------------------------------
# 虚时热通道：tti_imag_factors
# ---------------------------------------------------------------------------

"""
    tti_imag_factors(kernel, lattice; band) -> (; factors, audit)

热通道因子（τ/τ 核）：纯虚时格 = ExpDecaySumSum（上下三角级联 + 等时接触，
锚点列拟合池——现行 imag_factors 的极点内容）；混合格 τ 段 = Pairing（τ 站
全矩阵）。热配对审计保留（thermal_residual ∈ audit）。
"""
function tti_imag_factors(kernel::HybridizationKernel, lattice::AbstractGrassmannLattice;
                          band::Int, alg=ExactTTI())
    lattice isa Union{AbstractImagGrassmannLattice,AbstractMixedGrassmannLattice} ||
        return (; factors=GM.GrassmannFactor[], audit=NamedTuple())
    n = lattice.kτ
    M = _kernel_block(kernel, :τ, :τ, band)
    lower = ComplexF64[M[j+1, 1] for j in 1:n-1]
    upper = ComplexF64[M[1, j+1] for j in 1:n-1]
    thermal = norm(upper + reverse(lower)) / max(norm(lower), eps())
    thermal <= 1e-12 || throw(ArgumentError("imaginary kernel violates thermal pairing"))
    audit = (; thermal_residual=thermal)
    if lattice isa AbstractImagGrassmannLattice
        f = fit_exponentials(lower; rtol=alg.rtol, maxterms=alg.maxterms)
        audit = (; audit..., imag_expansion=(; f..., fit_interval=(1, n - 1),
                                               convention=:integrated_lag_one))
        # 环核的开放双线性分解（拟合约定 y[d] = Σₖ wₖ qₖ^{d−1}，环反周期
        # upper[m] = −lower[n−m] 已由热配对审计锚定）：
        #   上三角掩码 {ā_s, a_t}（链距 m = t−s）：系数 lower[m] = Σ (w/q)·q^m；
        #   下三角掩码 {a_s, ā_t}（链距 m）：系数 −upper[m] = lower[n−m]
        #     = Σ (w·q^{n−1})·(1/q)^m（环绕方向——反周期折叠）；
        #   等时对角：两取向级联的对角（Σα₁ 与 −Σα₂）相消后余 M[1,1] →
        #   λ=0 修正式携 M[1,1] − Σα₁ + Σα₂（每 τ 格同值——Toeplitz 环核等时元）。
        # 拟合节点含 ≈0 极点时（环绕级联 1/q 退化）→ Pairing 全矩阵精确兜底。
        B = lattice.bands
        stationτ = (cell=B, asite=band, bsite=band)
        if isempty(f.weights)
            return (; factors=GM.GrassmannFactor[], audit)
        elseif all(abs(q) > 1e-8 for q in f.nodes)
            α1 = f.weights ./ f.nodes
            α2 = f.weights .* f.nodes .^ (n - 1)
            sums = [GM.ExpDecaySum(α1, f.nodes; conja=true, conjb=false, stationτ...),
                    GM.ExpDecaySum(α2, 1 ./ f.nodes; conja=false, conjb=true, stationτ...)]
            contact = M[1, 1] - sum(α1) + sum(α2)
            iszero(contact) || push!(sums, GM.ExpDecaySum([contact], [0.0];
                conja=true, conjb=false, stationτ...))
            return (; factors=[GM.ExpDecaySumSum(sums)], audit)
        else
            stations = [_station(lattice, (:τ, j), band) for j in 1:n]
            C = Matrix{ComplexF64}(M)
            return (; factors=[GM.Pairing(C, stations, stations; conja=true, conjb=false)], audit)
        end
    end
    # 混合格：τ 段 Pairing（全矩阵，PolarExp SVD 承接压缩）
    stations = [_station(lattice, (:τ, j), band) for j in 1:n]
    C = Matrix{ComplexF64}(M)
    return (; factors=[GM.Pairing(C, stations, stations; conja=true, conjb=false)], audit)
end

# ---------------------------------------------------------------------------
# 实时双通道：tti_real_factors（greater/lesser 两侧列）
# ---------------------------------------------------------------------------

"Keldysh 尾审计（structured Lag 块近似：直读 kernel 的 Toeplitz 性；现行
real_factors 同款）+ (g, l) 拟合列。返回 (fg, fl, keldysh_residual)"
function _real_fits(kernel, lattice, band, alg)
    n = lattice.k
    g = ComplexF64[_kernel_block(kernel, :+, :+, band)[i+1, 1] for i in 1:n-1]
    l = ComplexF64[_kernel_block(kernel, :+, :-, band)[i+1, 1] for i in 1:n-1]
    locations = [(max(d, 0) + 1, max(-d, 0) + 1) for d in -(n - 1):n-1]
    error2 = 0.0; total2 = 0.0
    blk(a, b) = _kernel_block(kernel, a, b, band)
    for (i, j) in locations
        i == j && continue
        v = i > j ? (g[i-j], l[i-j], -g[i-j], -l[i-j]) :
                   (-conj(l[j-i]), conj(l[j-i]), -conj(g[j-i]), conj(g[j-i]))
        weight = n - abs(i - j)
        for (pair, x) in zip(((:+, :+), (:+, :-), (:-, :+), (:-, :-)), v)
            error2 += weight * abs2(x - blk(pair...)[i, j])
            total2 += weight * abs2(blk(pair...)[i, j])
        end
    end
    err = sqrt(error2) / max(sqrt(total2), eps())
    err <= max(alg.rtol, 1e-12) || throw(ArgumentError("ExactTTI real Keldysh relations fail: $err"))
    if lattice isa AbstractRealGrassmannLattice
        fg = fit_exponentials(g; rtol=alg.rtol, maxterms=alg.maxterms)
        fl = fit_exponentials(l; rtol=alg.rtol, maxterms=alg.maxterms)
        all(a -> a === :left, fg.anchors) && all(a -> a === :left, fl.anchors) || throw(ArgumentError(
            "growing real pole unsupported; use PartialIF"))
        return fg, fl, err
    end
    return nothing, nothing, err
end

"""
    tti_real_factors(kernel, lattice; band) -> (; factors, audit)

实时双通道因子（greater/lesser 两侧列）：Pairing × 4（(+,+)、(−,−)、(+,−)、
(−,+)），C = 拟合重构（Keldysh 关系审计通过的 (g,l) 列 ⟹ 与现行 pole 记忆
元胞因子同一算符）+ 等时接触（blocks[·][1,1]，现行 contact D 同源）；混合格
rt 块 = 原始块 C（无重构）。Keldysh 尾审计保留（keldysh_residual ∈ audit）。
"""
function tti_real_factors(kernel::HybridizationKernel, lattice::AbstractGrassmannLattice;
                          band::Int, alg=ExactTTI())
    lattice isa Union{AbstractRealGrassmannLattice,AbstractMixedGrassmannLattice} ||
        return (; factors=GM.GrassmannFactor[], audit=NamedTuple())
    n = _branch_len(lattice, :+)                    # 纯实时 = k；混合 = kt
    st_p = [_station(lattice, (:+, j), band) for j in 1:n]
    st_m = [_station(lattice, (:-, j), band) for j in 1:n]
    audit = NamedTuple()
    if lattice isa AbstractRealGrassmannLattice
        fg, fl, kerr = _real_fits(kernel, lattice, band, alg)
        audit = (keldysh_residual=kerr,
                 real_expansions=(greater=(; fg..., fit_interval=(1, n - 1), convention=:integrated_lag_one),
                                  lesser=(; fl..., fit_interval=(1, n - 1), convention=:integrated_lag_one)),
                 real_fit_residuals=[fg.relative_error, fl.relative_error])
        recg(m) = _fitcolumn(fg, m); recl(m) = _fitcolumn(fl, m)
        cpp(i, j) = i == j ? _kernel_block(kernel, :+, :+, band)[1, 1] :
                    i > j ? recg(i - j) : -conj(recl(j - i))
        cmm(i, j) = i == j ? _kernel_block(kernel, :-, :-, band)[1, 1] :
                    i > j ? -recl(i - j) : conj(recg(j - i))
        cpm(i, j) = i == j ? _kernel_block(kernel, :+, :-, band)[1, 1] :
                    i > j ? recl(i - j) : conj(recl(j - i))
        cmp(i, j) = i == j ? _kernel_block(kernel, :-, :+, band)[1, 1] :
                    i > j ? -recg(i - j) : -conj(recg(j - i))
        Cpp = [cpp(i, j) for i in 1:n, j in 1:n]
        Cmm = [cmm(i, j) for i in 1:n, j in 1:n]
        Cpm = [cpm(i, j) for i in 1:n, j in 1:n]
        Cmp = [cmp(i, j) for i in 1:n, j in 1:n]
        factors = [GM.Pairing(Cpp, st_p, st_p; conja=true, conjb=false),
                   GM.Pairing(Cmm, st_m, st_m; conja=true, conjb=false),
                   GM.Pairing(Cpm, st_p, st_m; conja=true, conjb=false),
                   GM.Pairing(Cmp, st_m, st_p; conja=true, conjb=false)]
        return (; factors, audit)
    end
    # 混合格：原始块 C（rt 块直读，无重构）
    factors = [GM.Pairing(Matrix{ComplexF64}(_kernel_block(kernel, a, b, band)),
                          a === :+ ? st_p : st_m, b === :+ ? st_p : st_m;
                          conja=true, conjb=false)
               for a in (:+, :-), b in (:+, :-)]
    return (; factors, audit)
end

# ---------------------------------------------------------------------------
# 混合交叉通道：tti_cross_factors
# ---------------------------------------------------------------------------

"""
    tti_cross_factors(kernel, lattice; band) -> (; factors, audit)

混合交叉通道因子（τ↔rt 六块）：Pairing × 2（τ→rt 行腿 = ā_τ、rt→τ 行腿 =
ā_rt；列腿恒 a）——块值直读（现行 rankone 发射的 GM 原生对应物；旧实现
rank-one 系数向量落在 a 型腿上 = 潜在味错位，新构造以 ā/a 腿直落核值、由
ExactTTI ≡ PartialIF 互换对拍裁决）。秩审计保留（cross_ranks/cross_residuals
∈ audit；SVD 秩一压缩由 PolarExp 按 cross_rtol 承接）。
"""
function tti_cross_factors(kernel::HybridizationKernel, lattice::AbstractGrassmannLattice;
                           band::Int, alg=ExactTTI())
    lattice isa AbstractMixedGrassmannLattice && lattice.Nt > 0 ||
        return (; factors=GM.GrassmannFactor[], audit=NamedTuple())
    st_τ = [_station(lattice, (:τ, j), band) for j in 1:lattice.Nτ]
    st_rt = vcat([_station(lattice, (:+, j), band) for j in 1:lattice.Nt],
                 [_station(lattice, (:-, j), band) for j in 1:lattice.Nt])
    factors = GM.GrassmannFactor[]
    ranks = Int[]; residuals = Float64[]
    for (rowbs, colbs, asites, bsites) in (((:τ,), (:+, :-), st_τ, st_rt),
                                           ((:+, :-), (:τ,), st_rt, st_τ))
        M = vcat([hcat([_kernel_block(kernel, a, b, band) for b in colbs]...) for a in rowbs]...)
        U, S, V = MAK.svd_compact(M); sv = diag(S); tot = sum(abs2, sv); r = length(sv)
        while r > 0 && sum(abs2, sv[r:end]) <= alg.cross_rtol^2 * tot; r -= 1; end
        push!(ranks, r); push!(residuals, iszero(tot) ? 0.0 : sqrt(sum(abs2, sv[r+1:end]) / tot))
        isempty(sv) || push!(factors, GM.Pairing(Matrix{ComplexF64}(M), asites, bsites;
                                                 conja=true, conjb=false))
    end
    audit = (cross_ranks=ranks, cross_residuals=residuals)
    return (; factors, audit)
end

# ---------------------------------------------------------------------------
# ExactTTI 装配：vacuumstate → 逐因子 expmult!×PolarExp
# ---------------------------------------------------------------------------

"ExactTTI 组装链：tti_real/imag/cross 因子逐个 expmult!×PolarExp（因子偶元
互相超交换 ⟹ 应用序无关）。PolarExp 对 Pairing/AnomalousPairing = SVD 秩一
精确指数（分量 Xₖ² = 0 ⟹ e^{Xₖ} = 1 + Xₖ）；对 ExpDecay 族（虚时 ESS）委托
TTIExp = ExpDecayExp 逐通道物化——**要求 GM ≥ expdecayexp-same-orbital-fix
（a7b39bb）**：同轨道（asite == bsite）触碰块已删，物化对朴素 exp 严格
（触碰项与 local × 过线配对 pfaffian 逐项相消 ⟹ 净零；修正前同轨道偏差
O(α²λ)，GTEMPO ESS 实测 3.2e-4）。截秩语义：real/imag 通道 rtol = 0（压缩
由拟合重构自带，不叠加截断），cross 通道 rtol = cross_rtol（口径同旧 cross
通道）。factor_rtol 映射到乘法压缩策略（旧逐因子 cleanup 语义）。"
function _exacttti_state(lattice::AbstractGrassmannLattice, kernel::HybridizationKernel,
                         alg::ExactTTI, band::Int, trunc::Truncation)
    diag = Dict{Symbol,Any}()
    mult = alg.factor_rtol === nothing ? SVDCompression(trunc=trunc) :
           SVDCompression(trunc=Truncation(D=typemax(Int), rtol=alg.factor_rtol))
    pol_exact = GM.PolarExp(rtol=0.0, alg_mult=mult)             # real/imag 通道：精确
    pol_cross = GM.PolarExp(rtol=alg.cross_rtol, alg_mult=mult)  # cross 通道：截秩
    state = vacuumstate(lattice)
    truncerrs = Float64[]
    for (tti, pol) in ((tti_real_factors, pol_exact), (tti_imag_factors, pol_exact),
                       (tti_cross_factors, pol_cross))
        (; factors, audit) = tti(kernel, lattice; band, alg)
        isempty(audit) || merge!(diag, pairs(audit))
        for f in factors
            state, err = GM.expmult!(state, f, pol)
            push!(truncerrs, err)
        end
    end
    return state, truncerrs, diag
end
