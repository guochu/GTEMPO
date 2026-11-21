# exacttti.jl：ExactTTI——任意耦合模式的 TTI 因子构造 + GM expmult!×PolarExp 装配
# （plan/if_expmult W4 + bcs_exacttti_plan Phase A′/B′（2026-11-19：is_tti 门控 +
#   band/味守卫退役 + 任意多带核适用面）；混合时间轴面退役——交叉通道
#   tti_cross_factors 与 cross_rtol 删除；审查二批：tti_factors 按轴 dispatch
#   （虚轴/实轴两方法，混合格无方法 = MethodError），通道函数退役）
#
# 因子构造按味对 (ba, ca, bb, cb) 分组遍历核键，**严格 TTI 门控**（2026-11-19
# 审查裁决：检查失败响抛，不降级直读）：
#   - `is_tti` 块级门（CorrelationMatrix 结构性成立，稠密块逐对角线 Toeplitz
#     检验）——非 TTI 块 ArgumentError（ExactTTI = 严格 TTI 算法，非 TTI 核是
#     调用方错误，PartialIF/PolarIF 承接任意块）；
#   - 物理恒等式审计响抛：虚时热配对（环反周期闭合）/ 实时 Keldysh 关系
#     （逐味对）；实时另要求四支块齐备；
#   - 唯一非 fit 出口 = 虚时退化拟合节点（|q| ≈ 0，环绕族 1/q 失格）→
#     Pairing 全矩阵兜底（拟合质量问题而非核形状问题，现行语义保留）。
#
# 虚时环族（任意味对，2026-11-19 广义化；band-diagonal (T,F) 现行公式 = ba==bb
# 特例，逐字退化；dense 算符对拍实证：单对 bitwise + 完整 BCS 四键 3.5e-18；
# 2026-11-20 v7 虚轴翻转——链站序与核胞序相同，通道按新方向重推，见
# plan/unified_chain_ordering_plan.md §7）：
# 链站序与核胞序相同（cell_block(:τ, j) = j）。拟合列 lower[m] = M[m+1, 1]
# （y[d] = Σ wₖ qₖ^{d−1}），三通道（取向读号：原取向 +系数、翻转取向 −系数）：
#   ch1 = ExpDecaySum(α₁ = −w·q^{n−1}, 1/q; ca, cb, asite=ba, bsite=bb)——链前向
#         {jr < jc} 条目，系数 upper[d] = −lower[n−d]（反周期折叠）；
#   ch2 = ExpDecaySum(α₂ = −w/q, q; cb, ca, asite=bb, bsite=ba)——后向族
#         （味/带双翻转，翻转读 −系数），{jr > jc} 条目系数 lower[m]；
#   ch3 = 等时接触修正 λ=0：Σα₃ 按 ba≷bb 镜像（同胞折叠取向随带序，带序在
#         翻转下不变；同味对 asite ≠ bsite 无 local 对角，同一公式）。
# 实时双通道（greater/lesser）：逐味对 (g,l) 列拟合 + Keldysh 关系四块重构，
# Pairing × 4（味旗 = 对的 (ca, cb)）；audit 失败降级原始块直读。
#
# 审计保留：热配对（虚时）/ Keldysh（实时）残差与拟合展开随 audit 返回，并入
# InfluenceFunctional.diagnostics（kernel diagnostics 字段已随 W2 删除）。

"""
    fit_exponentials(alg::ExpExp.ExponentialExpansionAlgorithm, y)
        -> (; weights, nodes, relative_error)

ExactTTI 的函数指数分解（alg_expan 面）：`ExpExp.exponential_expansion` 薄壳 +
GTEMPO 拟合约定对齐（y[d] = Σ wᵢ·qᵢ^(d−1) ⟹ w = α·β、q = β）+ IF 域策略层：
零序列 / `L == 1` 早退；分解不达 `alg.tol` 且 `L ≤ alg.n` 时以单位根
`q_j = 0.8·cis(2πj/L)` 做精确有界插值（合法有限窗展开，非回退到其它 IF 算
法——短窗对铅笔/总体最小二乘欠定时的无损失路）；仍不达则透传 ArgumentError
（消息带最优残差）。分解精度随 alg 配置（缺省 AnchoredMatrixPencil(40, 1e-10,
128)：精确单极点核机器精度；连续谱短核列受小样本极限）。
"""
function fit_exponentials(alg::ExpExp.ExponentialExpansionAlgorithm, y::Vector{<:Number})
    y = ComplexF64.(y)
    all(isfinite, y) || throw(ArgumentError("nonfinite fit samples"))
    (isempty(y) || all(iszero, y)) &&
        return (; weights=ComplexF64[], nodes=ComplexF64[], relative_error=0.0)
    length(y) == 1 && return (; weights=[y[1]], nodes=ComplexF64[0.0], relative_error=0.0)
    L = length(y)
    amps, bases = ExpExp.exponential_expansion(y, alg)
    w = amps .* bases
    recon = [sum(w[i] * bases[i]^(d - 1) for i in eachindex(w)) for d in 1:L]
    err = norm(recon - y) / max(norm(y), eps())
    if !(err <= alg.tol) && L <= alg.n
        # Short finite windows can be underdetermined for a pencil / total LS.
        # An exact bounded exponential interpolation is a valid finite-window
        # expansion, not a fit fallback to another IF algorithm.
        # !(err <= tol)（而非 err > tol）：铅笔在退化列上可返回 NaN 残差，
        # NaN 比较恒 false 会跳过兜底（BCS 反常列 [a,0,−a] 实证）。
        bases = ComplexF64[0.8 * cis(2pi * j / L) for j in 0:L-1]
        B = [q^(d - 1) for d in 1:L, q in bases]
        w = B \ y                      # B 已按 GTEMPO 约定 q^(d−1) 构建 ⟹ 权重 = amps
        recon = [sum(w[i] * bases[i]^(d - 1) for i in eachindex(w)) for d in 1:L]
        err = norm(recon - y) / max(norm(y), eps())
    end
    # 拟合质量契约：err ≤ alg.tol 或抛（含 L > alg.n 无回退资格的情形）
    err <= alg.tol || throw(ArgumentError(
        "exponential expansion unreachable: best residual $err exceeds tol $(alg.tol)"))
    return (; weights=ComplexF64.(w), nodes=ComplexF64.(bases), relative_error=err)
end

# ---------------------------------------------------------------------------
# is_tti：块的时间平移不变性（fit/direct 模式选择的块级门）
# ---------------------------------------------------------------------------

"""
    is_tti(M; rtol=1e-12) -> Bool

块的时间平移不变性（Toeplitz 性）。`CorrelationMatrix` 结构性成立恒 true；
一般矩阵逐对角线检验常数性（相对块幅值的 `rtol` 容差，零块恒 true；矩形块
按交叉块形状支持）。ExactTTI fit 模式的块级门——TTI 失败的味对降级 direct
（原始块 Pairing 直读，精确对任意块）。
"""
is_tti(::CorrelationMatrix) = true
function is_tti(M::AbstractMatrix; rtol::Real=1.0e-12)
    n1, n2 = size(M)
    amax = maximum(abs, M)
    iszero(amax) && return true
    for d in -(n2 - 1):(n1 - 1)
        i0 = max(d + 1, 1)
        j0 = i0 - d
        v0 = M[i0, j0]
        i = i0 + 1
        j = j0 + 1
        while i <= n1 && j <= n2
            isapprox(M[i, j], v0; rtol=rtol, atol=rtol * amax) || return false
            i += 1
            j += 1
        end
    end
    return true
end

# ---------------------------------------------------------------------------
# 站表与列读取助手
# ---------------------------------------------------------------------------

"链站号：格 (b,j) 的 band 带站 = (链块号 − 1)·bands + band"
_station(lattice::AbstractGrassmannLattice, cell::Tuple{Symbol,Int}, band::Int) =
    (cell_block(lattice, cell) - 1) * lattice.bands + band

_kernel_flavor_block(k::HybridizationKernel, b1::Symbol, ba::Int, ca::Bool,
                     b2::Symbol, bb::Int, cb::Bool) =
    k.blocks[(ContourBlockIndex(b1, ba, ca), ContourBlockIndex(b2, bb, cb))]

"拟合列的 lag 重建（ExactTTI 因子锚 = :left，:right 由增长极点守卫拒绝）"
_fitcolumn(f, m::Int) = sum(w * r^(m - 1) for (w, r) in zip(f.weights, f.nodes); init=0.0im)

"核键的味对标签（去重排序的 (ba, ca, bb, cb) 遍历表）"
_flavor_pairs(kernel::HybridizationKernel) =
    sort!(unique!([(key[1].band, key[1].conj, key[2].band, key[2].conj)
                   for key in keys(kernel.blocks)]))

# ---------------------------------------------------------------------------
# 虚轴因子构造：tti_factors(::AbstractImagGrassmannLattice)
# ---------------------------------------------------------------------------

"""
    tti_factors(kernel, lattice::AbstractImagGrassmannLattice; alg) -> (; factors, audit)

热通道因子（τ/τ 核，虚轴 dispatch）：逐味对环族（ExpDecaySumSum 三通道——前向 +
环绕 + 等时接触，2026-11-19 广义化：任意 (ba, ca, bb, cb) 味对，band-diagonal
normal 公式 = ba==bb 特例逐字退化）。**严格 TTI 门**：块非 Toeplitz（`is_tti`
失败）或热配对审计（环反周期闭合）失败 → ArgumentError（ExactTTI = 严格 TTI
算法，非 TTI 核走 PartialIF/PolarIF）；虚时退化拟合节点（|q| ≈ 0）→ Pairing
全矩阵兜底。audit 顶层 `thermal_residual` = 逐对最大值、`imag_expansions` = 逐对
(; pair, mode, thermal_residual, expansion)。
"""
function tti_factors(kernel::HybridizationKernel, lattice::AbstractImagGrassmannLattice;
                     alg=ExactTTI())
    factors = GM.GrassmannFactor[]
    entries = []
    thermal = 0.0
    for (ba, ca, bb, cb) in _flavor_pairs(kernel)
        fs, entry = _tti_imag_pair(kernel, lattice, ba, ca, bb, cb, alg)
        append!(factors, fs)
        push!(entries, entry)
        thermal = max(thermal, entry.thermal_residual)
    end
    return (; factors, audit=(; thermal_residual=thermal, imag_expansions=entries))
end

"单味对虚时环族（fit）；非 TTI / 热配对失败响抛，退化节点走 Pairing 兜底
（第二输出 = audit 条目）"
function _tti_imag_pair(kernel::HybridizationKernel, lattice::AbstractImagGrassmannLattice,
                        ba::Int, ca::Bool, bb::Int, cb::Bool, alg)
    n = lattice.kτ
    B = lattice.bands
    blk = _kernel_flavor_block(kernel, :τ, ba, ca, :τ, bb, cb)
    is_tti(blk) || throw(ArgumentError(
        "ExactTTI requires time-translation-invariant (Toeplitz) kernel blocks: " *
        "imag block for flavor pair (ba=$ba, ca=$ca) → (bb=$bb, cb=$cb) is not TTI; " *
        "use PartialIF/PolarIF for non-TTI kernels"))
    M = Matrix{ComplexF64}(blk)
    lower = ComplexF64[M[j+1, 1] for j in 1:n-1]
    upper = ComplexF64[M[1, j+1] for j in 1:n-1]
    thermal = norm(upper + reverse(lower)) / max(norm(lower), eps())
    thermal <= 1.0e-12 || throw(ArgumentError(
        "imaginary kernel violates thermal pairing (anti-periodic ring closure): " *
        "residual $thermal for flavor pair (ba=$ba, ca=$ca) → (bb=$bb, cb=$cb)"))
    rowst = [_station(lattice, (:τ, j), ba) for j in 1:n]
    colst = [_station(lattice, (:τ, j), bb) for j in 1:n]
    f = fit_exponentials(alg.alg_expan, lower)
    expansion = (; f..., fit_interval=(1, n - 1), convention=:integrated_lag_one)
    # 环核的开放双线性分解（拟合约定 y[d] = Σₖ wₖ qₖ^{d−1}，环反周期
    # upper[m] = −lower[n−m] 已由热配对审计锚定；v7 虚轴翻转：链站序与核胞序
    # **相同** ⟹ 核 s>t（下三角）条目在链上后向、s<t（上三角）条目链前向。
    # 取向读号（由旧降序家族的正确性钉死）：原取向（ca,cb / ba→bb）dense 读
    # = +系数，翻转取向（cb,ca / bb→ba）读 = −系数）：
    #   ch1 前向族（味对原取向，asite=ba/bsite=bb）：{jr < jc} 条目，距离 d
    #     系数 upper[d] = −lower[n−d] = Σ (−w·q^{n−1})·(1/q)^d；
    #   ch2 后向族（味/带双翻转，asite=bb/bsite=ba）：{jr > jc} 条目，距离 m
    #     系数 lower[m]：翻转读 −系数 ⟹ α₂ = −w/q、节点 q；
    #   ch3 等时接触 λ=0：同胞折叠取向随带序（WII 只物化胞内前向对；带序在
    #     翻转下不变 ⟹ 物化条件同旧降序家族）：
    #     ba < bb：ch1 同胞折 +T₁ 物化（T₁ = Σα₁）、ch2 同胞折被排除 ⟹
    #       Σα₃ = M[1,1] − T₁；
    #     ba > bb：ch2 同胞折 −T₂ 物化、ch1 被排除 ⟹ canonical（col@bb 早，
    #       row@ba 晚）= −T₂·(−读) 账 ⟹ Σα₃ = −M[1,1] − T₂（row-first = −canonical）；
    #     ba == bb：local 对角 (T,F)→−T₁、(F,T)→+T₂ ⟹ Σα₃ = M[1,1] − T₁ + T₂
    #     （符号组合经符号暴力搜索 + ET≡PartialIF 表逐位一致定谳，
    #     Δ = 1.6e-15；拟合节点含 ≈0 极点时 1/q 退化——节点门下走 direct）。
    if isempty(f.weights)
        # 零核列：仅等时接触可能非零（λ=0 通道单发）
        contact = M[1, 1]
        iszero(contact) && return GM.GrassmannFactor[],
            (; pair=(ba=ba, ca=ca, bb=bb, cb=cb), mode=:fit,
              thermal_residual=thermal, expansion)
        return [GM.ExpDecaySumSum([GM.ExpDecaySum([contact], [0.0]; conja=ca, conjb=cb,
                                                   cell=B, asite=ba, bsite=bb)])],
               (; pair=(ba=ba, ca=ca, bb=bb, cb=cb), mode=:fit,
                 thermal_residual=thermal, expansion)
    elseif all(abs(q) > 1.0e-8 for q in f.nodes)
        α1 = -(f.weights .* f.nodes .^ (n - 1))
        α2 = -(f.weights ./ f.nodes)
        T1 = sum(α1)
        T2 = sum(α2)
        sums = [GM.ExpDecaySum(α1, 1 ./ f.nodes; conja=ca, conjb=cb, cell=B, asite=ba, bsite=bb),
                GM.ExpDecaySum(α2, f.nodes; conja=cb, conjb=ca, cell=B, asite=bb, bsite=ba)]
        if ba < bb
            ch3 = GM.ExpDecaySum([M[1, 1] - T1], [0.0]; conja=ca, conjb=cb,
                                 cell=B, asite=ba, bsite=bb)
        elseif ba > bb
            ch3 = GM.ExpDecaySum([-M[1, 1] - T2], [0.0]; conja=cb, conjb=ca,
                                 cell=B, asite=bb, bsite=ba)
        else
            ch3 = GM.ExpDecaySum([M[1, 1] - T1 + T2], [0.0]; conja=ca, conjb=cb,
                                 cell=B, asite=ba, bsite=bb)
        end
        iszero(ch3.αs[1]) || push!(sums, ch3)
        return [GM.ExpDecaySumSum(sums)],
               (; pair=(ba=ba, ca=ca, bb=bb, cb=cb), mode=:fit,
                 thermal_residual=thermal, expansion)
    else
        # 退化节点（|q| ≈ 0）：环绕族 1/q 失格 → Pairing 全矩阵兜底（现行语义，
        # 拟合质量问题而非核形状问题）
        f2 = GM.Pairing(M, rowst, colst; conja=ca, conjb=cb)
        return [f2], (; pair=(ba=ba, ca=ca, bb=bb, cb=cb), mode=:direct,
                      thermal_residual=thermal, expansion)
    end
end

# ---------------------------------------------------------------------------
# 实轴因子构造：tti_factors(::AbstractRealGrassmannLattice)
# ---------------------------------------------------------------------------

"逐味对 (g,l) 拟合列（:± 支对块首列 lag 1..n−1；n = 块行数）"
function _real_columns(kernel, ba, ca, bb, cb)
    n = size(_kernel_flavor_block(kernel, :+, ba, ca, :+, bb, cb), 1)
    g = ComplexF64[_kernel_flavor_block(kernel, :+, ba, ca, :+, bb, cb)[i+1, 1] for i in 1:n-1]
    l = ComplexF64[_kernel_flavor_block(kernel, :+, ba, ca, :-, bb, cb)[i+1, 1] for i in 1:n-1]
    return g, l
end

"Keldysh 关系审计（structured Lag 块近似：直读 kernel 的 Toeplitz 性 + 四支
块恒等式；逐味对）。返回残差（相对加权范数）。"
function _keldysh_residual(kernel, ba, ca, bb, cb, n)
    g, l = _real_columns(kernel, ba, ca, bb, cb)
    locations = [(max(d, 0) + 1, max(-d, 0) + 1) for d in -(n - 1):n-1]
    error2 = 0.0
    total2 = 0.0
    blk(a, b) = _kernel_flavor_block(kernel, a, ba, ca, b, bb, cb)
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
    return sqrt(error2) / max(sqrt(total2), eps())
end

"""
    tti_factors(kernel, lattice::AbstractRealGrassmannLattice; alg) -> (; factors, audit)

实时双通道因子（greater/lesser，实轴 dispatch）：逐味对 (g,l) 列拟合 + Keldysh
关系四块重构，Pairing × 4（味旗 = 对的 (ca, cb)——normal (ā,a) 与 BCS 同味
(ā,ā)/(a,a) 统一）。**严格门**：四支块齐备、`is_tti` 通过、Keldysh 审计（逐味
对）通过，任一失败 → ArgumentError（ExactTTI = 严格 TTI 算法，非 TTI 核走
PartialIF/PolarIF）。audit 顶层 `keldysh_residual` = 逐对最大值、`real_expansions`
= 逐对 (; pair, mode, keldysh_residual, greater, lesser)。
"""
function tti_factors(kernel::HybridizationKernel, lattice::AbstractRealGrassmannLattice;
                     alg=ExactTTI())
    factors = GM.GrassmannFactor[]
    entries = []
    kerr = 0.0
    for (ba, ca, bb, cb) in _flavor_pairs(kernel)
        fs, entry = _tti_real_pair(kernel, lattice, ba, ca, bb, cb, alg)
        append!(factors, fs)
        push!(entries, entry)
        kerr = max(kerr, entry.keldysh_residual)
    end
    return (; factors, audit=(; keldysh_residual=kerr, real_expansions=entries))
end

"单味对实时双通道（fit 重构）；缺块 / 非 TTI / Keldysh 审计失败响抛
（第二输出 = audit 条目）"
function _tti_real_pair(kernel::HybridizationKernel, lattice::AbstractRealGrassmannLattice,
                        ba::Int, ca::Bool, bb::Int, cb::Bool, alg)
    n = lattice.k                       # 纯实时分支格点数
    st = br -> [_station(lattice, (br, j), ba) for j in 1:n]
    stcol = br -> [_station(lattice, (br, j), bb) for j in 1:n]
    key(a, b) = (ContourBlockIndex(a, ba, ca), ContourBlockIndex(b, bb, cb))
    missing_keys = [(a, b) for (a, b) in ((:+, :+), (:+, :-), (:-, :+), (:-, :-))
                    if !haskey(kernel.blocks, key(a, b))]
    isempty(missing_keys) || throw(ArgumentError(
        "ExactTTI real kernel is missing branch blocks $missing_keys for flavor pair " *
        "(ba=$ba, ca=$ca) → (bb=$bb, cb=$cb)"))
    all(is_tti(kernel.blocks[key(a, b)]) for (a, b) in ((:+, :+), (:+, :-), (:-, :+), (:-, :-))) ||
        throw(ArgumentError(
            "ExactTTI requires time-translation-invariant (Toeplitz) kernel blocks: " *
            "real blocks for flavor pair (ba=$ba, ca=$ca) → (bb=$bb, cb=$cb) are not TTI; " *
            "use PartialIF/PolarIF for non-TTI kernels"))
    kerr = _keldysh_residual(kernel, ba, ca, bb, cb, n)
    kerr <= 1.0e-12 || throw(ArgumentError(
        "ExactTTI real Keldysh relations fail: $kerr for flavor pair " *
        "(ba=$ba, ca=$ca) → (bb=$bb, cb=$cb)"))
    # fit 模式：(g,l) 列拟合 + Keldysh 关系四块重构（接触直读 blocks[·][1,1]）
    fg = fit_exponentials(alg.alg_expan, first(_real_columns(kernel, ba, ca, bb, cb)))
    fl = fit_exponentials(alg.alg_expan, last(_real_columns(kernel, ba, ca, bb, cb)))
    expansion = e -> (; e..., fit_interval=(1, n - 1), convention=:integrated_lag_one)
    recg(m) = _fitcolumn(fg, m)
    recl(m) = _fitcolumn(fl, m)
    cpp(i, j) = i == j ? _kernel_flavor_block(kernel, :+, ba, ca, :+, bb, cb)[1, 1] :
                i > j ? recg(i - j) : -conj(recl(j - i))
    cmm(i, j) = i == j ? _kernel_flavor_block(kernel, :-, ba, ca, :-, bb, cb)[1, 1] :
                i > j ? -recl(i - j) : conj(recg(j - i))
    cpm(i, j) = i == j ? _kernel_flavor_block(kernel, :+, ba, ca, :-, bb, cb)[1, 1] :
                i > j ? recl(i - j) : conj(recl(j - i))
    cmp(i, j) = i == j ? _kernel_flavor_block(kernel, :-, ba, ca, :+, bb, cb)[1, 1] :
                i > j ? -recg(i - j) : -conj(recg(j - i))
    factors = [GM.Pairing(Matrix{ComplexF64}([cpp(i, j) for i in 1:n, j in 1:n]),
                          st(:+), stcol(:+); conja=ca, conjb=cb),
               GM.Pairing(Matrix{ComplexF64}([cmm(i, j) for i in 1:n, j in 1:n]),
                          st(:-), stcol(:-); conja=ca, conjb=cb),
               GM.Pairing(Matrix{ComplexF64}([cpm(i, j) for i in 1:n, j in 1:n]),
                          st(:+), stcol(:-); conja=ca, conjb=cb),
               GM.Pairing(Matrix{ComplexF64}([cmp(i, j) for i in 1:n, j in 1:n]),
                          st(:-), stcol(:+); conja=ca, conjb=cb)]
    entry = (; pair=(ba=ba, ca=ca, bb=bb, cb=cb), mode=:fit, keldysh_residual=kerr,
             greater=expansion(fg), lesser=expansion(fl))
    return factors, entry
end

# ---------------------------------------------------------------------------
# ExactTTI 装配
# ---------------------------------------------------------------------------

"ExactTTI 组装链：tti_factors 逐因子 expmult!×PolarExp（因子偶元互相超交换
⟹ 应用序无关）。PolarExp 对 Pairing/AnomalousPairing = SVD 秩一精确指数（分
量 Xₖ² = 0 ⟹ e^{Xₖ} = 1 + Xₖ）；对 ExpDecay 族（虚时 ESS）委托 TTIExp =
ExpDecayExp 逐通道物化——**要求 GM ≥ expdecayexp-same-orbital-fix（a7b39bb）**：
同轨道（asite == bsite）触碰块已删，物化对朴素 exp 严格（触碰项与 local ×
过线配对 pfaffian 逐项相消 ⟹ 净零；修正前同轨道偏差 O(α²λ)，GTEMPO ESS 实测
3.2e-4）。截秩语义：rtol = 0（压缩由拟合重构/direct 直读自带精确性，不叠加
截断——混合 cross 通道的截秩口径随其退役删除）。标量型：拟合层因子恒
ComplexF64（fit_exponentials 契约）——实链（虚时格 × 实谱核）在装配入口复化
（GM mult! 不做就地跨型提升 ⟹ ExactTTI 返回态恒 ComplexF64、≠ 输入别名）。

装配路由（plan bcs_exacttti §1.3）：**band 对角 normal 核 → 逐带独立装配，
返回逐带元组**（带间不耦合 ⟹ 各带链键维独立；单带核恒就地单链——golden 锚）；
**含跨带/反常键 → 恒单链**（PartialIF 第九批先例：味对因子单序列覆盖全部
变量，`gmps` 就地改写）。"
function _exacttti_state!(gmps, err::Float64, lattice::AbstractGrassmannLattice,
                          kernel::HybridizationKernel, alg::ExactTTI)
    diag = Dict{Symbol,Any}()
    mult = _mult(alg)
    pol = GM.PolarExp(rtol=0.0, alg_mult=mult)     # 因子级精确（rtol=0）
    TK.scalartype(gmps) <: Real && (gmps = complex(gmps))   # 实链复化（GM complex：
    # 拟合层因子恒 ComplexF64，GM mult! 不做就地跨型提升 ⟹ ExactTTI 返回态恒
    # ComplexF64、≠ 输入别名；已复链 complex 原样返回）
    if !kernel.anomalous && !_has_cross_band_keys(kernel)
        # band 对角 normal：逐带元组路由（golden 锚；单带核恒就地单链）
        single = length(kernel.bands) == 1
        band_states = Any[single ? gmps : copy(gmps) for _ in kernel.bands]
        for (bi, b) in enumerate(kernel.bands)
            (; factors, audit) = tti_factors(_band_kernel(kernel, b), lattice; alg)
            _merge_audit!(diag, audit)
            state = band_states[bi]
            for f in factors
                state, e = GM.expmult!(state, f, pol)
                err += e
            end
            band_states[bi] = state
        end
        return (single ? band_states[1] : tuple(band_states...)), err, diag
    end
    # 跨带/反常核：恒单链（味对因子单序列；就地装配）
    (; factors, audit) = tti_factors(kernel, lattice; alg)
    _merge_audit!(diag, audit)
    for f in factors
        gmps, e = GM.expmult!(gmps, f, pol)
        err += e
    end
    return gmps, err, diag
end

"audit 合并（跨带/逐因子多次 merge：残差标量取 max、逐对条目表 vcat）"
function _merge_audit!(diag::Dict{Symbol,Any}, audit::NamedTuple)
    for p in (:thermal_residual, :keldysh_residual)
        haskey(audit, p) || continue
        diag[p] = haskey(diag, p) ? max(diag[p], getproperty(audit, p)) :
                  getproperty(audit, p)
    end
    for p in (:imag_expansions, :real_expansions)
        haskey(audit, p) || continue
        diag[p] = haskey(diag, p) ? vcat(diag[p], getproperty(audit, p)) :
                  getproperty(audit, p)
    end
    return diag
end
