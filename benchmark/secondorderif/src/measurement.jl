# ============================================================================
# measurement.jl — 观测量（occupation / greater / lesser）
#
# 全部通过"带插入重建 K"计算（U1GTEMPO2 的 scaled_raw_expectation 路线）：
# 分子 = integrate(K_with_insertions, IF)，归一化用无插入的 Z。
#
# 插入约定（ED 校准，validation/unit_weyl_step.jl 逐一验证）：
#   occupation(t_j)   j ≥ 1：n̂ 在 (+, j) 边界 :left 折叠（A_j ← n̂·A_j）
#   greater(t,0)      d 在 (+, t) :left；ddag 在 (+, 0) :left（ρ₀ ← d†ρ₀）
#   lesser(t,0)       d 在 (+, t) :left；ddag 在 (−, 0) :right（ρ₀ ← ρ₀d†）
# 两个时刻均支持（j ≠ 0 时 ddag 插入对应分支的对应边界）。
# ============================================================================

_fock_ops(bands::Int) = _fock_annihilators(bands)

"""
    occupation(sol::Solution, j; branch=default, band=nothing) -> ComplexF64

⟨n̂(t_j)⟩；默认分支为轮廓的第一支（虚时 :τ、实时 :+）。
`band=b` 时插入逐带数算符 n̂_b = d†_b d_b（其余带不折叠）——
U1GTEMPO2 的 occupation 本身就是逐自旋的；两带金标准 per-band = ED 逐带值
（对称双带下 ⟨n₁⟩ = ⟨n₂⟩）。缺省 = 总占据 ⟨n₁+n₂+…⟩。
"""
function occupation(s::Solution, j::Integer; branch::Symbol = first(branches(s.prob.contour)),
                    band::Union{Nothing,Int} = nothing)
    prob = s.prob
    c = prob.contour
    # 逐支边界检查（混合轮廓各支边界数不同）
    jmax = c isa WeylMixedTime ?
           (branch === :τ ? c.Nτ : c.Nt) : ncells(c)
    (0 <= j <= jmax) || throw(BoundsError(0:jmax, j))
    M = c.bands
    d = _fock_ops(M)
    n̂ = if isnothing(band)
        acc = d[1]' * d[1]
        for b in 2:M
            acc += d[b]' * d[b]
        end
        acc
    else
        (1 <= band <= M) || throw(BoundsError(1:M, band))
        d[band]' * d[band]
    end
    # j=0：n̂ 折进 ρ₀（tau：R = n̂·R；实时：cutop(0)·ρ₀）——平衡态下与任意
    # 边界一致（τ 热态循环）；旧实现 j=0 返回空插入（恒 1，占位行为）。
    ins = [WeylInsertion(branch, j, n̂)]
    K = system_factor(c, prob.model; order = prob.order,
                      coupling_norm = coupling_norm(prob.kernel), rho0 = prob.rho0,
                      insertions = ins)
    return expectation(prob.lattice, K, s.mpsIs...; Z = s.Z)
end

"""
    raw_correlation(sol, i, j; kind=:greater) -> ComplexF64

无前因子的两点函数 ⟨d(t_i) ddag(t_j)⟩（greater）/ ⟨ddag(t_j) d(t_i)⟩（lesser）。

**历史备注（2026-11-06）**：M≥2 且 i ≠ j 的跨格奇插入两点曾系统性错误，
根因 = K 装配的 ysv↔site 带映射互换 + 补偿缺双向跨带交叉项
⊕_{b≠b'} c_b·p_{b'}（详见 HISTORY 第十轮），均已修复并经 ED 全网格验证
（`validation/mixed_twopoint_check.jl`）；API 承诺区间 i ≥ j
（i<j 返回轮廓时序积类量，不在本文档语义内）。
"""
function raw_correlation(s::Solution, i::Integer, j::Integer; kind::Symbol = :greater)
    prob = s.prob
    c = prob.contour
    c isa WeylRealTime || c isa WeylMixedTime ||
        throw(ArgumentError("实时/混合轮廓专用；纯虚时用 matsubara"))

    M = c.bands
    d = _fock_ops(M)[1]
    ddag = d'
    ins = WeylInsertion[]
    if kind === :greater
        # ⟨d(t_i) ddag(t_j)⟩（U1GTEMPO2 约定）：d 在 (+, i) 折叠（B 线程）；
        # ddag 在 (+, j)（j=0 折进 ρ₀ 左侧）。
        # 0⁺ 约定（2026-11-05 修正）：同切位（i==j，含 j=0）时 cutop 把
        # 后 push 的矩阵乘在左 = 轮廓更晚。greater 取 d(t⁺) 更晚 ⟹
        # ddag 先 push、d 后 push（等时核 = ⟨d d†⟩ = 1−⟨n⟩；旧 push 序
        # 给出 d†·d = ⟨n⟩，等时两点恒错）。
        push!(ins, WeylInsertion(:+, j, ddag; side = :left))
        push!(ins, WeylInsertion(:+, i, d))
    elseif kind === :lesser
        # ⟨ddag(t_j) d(t_i)⟩（U1GTEMPO2 约定）：d 在 (+, i) 折叠；
        # ddag 在 (−, j) 右侧折叠（j=0 折进 ρ₀ 右侧）
        push!(ins, WeylInsertion(:+, i, d))
        push!(ins, WeylInsertion(:-, j, ddag; side = :right))
    else
        throw(ArgumentError("kind 必须为 :greater 或 :lesser"))
    end
    K = system_factor(c, prob.model; order = prob.order,
                      coupling_norm = coupling_norm(prob.kernel), rho0 = prob.rho0,
                      insertions = ins)
    return expectation(prob.lattice, K, s.mpsIs...; Z = s.Z)
end

"greater(t_i, t_j) = −i⟨d(t_i)ddag(t_j)⟩"
greater(s::Solution, i::Integer, j::Integer = 0) =
    -1.0im * raw_correlation(s, i, j; kind = :greater)

"lesser(t_i, t_j) = +i⟨ddag(t_j)d(t_i)⟩"
lesser(s::Solution, i::Integer, j::Integer = 0) =
    1.0im * raw_correlation(s, i, j; kind = :lesser)

"""
    matsubara(sol, i, j=0) -> ComplexF64

虚时格林函数（Matsubara）：τ 环上 d 在边界 i、ddag 在边界 j 处插入（均 :left，
沿 A 链的算符乘积位置即轮廓时间序）。
"""
function matsubara(s::Solution, i::Integer, j::Integer = 0; force::Bool = false)
    prob = s.prob
    c = prob.contour
    c isa WeylImagTime || throw(ArgumentError("虚时轮廓专用"))
    # 修复边界（2026-11-06 第十轮补遗）：τ 装配的 F 区门控跨带补偿闭合了
    # 等时 + 全部 j≥1 的两点（含 close 格 d）；唯一残留 = j==0（ddag 折入
    # fold 的 R——R 线 s0→tL 的逐项修正未生效，~2e-2），显式抛错
    # （force=true 绕过）。M=1 任意、occupation（偶插入）不受影响。
    if !force && c.bands >= 2 && j == 0
        throw(ArgumentError(
            "matsubara: M≥2 的 j==0（ddag 折入 fold-R 侧）两点当前不可靠" *
            "（fold 元素的跨带修正未生效，HISTORY 第十轮补遗）。" *
            "可用：j≥1 的两点、等时、occupation、M=1 任意。传 force=true 绕过。"))
    end
    M = c.bands
    d = _fock_ops(M)[1]
    K = system_factor(c, prob.model; order = prob.order,
                      coupling_norm = coupling_norm(prob.kernel),
                      insertions = [WeylInsertion(:τ, i, d),
                                    WeylInsertion(:τ, j, d')])
    return expectation(prob.lattice, K, s.mpsIs...; Z = s.Z)
end
