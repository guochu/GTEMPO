# ============================================================================
# kernels.jl — Weyl N×N 浴影响核（单元级）与耦合范数
#
# IF_raw = Z_bath · exp(Φ), Φ = Σ_{iα,jγ} ā_{iα} Δ_{iα,jγ} a_{jγ}，
# Δ = 原浴二点核在单元上的精确双重积分：跨单元为矩形区域，同单元保留两个
# 有序三角区域（推导：benchmark/weyl-secondorder.html）。
#
# 实时核（关键）：QuAPI `Δt` 的原始值是**顶点约定**（_vs_rt_acc 映射下的
# 节点读数），不是 Weyl 全宽单元的积分值。这里对单模 Dirac-δ 浴直接从
# **精确浴迹多项式**提取单元级核：单元演化算符
#     E_j^± = e^{∓iδt·H_b}·[I + G ± G²/2]，G = ∓iδt(ā F' + F'† a)
# （F' 为海森堡共轭后的浴湮灭矩阵，G³ = 0 使级数精确闭合），按轮廓时序
# (+ 支 1..N，再 − 支 N..1) 乘积并取浴迹，其线性单项式系数恰为
# exp(Φ) 的 Φ 核。交叉验证：τ 核直接用 QuAPI `Δτ`（本身按单元给出，
# 已与 ED 验证一致）。
#
# 块值沿用 QuAPI 的存储约定（跨支块内嵌 Keldysh 相对负号），因此可以逐块
# 交给 GTEMPO 的 partialif_hybrid/hybriddynamics 消费（它在跨支读数处补
# 负号）；同时除以浴迹常数（QuAPI 核的归一化约定，归一化观测量下抵消）。
# ============================================================================
using GTEMPO: Δt, Δτ, Δm, fermionicbath
using ImpurityModelBase: DiscreteSpectrum, spectrumvalues

"""
    HybridizationKernel

单元级杂化核：`blocks[(b1,b2)]` 的 (i,j) 元 = ā_i^{b1} 与 a_j^{b2} 之间的
IF 指数系数（含 QuAPI 存储约定的跨支 Keldysh 负号），单元下标 1..N。
`coupling_norm = Σ_k |V_k|²` 是 Weyl 标量补偿所需的单带总耦合权重
（U1GTEMPO2：order=2 必须提供 coupling_norm；有限采样核无法可靠推断）。
"""
struct HybridizationKernel
    contour::AbstractWeylContour
    blocks::Dict{Tuple{Symbol,Symbol},Matrix{ComplexF64}}
    coupling_norm::Float64
end

coupling_norm(k::HybridizationKernel) = k.coupling_norm

# DiracDelta 谱的积分权重 α（精确）；一般有界谱数值积分。单模基准只需前者。
α_weight(spec::DiracDelta) = spec.α
# 离散（多 δ）谱：总耦合 = Σ αₙ（quadgk 对 ±Inf 上的离散谱无意义）
α_weight(spec::DiscreteSpectrum) = sum(spectrumvalues(spec))
α_weight(spec) = quadgk(spec, lowerbound(spec), upperbound(spec))[1]
α_weight(bath::AbstractBath) = α_weight(bath.spectrum)

"""
    hybridization(bath, c::AbstractWeylContour; band=1) -> HybridizationKernel

单元级浴核：
- 虚时：QuAPI `Δτ`（本身按 N 个单元给出，一般有界谱浴，ED 已验证）；
- 实时：按谱分派——
  * Dirac-δ 单模浴：**精确浴迹提取**（单元演化算符解析精确，见下）；
  * 一般有界谱浴：**QuAPI `Δt` 谱路径**——其块值即单元级全宽积分
    （`_gⱼₖ` 矩形 / `_gⱼⱼₖₖ` 三角，2026-11-05 再定性实证），Weyl 单元核 =
    `[1:N, 1:N]` 直接截取，与精确路径在 Dirac-δ 上逐位一致（≤1e-12，
    validation/kernel_align_probe.jl）。两条路径互为交叉验证。
"""
function hybridization(bath, c::AbstractWeylContour; band::Int = 1)
    band == 1 || throw(ArgumentError("单元级核按带构造；多带用 influence_functional 复制"))
    blocks = Dict{Tuple{Symbol,Symbol},Matrix{ComplexF64}}()
    if c isa WeylImagTime
        corr = Δτ(bath; N = c.Nτ, δτ = c.dτ)
        D = corr.data                             # ImagCorrelationFunction 内的 CorrelationMatrix，N 个单元
        blocks[(:τ, :τ)] = Matrix{ComplexF64}(D[1:c.Nτ, 1:c.Nτ])
    elseif c isa WeylMixedTime
        blocks = mixed_kernel_quapi(bath, c)
    elseif bath.spectrum isa DiracDelta
        blocks = realtime_kernel_trace(bath, c)
    else
        blocks = realtime_kernel_quapi(bath, c)
    end
    return HybridizationKernel(c, blocks, Float64(α_weight(bath)))
end

"""
    mixed_kernel_quapi(bath, c::WeylMixedTime) -> Dict{Tuple{Symbol,Symbol},Matrix{ComplexF64}}

混合（Kadanoff–Baym L 形）轮廓的单元级核：QuAPI `Δm` 谱积分（任意有界谱），
九个分支块（(+,+)、(+,−)、(+,τ)、(−,+)、(−,−)、(−,τ)、(τ,+)、(τ,−)、(τ,τ)），
实时索引截取 [1:Nt]、虚时索引 [1:Nτ]——与实时核的单元截取同一语义。
"""
function mixed_kernel_quapi(bath, c::WeylMixedTime)
    qm = Δm(bath; Nτ = c.Nτ, t = c.Nt * c.dt, Nt = c.Nt, δτ = c.dτ)
    blocks = Dict{Tuple{Symbol,Symbol},Matrix{ComplexF64}}()
    for b1 in (:+, :-, :τ), b2 in (:+, :-, :τ)
        Mq = GTEMPO.QuAPI.branch(qm; b1 = b1, b2 = b2)
        ri(c2) = c2 === :τ ? (1:c.Nτ) : (1:c.Nt)
        blocks[(b1, b2)] = Matrix{ComplexF64}(Mq[ri(b1), ri(b2)])
    end
    return blocks
end


"""
    realtime_kernel_quapi(bath, c::WeylRealTime) -> Dict{Tuple{Symbol,Symbol},Matrix{ComplexF64}}

一般有界谱浴的单元级实时核：QuAPI `Δt` 的谱积分（任意 `AbstractBoundedFunction`
谱，内部 quadgk）按单元索引 `[1:N, 1:N]` 截取。块值语义与
`realtime_kernel_trace`（Dirac-δ 精确路径）完全一致——同一单元积分核的两条
独立计算路径，Dirac-δ 上互验。
"""
function realtime_kernel_quapi(bath, c::WeylRealTime)
    qt = Δt(bath; N = c.Nt, t = c.Nt * c.dt)
    blocks = Dict{Tuple{Symbol,Symbol},Matrix{ComplexF64}}()
    for (b1, b2) in ((:+, :+), (:+, :-), (:-, :+), (:-, :-))
        Mq = GTEMPO.QuAPI.branch(qt; b1 = b1, b2 = b2)
        size(Mq, 1) == size(Mq, 2) == c.Nt + 1 ||
            throw(ArgumentError("QuAPI Δt 块尺寸 $(size(Mq)) 与单元数不符（应为 $(c.Nt+1)² 节点阵）"))
        blocks[(b1, b2)] = Matrix{ComplexF64}(Mq[1:c.Nt, 1:c.Nt])
    end
    return blocks
end

# ---------------------------------------------------------------------------
# 实时核的精确浴迹提取（单模 Dirac-δ 浴）
# ---------------------------------------------------------------------------
# 浴 Fock 2×2；F = √α·b（行=末态，列=初态：F[empty, occ]）
function _singlemode_ops(α::Real, ω::Real)
    F = zeros(ComplexF64, 2, 2)
    F[1, 2] = sqrt(α)                 # b：|occ⟩ → |empty⟩
    return F, F', ω
end

"单元演化算符的精确多项式（局部位 1=a, 2=ā；Dict: mask → 2×2 浴矩阵）"
function _cell_E(dt::Real, sgn::Int, F, Fdg, ω::Real)
    hb = ComplexF64[0 0; 0 ω]
    lam = sgn * im * dt
    Ebath = exp(-lam * hb)   # e^{-λH_b}（U1 约定：+支 λ=+i dt ⟹ 正向演化）
    # E = e^{-λH_b}·Texp[−λ∫₀¹ds V(s)]，V(s) = ā·F(s) + F†(s)·a，
    # F(s) = F·e^{−sλω}（浴侧共轭积分在 s 内完成）。Grassmann 每支每格仅
    # 2 个槽 ⟹ 三阶及以上单项式恒零，展开精确闭合。
    # 一阶：−λ∫F(s)ds = −F·(1−e^{−λω})/(λω)（固定共轭 F(1) 会差 O(λω)，
    #   经线性×线性积污染所有跨格核对——2026-11-04 定位）
    # 二次有序三角（2026-11-04 IF 系数交叉核对定位）：宗量 = λω（U1 的 e·li，
    #   ∫₀¹du(1−u)e^{±λωu}），此前误用 λω·dt（虚部差 dt 倍 → 实时残差随时间
    #   累积的根因）；f/g 通道与 tri(±) 的配对经浴迹/U1 双重标定：
    #   aā 系数 = f·tri(+λω) − g·tri(−λω) 经 flip 链 = U1 mode_block 对角。
    lw = lam * ω
    if abs(lw) < 1e-12
        Fbar = F; Fbardg = Fdg
        tm = dt^2 / 2; tp = dt^2 / 2
    else
        Fbar = F * (1 - exp(-lw)) / lw
        Fbardg = Fdg * (exp(lw) - 1) / lw
        tm = dt^2 * (exp(lw) - 1 - lw) / lw^2       # dt²·tri(+λω)（f 通道）
        tp = dt^2 * (exp(-lw) - 1 + lw) / lw^2      # dt²·tri(−λω)（g 通道）
    end
    pol = Dict{Int,Matrix{ComplexF64}}(
        0 => Ebath,
        1 => Ebath * (-lam * Fbardg),
        2 => Ebath * (-lam * Fbar),
        3 => Ebath * (Fdg * F .* tm - F * Fdg .* tp),
    )
    return pol
end

"Grassmann 多项式乘法（规范序 = 位升序；存储符号 = 合并换位宇称）"
function _gmul(X::Dict{Int,Matrix{ComplexF64}}, Y::Dict{Int,Matrix{ComplexF64}}, n::Int)
    out = Dict{Int,Matrix{ComplexF64}}()
    for (mx, Mx) in X, (my, My) in Y
        (mx & my) != 0 && continue
        sgn = 1
        for i in 1:n
            ((mx >> (i - 1)) & 1) == 1 || continue
            for j2 in (i+1):n
                ((my >> (j2 - 1)) & 1) == 1 && (sgn = -sgn)
            end
        end
        m = mx | my
        out[m] = get(out, m, 0.0im) .+ sgn .* (Mx * My)
    end
    return out
end

"""
    realtime_kernel_trace(bath, c::WeylRealTime; storage=:u1) -> Dict{...}

单模 Dirac-δ 浴的单元级实时核：精确浴迹多项式的线性单项式系数。
浴迹：Tr_b[E₁⁻···E_N⁻ · E_N⁺···E₁⁺]（右乘在先 = 轮廓时序）；
核提取：线性项 ā_i^{b₁} a_j^{b₂} 的系数（规范序修正后），并除以浴迹常数。

`storage`（2026-11-04 第二轮定标）：
- `:u1`（默认）= U1GTEMPO2 mode_block 原始存储（与 U1 核逐位一致）——
  保持单模浴核的秩 1/Toeplitz 结构（B′ 三角翻号会在 Pfaffian 组合层面
  破坏它：IF 四阶扇区多出 1429 个 U1 精确为零的实质单项式）；
- `:B′` = 旧补偿存储，仅供对照。
"""
function realtime_kernel_trace(bath, c::WeylRealTime; storage::Symbol = :u1)
    spec = bath.spectrum
    spec isa DiracDelta || throw(ArgumentError("实时核精确提取仅支持单模 Dirac-δ 浴"))
    ω, α = spec.ω, spec.α
    F, Fdg, ωv = _singlemode_ops(α, ω)
    ω == ωv || throw(ArgumentError())
    n = c.Nt
    hb = ComplexF64[0 0; 0 ω]
    # 浴初始态：热 ρ_b = e^{−βh_b}/Z_b（与 τ 核 Δτ 的 β 权重同一约定）。
    # 无权重闭合迹 = β=0 浴，与 ED 参考的热初态不符（2026-11-04 定位的
    # α 线性反向漂移根因）。
    β = bath.β

    # cell 内位：bit0 = ā, bit3 = a（正支）；负支：bit1 = ā, bit2 = a
    function gpol_branch(j::Int, branch::Symbol)
        pol = _cell_E(c.dt, branch === :+ ? 1 : -1, F, Fdg, ω)
        off = 4 * (j - 1)
        bā, ba = branch === :+ ? (0, 3) : (1, 2)
        out = Dict{Int,Matrix{ComplexF64}}()
        for (m, M) in pol
            g = 0
            (m & 1) != 0 && (g |= 1 << (off + ba))
            (m & 2) != 0 && (g |= 1 << (off + bā))
            out[g] = M
        end
        return out
    end

    ρb = exp(-β * hb)
    Zb = tr(ρb)                                               # 归一化 = tr e^{−βh_b}
    # cell-内变量位（全局）：+支：ā = 4(j−1)+1，a = 4(j−1)+4；−支：ā = +2，a = +3
    function bits(b::Symbol, j::Int, conj::Bool)
        off = 4 * (j - 1)
        bā, ba = b === :+ ? (1, 4) : (2, 3)
        return conj ? off + bā : off + ba
    end
    # 单因子项矩阵：E^(0)（真空）、E^(ā)、E^(a)、E^(aā)（同单元同支）。
    # E = e^{−λH_b}·Texp[−λ∫₀¹ds V(s)] 的精确展开：线性 −λ∫F(s)ds、
    # 二次 T 序 (F†F−FF†)·S（见 _cell_E 注）。规范序 a 先 ā 后。
    # 精确单元多项式（U1 mode_block 同构）：线性 = −λ∫₀¹F(s)ds；
    # 二次 T 序 = ∫_{s>s'}V(s)V(s')（普通算符时序，无费米子反对称化）：
    #   ā·a 系数 = FF†·tri(−z) − F†F·tri(z)，tri(z) = (e^z−1−z)/z²（U1 tri）。
    function factor(fj::Int, fbr::Symbol, kind::Symbol)
        lam = fbr === :+ ? im * c.dt : -im * c.dt
        Ebath = exp(-lam * hb)   # e^{-λH_b}（U1 约定：+支 λ=+i dt ⟹ 正向演化）
        lw = lam * ω
        if abs(lw) < 1e-12
            Fbar = F; Fbardg = Fdg
            tm, tp = c.dt^2 / 2, c.dt^2 / 2
        else
            Fbar = F * (1 - exp(-lw)) / lw
            Fbardg = Fdg * (exp(lw) - 1) / lw
            # 有序三角（2026-11-04 IF 系数交叉核对定位）：宗量 λω = e·li
            #（U1 mode_block 同构；此前误用 λω·dt，虚部差 dt 倍且 f/g
            # 通道配对相反——实时残差随时间累积的根因）
            tm = c.dt^2 * (exp(lw) - 1 - lw) / lw^2     # tri(+λω)，f 通道
            tp = c.dt^2 * (exp(-lw) - 1 + lw) / lw^2    # tri(−λω)，g 通道
        end
        if kind === :vac
            return Ebath
        elseif kind === :abar
            return Ebath * (-lam * Fbar)
        elseif kind === :a
            return Ebath * (-lam * Fbardg)
        else # :aab 同单元同支的 ā·a 二次项（a·ā 存储序）
            return Ebath * ((Fdg * F .* tm - F * Fdg .* tp))
        end
    end
    hb = ComplexF64[0 0; 0 ω]

    blocks = Dict{Tuple{Symbol,Symbol},Matrix{ComplexF64}}()
    for b1 in (:+, :-), b2 in (:+, :-)
        M = zeros(ComplexF64, n, n)
        for i in 1:n, j in 1:n
            # 轮廓序（左→右）：− 支 1..N，再 + 支 N..1（右端最先作用）
            # P = E₁⁻⋯E_N⁻ · E_N⁺⋯E₁⁺：左 = E₁⁻
            factors = Tuple{Int,Symbol}[]
            for jj in 1:n
                push!(factors, (jj, :-))
            end
            for jj in n:-1:1
                push!(factors, (jj, :+))
            end
            acc = Matrix{ComplexF64}(ρb)                       # 热初态权重入迹
            for (fj, fbr) in factors
                if fj == i && fbr == b1
                    term = (fj == j && fbr == b2) ? factor(fj, fbr, :aab) : factor(fj, fbr, :abar)
                elseif fj == j && fbr == b2
                    term = factor(fj, fbr, :a)
                else
                    term = factor(fj, fbr, :vac)
                end
                acc = acc * term
            end
            coef = tr(acc) / Zb
            bā_ = bits(b1, i, true)
            ba_ = bits(b2, j, false)
            Δ = ba_ > bā_ ? coef : -coef
            Δ = -Δ
            # 基础 B′ 三角补偿（2026-11-04 定标）
            if b1 === :+ && b2 === :+
                i <= j && (Δ = -Δ)
            elseif b1 === :- && b2 === :-
                i >= j && (Δ = -Δ)
            elseif b1 === :-
                Δ = -Δ
            end
            if storage === :u1
                # 再复合自逆映射 S（IF 系数交叉核对实证：同支 i<j 翻、
                # 跨支 i>j 翻）⟹ 与 U1GTEMPO2 mode_block 核逐位一致
                if b1 === b2
                    i < j && (Δ = -Δ)
                else
                    i > j && (Δ = -Δ)
                end
            end
            M[i, j] = abs(Δ) < 1e-13 ? zero(Δ) : Δ
        end
        blocks[(b1, b2)] = M
    end
    return blocks
end

# ---------------------------------------------------------------------------
# GTEMPO hybriddynamics 的核视图：partialif_hybrid 通过
#   corr[conj(row ā), conj(col a)] 读取（"a-ā 约定 x[a,ā]=Δ"），
# 普通费米浴：a-ā → +Δ，ā-a → −Δ，同 conj → 0；跨支负号已内嵌在块值里，
# 由 hybriddynamics 的读数规则补回物理块 —— 与 QuAPI RealCorrelationFunction
# 的 getindex 语义完全一致，只是下标直接用单元 1..N。
# ---------------------------------------------------------------------------
struct KernelCorrView <: QuAPI.AbstractCorrelationFunction
    kernel::HybridizationKernel
    band::Int
end

function Base.getindex(v::KernelCorrView, i::GrassmannContourIndex, j::GrassmannContourIndex)
    (i.band == v.band && j.band == v.band) || return zero(ComplexF64)
    (i.conj != j.conj) || return zero(ComplexF64)
    Δ = v.kernel.blocks[(branch(i), branch(j))]
    # 混合轮廓的端点顶点（rt Nt+1）不进 Weyl 单元：越界读数 = 零耦合
    (1 <= i.j <= size(Δ, 1) && 1 <= j.j <= size(Δ, 2)) || return zero(ComplexF64)
    val = Δ[i.j, j.j]
    # 跨支读数符号约定（与 U1GTEMPO2 核值逐位核对，2026-11-05）：
    # ① 实时 (:+,:)：partialif_hybrid 消费端内嵌跨支负号（keldysh gate），
    #   视图预翻与之相消 ⟹ 净读数 = 核原始值；
    # ② 混合轮廓：消费端 gate 关闭；QuAPI Δm 的 τ↔rt 交叉块相对物理核差
    #   一个整体负号（ττ/±± 块无差——U1 对照实证）⟹ 视图仅翻 τ↔rt 读数；
    #   +↔− 读数沿用 ① 的双翻相消。
    if v.kernel.contour isa WeylRealTime
        # 实时：消费端 keldysh gate 内嵌跨支负号，视图预翻相消 ⟹ 净 = 原始值
        branch(i) !== branch(j) && (val = -val)
    elseif v.kernel.contour isa WeylMixedTime
        # 混合：消费端 gate 关闭。Δm 的 τ↔rt 交叉块相对物理核差一个整体负号
        # （ττ / ±± 块无差——U1 对照实证）⟹ 仅翻 τ↔rt；±↔∓ 块直读原始值
        # （单翻会给出负的物理值——Anderson mixed 斜率 0.64→2.0 实证）。
        (branch(i) === :τ) !== (branch(j) === :τ) && (val = -val)
    end
    return i.conj ? -val : val
end

Z2Tensors.scalartype(::Type{KernelCorrView}) = ComplexF64

# KernelCorrView 的单元级块访问（扩展 QuAPI/GTEMPO 的 branch；勿新建同名函数，
# 否则会遮蔽 branch(::GrassmannContourIndex)）
import GTEMPO: branch
branch(x::KernelCorrView; b1::Symbol, b2::Symbol) = x.kernel.blocks[(b1, b2)]

