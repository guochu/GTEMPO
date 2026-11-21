# correlationfunction.jl：浴关联 → 单元级杂化核（Abstract/Frozen/Exact
# HybridizationKernel 类型层级，exactif_port_plan E6）
#
# plan/if_expmult_backend_plan.md W2（2026-11-17 用户五裁决）：旧
# HybridizationKernel{L,C}（lattice + QuAPI corr 本体 + (Symbol,Symbol) 视图 +
# diagnostics）重定义为纯数据核——blocks 以 ContourBlockIndex 对键入、值 =
# Union{Matrix,CorrelationMatrix}（QuAPI Toeplitz O(N) 存储，|i−j| 查表）、
# coupling_norm::Float64、bands::Vector{Int}。normal 与 BCS bath 统一同一类型
# （双构造方法 + changebands + Base.:+ 多浴耦合）。
#
# 2026-11-19：ContourBlockIndex 增 conj::Bool 腿味（true = ā、false = a）——
# 核键 = (branch, band, conj) 完备腿信息构造期写定，装配端（fullfactor/
# hybriddynamics）按键直读、bcs 味旗标退役。Kernel 从此对任意 bath 自完备。
#
# 2026-11-21（exactif_port_plan Phase A，E6 用户裁决）：HybridizationKernel
# 类型层级化——AbstractHybridizationKernel + FrozenHybridizationKernel（现行
# QuAPI 冻结核改名）/ ExactHybridizationKernel（IFD 严格核对接，带号
# normdelta 单一存储）；
# `HybridizationKernel` 名退役（统一类型回收与否 = Phase E 合并接口裁决）；
# 常数消费 = 类型分派 `_constant_lscale`（sysdynamics/fockmpo.jl）。
# 存储约定（构造期符号收编——读点直读）：
#   blocks[(CBI(b1,band,ca), CBI(b2,band,cb))][i, j] = IF 指数系数（i, j =
#   分支局部格点下标；行腿味 ca、列腿味 cb——normal 键 = (ā, a) 即
#   (true, false)，BCS 反常键 = (ā, ā) / (a, a)）。
#   - 纯虚时/纯实时：QuAPI branch 块原值（实时跨支块的 Keldysh 相对号内嵌于
#     块值 = 指数系数的一部分；U1GTEMPO2 mode_block 口径逐位一致）；
#   - 混合轮廓：τ↔rt 交叉块构造期整体取负（QuAPI Δm 的 τ↔rt 块相对物理核
#     差一个负号——旧 _NegView 惰性预翻的物化版本）；
#   - BCS 三轮廓统一 = 味块 branch 视图原值（raw）。2026-11-19 v3-① 裁决：
#     旧「实时访问器物化（顶点移位内嵌，bcsweyl 锚定口径）」整体移除——移位
#     与标准网格装配错位（Δ=0 对照 normal 实证：首行接触缺失/条目移位一格），
#     Weyl2×BCS/ExactTTI fit 全消费端统一吃 raw；QuAPI 包内 _vs_rt_acc 语义
#     不动（GTEMPO 不再消费任何移位读取）。T6 探针进一步否决了计划预留的
#     「跨支构造期取负」：bcs_Δt 支视图自带与 normal 相同的指数系数 Keldysh
#     号（Δ=0 四支对逐位 0.0），取反即错。
#
# 电流公式注记：η⁺⁺ = blocks 直读；η⁺⁻ = −blocks[(:+,:-)]——电流算符的
# η 系数 = 无号物理块（backup/GTEMPO currents.jl 同款），与指数系数差
# Keldysh 基的一个相对号，公式层常数（2026-11-18 裁决保留，见 inserted_k.jl）。

# ---------------------------------------------------------------------------
# ContourBlockIndex + 核类型层级（Abstract/Frozen/Exact——exactif_port_plan E6）
# ---------------------------------------------------------------------------

"""
    ContourBlockIndex(branch::Symbol, band::Int, conj::Bool)

核块的轮廓索引：(branch ∈ (:τ, :+, :-), band = 带号, conj = 腿味)。核块键 =
`(ContourBlockIndex, ContourBlockIndex)`（行腿 × 列腿）。`conj = true` ⟺ 腿为
ā（GM 变量基的 η̄，变量号偶）、`false` ⟺ a（η，变量号奇）。normal 块键 =
`(ā, a)` 即 `(true, false)`；BCS 反常块（ā↑ā↓ / a↓a↑）两腿同味。腿味由
correlationfunction 构造期写定，装配端按键直读（kernel 自完备，无调用端
味旗标）。
"""
struct ContourBlockIndex
    branch::Symbol     # :τ / :+ / :-
    band::Int
    conj::Bool         # true = ā（η̄）/ false = a（η）
end

"""
    AbstractHybridizationKernel

单元级杂化核的抽象基类（exactif_port_plan E6，2026-11-20 用户裁决：类型层级
取代 scheme 字段）。两具体型共享同一数据面：`blocks`（ContourBlockIndex 对
键）、`bands::Vector{Int}`、`anomalous::Bool`（构造器从键派生）；常数面按型
分派——`FrozenHybridizationKernel.coupling_norm`（Σα，Weyl 标量）/
`ExactHybridizationKernel.normdelta`（带号 Z 常数）。下游
（`correlationfunction`/`fullfactor`/`hybriddynamics`/PartialIF 系/
measurements）一律消费抽象型。
"""
abstract type AbstractHybridizationKernel end

"""
    FrozenHybridizationKernel{T<:Number} <: AbstractHybridizationKernel

单元级杂化核（canonical 指数系数约定，见文件头存储约定）——QuAPI 冻结核
离散化（原 `HybridizationKernel` 改名；`HybridizationKernel` 名退役，统一
类型回收与否 = exactif_port_plan Phase E 合并接口裁决点）：

- `blocks::Dict{Tuple{ContourBlockIndex,ContourBlockIndex}, Union{Matrix{T},CorrelationMatrix{T}}}`
  ——块键 = (行腿, 列腿) 的 ContourBlockIndex 对；`CorrelationMatrix` = QuAPI
  Toeplitz 向量存储（O(N)、|i−j| 查表、等时 = ηⱼₖ[1]+ηₖⱼ[1] 接触合并），
  `Matrix` = 已物化稠密（混合 τ↔rt 折叠块 / BCS 实时移位块）。两表示的
  getindex 同值，物化仅在构造端需要时发生。
- `coupling_norm::Float64` = 全带总耦合 Σ_b α_b（Weyl 标量补偿用；order=2
  必需，单带 = coupling_norm(bath)，全带 = bands·倍，BCS 轨道 = 2×）。
- `bands::Vector{Int}` = 本核覆盖的带号有序表（单带 → [b]；全带 → [1..B]；
  BCS 轨道 → [band1, band2]）。
- `anomalous::Bool` = 是否含异味键（a·a / ā·ā——不满足荷守恒的块，BCS 反常
  通道）。**构造器从核键派生**（`any(key -> key[1].conj == key[2].conj)`）——
  单一事实来源，`changebands`/`Base.:+` 等重装路径自动一致（normal + 反常
  相加 = 反常）。

**canonical 化（构造期）**：跨味键一律等效转换为 (ā, a) 取向——(a, ā) 键在
构造期换腿 + 转置 + 换序号取负（`v·a_i·ā_j = (−vᵀ)·ā_j·a_i`）；同味键
（(ā,ā)/(a,a)）无 normal-conj 等效形、原样保留。核存储恒 canonical，消费端
（`fullfactor` 的 Pairing 路由、ExactTTI 味对分组）按键直读、无取向检查。

`scalartype` 仅 Type 上定义（TensorKit 泛型回退 scalartype(x) =
scalartype(typeof(x)) 自动覆盖实例）；T 由构造端 = `scalartype(corr)` 动态定
（实谱虚时核 Float64、复谱/实时/混合核 ComplexF64——见下方 QuAPI 容器
scalartype 表）。
"""
struct FrozenHybridizationKernel{T<:Number} <: AbstractHybridizationKernel
    blocks::Dict{Tuple{ContourBlockIndex,ContourBlockIndex},Union{Matrix{T},CorrelationMatrix{T}}}
    coupling_norm::Float64
    bands::Vector{Int}
    anomalous::Bool     # 异味键存在（构造器从键派生）

    function FrozenHybridizationKernel(blocks::Dict{Tuple{ContourBlockIndex,ContourBlockIndex},
                                                   Union{Matrix{T},CorrelationMatrix{T}}},
                                       coupling_norm::Real, bands::Vector{Int}) where {T<:Number}
        isfinite(coupling_norm) && coupling_norm >= 0 ||
            throw(ArgumentError("invalid coupling norm $coupling_norm"))
        isempty(bands) && throw(ArgumentError("empty band list"))
        allunique(bands) && issorted(bands) || throw(ArgumentError(
            "bands must be a sorted list of unique band indices, got $bands"))
        all(b -> b >= 1, bands) || throw(ArgumentError("band indices must be ≥ 1"))
        canon, anomalous = _canon_kernel_blocks(blocks, bands, T)
        new{T}(canon, Float64(coupling_norm), bands, anomalous)
    end
end

"""
    ExactHybridizationKernel{T<:Number} <: AbstractHybridizationKernel

严格核（IFD ExactKernel 族）的核容器（exactif_port_plan E4′/E6）：blocks 面
与 Frozen 完全同构（imag + Exact = τ 环 Toeplitz 点值块，CorrelationMatrix；
real/mixed + Exact = 精确消元核，非平移不变，dense `Matrix`），整体 Z 常数 =
`normdelta::Float64`（**带号**，单一存储，不折进 blocks）：

- Δτ 无耦合极限 `(−1)^{NK}·Z_bath`、Δt = `Π n̄ₖ·det(BB)`、Δm =
  `(−1)^{n(K+1)}·det(BB)/Z_bath`（号约定锚 = IFD 金标准
  `test/exact/fermionic_{imagtime,realtime,mixedtime}.jl`）；
- 构造守卫 = `isfinite` + 实值尾取实（IFD 元组第二元随容器标量型带 ~1e-15
  虚部尾，det 实值）；
- 下游常数消费 = 类型分派（`_constant_lscale`：log|normdelta| 折进链幅
  logscale、号经 `lmul!` 折进首站）⟹ Z 与分子同携带，比值观测量严格相消、
  绝对 Z 与真空形式 n̄_vac 需要（IFD normdelta 语义定谳一致）。
"""
struct ExactHybridizationKernel{T<:Number} <: AbstractHybridizationKernel
    blocks::Dict{Tuple{ContourBlockIndex,ContourBlockIndex},Union{Matrix{T},CorrelationMatrix{T}}}
    normdelta::Float64
    bands::Vector{Int}
    anomalous::Bool     # 异味键存在（构造器从键派生）

    function ExactHybridizationKernel(blocks::Dict{Tuple{ContourBlockIndex,ContourBlockIndex},
                                                  Union{Matrix{T},CorrelationMatrix{T}}},
                                      normdelta::Number, bands::Vector{Int}) where {T<:Number}
        isempty(bands) && throw(ArgumentError("empty band list"))
        allunique(bands) && issorted(bands) || throw(ArgumentError(
            "bands must be a sorted list of unique band indices, got $bands"))
        all(b -> b >= 1, bands) || throw(ArgumentError("band indices must be ≥ 1"))
        canon, anomalous = _canon_kernel_blocks(blocks, bands, T)
        new{T}(canon, _real_normdelta(normdelta), bands, anomalous)
    end
end

"normdelta 实值化守卫：isfinite + 虚部尾 ≤ 1e-12·max(1,|Re|) 取实（det 实值
契约； genuinely 复值响抛）"
function _real_normdelta(nd::Number)
    abs(imag(nd)) <= 1.0e-12 * max(1.0, abs(real(nd))) || throw(ArgumentError(
        "normdelta must be real (Z constant contract), got imag tail $(imag(nd))"))
    v = Float64(real(nd))
    isfinite(v) || throw(ArgumentError("invalid normdelta $nd"))
    return v
end

"核 blocks 的 canonical 化 + 键/带一致性守卫（两具体型构造器共用）：
跨味键 (a, ā) 等效转换为 (ā, a)——换腿 + 转置 + 换序号取负（费米换序号 −1）；
同味键原样保留。碰撞（同一 canonical 键多来源）= 加法合并（双线性项求和的
代数事实；well-formed 核无碰撞）。返回 (canon, anomalous)。"
function _canon_kernel_blocks(blocks::Dict{Tuple{ContourBlockIndex,ContourBlockIndex},
                                          Union{Matrix{T},CorrelationMatrix{T}}},
                              bands::Vector{Int}, ::Type{T}) where {T<:Number}
    canon = Dict{Tuple{ContourBlockIndex,ContourBlockIndex},
                 Union{Matrix{T},CorrelationMatrix{T}}}()
    for (key, value) in blocks
        if !key[1].conj && key[2].conj
            newkey = (ContourBlockIndex(key[2].branch, key[2].band, true),
                      ContourBlockIndex(key[1].branch, key[1].band, false))
            # CM = IFD transpose（两 Toeplitz 向量互换）+ 取负，容器原样；
            # 稠密块物化（-transpose 为惰性 Transpose 包装器，须落 Matrix）
            v = value isa CorrelationMatrix ? -transpose(value) :
                Matrix(-transpose(value))
        else
            newkey = key
            v = value
        end
        canon[newkey] = haskey(canon, newkey) ? canon[newkey] + v : v
    end
    for (key, value) in canon
        key[1].band in bands && key[2].band in bands || throw(ArgumentError(
            "block key $(key) references bands outside the kernel band list $bands"))
        all(isfinite, value) || throw(ArgumentError(
            "nonfinite kernel block $(key[1].branch)/$(key[2].branch) band $(key[1].band)"))
    end
    return canon, any(key -> key[1].conj == key[2].conj, keys(canon))
end

TK.scalartype(::Type{FrozenHybridizationKernel{T}}) where {T} = T
TK.scalartype(::Type{ExactHybridizationKernel{T}}) where {T} = T

# QuAPI 关联函数容器的标量型（kernel 构造端 T = scalartype(corr)——实谱虚时
# 核全程 Float64、复谱/实时/混合核 ComplexF64；Real/Mixed 在 QuAPI 结构上
# 即 ComplexF64 约束）。BCS 容器 = 味块类型的委托。
TK.scalartype(::Type{ImagCorrelationFunction{T}}) where {T} = T
TK.scalartype(::Type{<:RealCorrelationFunction}) = ComplexF64
TK.scalartype(::Type{<:MixedCorrelationFunction}) = ComplexF64
TK.scalartype(::Type{<:BCSCorrelationFunction{C}}) where {C} = TK.scalartype(C)

"分支格点数（窗口化基准）"
_branch_len(lattice::AbstractGrassmannLattice, b::Symbol) =
    b === :τ ? lattice.kτ :
    lattice isa AbstractRealGrassmannLattice ? lattice.k : lattice.Nt

"""
    _window_block(B, k1, k2; T=ComplexF64) -> Union{CorrelationMatrix{T},Matrix{T}}

branch 块的 lattice 网格窗口化 + 标量型归一：`CorrelationMatrix` = Toeplitz
向量前缀切片（零拷贝语义、数据不物化；换型仅当 T ≠ 块元素类型——动态
scalartype 缺省恒等）；一般矩阵 = 视图（需要窗口化/换型时物化）。窗口不足
即形状守卫响抛。
"""
function _window_block(B::CorrelationMatrix, k1::Int, k2::Int, ::Type{T}) where {T<:Number}
    n = size(B, 1)
    (size(B, 2) == n) || throw(DimensionMismatch("non-square CorrelationMatrix block"))
    (n >= k1 && n >= k2) || throw(DimensionMismatch(
        "kernel shape: corr window $n too small for lattice block ($k1,$k2)"))
    v1 = eltype(B) == T ? B.ηⱼₖ : T.(B.ηⱼₖ)
    v2 = eltype(B) == T ? B.ηₖⱼ : T.(B.ηₖⱼ)
    (k1 == n && k2 == n) && return CorrelationMatrix(v1, v2)
    return CorrelationMatrix(v1[1:k1], v2[1:k2])
end

function _window_block(B::AbstractMatrix, k1::Int, k2::Int, ::Type{T}) where {T<:Number}
    (size(B, 1) >= k1 && size(B, 2) >= k2) || throw(DimensionMismatch(
        "kernel shape: corr window $(size(B)) too small for lattice block ($k1,$k2)"))
    (size(B, 1) == k1 && size(B, 2) == k2 && eltype(B) == T) && return B
    return T.(view(B, 1:k1, 1:k2))
end

# 构造期符号折叠（混合 τ↔rt 交叉块整体取负；_NegView 惰性视图的物化替代）：
# -W——CorrelationMatrix 走 QuAPI Base.:- 向量取负、一般矩阵 = view 取负物化。

# ---------------------------------------------------------------------------
# 构造：normal bath（QuAPI Δτ/Δt/Δm 直提取）
# ---------------------------------------------------------------------------

"""
    correlationfunction(bath::AbstractNormalBath, lattice; band=nothing)
        -> FrozenHybridizationKernel
    correlationfunction(bath::AbstractNormalBath, lattice::AbstractImagGrassmannLattice,
                        scheme::IFD.ExactEqualDistantKernel; band=nothing)
        -> ExactHybridizationKernel

从 QuAPI 离散化函数直接提取 canonical 单元级核（IF = exp(Σ ā Δ_ij a_j)）——
secondorderif kernels.jl 同款三路径（ED 验证过）：

- 虚时：Δτ 的单元矩阵窗口 [1:Nτ, 1:Nτ]；
- 实时：Δt 谱路径（离散谱精确极点装配），四块 [1:Nt, 1:Nt]（块值 = 单元级
  全宽积分；实时跨支块按存储约定保留 QuAPI 内嵌 Keldysh 号 = 指数系数）；
- 混合：Δm 九块，rt [1:Nt]/τ [1:Nτ]；τ↔rt 交叉块构造期整体取负（QuAPI
  Δm 相对物理核差一个负号）。

`scheme`（缺省 `IFD.QuAPI()`）= 冻结核路线（Frozen 容器）。**虚时 + Exact**
（`IFD.ExactEqualDistantKernel`）= 严格核路线：`exact_fermionic_Δτ` 返回
`(ImagCorrelationFunction, normdelta)` 元组——blocks 照常窗口化（τ 环点值
核，CorrelationMatrix 存储）、normdelta 透传 `ExactHybridizationKernel`
（E4′：严格性全浓缩在带号标量 Z 常数；real/mixed + Exact 的装配级容器
（骨架+踢契约）不走本路线，见 `exactif_kernel`）。

`band` = 单带构造（blocks 仅该带的键、bands = [band]）；缺省全带（各带键
值共享同一底层块对象，零拷贝）。`coupling_norm` 自动 = 带数 ×
`coupling_norm(bath)`（全带总量 Σ_b α_b）。β 不一致告警保留。
"""
function correlationfunction(bath::AbstractNormalBath, lattice::AbstractGrassmannLattice,
                             scheme::IFD.DiscretizationScheme=IFD.QuAPI();
                             band::Union{Nothing,Int}=nothing)
    # scheme 路由（exactif_port_plan Phase A 构造路线归位）：QuAPI = 冻结核
    # （Frozen）；ExactKernel 仅虚时格放行（严格 τ 环点值核 + normdelta 元组，
    # Exact 容器）；real/mixed + Exact 的装配级骨架+踢契约待 exactif_kernel
    # 路由（Phase B）
    scheme isa IFD.QuAPI || scheme isa IFD.ExactKernel || throw(ArgumentError(
        "unsupported discretization scheme $(typeof(scheme))"))
    exact = scheme isa IFD.ExactKernel
    exact && !(lattice isa AbstractImagGrassmannLattice) && throw(ArgumentError(
        "real/mixed lattices do not support the Exact scheme in correlationfunction: " *
        "the strict containers carry assembly-level contracts (skeleton + kicks, " *
        "Z = normdelta · Pf(B)) — use exactif_kernel (see plan/exactif_port_plan.md)"))
    bands = band === nothing ? collect(1:lattice.bands) :
            (1 <= band <= lattice.bands || throw(BoundsError(1:lattice.bands, band));
             [band])
    if lattice isa Union{AbstractImagGrassmannLattice,AbstractMixedGrassmannLattice}
        isapprox(lattice.β, bath.β; rtol=1.0e-12, atol=1.0e-14) ||
            @warn "lattice.β=$(lattice.β), but bath.β=$(bath.β)"
    end
    corr = if lattice isa AbstractImagGrassmannLattice
        res = IFD.Δτ(bath; N=lattice.Nτ, δτ=lattice.δτ, scheme)
        if exact
            normdelta = res[2]                     # exact 元组解包（E4′）
            res[1]
        else
            res
        end
    elseif lattice isa AbstractRealGrassmannLattice
        IFD.Δt(bath; N=lattice.Nt, t=lattice.Nt * lattice.δt, scheme)
    else
        IFD.Δm(bath; Nτ=lattice.Nτ, t=lattice.Nt * lattice.δt, Nt=lattice.Nt, δτ=lattice.δτ, scheme)
    end
    T = TK.scalartype(corr)
    blocks = Dict{Tuple{ContourBlockIndex,ContourBlockIndex},Union{Matrix{T},CorrelationMatrix{T}}}()
    mixed = lattice isa AbstractMixedGrassmannLattice
    for b1 in branches(lattice), b2 in branches(lattice)
        fold = mixed && ((b1 === :τ) !== (b2 === :τ))
        W = _window_block(branch(corr; b1=b1, b2=b2), _branch_len(lattice, b1),
                          _branch_len(lattice, b2), T)
        blocks[(ContourBlockIndex(b1, first(bands), true),
                ContourBlockIndex(b2, first(bands), false))] = fold ? -W : W
        for bd in bands[2:end]                     # 多带键共享同一块对象（零拷贝）
            blocks[(ContourBlockIndex(b1, bd, true), ContourBlockIndex(b2, bd, false))] =
                blocks[(ContourBlockIndex(b1, first(bands), true),
                        ContourBlockIndex(b2, first(bands), false))]
        end
    end
    if exact
        return ExactHybridizationKernel(blocks, normdelta, bands)
    end
    return FrozenHybridizationKernel(blocks, Float64(length(bands) * coupling_norm(bath)), bands)
end

# ---------------------------------------------------------------------------
# 构造：BCS bath（四味块按 (branch, band) 键入同一 Dict）
# ---------------------------------------------------------------------------

# BCS 味键表：同带键 = 正常通道（行 ā、列 a，uu 味块）；band1→band2 = 反常
# ā↑ā↓（ud 味块，双腿 ā）；band2→band1 = 反常 a↓a↑（du 味块，双腿 a）。
# band1 = 2·orbital−1、band2 = 2·orbital（QuAPI GrassmannContourIndex 的
# 奇偶带味约定）。味表构造期写进 ContourBlockIndex.conj（键自完备）。

"BCS 块键的腿 conj 旗标 (行, 列)（构造期写进 ContourBlockIndex.conj）：
同带 = (ā, a) = (true, false)；band1→band2 = (ā, ā) = (true, true)；
band2→band1 = (a, a) = (false, false)"
function _bcs_leg_conjs(ba::Int, bb::Int)
    ba == bb && return (true, false)
    return isodd(ba) ? (true, true) : (false, false)
end

function _bcs_flavor(corr, ba::Int, bb::Int)
    ba == bb && return corr.uu
    return isodd(ba) ? corr.ud : corr.du
end

"""
    correlationfunction(bath::AbstractBCSBath, lattice; orbital=nothing)
        -> FrozenHybridizationKernel

BCS bath 的核（Frozen 容器；BCS + Exact 无 normdelta 契约——IFD 刻意留白，
待 BCS 后端立项）：`orbital` 选 Nambu 对（band1 = 2·orbital−1、
band2 = 2·orbital；lattice 需 bands ≥ band2）；`orbital = nothing`（缺省，行为
同 AbstractNormalBath 的全带模式）= 全部 Nambu 对（1:lattice.bands÷2）——每个
orbital 共享同一套杂化函数（四味块每支对只装配一次，多对键零拷贝别名）。
对称性守卫 `bcs_symmetry_guard`（仅 :z2，D7）。`coupling_norm` = bands 数 ×
`coupling_norm(bath)`（带对各 α：单对 = 2×，全轨道 = 2·轨道数×）。轮廓约定：
三轮廓统一 = 味块 branch 视图原值（raw，CorrelationMatrix/Toeplitz，与
normal 同构——Δ=0 时逐位等于 normal 核）；旧顶点移位口径与跨支取负均已
移除（v3-① + T6 探针否决）。
"""
function correlationfunction(bath::AbstractBCSBath, lattice::AbstractGrassmannLattice,
                             scheme::IFD.DiscretizationScheme=IFD.QuAPI();
                             orbital::Union{Nothing,Int}=nothing)
    scheme isa IFD.QuAPI || throw(ArgumentError(
        "BCS correlationfunction is restricted to the QuAPI scheme (frozen kernel), " *
        "got $(typeof(scheme)) — strict BCS normdelta is an IFD leave-blank (BCS backend project)"))
    bcs_symmetry_guard(symmetry(lattice))
    orbitals = orbital === nothing ? collect(1:lattice.bands ÷ 2) :
               (2 * orbital <= lattice.bands || throw(ArgumentError(
                    "BCS orbital $orbital requires bands ≥ $(2 * orbital), got $(lattice.bands)"));
                [orbital])
    isempty(orbitals) && throw(ArgumentError(
        "BCS bath requires at least one Nambu pair, got bands = $(lattice.bands)"))
    if lattice isa Union{AbstractImagGrassmannLattice,AbstractMixedGrassmannLattice}
        isapprox(lattice.β, bath.β; rtol=1.0e-12, atol=1.0e-14) ||
            @warn "lattice.β=$(lattice.β), but bath.β=$(bath.β)"
    end
    corr = if lattice isa AbstractImagGrassmannLattice
        IFD.Δτ(bath; N=lattice.Nτ, δτ=lattice.δτ, scheme)
    elseif lattice isa AbstractRealGrassmannLattice
        IFD.Δt(bath; N=lattice.Nt, t=lattice.Nt * lattice.δt, scheme)
    else
        IFD.Δm(bath; Nτ=lattice.Nτ, t=lattice.Nt * lattice.δt, Nt=lattice.Nt, δτ=lattice.δτ, scheme)
    end
    T = TK.scalartype(corr)
    blocks = Dict{Tuple{ContourBlockIndex,ContourBlockIndex},Union{Matrix{T},CorrelationMatrix{T}}}()
    for b1 in branches(lattice), b2 in branches(lattice)
        # 四味块每支对只装配一次（多 Nambu 对共享同一套块对象——零拷贝别名）；
        # 味选择按 QuAPI 奇偶带约定（代表带 (1, 2)，带号超出奇偶无关），
        # orbital o 的实际带 = 2o−2+dba / 2o−2+dbb
        for (dba, dbb) in ((1, 1), (1, 2), (2, 2), (2, 1))
            ca, cb = _bcs_leg_conjs(dba, dbb)
            # raw 味块直存——bcs_Δt 味块 = 标准 fermionic_Δt 产物，支视图自带与
            # normal Δt 相同的指数系数 Keldysh 号（Δ=0 对照 normal 逐位实证 0.0，
            # 2026-11-19 T6 探针：计划预留的「跨支构造期取负」被否决——取反即错）
            W = _window_block(branch(_bcs_flavor(corr, dba, dbb); b1=b1, b2=b2),
                              _branch_len(lattice, b1), _branch_len(lattice, b2), T)
            for o in orbitals
                blocks[(ContourBlockIndex(b1, 2o - 2 + dba, ca),
                        ContourBlockIndex(b2, 2o - 2 + dbb, cb))] = W
            end
        end
    end
    bands = reduce(vcat, ([2o - 1, 2o] for o in orbitals))   # orbitals 升序 ⟹ 有序
    return FrozenHybridizationKernel(blocks, length(bands) * coupling_norm(bath), bands)
end

# ---------------------------------------------------------------------------
# changebands：带重排（单带核免重建 corr 复用出多带键）
# ---------------------------------------------------------------------------

"""
    changebands(H::AbstractHybridizationKernel, x::Vector{Pair{Int,Int}})
        -> 同型 HybridizationKernel

带重排：按 `x = [old => new, …]` 映射 blocks 键的 `ContourBlockIndex.band`
（同一 old 可映射到多个 new——块对象别名共享，标量/矩阵数据零改动），
`bands` = 重排后的有序并集。常数面随型保留（Frozen `coupling_norm` /
Exact `normdelta`）。用途：单带核建一份、重排出多带装配所需带号。
"""
function changebands(H::AbstractHybridizationKernel, x::Vector{Pair{Int,Int}})
    for old in first.(x)
        old in H.bands || throw(ArgumentError(
            "changebands: source band $old not in kernel bands $(H.bands)"))
    end
    news = last.(x)
    allunique(news) && all(b -> b >= 1, news) || throw(ArgumentError(
        "changebands: target bands must be distinct positive integers, got $news"))
    blocks = Dict{Tuple{ContourBlockIndex,ContourBlockIndex},eltype(values(H.blocks))}()
    for (old, new) in x                          # 逐对迭代——同一 old 可映射多个 new
        for (key, value) in H.blocks             # （单带核扩多带：块对象别名共享）
            ba = key[1].band == old ? new : key[1].band
            bb = key[2].band == old ? new : key[2].band
            blocks[(ContourBlockIndex(key[1].branch, ba, key[1].conj),
                    ContourBlockIndex(key[2].branch, bb, key[2].conj))] = value
        end
    end
    bands = sort(union!(collect(Set(last.(x))), collect(Set(
        b for b in H.bands if b ∉ first.(x)))))
    return _retype_kernel(H, blocks, bands)
end

"同型核重建（changebands/_band_kernel 共用；常数面按型透传）"
function _retype_kernel(H::AbstractHybridizationKernel, blocks, bands::Vector{Int})
    H isa FrozenHybridizationKernel &&
        return FrozenHybridizationKernel(blocks, H.coupling_norm, bands)
    H isa ExactHybridizationKernel &&
        return ExactHybridizationKernel(blocks, H.normdelta, bands)
    throw(ArgumentError("unsupported kernel type $(typeof(H))"))
end

# ---------------------------------------------------------------------------
# Base.:+：多 bath 耦合同一杂质（双 lead 逐 lead 建核后相加）
# ---------------------------------------------------------------------------

"""
    Base.:+(A::FrozenHybridizationKernel, B::FrozenHybridizationKernel)
        -> FrozenHybridizationKernel

多 bath 耦合的核相加：blocks 键 = 并集（同键块相加——`CorrelationMatrix` +
`CorrelationMatrix` = QuAPI Toeplitz 向量相加，混合表示 = 泛型逐元素（任一方
dense 即 dense 结果）；单侧缺键 = 透传并标量型对齐 T）；`coupling_norm` 相加；
`bands` = 有序并集；标量型 = 两者 `scalartype` 的 promote。用途：双 lead（异
μ/β 无法合成单 bath）逐 lead 建 kernel 后相加——总耦合标量自动正确。

Exact 核不定义加法（`+(::Exact, ::Exact)` 与混合 `+(::Frozen, ::Exact)` 均
`MethodError`/ArgumentError）：normdelta 是行列式型常数，多浴联合装配是否
恰好因子化（blocks 加 + normdelta 乘）= exactif_port_plan Phase B 数值裁决门，
裁决落地前不提供（恒错禁止：normdelta 相加 / 绝对 Z 静默降级）。
"""
function Base.:+(A::FrozenHybridizationKernel, B::FrozenHybridizationKernel)
    T = promote_type(TK.scalartype(A), TK.scalartype(B))
    blocks = Dict{Tuple{ContourBlockIndex,ContourBlockIndex},Union{Matrix{T},CorrelationMatrix{T}}}()
    for key in union(keys(A.blocks), keys(B.blocks))
        inA = haskey(A.blocks, key); inB = haskey(B.blocks, key)
        blocks[key] = inA && inB ? A.blocks[key] + B.blocks[key] :
                      _convert_block(T, inA ? A.blocks[key] : B.blocks[key])
    end
    return FrozenHybridizationKernel(blocks, A.coupling_norm + B.coupling_norm,
                                     sort(union(A.bands, B.bands)))
end

Base.:+(A::ExactHybridizationKernel, B::ExactHybridizationKernel) =
    throw(ArgumentError(
        "Exact kernel addition is pending the factorization verdict " *
        "(exactif_port_plan Phase B: joint two-bath normdelta · blocks+ must be " *
        "verified against the IFD gold standard first)"))
Base.:+(A::AbstractHybridizationKernel, B::AbstractHybridizationKernel) =
    throw(ArgumentError(
        "mixed kernel addition Frozen+Exact is not defined (different Z-constant " *
        "contracts: coupling_norm Σα vs signed normdelta)"))

"块对象标量型对齐（跨标量核相加的缺键透传；同型零拷贝）"
_convert_block(::Type{T}, v::CorrelationMatrix) where {T} =
    eltype(v) == T ? v : CorrelationMatrix(T.(v.ηⱼₖ), T.(v.ηₖⱼ))
_convert_block(::Type{T}, v::AbstractMatrix) where {T} = eltype(v) == T ? v : T.(v)


"逐带 kernel 视图：只保留 band 的键、bands = [band]（块对象别名共享，零拷贝）。
多带 normal 核的逐带装配内部切片用。"
function _band_kernel(kernel::AbstractHybridizationKernel, band::Int)
    band in kernel.bands || throw(BoundsError(kernel.bands, band))
    blocks = Dict{Tuple{ContourBlockIndex,ContourBlockIndex},eltype(values(kernel.blocks))}()
    for (key, value) in kernel.blocks
        (key[1].band == band && key[2].band == band) || continue
        blocks[key] = value
    end
    return _retype_kernel(kernel, blocks, [band])
end

# ---------------------------------------------------------------------------
# coupling_norm：浴谱总强度（kernel.coupling_norm 的单浴因子；多带总量 =
# Σ_b coupling_norm(bath_b)。μ/β 不进谱权重——它们只进 Δ 装配的核因子
# （spectrumshift(f, μ) 等），离散谱求和天然 μ 不变，bath 层直取谱）
# ---------------------------------------------------------------------------

"""
    coupling_norm(f) -> Real
    coupling_norm(bath::AbstractNormalBath) -> Real

浴谱总强度 α（= 单带核的 `coupling_norm` 因子；全带总量 = 带数 × 该值，
多浴 = 各浴核相加自动累加）：

- `DiracDelta`：α；
- `DiscreteSpectrum`：Σₙ αₙ（fs 为 δ-峰权重，vₙ = √αₙ = `spectrumcouplings`）；
- 连续 `AbstractBoundedFunction`：∫ f(ε)dε over [lb, ub]（quadgk 数值积分，
  容差同 IMB QuadGKTolerance）；
- δ×函数乘积（`DeltaMultF`/`DeltasMultF`，如热流的 ω 加权谱）：走各自
  quadgkwrapper 求和/取值规则（g ≡ 1）——同为 AbstractBoundedFunction，
  由 `quadgkwrapper(f)` 一并覆盖；
- `AbstractNormalBath`：直接取其谱——**与 μ、β 无关**；
- `AbstractBCSBath`：裸谱强度 ∫f dε（每带，`bath.spectrum` 直取——与 normal
  同款）。**Δ、μ、β 都不进**：Weyl 标量只数对角接触（āσaσ 槽）的宽度，
  Bogoliubov 正则性（u²+v²=1）保证反常通道（uv 块，ā↑ā↓/a↓a↑ 槽）对标量
  零贡献、总宽度 = 裸谱 normal 范数（推导与数值裁决见 `bcsweyl/README.md`：
  平凡系统接触缺口提取 cn*(Δ)/Σα ≡ 1（≤4e-4，Δ∈[0,2]）；generic H_S 下
  Z 收敛斜率 1.97–1.98）。一个 BCS 浴 = 一个带对（两带各 α）⟹ K 侧
  全带总量对该对按 2× 计入。

扩展点：自定义谱子类实现 `quadgkwrapper` 乘积路径或直接补 `coupling_norm`
方法。
"""
coupling_norm(f::DiracDelta) = f.α
coupling_norm(f::DiscreteSpectrum) = sum(spectrumvalues(f))
coupling_norm(f::AbstractBoundedFunction) = quadgkwrapper(f)
coupling_norm(bath::AbstractNormalBath) = coupling_norm(bath.spectrum)
coupling_norm(bath::AbstractBCSBath) = coupling_norm(bath.spectrum)
