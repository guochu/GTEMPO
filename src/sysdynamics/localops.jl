# localops.jl：符号化局域算符（LocalOp/FockOp）——插入算符的纯符号容器层
# （plan/fockmpo_dress_selfread_plan.md §2.1；2026-11-20 合并裁决 V1–V7）
#
# 设计钉定：
#   - V2/V3：插入全程符号态；荷从符号按位派生，无手工标签——"prefix 真凶"类
#     标签错误从源头消灭；
#   - V4：FockOp 无乘法（被动容器；LocalOp↔FockOp 混合乘法亦不定义）；
#   - V5：LocalOp * LocalOp → LocalOp 符号级乘法表封闭（M₂ 代数 + 零吸收）；
#   - V6：LocalOp 无带信息（纯 Symbol）；FockOp.data[b] = 第 b 带因子
#     （1-based，长度 = bands 在物化/荷派生点校验）。跨带积算符序 = 带序
#     （band 1 最左，`local_weights` fold 同款）；
#   - V7：物化只发生在生成 FockMPO 的装配点（`_fockop_matrix`），不保留
#     raw-matrix 构造器。物化 = 逐带矩阵 fold，每带矩阵自带 JW 符号，普通
#     矩阵积符号天然正确——不加任何 graded⊗ 符号（那是槽编码装置、
#     `graded_sign_factor = 1` 关闭，与算符代数无关）。

# ---------------------------------------------------------------------------
# 类型（被动容器；结构化 ==/hash 是插入缓存键的先决）
# ---------------------------------------------------------------------------

"""
    LocalOp(s::Symbol)

单带局域算符符号（无带信息——V6；带定位由 FockOp 的位置承载）。六符号张成
单带 M₂ 代数 + 零吸收（`:a·:n = :a` 而 `:n·:a = :empty` 等非对称积见 `*`）：

| s | 单带矩阵 | 荷 |
|---|---|---|
| `:a` | d | −q |
| `:adag` | d† | +q |
| `:n` | d†d | 0 |
| `:aadg` | d d†（= I−n，物化为积保留语义） | 0 |
| `:I` | I | 0 |
| `:empty` | 0（零算符） | 0 |
"""
struct LocalOp
    s::Symbol
    function LocalOp(s::Symbol)
        s in (:a, :adag, :n, :aadg, :I, :empty) ||
            throw(ArgumentError("unknown LocalOp symbol $s (valid: :a :adag :n :aadg :I :empty)"))
        new(s)
    end
end

Base.hash(op::LocalOp, h::UInt) = hash(op.s, h)

"""
    FockOp(data::Vector{LocalOp})

带序局域算符元组（被动容器，无乘法——V4）：`data[b]` = 第 b 带因子
（1-based；长度 = bands 在物化/荷派生点校验，构造不自校验）。物化（V7）=
逐带矩阵 fold，band 1 最左——跨带积算符序 = 带序（设计属性；其它序 = 未来
扩展点）。`data` 含 `:empty` ⟹ 整积零矩阵（物化短路、荷 0）。
"""
struct FockOp
    data::Vector{LocalOp}
end

Base.:(==)(a::FockOp, b::FockOp) = a.data == b.data
Base.hash(op::FockOp, h::UInt) = hash(op.data, h)

# ---------------------------------------------------------------------------
# LocalOp 符号级乘法表（V5）：行 = 左因子、列 = 右因子，右先作用（矩阵积序）。
# 单带 span{I, d, d†, n} = M₂ 代数 + 零吸收；36 项由 2×2 矩阵乘逐位验证
# （test/api/localops.jl）。`:empty` = 零算符吸收一切（a·a = 0 等封闭性要求）。
# ⚠️ 计划 §2.1 表中 (adag,n)/(adag,aadg)/(aadg,n)/(aadg,aadg) 四项与矩阵代数
# 矛盾（两行的 n↔aadg 列互换）——本表按矩阵代数落地（d†·n = 0、
# d†·(I−n) = d†、(I−n)·n = 0、(I−n)² = I−n），见计划 §8 实施偏差。
# ---------------------------------------------------------------------------

const _LOP_SYMBOLS = (:I, :a, :adag, :n, :aadg, :empty)

const _LOP_MULT = [:I     :a      :adag   :n      :aadg   :empty;
                   :a     :empty  :aadg   :a      :empty  :empty;
                   :adag  :n      :empty  :empty  :adag   :empty;
                   :n     :empty  :adag   :n      :empty  :empty;
                   :aadg  :a      :empty  :empty  :aadg   :empty;
                   :empty :empty  :empty  :empty  :empty  :empty]

_lop_index(s::Symbol) = findfirst(==(s), _LOP_SYMBOLS)::Int

"""
    LocalOp(x) * LocalOp(y) -> LocalOp

符号级乘法（V5；唯一乘法）：矩阵积序（右因子先作用）。非对称注记：
`:a·:n = :a` 而 `:n·:a = :empty`（d·d†·d = d vs d†·d·d = 0）——按行列序查。
自动化简不做；跨带组合走 FockOp 位置排布（不是乘法）。
"""
function Base.:(*)(x::LocalOp, y::LocalOp)
    return LocalOp(_LOP_MULT[_lop_index(x.s), _lop_index(y.s)])
end

# ---------------------------------------------------------------------------
# 构造 / 物化 / 荷派生
# ---------------------------------------------------------------------------

"单带因子 helper（测量 API 主构造路径）：目标带放符号、其余 `:I` 填充。"
function _single_band_op(s::Symbol, band::Int, bands::Int)
    s in (:a, :adag, :n, :aadg, :I, :empty) ||
        throw(ArgumentError("unknown local operator symbol $s"))
    (1 <= band <= bands) || throw(BoundsError(1:bands, band))
    return FockOp([LocalOp(b == band ? s : :I) for b in 1:bands])
end

"单带符号 → 该带矩阵（band b 的 JW 穿越号由 _fock_annihilators 自带）"
function _localop_matrix(op::LocalOp, d::Vector{Matrix{ComplexF64}}, b::Int)
    s = op.s
    s === :a && return d[b]
    s === :adag && return d[b]'
    s === :n && return d[b]' * d[b]
    s === :aadg && return d[b] * d[b]'
    s === :I && return Matrix{ComplexF64}(I, size(d[b]))
    return zeros(ComplexF64, size(d[b]))   # :empty（构造器守卫后不可达）
end

"""
    _fockop_matrix(op::FockOp, bands) -> Matrix{ComplexF64}

V7 装配点物化：逐带矩阵 fold，band 1 最左（`O = O * mat_b, b = 1..bands`，
`local_weights` fold 同款）。每带全空间矩阵自带 JW 符号，普通矩阵积符号天然
正确——不加任何 graded⊗ 符号。`data` 含 `:empty` ⟹ 整积零矩阵（短路）。
"""
function _fockop_matrix(op::FockOp, bands::Int)
    length(op.data) == bands ||
        throw(ArgumentError("FockOp has $(length(op.data)) band factors, expected bands = $bands"))
    any(x -> x.s === :empty, op.data) && return zeros(ComplexF64, 2^bands, 2^bands)
    d = _fock_annihilators(bands)
    O = Matrix{ComplexF64}(I, 2^bands, 2^bands)
    for b in 1:bands
        O = O * _localop_matrix(op.data[b], d, b)
    end
    return O
end

"""
    _fockop_charge(op::FockOp, G, bands) -> Sector

V3 荷派生：`Σ_b charge(data[b], b)`（:a = −q_b、:adag = +q_b、:n/:aadg/:I
= 0）；`data` 含 `:empty` ⟹ 零算符荷 0。奇偶 = `GM.fermionparity(荷)`。
"""
function _fockop_charge(op::FockOp, ::Type{G}, bands::Int) where {G<:Sector}
    length(op.data) == bands ||
        throw(ArgumentError("FockOp has $(length(op.data)) band factors, expected bands = $bands"))
    any(x -> x.s === :empty, op.data) && return one(G)
    q = one(G)
    for b in 1:bands
        s = op.data[b].s
        s === :adag && (q = qadd(q, var_charge(G, b, bands)))
        s === :a && (q = qadd(q, qneg(var_charge(G, b, bands))))
    end
    return q
end
