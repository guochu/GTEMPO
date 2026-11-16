# operators.jl：算符路线观测量（GM ExpectationCache 快照 + SparseGMPS 单项式，
# D5 两轨之二）——旧栈 cached_gf/TwosideExpectationCache 同构：缓存一次、
# 任意多观测量只扫算符支撑；与 inserted_k.jl 的 inserted-K 重建路线在
# golden 对拍下等价（D5 两轨）。
#
# ⚠️ 本轨语义 = **胞坐标 0⁺ 点**（见 cached_gf 的 D5 注原样搬运）：Chrono 格
# 的胞变量位于 Weyl 胞两端（a = 胞左端 t_{c−1}、ā = 胞右端 t_c），单项式
# 关联是胞坐标 0⁺ 点的量，与 occupation/greater/lesser 的边界矩阵插入点差
# O(胞内接触项)（δt=0.15 时 ~1e-3 档，随 δt 收敛）——边界观测量请用
# inserted-K 路线（inserted_k.jl）。

"""
    lattice_monomial(lattice, indices::NTuple{N,GrassmannContourIndex}; coeff=1) -> GM.SparseGMPS

轮廓格点处的 Grassmann 单项式算符（旧 GTerm 同语义）：变量位 = lattice 格点
映射，表系数 = coeff·(−1)^{给定序逆序数}（升序变量基）；站张量键恒平凡、
每站一个 digit 槽置位（vac/pair/η/η̄），支撑 = 触及站的闭包。
"""
function lattice_monomial(lattice::AbstractGrassmannLattice, indices::NTuple{N,GrassmannContourIndex};
                          coeff::Number=1) where {N}
    spec = TK.sectortype(lattice)
    positions = [lattice[ix] for ix in indices]   # 变量位（给定序）
    inv = 0
    for i in eachindex(positions), j in (i + 1):length(positions)
        positions[i] > positions[j] && (inv += 1)
    end
    c0 = coeff * (-1)^inv
    touched = sort!(unique((p + 1) ÷ 2 for p in positions))
    first_s, last_s = first(touched), last(touched)
    spaces = site_spaces(lattice)
    charges = _site_charges(spec, lattice)
    z = one(spec)
    data = Any[]; q = z   # q = 已吞并变量的累积荷（键穿线；末站须回到平凡）
    for s in first_s:last_s
        o1 = (2s - 1) in positions; o2 = (2s) in positions
        E = spaces[s]
        T = TK.zeros(ComplexF64, qspace(spec, [q]) ⊗ E ← qspace(spec, [q]))
        T[(q, z, TK.dual(q))][:, 1, :] .= 1
        if o1 || o2
            basis = GM.grassmann_basis(E)
            tag = o1 ? (o2 ? basis.pair : basis.eta) : basis.etabar
            qout = qadd(qadd(q, o1 ? charges[2s - 1] : z), o2 ? charges[2s] : z)
            T = TK.zeros(ComplexF64, qspace(spec, [q]) ⊗ E ← qspace(spec, [qout]))
            T[(q, tag.sector, TK.dual(qout))][:, tag.slot, :] .=
                s == first_s ? c0 : 1
            q = qout
        end
        push!(data, T)
    end
    q == z || throw(ArgumentError("operator monomial is not charge-neutral (open legs)"))
    return GM.SparseGMPS(convert(Vector{typeof(first(data))}, data), collect(first_s:last_s))
end

"""
    environments_cache(K, Is...) -> GM.ExpectationCache

多链联合收缩缓存（K 为第一条链 = 算符目标链，GM A/Bs 约定）；Z 与全部
cached_* 观测量共用，构建一次 O(L)，快照不可变、无限复用。
"""
environments_cache(K::GM.GMPS, Is::GM.GMPS...) = GM.environments(K, Is...)

"缓存 Z（无插入分母；与 integrate(K, Is...) 同值、读左端 O(1)）"
cached_Z(cache::GM.ExpectationCache) = value(GM.Zvalue(cache))

"""
    cached_gf(lattice, indices, K, Is...) -> Complex

轮廓格点单项式的原始 Grassmann 关联函数 ⟨Π g⟩（已归一化）。

实现 = `op ∧ K` 稀疏楔积（只触碰算符支撑站）+ 引擎单发全迹
（`GM.Zvalue`，O(L) 左扫、无右扫补偿通道）——相对 `raw_correlation`
省掉逐点 K 重建（M=2 实测单点 ~88% 成本）。

⚠️ 语义（D5 裁决）：Chrono 格的胞变量位于 Weyl 胞两端（a = 胞左端 t_{c−1}、
ā = 胞右端 t_c），单项式关联是**胞坐标 0⁺ 点**的量，与 `occupation`/
`greater`/`lesser` 的边界矩阵插入点差 O(胞内接触项)（δt=0.15 时 ~1e-3 档，
随 δt 收敛）——边界观测量请用 inserted-K 路线（inserted_k.jl）。

GM 快照扫掠路线（`GM.expectation(op, env)`，只扫算符支撑）2026-11-16 复验：
与本精确路径 12 位吻合（cprime-logscale 守卫已闭合旧右扫缺口）——本包装保留
精确路径为缺省（每次查询 O(L) 且无快照生命周期约束）；多算符高频查询可换
快照路线，值等价。
"""
function cached_gf(lattice::AbstractGrassmannLattice, indices::NTuple{N,GrassmannContourIndex},
                   K::GM.GMPS, Is::GM.GMPS...; cache=nothing) where {N}
    # cache 仅为与 cached_correlation 的透传接口兼容（D5 裁决：本路线走
    # op∧K + 单发 Zvalue 精确路径，不用右扫快照——缺口见文件头 D5 注）
    op = lattice_monomial(lattice, indices)
    W = op * K
    return value(GM.Zvalue([W, Is...]))
end

"""
    cached_correlation(lattice, K, Is...; i, j=0, kind=:greater, band=1, cache) -> Complex

原始两点单项式关联（胞坐标 0⁺ 点，语义与容差见 cached_gf 的 D5 说明）：
  greater：⟨a(+,i) ā(+,j)⟩  lesser：⟨ā(−,i) a(+,j)⟩  matsubara：⟨a(τ,i) ā(τ,j)⟩
"""
function cached_correlation(lattice::AbstractGrassmannLattice, K::GM.GMPS, Is::GM.GMPS...;
                            i::Integer, j::Integer=0, kind::Symbol=:greater, band::Int=1,
                            cache::GM.ExpectationCache=GM.environments(K, Is...))
    if kind === :matsubara
        lattice isa AbstractImagGrassmannLattice || throw(ArgumentError("matsubara 需虚时轮廓"))
        idx = (GrassmannContourIndex(Int(i), :τ, band, false),
               GrassmannContourIndex(Int(j), :τ, band, true))
    elseif kind === :greater
        idx = (GrassmannContourIndex(Int(i), :+, band, false),
               GrassmannContourIndex(Int(j), :+, band, true))
    elseif kind === :lesser
        idx = (GrassmannContourIndex(Int(i), :-, band, true),
               GrassmannContourIndex(Int(j), :+, band, false))
    else
        throw(ArgumentError("kind must be :greater, :lesser or :matsubara"))
    end
    return cached_gf(lattice, idx, K, Is...; cache=cache)
end
