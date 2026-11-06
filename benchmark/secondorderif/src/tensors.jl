# ============================================================================
# tensors.jl — 单元核 → Z2 site tensors（K 链的局域张量装配）
#
# Weyl 系统因子是一条 automaton 链：虚拟键携带 (第一成员, 第二成员) 两个 Fock
# 指标（对齐 U1GTEMPO2 system_factor_reference 的 pair bond；虚时分支的第二
# 成员恒等）。每个单元核
#
#     C[(s_L,t_L), y₁..y_V, (s_R,t_R)]
#
# 沿 site 顺序做逐变量的 Z2 分级 SVD 拆分（U1GTEMPO2 `split_cell` 的
# Z2Tensors 等价物；尺寸 ≤ 16×2^V×16，精确、无截断）。拆出的 site tensors
# 直接铺入 vacuum GrassmannMPS —— 各单元的窗口互不相交，平铺即精确算符乘积
# （同 GTEMPO sysdynamics_fast 的思路）。
#
# 约定：偶/奇 Grassmann 变量携带 Z2Irrep(0)/Z2Irrep(1)，键宇称 = 已固定变量
# 宇称之和（mod 2）；键空间 Z2Space(0=>n₀, 1=>n₁)，扇区 0 状态排在前，扇区内
# 顺序由拆分过程显式记录并在相邻张量间保持一致。
# ============================================================================

⊕(a::Integer, b::Integer) = Int((a + b) % 2)
paritybits(x::Integer) = Int(count_ones(x) % 2)

_z2space(d0::Int, d1::Int) = (d1 == 0) ? Z2Space(0 => d0) :
                             (d0 == 0) ? Z2Space(1 => d1) : Z2Space(0 => d0, 1 => d1)

"""
    CellCore

单元核：`data[(in, y₁, …, y_V, out)]`（y 按 site 顺序的变量占据位 0/1），
`inpar`/`outpar` 为边界键各状态的宇称。
"""
struct CellCore
    data::Array{ComplexF64}
    inpar::Vector{Int}
    outpar::Vector{Int}
    function CellCore(data::Array{ComplexF64}, inpar::Vector{Int}, outpar::Vector{Int})
        (ndims(data) >= 3) || throw(ArgumentError("单元核至少要有一个变量"))
        iseven(ndims(data) - 2) || throw(ArgumentError("变量数必须为偶数"))
        length(inpar) == size(data, 1) || throw(DimensionMismatch("inpar"))
        length(outpar) == size(data, ndims(data)) || throw(DimensionMismatch("outpar"))
        new(data, inpar, outpar)
    end
end

"""
    split_core(core::CellCore) -> Vector{TensorMap}

沿变量顺序逐个拆分，返回 V 个 site tensor，第 k 个形如
`(bond_{k−1} ⊗ phys) ← bond_k`，bond₀ / bond_V 为平凡真空空间。
"""
function split_core(core::CellCore)
    (; data, inpar, outpar) = core
    T = ComplexF64
    V = ndims(data) - 2
    outsize = size(data, ndims(data))
    tensors = Vector{TensorMap}(undef, V)

    work = data                        # 轴 (bond, v_k, …, v_V, out)
    bondpar = copy(inpar)
    for k in 1:V-1
        B = size(work, 1)
        tail = collect(size(work)[3:end])       # (v_{k+1}…v_V, out)
        insec = [1 + count(==(bondpar[b]), bondpar[1:b-1]) for b in 1:B]  # 扇区内序号
        tailsize = prod(tail)
        mat = reshape(work, B * 2, tailsize)    # 行 = (bond × v_k)，v_k 最快
        remvars = 2^(V - k)                     # 剩余变量的占据位串宽度
        rowpar = [bondpar[(r - 1) % B + 1] ⊕ ((r - 1) ÷ B) for r in 1:B*2]
        colpar = [paritybits((j - 1) % remvars) ⊕ outpar[(j - 1) ÷ remvars + 1]
                  for j in 1:tailsize]
        # 按中间键宇称 σ 分块做精确 SVD；U 的列序 = 新键扇区内状态序
        Us, Ss, Vhs, colsσ = Dict{Int,Matrix{ComplexF64}}(), Dict{Int,Vector{Float64}}(),
                            Dict{Int,Matrix{ComplexF64}}(), Dict{Int,Vector{Int}}()
        for σ in (0, 1)
            rows = findall(==(σ), rowpar)
            cols = findall(==(σ), colpar)
            (isempty(rows) || isempty(cols)) && continue
            U, S, Vh = svd_thin(mat[rows, cols])
            isempty(S) && continue
            Us[σ] = U
            Ss[σ] = S
            Vhs[σ] = Vh
            colsσ[σ] = cols
        end
        n0 = haskey(Us, 0) ? size(Us[0], 2) : 0
        n1 = haskey(Us, 1) ? size(Us[1], 2) : 0
        midpar = vcat(fill(0, n0), fill(1, n1))

        # site tensor k：(bond ⊗ v_k) ← mid；通道 (c_bond, c_y, σ=c_b⊕c_y)。
        # U 的行 = rows_σ（r 序：y 大股、b 细股），通道行须按位过滤
        lspace = _z2space(count(==(0), bondpar), count(==(1), bondpar))
        t = zeros(T, lspace ⊗ _PH ← _z2space(n0, n1))
        rowsof = [(b, y) for y in 0:1 for b in 1:B]
        for (c_b, c_y) in ((0, 0), (0, 1), (1, 0), (1, 1))
            σ = c_b ⊕ c_y
            haskey(Us, σ) || continue
            U = Us[σ]
            nσ = σ == 0 ? n0 : n1
            rows_σ = [by for by in rowsof if bondpar[by[1]] ⊕ by[2] == σ]
            chanpos = [pos for (pos, (b, y)) in enumerate(rows_σ)
                       if bondpar[b] == c_b && y == c_y]
            isempty(chanpos) && continue
            block = t[(Z2Irrep(c_b), Z2Irrep(c_y), Z2Irrep(σ))]   # (d_l, d_ph, d_r) 视图
            for (ri, pos) in enumerate(chanpos)
                b = rows_σ[pos][1]
                for j in 1:nσ
                    block[insec[b], 1, j] = U[pos, j]
                end
            end
        end
        tensors[k] = t

        # 新 work = S·Vᵀ（行 = mid 状态：先扇区 0 后扇区 1，与 midpar 一致）
        newwork = zeros(T, n0 + n1, tailsize)
        for σ in keys(Us)
            S, Vh, cols = Ss[σ], Vhs[σ], colsσ[σ]
            offset = σ == 0 ? 0 : n0
            for ci in eachindex(S), cs in eachindex(cols)
                newwork[offset+ci, cols[cs]] = S[ci] * Vh[ci, cs]
            end
        end
        work = reshape(newwork, n0 + n1, tail...)
        bondpar = midpar
    end

    # 最后一个变量：work 轴 (bond, v_V, out)：(bond ⊗ v_V) ← out 键。
    # out 空间按宇比分扇区（扇区 0 状态在前），与相邻单元的 in 空间逐一对应
    B = size(work, 1)
    n0 = count(==(0), bondpar)
    outn0 = count(==(0), outpar)
    t = zeros(T, _z2space(n0, B - n0) ⊗ _PH ←
                 _z2space(outn0, outsize - outn0))
    for b in 1:B, o in 1:outsize
        y = bondpar[b] ⊕ outpar[o]              # 宇称守恒决定 v_V 的占据
        idx = 1 + count(==(bondpar[b]), bondpar[1:b-1])
        oidx = 1 + count(==(outpar[o]), outpar[1:o-1])
        block = t[(Z2Irrep(bondpar[b]), Z2Irrep(y), Z2Irrep(outpar[o]))]
        block[idx, 1, oidx] = work[b, y + 1, o]
    end
    tensors[V] = t
    return tensors
end

"瘦 SVD：丢弃相对 1e-14 的数值零奇异值；返回 (U, S, Vᴴ)，M ≈ U·diag(S)·Vᴴ"
function svd_thin(M::AbstractMatrix{ComplexF64})
    U, S, V = svd(M)
    smax = isempty(S) ? 0.0 : maximum(S)
    r = count(>(1e-14 * smax), S)
    return U[:, 1:r], S[1:r], Matrix{ComplexF64}(V[:, 1:r]')
end

"""
    tile_cells!(mps::GrassmannMPS, positions, tensors) -> GrassmannMPS

把互不相交窗口的 site tensors 铺入 vacuum GMPS 的对应位置。
"""
function tile_cells!(mps::GrassmannMPS, positions::Vector{Int},
                     tensors::Vector{<:TensorMap})
    length(positions) == length(tensors) ||
        throw(DimensionMismatch("positions 与 tensors 数目不一致"))
    allunique(positions) || throw(ArgumentError("窗口重叠"))
    for (p, t) in zip(positions, tensors)
        mps[p] = t
    end
    GTEMPO.unset_svectors!(mps)
    return mps
end
