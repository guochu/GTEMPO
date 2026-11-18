# partialif.jl：PartialIF——整形式幂零指数装配（plan/if_expmult W3 语义换代）
#
# 内核从「逐行秩一因子积 + 乘间压缩」（旧 influence.jl 逐行路径，已随
# W1/W3 换代 + 2026-11-18 退役删除）换代为「AnomalousPairing 整形式 × GM
# expmult! PartialExp 幂零指数 + 压缩」。两者是同一 IF 算符 e^{Σ ā Δ a} 的
# 两种精确表示（偶元互相超交换 ⟹ 逐行因子积 ≡ 整形式单指数）；截断 =
# 压缩（SVDCompression）。
#
# fullfactor 的 C 矩阵约定（GM AnomalousPairing，变量 α ↔ 腿 α：奇 = η（a）、
# 偶 = η̄（ā），站点 = (α+1)÷2）：**有向单一取向存储、不存镜像**——
#   normal（band 键）：C[2·st(ā_row), 2·st(a_col) − 1] = blocks 值（ηη̄ 扇区；
#   deta 折叠 C[2i−1,2j] − C[2j,2i−1] 把唯一取向展开为全单项式系数）；
#   BCS（轨道带对键）：同带键 = 正常象限（同上式）；band1→band2 键 =
#   C[2·st(ā_b1), 2·st(ā_b2)] = ud（η̄η̄ 扇区）；band2→band1 键 =
#   C[2·st(a_b2), 2·st(a_b1)] = du（ηη 扇区）。
# U1/U1U1：normal 仅 ηη̄ 条目（荷中性）⟹ 可物化；BCS 同味扇区带净荷 ⟹
# 仅 Z2/Z2Z2（D7 守卫在 kernel 构造端已断言）。

"IF 变量位：格 (b,j) 的 band 带变量（block_offset + 2(band−1) + 1 + conj；
conj=false → a、true → ā）"
_cell_var_index(lattice::AbstractGrassmannLattice, cell::Tuple{Symbol,Int}, band::Int, conj::Bool) =
    block_offset(lattice, cell_block(lattice, cell)) + 2 * (band - 1) + Int(conj)

"kernel 是否 BCS 轨道核（存在跨带键）"
_is_orbital_kernel(kernel::HybridizationKernel) =
    any(key[1].band != key[2].band for key in keys(kernel.blocks))

"""
    fullfactor(kernel, lattice) -> GM.AnomalousPairing

IF 指数的 2L×2L 一般二次型系数矩阵（L = 链站数；块布局与变量序规范表见
文件头）。装配范围由 kernel 决定（band 信息在 `kernel.bands`，不再重复输入）：
BCS 轨道 kernel = 轨道全域因子（同带键按正常象限、跨带键按反常象限 ud/du，
腿 conj 由 `_bcs_leg_conjs` 味表回推）；normal 单带 kernel = 该带因子。
多带 normal kernel 请先切 per-band kernel（`correlationfunction(...; band=)`
或 `changebands`）或走 `influence_functional`/`hybriddynamics`（内部逐带装配）。
块值 = 构造期符号收编的指数系数（correlationfunction.jl 存储约定），直读不补号。
"""
function fullfactor(kernel::HybridizationKernel, lattice::AbstractGrassmannLattice)
    if !_is_orbital_kernel(kernel) && length(kernel.bands) != 1
        throw(ArgumentError(
            "fullfactor requires a single-band normal kernel or a BCS orbital kernel, " *
            "got bands $(kernel.bands); use influence_functional/hybriddynamics for " *
            "multi-band assembly"))
    end
    return _fullfactor(kernel, lattice)
end

function _fullfactor(kernel::HybridizationKernel, lattice::AbstractGrassmannLattice)
    L = length(lattice) ÷ 2
    C = zeros(ComplexF64, 2L, 2L)
    bcs = _is_orbital_kernel(kernel)
    for (key, blk) in kernel.blocks
        ca, cb = bcs ? _bcs_leg_conjs(key[1].band, key[2].band) : (true, false)
        b1, b2 = key[1].branch, key[2].branch
        for j1 in 1:_branch_len(lattice, b1), j2 in 1:_branch_len(lattice, b2)
            v = blk[j1, j2]
            iszero(v) && continue
            C[_cell_var_index(lattice, (b1, j1), key[1].band, ca),
              _cell_var_index(lattice, (b2, j2), key[2].band, cb)] += v
        end
    end
    return GM.AnomalousPairing(C)
end

"PartialIF 组装（原地）：expmult!(gmps, fullfactor, PartialExp)，逐 kernel.bands
装配（单带/BCS = 就地装配 gmps；多带 normal = 逐带从 copy(gmps) 起步）；乘法
压缩截断 = alg.trunc"
function _partialif_state(gmps, lattice::AbstractGrassmannLattice, kernel::HybridizationKernel,
                          alg::PartialIF)
    exp = GM.PartialExp(alg_mult=SVDCompression(trunc=alg.trunc))
    if length(kernel.bands) == 1 || _is_orbital_kernel(kernel)
        state, err = GM.expmult!(gmps, _fullfactor(kernel, lattice), exp)
        return state, [err]
    end
    states = [GM.expmult!(_copy_mps(gmps),
                          _fullfactor(_band_kernel(kernel, b), lattice), exp)
              for b in kernel.bands]
    return tuple(first.(states)...), map(last, states)
end

"""
    bcs_hybriddynamics(lattice, bath::AbstractBCSBath, alg=PartialIF(); orbital=1, trunc) -> GMPS

BCS bath 的 IF 链便捷入口：`orbital` 选 Nambu 对（band1 = 2·orbital−1、
band2 = 2·orbital），对称性守卫 :z2（D7）。kernel 构造 =
`correlationfunction(bath, lattice; orbital)`，装配走 fullfactor 轨道路径。
"""
function bcs_hybriddynamics(lattice::AbstractGrassmannLattice, bath::AbstractBCSBath,
                            alg::InfluenceFunctionalAlgorithm=PartialIF();
                            orbital::Int=1)
    alg isa PartialIF || throw(ArgumentError(
        "BCS IF supports only PartialIF on the GM stack (D7), got $(typeof(alg))"))
    kernel = correlationfunction(bath, lattice; orbital)
    return hybriddynamics(lattice, kernel, alg)
end
