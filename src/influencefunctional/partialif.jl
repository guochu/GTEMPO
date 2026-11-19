# partialif.jl：PartialIF——整形式幂零指数装配（plan/if_expmult W3 语义换代）
#
# 内核从「逐行秩一因子积 + 乘间压缩」（旧 influence.jl 逐行路径，已随
# W1/W3 换代 + 2026-11-18 退役删除）换代为「AnomalousPairing 整形式 × GM
# expmult! PartialExp 幂零指数 + 压缩」。两者是同一 IF 算符 e^{Σ ā Δ a} 的
# 两种精确表示（偶元互相超交换 ⟹ 逐行因子积 ≡ 整形式单指数）；截断 =
# 压缩（SVDCompression）。
#
# fullfactor 的 C 矩阵约定（GM AnomalousPairing，变量 α ↔ 腿 α：奇 = η（a）、
# 偶 = η̄（ā），站点 = (α+1)÷2）：**有向单一取向存储、不存镜像**——腿味 =
# ContourBlockIndex.conj（构造期写定，键自完备——normal 同带键 = (ā, a)、
# BCS 反常键 = (ā↑, ā↓) / (a↓, a↑)，装配端按键直读、无味旗标）：
#   normal（band 键）：C[2·st(ā_row), 2·st(a_col) − 1] = blocks 值（ηη̄ 扇区；
#   deta 折叠 C[2i−1,2j] − C[2j,2i−1] 把唯一取向展开为全单项式系数）；
#   BCS（轨道带对键）：同带键 = 正常象限（同上式）；band1→band2 键 =
#   C[2·st(ā_b1), 2·st(ā_b2)] = ud（η̄η̄ 扇区）；band2→band1 键 =
#   C[2·st(a_b2), 2·st(a_b1)] = du（ηη 扇区）。
# U1/U1U1：normal 仅 ηη̄ 条目（荷中性）⟹ 可物化；BCS 同味扇区带净荷 ⟹
# 仅 Z2/Z2Z2（D7 守卫在 kernel 构造端已断言）。

"kernel 是否含跨带键（band_i ≠ band_j 的键——跨带耦合核或 BCS 轨道核）"
_has_cross_band_keys(kernel::HybridizationKernel) =
    any(key[1].band != key[2].band for key in keys(kernel.blocks))

"kernel 是否纯 normal 取向（所有键 = 行 ľ conj=true、列 conj=false——ExactTTI
三通道因子构造的核形状契约；反常键 = BCS 专属，D7 PartialIF-only）"
_is_normal_conj(kernel::HybridizationKernel) =
    all(key -> key[1].conj && !key[2].conj, keys(kernel.blocks))

"""
    fullfactor(kernel, lattice) -> GM.AnomalousPairing

IF 指数的一般二次型系数（**紧凑站点表**）：C 覆盖 kernel 作用域的全部变量对
（2M×2M，M = 作用域站数），`sites[k]` = 第 k 对变量的链站（严格递增，GM
AnomalousPairing 契约）——lattice 可包含比 kernel 更多的 bands（核外带 =
观测站，不进因子）。变量 α ↔ (sites[(α+1)÷2], conj = iseven(α))。

腿味 = 核键的 `ContourBlockIndex.conj`（correlationfunction 构造期写定）：
同带键 = (ā, a) 法向取向 `C[2·st(ā), 2·st(a)−1]`；跨带键按其键上味落象限
（法向有向耦合 (ā, a) 或 BCS 反常 (ā, ā) / (a, a)）。块值 = 构造期符号收编
的指数系数（correlationfunction.jl 存储约定），直读不补号。
"""
function fullfactor(kernel::HybridizationKernel, lattice::AbstractGrassmannLattice)
    return _fullfactor(kernel, lattice)
end

function _fullfactor(kernel::HybridizationKernel, lattice::AbstractGrassmannLattice)
    B = lattice.bands
    # kernel 作用域变量槽：(链站, band)，站升序（GM sites 严格递增契约）
    slots = sort!(collect(((cell_block(lattice, cell) - 1) * B + b, b)
                           for cell in _storage_cells(lattice) for b in kernel.bands))
    M = length(slots)
    slotidx = Dict(slots[k] => k for k in 1:M)      # (站, band) → 紧凑变量对号
    C = zeros(ComplexF64, 2M, 2M)
    for (key, blk) in kernel.blocks
        ca, cb = key[1].conj, key[2].conj           # 腿味在键上（构造期写定）
        br1, br2 = key[1].branch, key[2].branch
        ba, bb = key[1].band, key[2].band
        for j1 in 1:_branch_len(lattice, br1), j2 in 1:_branch_len(lattice, br2)
            v = blk[j1, j2]
            iszero(v) && continue
            st1 = (cell_block(lattice, (br1, j1)) - 1) * B + ba
            st2 = (cell_block(lattice, (br2, j2)) - 1) * B + bb
            # 变量号 = 2·slot − 1 + Int(conj)（奇 = a/η、偶 = ā/η̄）
            C[2 * slotidx[(st1, ba)] - 1 + Int(ca),
              2 * slotidx[(st2, bb)] - 1 + Int(cb)] += v
        end
    end
    return GM.AnomalousPairing(C, [first(s) for s in slots])
end

"PartialIF 组装（原地）：expmult!(gmps, fullfactor, PartialExp)，逐 kernel.bands
装配（单带/BCS = 就地装配 gmps；多带 normal = 逐带从 copy(gmps) 起步）；乘法
压缩截断 = alg.trunc"
function _partialif_state(gmps, lattice::AbstractGrassmannLattice, kernel::HybridizationKernel,
                          alg::PartialIF)
    exp = GM.PartialExp(alg_mult=SVDCompression(trunc=alg.trunc))
    if length(kernel.bands) == 1 || _has_cross_band_keys(kernel)
        # 单带 / BCS 轨道 / 一般跨带耦合核：二次型不可按带分解 ⟹ 单链整体装配
        state, err = GM.expmult!(gmps, _fullfactor(kernel, lattice), exp)
        return state, [err]
    end
    states = [GM.expmult!(_copy_mps(gmps),
                          _fullfactor(_band_kernel(kernel, b), lattice), exp)
              for b in kernel.bands]
    return tuple(first.(states)...), map(last, states)
end

"""
    bcs_hybriddynamics(lattice, bath::AbstractBCSBath, alg=PartialIF(); orbital=1) -> GMPS

BCS bath 的 IF 链便捷入口：`orbital` 选 Nambu 对（band1 = 2·orbital−1、
band2 = 2·orbital），对称性守卫 :z2（D7）。kernel 构造 =
`correlationfunction(bath, lattice; orbital)`，装配走 fullfactor 轨道路径
（反常腿味在核键上，与 normal 同一装配函数）。
"""
function bcs_hybriddynamics(lattice::AbstractGrassmannLattice, bath::AbstractBCSBath,
                            alg::InfluenceFunctionalAlgorithm=PartialIF();
                            orbital::Int=1)
    alg isa PartialIF || throw(ArgumentError(
        "BCS IF supports only PartialIF on the GM stack (D7), got $(typeof(alg))"))
    kernel = correlationfunction(bath, lattice; orbital)
    return first(hybriddynamics(lattice, kernel, alg))
end
