# partialif.jl：PartialIF——整形式幂零指数装配（plan/if_expmult W3 语义换代）
#
# 内核从「逐行秩一因子积 + 乘间压缩」（旧 partial_row 路径，influence.jl
# 冻结档案）换代为「AnomalousPairing 整形式 × GM expmult! PartialExp 幂零
# 指数 + 压缩」。两者是同一 IF 算符 e^{Σ ā Δ a} 的两种精确表示（偶元互相
# 超交换 ⟹ 逐行因子积 ≡ 整形式单指数）；截断 = 压缩（SVDCompression），
# 与旧乘间压缩同档。验收闸门 = 每因子 dense 对拍 vs 旧链 + ExactTTI 互换
# 对拍（tests/gates/）。
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

"""
    fullfactor(kernel, lattice; band=nothing) -> GM.AnomalousPairing

IF 指数的 2L×2L 一般二次型系数矩阵（L = 链站数；块布局与变量序规范表见
文件头）。`band` = normal 单带装配（blocks 按 `(b1,band),(b2,band)` 键取值，
腿味 (ā, a)）；`band = nothing` = 轨道全域装配（BCS kernel——同带键按正常
象限、跨带键按反常象限（ud/du）落入 C 的对应扇区；腿 conj 由 `_bcs_leg_conjs`
味表回推）。块值 = 构造期符号收编的指数系数（correlationfunction.jl 存储
约定），直读不补号。
"""
function fullfactor(kernel::HybridizationKernel, lattice::AbstractGrassmannLattice;
                    band::Union{Nothing,Int}=nothing)
    L = length(lattice) ÷ 2
    C = zeros(ComplexF64, 2L, 2L)
    for (key, blk) in kernel.blocks
        if band === nothing
            (1 <= key[1].band <= lattice.bands && 1 <= key[2].band <= lattice.bands) ||
                throw(BoundsError(1:lattice.bands, (key[1].band, key[2].band)))
            ca, cb = _bcs_leg_conjs(key[1].band, key[2].band)
        else
            (key[1].band == band && key[2].band == band) || continue
            ca, cb = true, false                       # normal 键 = ā 行 × a 列
        end
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

"PartialIF 组装链：vacuumstate → expmult!(fullfactor, PartialExp)"
function _partialif_state(lattice::AbstractGrassmannLattice, kernel::HybridizationKernel;
                          band::Union{Nothing,Int}, trunc::Truncation)
    factor = fullfactor(kernel, lattice; band)
    state, err = GM.expmult!(vacuumstate(lattice), factor,
                             GM.PartialExp(alg_mult=SVDCompression(trunc=trunc)))
    return state, [err]
end
