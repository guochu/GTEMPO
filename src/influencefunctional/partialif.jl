# partialif.jl：PartialIF——整形式幂零指数装配（plan/if_expmult W3 语义换代）
#
# 内核从「逐行秩一因子积 + 乘间压缩」（旧 influence.jl 逐行路径，已随
# W1/W3 换代 + 2026-11-18 退役删除）换代为「整形式因子 × GM expmult! PartialExp
# 幂零指数 + 压缩」。两者是同一 IF 算符 e^{Σ ā Δ a} 的两种精确表示（偶元互相
# 超交换 ⟹ 逐行因子积 ≡ 整形式单指数）；截断 = 压缩（SVDCompression）。
#
# fullfactor 的因子路由（2026-11-19 用户裁决，按 kernel.anomalous 派生旗标）：
# - **normal（荷守恒）→ GM.Pairing**：矩形站点表（变量槽站号 = 带解析链站，
#   天然互异），C = M×M 全矩阵（行 = ā 槽、列 = a 槽；前后向/同站三路条目
#   由 Pairing 构造器分派，接触项折 hloc 角）——含多带带对角与跨带法向耦合
#   的一般 Δ_b1b2 全部条目，**单因子覆盖全部变量**（多带不再逐带拆分）；
# - **anomalous（异味通道）→ GM.AnomalousPairing**：紧凑变量表（2M×2M，
#   变量 α ↔ (sites[(α+1)÷2], conj = iseven(α))）——四种味组合全表达：
#   normal 键 = (ā, a) 象限（deta 折叠展开全单项式系数）；BCS 反常键 =
#   (ā↑, ā↓) ud / (a↓, a↑) du 同味象限。
# 腿味 = ContourBlockIndex.conj（构造期写定，键自完备），装配端按键直读。
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
    fullfactor(kernel, lattice) -> GM.Pairing | GM.AnomalousPairing

IF 指数的一般二次型因子（**单因子覆盖 kernel.bands 全部变量**，多带不拆分）。
lattice 可包含比 kernel 更多的 bands（核外带 = 观测站，不进因子）。路由按
`kernel.anomalous`（构造器从键派生）：

- **normal（荷守恒）→ `GM.Pairing(C, sites, sites; conja=true, conjb=false)`**：
  站点表 = 变量槽链站（带解析站号，严格递增互异），C = M×M 全矩阵——
  C[st(ā_row), st(a_col)] = blocks 值（行 = ā 腿、列 = a 腿；含多带带对角与
  一般跨带 Δ_b1b2 全部条目）；
- **anomalous（异味）→ `GM.AnomalousPairing(C, sites)`**：紧凑变量表
  （2M×2M），变量 α ↔ (sites[(α+1)÷2], conj = iseven(α))——BCS 反常键
  (ā↑, ā↓) / (a↓, a↑) 落同味象限。

块值 = 构造期符号收编的指数系数（correlationfunction.jl 存储约定），直读
不补号。non-anomalous 核若带非 canonical 键（行腿非 ā）响抛（Pairing 单取向
全局味旗无法混合表达）。
"""
function fullfactor(kernel::HybridizationKernel, lattice::AbstractGrassmannLattice)
    return _fullfactor(kernel, lattice)
end

function _fullfactor(kernel::HybridizationKernel, lattice::AbstractGrassmannLattice)
    B = lattice.bands
    # kernel 作用域变量槽：(链站, band)，站升序（带解析站号 = (块号−1)·B + band
    # 天然互异——Pairing 站点表唯一性契约自动满足）
    slots = sort!(collect(((cell_block(lattice, cell) - 1) * B + b, b)
                           for cell in _storage_cells(lattice) for b in kernel.bands))
    M = length(slots)
    slotidx = Dict(slots[k] => k for k in 1:M)      # (站, band) → 紧凑变量槽号
    stations = [first(s) for s in slots]
    T = TK.scalartype(kernel)                       # 核标量型动态（实谱核 = Float64）
    if kernel.anomalous
        # 异味通道：紧凑变量表 AnomalousPairing（变量号 = 2·槽 − 1 + Int(conj)）
        C = zeros(T, 2M, 2M)
        for (key, blk) in kernel.blocks
            ca, cb = key[1].conj, key[2].conj       # 腿味在键上（构造期写定）
            br1, br2 = key[1].branch, key[2].branch
            ba, bb = key[1].band, key[2].band
            for j1 in 1:_branch_len(lattice, br1), j2 in 1:_branch_len(lattice, br2)
                v = blk[j1, j2]
                iszero(v) && continue
                st1 = (cell_block(lattice, (br1, j1)) - 1) * B + ba
                st2 = (cell_block(lattice, (br2, j2)) - 1) * B + bb
                C[2 * slotidx[(st1, ba)] - 1 + Int(ca),
                  2 * slotidx[(st2, bb)] - 1 + Int(cb)] += v
            end
        end
        return GM.AnomalousPairing(C, stations)
    end
    # normal（荷守恒）：矩形站点表 Pairing——C = 槽全矩阵（行 ā、列 a）
    _is_normal_conj(kernel) || throw(ArgumentError(
        "non-anomalous kernel must carry canonical (ā, a) keys for the Pairing route"))
    C = zeros(T, M, M)
    for (key, blk) in kernel.blocks
        br1, br2 = key[1].branch, key[2].branch
        ba, bb = key[1].band, key[2].band
        for j1 in 1:_branch_len(lattice, br1), j2 in 1:_branch_len(lattice, br2)
            v = blk[j1, j2]
            iszero(v) && continue
            st1 = (cell_block(lattice, (br1, j1)) - 1) * B + ba
            st2 = (cell_block(lattice, (br2, j2)) - 1) * B + bb
            C[slotidx[(st1, ba)], slotidx[(st2, bb)]] += v
        end
    end
    return GM.Pairing(C, stations, stations; conja=true, conjb=false)
end

"PartialIF 组装（原地）：expmult!(gmps, fullfactor, PartialExp)——单因子覆盖
kernel.bands 全部变量（恒就地单链装配，多带逐带拆分契约退役）；乘法压缩
截断 = alg.trunc"
function _partialif_state(gmps, lattice::AbstractGrassmannLattice, kernel::HybridizationKernel,
                          alg::PartialIF)
    exp = GM.PartialExp(alg_mult=SVDCompression(trunc=alg.trunc))
    state, err = GM.expmult!(gmps, _fullfactor(kernel, lattice), exp)
    return state, [err]
end

"PolarIF 组装（原地）：expmult!(gmps, fullfactor, PolarExp)——对角剥离/味扇区
折叠后逐秩一分量 mult!（精确因子化），乘法压缩截断 = alg.trunc；返回
(state, [truncerr])（与 _partialif_state 同形）"
function _polarif_state(gmps, lattice::AbstractGrassmannLattice, kernel::HybridizationKernel,
                        alg::PolarIF)
    exp = GM.PolarExp(rtol=alg.rtol, alg_mult=SVDCompression(trunc=alg.trunc))
    state, err = GM.expmult!(gmps, _fullfactor(kernel, lattice), exp)
    return state, [err]
end

"""
    bcs_hybriddynamics(lattice, bath::AbstractBCSBath, alg=PartialIF(); orbital=nothing) -> GMPS

BCS bath 的 IF 链便捷入口：`alg` ∈ {PartialIF, PolarIF}（D7：GM 栈仅此二），
`orbital` 选 Nambu 对（band1 = 2·orbital−1、
band2 = 2·orbital），缺省 nothing = 全部 Nambu 对共享同一套杂化函数
（correlationfunction 的轨道模式），对称性守卫 :z2（D7）。kernel 构造 =
`correlationfunction(bath, lattice; orbital)`，装配走 fullfactor 轨道路径
（反常腿味在核键上，与 normal 同一装配函数）。
"""
function bcs_hybriddynamics(lattice::AbstractGrassmannLattice, bath::AbstractBCSBath,
                            alg::InfluenceFunctionalAlgorithm=PartialIF();
                            orbital::Union{Nothing,Int}=nothing)
    alg isa Union{PartialIF,PolarIF} || throw(ArgumentError(
        "BCS IF supports only PartialIF/PolarIF on the GM stack (D7), got $(typeof(alg))"))
    kernel = correlationfunction(bath, lattice; orbital)
    return first(hybriddynamics(lattice, kernel, alg))
end
