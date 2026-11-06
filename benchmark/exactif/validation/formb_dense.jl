# =============================================================================
# form_b/formb_dense.jl —— 形态 B（消元式单扇区）稠密系数级验证
#
# 把 IFGTEMPO 两扇区 lattice 中的 IF Grassmann 变量从多项式中精确积掉：
#
#   两扇区多项式   P(a,x) = ν·exp(½aᵀGa) · ∏_t exp(−ā_t x_t − x̄_t a_t) · K[x] · M₀
#   --∫Da-->       P_elim(x) = ν·Pf(G) · exp(Σ_{λ<λ'} D[λ,λ'] θ_λ θ_λ') · K[x] · M₀
#
# 其中 D（Δ_eff 的反对称多项式形式）只作用在接口轨迹腿 x_t± / x̄_t± 上，
# 层 0 与开放腿无源（README"边界层全零"）。接口测度被消元吸收后，收缩中
# 不得再出现任何内部自耦合测度（README"测度簿记"③）；layer 0 的 M₀ 保留。
#
# 验证协议（全部在相干态表示内端到端，勿用占据基 path values）：
#   Step A  稠密两扇区副本：Z == 1 且 ρ(T) == Fock 离散 Trotter 链
#   Step B  消元代数：暴力 ∫Da == ν·Pf(G)·exp(D)（系数级）
#   Step C  消元式收缩：Z == 1、ρ(T)/⟨n(T)⟩ == Fock（N=1..4 × order 1/2/4）
#   Step D  两自旋 + IFGTEMPO GMPS 交叉验证（gtempo_density_matrix_direct）
#
# 约定来源（全部照搬 IFGTEMPO 工作代码，不新造）：
#   · lattice 布局 = timelocallattice（单自旋取其单 flavor 投影，Step D 断言）
#   · kernel 单项式符号 = impurity_kernel_gmps 配方：
#       coeff(bra,ket) = (−1)^{C(count_ones(ket),2)}·K[bra+1,ket+1]，local 序 (bar…,eta…)
#   · 接口测度 = IF_DERIVATION_AND_PACKAGE_LOGIC.md §9：
#       M_t = exp[−ā_t⁺x_t⁺ − x̄_t⁺a_t⁺ − ā_t⁻x_t⁻ − x̄_t⁻a_t⁻]（逐相邻对）
#   · 算符链 = §2：U = A(link_{Nc+1})·Q_{Nc}·A(link_{Nc})···Q_1·A(link_1)
#   · 高斯积分恒等式（n=2 直算核验）：
#       ∫exp(½aᵀGa + jᵀa) = Pf(G)·exp(−½jᵀG⁻¹j)
#     源 j_i = −partner_i 给出递增 monomial θ_λθ_λ' 的系数
#       D[λ,λ'] = −½(G⁻¹[i,k] − G⁻¹[k,i]) = −Ā[i,k]（i,k = λ,λ' 的 IF 槽位）
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/form_b/formb_dense.jl
# =============================================================================
using LinearAlgebra, Printf
using GTEMPO

const IFT_ROOT = joinpath(dirname(dirname(dirname(pathof(GTEMPO)))), "IFGTEMPO")
isdir(IFT_ROOT) || error("IFGTEMPO 仓库未找到：$IFT_ROOT")
push!(LOAD_PATH, IFT_ROOT)
using IFGTEMPO

# ---------------------------------------------------------------------------
# 稠密 Grassmann 多项式机器
#
# 表示：Vector{ComplexF64}(2^L)，index 的 bit s ↔ 生成元 θ_s（site s），
# 单项式恒按 site 递增序存放。因子以"局部顺序单项式表"给出，嵌入时自动
# 计入排列符号。
# ---------------------------------------------------------------------------

const Monomials = Vector{Tuple{Vector{Int},ComplexF64}}

"""局部顺序单项式表的逆序数（嵌入时的排列符号）。"""
function mono_inversions(sites::Vector{Int})
    c = 0
    for i in eachindex(sites), j in (i + 1):length(sites)
        sites[i] > sites[j] && (c += 1)
    end
    return c
end

"""把局部顺序单项式表展开为稠密多项式（site 递增约定）。空 site 表 = 标量。"""
function embed(L::Int, monos::Monomials)
    out = zeros(ComplexF64, 1 << L)
    for (sites, coeff) in monos
        all(0 .<= sites .< L) || error("site 越界：$sites (L=$L)")
        length(unique(sites)) == length(sites) || error("重复 site：$sites")
        sgn = isodd(mono_inversions(sites)) ? -1 : 1
        mask = 0
        for s in sites
            mask |= 1 << s
        end
        out[mask + 1] += sgn * coeff
    end
    return out
end

"""合并符号：递增序 monomial(smasks) × 递增序 monomial(sub) → 递增序乘积。"""
function merge_sign(smask::UInt64, sub::UInt64)
    s = 0
    m = smask
    while m != 0
        b = m & -m
        s += count_ones(sub & (b - one(UInt64)))
        m &= m - one(UInt64)
    end
    return isodd(s) ? -1 : 1
end

"""稠密多项式 × 单项式因子（稀疏 × 稠密，就地写入 out，O(#monos·2^L)）。"""
function mul_monomials!(out::Vector{ComplexF64}, dense::Vector{ComplexF64},
                        L::Int, monos::Monomials)
    fill!(out, zero(ComplexF64))
    total = (one(UInt64) << L) - one(UInt64)
    for (sites, coeff) in monos
        all(0 .<= sites .< L) || error("site 越界：$sites (L=$L)")
        length(unique(sites)) == length(sites) || error("重复 site：$sites")
        smask = UInt64(0)
        for s in sites
            smask |= one(UInt64) << s
        end
        base = coeff * (isodd(mono_inversions(sites)) ? -1 : 1)
        comp = ~smask & total
        sub = comp
        while true
            @inbounds out[(sub | smask) + 1] +=
                base * merge_sign(smask, sub) * dense[sub + 1]
            sub == 0 && break
            sub = (sub - one(UInt64)) & comp
        end
    end
    return out
end

"""顺序 Berezin 积分：按 site 序逐个积掉（∫dθ θ·f → step_sign·f）。

`order` 取 :inc（site 从小到大）或 :dec。这是整体定向的约定旋钮，
由 Step A 对 Z=1/ρ=Fock 的匹配锁定后全程使用。"""
function integrate_out(dense::Vector{ComplexF64}, L::Int, sites::Vector{Int};
                       step_sign::Int=1, order::Symbol=:inc)
    cur = copy(dense)
    seq = order === :dec ? reverse(sites) : sites
    for s in seq
        nxt = zeros(ComplexF64, 1 << L)
        bit = one(UInt64) << s
        @inbounds for idx in 1:length(cur)
            m = UInt64(idx - 1)
            m & bit == 0 && continue
            nxt[m ⊻ bit + 1] += step_sign * cur[idx]
        end
        cur = nxt
    end
    return cur
end

"""多模式 kernel（复矩阵 K[bra+1,ket+1]）→ 单项式表（impurity_kernel_gmps 配方）。

local 顺序：(bar(mode1..), eta(mode1..))；
coeff(bra,ket) = (−1)^{C(count_ones(ket),2)}·K[bra+1,ket+1]。"""
function kernel_monomials(barsites::Vector{Int}, etasites::Vector{Int},
                          K::AbstractMatrix)
    n = length(barsites)
    (n == length(etasites) && size(K) == (1 << n, 1 << n)) ||
        throw(DimensionMismatch("kernel 矩阵与腿数不匹配：$n / $(size(K))"))
    monos = Monomials()
    for bra in 0:((1 << n) - 1), ket in 0:((1 << n) - 1)
        p = count_ones(ket)
        coeff = (isodd(p * (p - 1) ÷ 2) ? -1 : 1) * K[bra + 1, ket + 1]
        iszero(coeff) && continue
        loc = Int[]
        for i in 0:(n - 1)
            (bra >> i) & 1 == 1 && push!(loc, barsites[i + 1])
        end
        for i in 0:(n - 1)
            (ket >> i) & 1 == 1 && push!(loc, etasites[i + 1])
        end
        push!(monos, (loc, ComplexF64(coeff)))
    end
    return monos
end

# ---------------------------------------------------------------------------
# lattice 注册表：与 timelocallattice 完全一致的 site 布局
# （nflav=2 逐位一致，Step D 中断言；nflav=1 为其单 flavor 投影）
# ---------------------------------------------------------------------------

struct Var
    typ::Symbol   # :imp / :IF
    t::Int
    br::Symbol    # :+ / :-
    stat::Symbol  # :bar / :eta
    f::Int
end

function build_sites(N::Int, nflav::Int)
    vars = Var[]
    for f in 1:nflav
        push!(vars, Var(:imp, 0, :+, :bar, f)); push!(vars, Var(:imp, 0, :+, :eta, f))
    end
    for f in 1:nflav
        push!(vars, Var(:imp, 0, :-, :bar, f)); push!(vars, Var(:imp, 0, :-, :eta, f))
    end
    for t in 1:N, f in 1:nflav
        push!(vars, Var(:IF, t, :+, :bar, f)); push!(vars, Var(:imp, t, :+, :eta, f))
        push!(vars, Var(:imp, t, :+, :bar, f)); push!(vars, Var(:IF, t, :+, :eta, f))
        push!(vars, Var(:IF, t, :-, :bar, f)); push!(vars, Var(:imp, t, :-, :eta, f))
        push!(vars, Var(:imp, t, :-, :bar, f)); push!(vars, Var(:IF, t, :-, :eta, f))
    end
    for f in 1:nflav
        push!(vars, Var(:imp, N + 1, :+, :bar, f)); push!(vars, Var(:imp, N + 1, :-, :eta, f))
    end
    smap = Dict(v => i - 1 for (i, v) in enumerate(vars))   # 0-based site
    return smap, vars, length(vars)
end

"""压缩映射：只保留给定类型的变量，按原相对顺序重编号 0..n−1。"""
function compact_map(smap::Dict{Var,Int}, vars::Vector{Var}, keep::Function)
    cvars = filter(keep, vars)
    cmap = Dict(v => i - 1 for (i, v) in enumerate(cvars))
    return cmap, length(cvars)
end

"""IF 变量的 G 矩阵槽位（1-based）：flat = 4(t−1)+slot+1，slot 0..3 =
(+bar, +eta, −bar, −eta) —— influence_matrix_G 的约定。"""
gslot(t::Int, br::Symbol, stat::Symbol) =
    4 * (t - 1) + (br === :+ ? 1 : 3) + (stat === :bar ? 0 : 1)

if_var(f::Int, i::Int) =   # G 槽位 i（1-based）→ IF 变量
    Var(:IF, div(i - 1, 4) + 1, mod(i - 1, 4) < 2 ? :+ : :-,
        iseven(mod(i - 1, 4)) ? :bar : :eta, f)

"""IF slot (t,br,stat) 的测度 partner（IMP 腿，bar↔eta 交叉）。"""
partner_var(t::Int, br::Symbol, stat::Symbol, f::Int) =
    Var(:imp, t, br, stat === :bar ? :eta : :bar, f)

# ---------------------------------------------------------------------------
# 因子构造（全部只依赖 site 映射）
# ---------------------------------------------------------------------------

"""把共享首变量的双线性对 (site_i < site_j, coeff) 按"行"分组成独立因子。

exp(Σ 双线性) = ∏_i (1 + θ_i Σ_{j>i} c_ij θ_j)（行首唯一 + 幂零性 ⇒ 精确，
get_if_gmps 的逐行 partialmpo 同款技巧）。每个行因子含常数 1。"""
function bilinear_row_factors(pairs::Vector{Tuple{Int,Int,ComplexF64}})
    byrow = Dict{Int,Vector{Tuple{Int,ComplexF64}}}()
    for (i, j, c) in pairs
        push!(get!(byrow, i, Tuple{Int,ComplexF64}[]), (j, c))
    end
    factors = Monomials[]
    for i in sort!(collect(keys(byrow)))
        monos = Monomials([([], 1.0 + 0im)])
        for (j, c) in byrow[i]
            push!(monos, ([i, j], c))
        end
        push!(factors, monos)
    end
    return factors
end

"""IF Gaussian：exp(Σ_{i<j} G[i,j] θ_i θ_j)（某 flavor），逐行因子化。
返回 (因子列表, 归一化 ν = ZN/Z)。"""
function if_gaussian_factors(smap, data, N::Int, f::Int)
    pairs = Tuple{Int,Int,ComplexF64}[]
    for i in 1:(4N), j in (i + 1):(4N)
        v = data.G[i, j]
        iszero(v) && continue
        push!(pairs, (smap[if_var(f, i)], smap[if_var(f, j)], v))
    end
    return bilinear_row_factors(pairs), data.normalization
end

"""接口交叉测度（§9）：每 (接口, flavor, identity) 一个独立因子 exp(−·) = (1 − ·)。

M_t = exp[−ā_t⁺x_t⁺ − x̄_t⁺a_t⁺ − ā_t⁻x_t⁻ − x̄_t⁻a_t⁻]——四个恒等互相独立，
不同恒等的多项式必须能相乘，故逐恒等拆分（不可合并成单因子表）。"""
function interface_measure_factors(smap, N::Int, nflav::Int)
    factors = Monomials[]
    for t in 1:N, f in 1:nflav, br in (:+, :-)
        push!(factors, Monomials([([], 1.0 + 0im),
                                  ([smap[Var(:IF, t, br, :bar, f)],
                                    smap[Var(:imp, t, br, :eta, f)]], -1.0 + 0im)]))
        push!(factors, Monomials([([], 1.0 + 0im),
                                  ([smap[Var(:imp, t, br, :bar, f)],
                                    smap[Var(:IF, t, br, :eta, f)]], -1.0 + 0im)]))
    end
    return factors
end

"""layer 0 测度 M₀：每 (flavor, branch) 的 resolution 独立因子 exp[−x̄₀±x₀±]。"""
function layer0_measure_factors(smap, nflav::Int)
    factors = Monomials[]
    for f in 1:nflav, br in (:+, :-)
        push!(factors, Monomials([([], 1.0 + 0im),
                                  ([smap[Var(:imp, 0, br, :bar, f)],
                                    smap[Var(:imp, 0, br, :eta, f)]], -1.0 + 0im)]))
    end
    return factors
end

"""K 因子：初态 ρ₀ + 前向/后向链 kernel，每个 kernel 一个独立因子。

单 flavor：kernel 矩阵 2×2；两 flavor：4×4（mode 序 f=1,2）。"""
function k_kernel_factors(smap, N::Int, nflav::Int, links::Vector{Float64},
                          Him::AbstractMatrix, rho0::AbstractMatrix)
    factors = Monomials[]
    bars(t, br) = [smap[Var(:imp, t, br, :bar, f)] for f in 1:nflav]
    etas(t, br) = [smap[Var(:imp, t, br, :eta, f)] for f in 1:nflav]
    nc = length(links) - 1
    # 初态密度：⟨x̄₀⁺|ρ₀|x₀⁻⟩（fwd bar × bwd eta）
    push!(factors, kernel_monomials(bars(0, :+), etas(0, :-), rho0))
    # 前向链 link t：⟨x̄_t⁺|A|x_{t−1}⁺⟩
    for t in 1:(nc + 1)
        push!(factors, kernel_monomials(bars(t, :+), etas(t - 1, :+),
                                        exp((-1im * links[t]) .* Him)))
    end
    # 后向链 link t：⟨x̄_{t−1}⁻|A†|x_t⁻⟩（开放端往回收缩）
    for t in (nc + 1):-1:1
        push!(factors, kernel_monomials(bars(t - 1, :-), etas(t, :-),
                                        exp((+1im * links[t]) .* Him)))
    end
    return factors
end

"""迹闭合因子：每 flavor 一个独立 resolution (1 + closure_k·x̄⁺x⁻)。

closure_k = −1 是恒等 resolution e^{−x̄x}（配合 ⟨−x̄| 翻转）；+1 是把翻转
并入闭合（⟨x̄|→⟨−x̄| 等价于 measure 取 +x̄x）。由 Step A 校准锁定。
两 flavor 的 resolution 相互独立，必须逐 flavor 拆分（与测度同理）。"""
function closure_factors(smap, N::Int, nflav::Int, k::Int)
    factors = Monomials[]
    for f in 1:nflav
        push!(factors, Monomials([([], 1.0 + 0im),
                                  ([smap[Var(:imp, N + 1, :+, :bar, f)],
                                    smap[Var(:imp, N + 1, :-, :eta, f)]],
                                   ComplexF64(k))]))
    end
    return factors
end

"""开放腿的 ρ(T) 提取（impurity_kernel_gmps 配方的逆映射）。

bra 位 = bar 腿（mode 序），ket 位 = eta 腿；系数逆映射含 (−1)^{C(ket,2)}
与 local→global 排列符号（±1 自逆）。"""
function open_rho(integrated::Vector{ComplexF64}, smap, N::Int, nflav::Int)
    barsites = [smap[Var(:imp, N + 1, :+, :bar, f)] for f in 1:nflav]
    etasites = [smap[Var(:imp, N + 1, :-, :eta, f)] for f in 1:nflav]
    n = nflav
    dims = 1 << n
    rho = zeros(ComplexF64, dims, dims)
    for bra in 0:(dims - 1), ket in 0:(dims - 1)
        loc = Int[]
        for i in 0:(n - 1)
            (bra >> i) & 1 == 1 && push!(loc, barsites[i + 1])
        end
        for i in 0:(n - 1)
            (ket >> i) & 1 == 1 && push!(loc, etasites[i + 1])
        end
        mask = 0
        for s in loc
            mask |= 1 << s
        end
        p = count_ones(ket)
        sgn = (isodd(mono_inversions(loc)) ? -1 : 1) *
              (isodd(p * (p - 1) ÷ 2) ? -1 : 1)
        rho[bra + 1, ket + 1] = sgn * integrated[mask + 1]
    end
    return rho
end

"""消元后的有效双线性（Δ_eff 的反对称多项式形式），逐行因子化。

源嵌入 j_i = −partner_i；递增 monomial θ_λθ_λ' 的系数
    D[λ,λ'] = −½(G⁻¹[i,k] − G⁻¹[k,i])，
i = slot(λ), k = slot(λ')（λ,λ' 的 IF 槽位）。`elim_sign` 为校准旋钮。"""
function delta_eff_factors(smap, Ginv::AbstractMatrix, N::Int, f::Int,
                           elim_sign::Int)
    GA = (Ginv .- transpose(Ginv)) ./ 2
    legs = Tuple{Int,Int}[]   # (partner site, G 槽位)
    for t in 1:N, br in (:+, :-), stat in (:bar, :eta)
        push!(legs, (smap[partner_var(t, br, stat, f)], gslot(t, br, stat)))
    end
    sort!(legs, by = x -> x[1])   # site 递增 → monomial 序
    pairs = Tuple{Int,Int,ComplexF64}[]
    for i in eachindex(legs), k in (i + 1):length(legs)
        v = elim_sign * GA[legs[i][2], legs[k][2]]
        iszero(v) && continue
        push!(pairs, (legs[i][1], legs[k][1], v))
    end
    return bilinear_row_factors(pairs)
end

# ---------------------------------------------------------------------------
# Fock 离散链参照（同一 schedule 的精确离散演化；§2 算符链）
# ---------------------------------------------------------------------------

function annih_op(dim::Int, mode::Int)
    op = zeros(ComplexF64, dim, dim)
    for n in 0:(dim - 1)
        (n >> mode) & 1 == 1 || continue
        sgn = isodd(count_ones(n & ((1 << mode) - 1))) ? -1 : 1
        op[(n ⊻ (1 << mode)) + 1, n + 1] = sgn
    end
    return op
end

"""U = A(link_{Nc+1})·Q_{Nc}·A(link_{Nc})···Q_1·A(link_1)（右端先作用）。"""
function chain_operator(Hq::AbstractMatrix, Him::AbstractMatrix,
                        cells::Vector{Float64}, links::Vector{Float64})
    U = exp((-1im * links[end]) .* Him)
    for t in length(cells):-1:1
        U = U * exp(-1im * cells[t] .* Hq) * exp((-1im * links[t]) .* Him)
    end
    return U
end

"""单自旋离散链：ρ_disc = Tr_bath[U (ρ₀⊗ρ_b) U†]（mode 序：imp=0, bath≥1）。"""
function fock_chain_single(up::FiniteBath, eps_d::Real,
                           sched::TrotterSchedule, rho0::AbstractMatrix)
    m = bath_modes(up)
    nt = 1 + m
    dim = 1 << nt
    ops = [annih_op(dim, k) for k in 0:(nt - 1)]
    hQ = bath_one_particle_hamiltonian(up)          # impurity level 0
    Hq = zeros(ComplexF64, dim, dim)
    for i in 1:(m + 1), j in 1:(m + 1)
        iszero(hQ[i, j]) && continue
        Hq .+= hQ[i, j] .* (ops[i]' * ops[j])
    end
    Him = eps_d .* (ops[1]' * ops[1])
    Hb = zeros(ComplexF64, dim, dim)
    for i in 2:(m + 1), j in 2:(m + 1)
        iszero(up.bath_hamiltonian[i - 1, j - 1]) && continue
        Hb .+= up.bath_hamiltonian[i - 1, j - 1] .* (ops[i]' * ops[j])
    end
    # ρ_b 嵌入全空间：exp(−βH_bath) 的迹含杂质维度，须按 bath 配分函数归一
    Zb = tr(exp(-up.beta .* Hb)) / 2
    rb = exp(-up.beta .* Hb) ./ Zb
    U = chain_operator(Hq, Him, sched.cells, sched.links)
    r0 = zeros(ComplexF64, dim, dim)
    for b in 0:((1 << m) - 1), a in 0:1, a2 in 0:1
        r0[(b << 1 | a) + 1, (b << 1 | a2) + 1] = rho0[a + 1, a2 + 1]
    end
    rtot = U * (r0 * rb) * adjoint(U)
    rho = zeros(ComplexF64, 2, 2)
    for b in 0:((1 << m) - 1), a in 0:1, a2 in 0:1
        rho[a + 1, a2 + 1] += rtot[(b << 1 | a) + 1, (b << 1 | a2) + 1]
    end
    return rho
end

"""两自旋离散链（mode 序：imp_up=0, imp_dn=1, up-bath, dn-bath）。"""
function fock_chain_two_spin(up::FiniteBath, dn::FiniteBath,
                             imp::TwoSpinImpurityModel, sched::TrotterSchedule)
    mu, md = bath_modes(up), bath_modes(dn)
    nt = 2 + mu + md
    dim = 1 << nt
    ops = [annih_op(dim, k) for k in 0:(nt - 1)]
    Hq = zeros(ComplexF64, dim, dim)
    Hb = zeros(ComplexF64, dim, dim)
    for (model, off, impmode) in ((up, 2, 0), (dn, 2 + mu, 1))
        for i in 1:bath_modes(model), j in 1:bath_modes(model)
            Hq .+= model.bath_hamiltonian[i, j] .* (ops[off + i]' * ops[off + j])
            Hb .+= model.bath_hamiltonian[i, j] .* (ops[off + i]' * ops[off + j])
        end
        for i in 1:bath_modes(model)
            Hq .+= conj(model.hybridization[i]) .* (ops[impmode + 1]' * ops[off + i])
            Hq .+= model.hybridization[i] .* (ops[off + i]' * ops[impmode + 1])
        end
    end
    Him = imp.epsilon_up .* (ops[1]' * ops[1]) +
          imp.epsilon_down .* (ops[2]' * ops[2]) +
          imp.interaction_u .* (ops[1]' * ops[1]) * (ops[2]' * ops[2])
    # ρ_b 嵌入全空间，按 bath 配分函数归一（杂质维度 4）
    Zb = tr(exp(-up.beta .* Hb)) / 4
    rb = exp(-up.beta .* Hb) ./ Zb
    U = chain_operator(Hq, Him, sched.cells, sched.links)
    r0 = zeros(ComplexF64, dim, dim)
    for b in 0:((1 << (nt - 2)) - 1), a in 0:3, a2 in 0:3
        r0[(b << 2 | a) + 1, (b << 2 | a2) + 1] = imp.initial_density[a + 1, a2 + 1]
    end
    rtot = U * (r0 * rb) * adjoint(U)
    rho = zeros(ComplexF64, 4, 4)
    for b in 0:((1 << (nt - 2)) - 1), a in 0:3, a2 in 0:3
        rho[a + 1, a2 + 1] += rtot[(b << 2 | a) + 1, (b << 2 | a2) + 1]
    end
    return rho
end

# ---------------------------------------------------------------------------
# 单例收缩（两扇区 / 消元式），返回 (Z, ρ, T_tot, dt)
# ---------------------------------------------------------------------------

"""组装多项式并收缩。

mode = :twosector：IF Gaussian × 接口测度 × M₀ × K × closure（完整两扇区，
暴力稠密；只用于 N ≤ 2 的校准/对照）。
mode = :elim：ν·Pf(G)·exp(D) × M₀ × K × closure（消元式单扇区）。

约定旋钮 conv = (step_sign, order, closure_k)；ρ 一律从开放核提取
（不乘 closure、保留开放腿），Z 另经闭合收缩独立计算。"""
function dense_case(; order::Int, N::Int, nflav::Int, models, imp,
                    mode::Symbol, conv, elim_sign::Int=-1, pf_scale=1.0 + 0im)
    sched = trotter_schedule(order, N, T_TOT / N)
    # lattice 层数 = cells 数（order 4 每宏步 3 个 subcell ⇒ 3N 层）
    nl = length(sched.cells)
    smap, vars, Lfull = build_sites(nl, nflav)
    rho0 = nflav == 1 ? imp.rho0 : imp.initial_density
    Him = nflav == 1 ? ComplexF64[0 0; 0 imp.eps] : impurity_hamiltonian(imp)

    if mode === :twosector
        L = Lfull
        stack = Monomials[]
        scale = 1.0 + 0im
        for f in 1:nflav
            factors, ν = if_gaussian_factors(smap, models[f].ifdata, N, f)
            append!(stack, factors)
            scale *= ν
        end
    elseif mode === :elim
        # 压缩到 x 腿（跳过 IF 变量）
        cmap, L = compact_map(smap, vars, v -> v.typ === :imp)
        smap = cmap
        stack = Monomials[]
        scale = pf_scale
        for f in 1:nflav
            append!(stack, delta_eff_factors(smap, models[f].Ginv, nl, f, elim_sign))
        end
    else
        error("mode 必须是 :twosector 或 :elim")
    end
    if mode === :twosector
        append!(stack, interface_measure_factors(smap, nl, nflav))
    end
    append!(stack, layer0_measure_factors(smap, nflav))
    append!(stack, k_kernel_factors(smap, nl, nflav, sched.links, Him, rho0))

    dense = embed(L, stack[1])
    dense .*= scale
    buf = similar(dense)
    for k in 2:length(stack)
        mul_monomials!(buf, dense, L, stack[k])
        dense, buf = buf, dense
    end

    # ρ(T)：保留开放腿，积掉其余
    openvars = vcat([Var(:imp, nl + 1, :+, :bar, f) for f in 1:nflav],
                    [Var(:imp, nl + 1, :-, :eta, f) for f in 1:nflav])
    keep = [smap[v] for v in openvars]
    intsites = [s for s in 0:(L - 1) if all(s != k for k in keep)]
    poly = integrate_out(dense, L, intsites; step_sign = conv.step_sign,
                         order = conv.order)
    rho = open_rho(poly, smap, nl, nflav)

    # Z：加闭合因子（逐 flavor resolution）后全收缩
    zclosures = closure_factors(smap, nl, nflav, conv.closure_k)
    Zbuf = similar(dense)
    Zpoly = dense
    for zc in zclosures
        mul_monomials!(Zbuf, Zpoly, L, zc)
        Zpoly, Zbuf = Zbuf, Zpoly
    end
    Z = integrate_out(Zpoly, L, collect(0:(L - 1)); step_sign = conv.step_sign,
                      order = conv.order)[1]
    return (Z = Z, rho = rho, T = T_TOT, dt = T_TOT / N)
end

# ---------------------------------------------------------------------------
# 消元代数的暴力对照（Step B）与 Pf_mine
# ---------------------------------------------------------------------------

"""Pf_mine(G)：exp(Σ_{i<j} G[i,j] a_i a_j) 在给定定向下的积分（数值 Pfaffian），
含归一化 ν 乘子。`nl` = 接口层数（= length(cells)，order 4 为 3N）；
定向旋钮与消元/收缩一致。"""
function pf_mine(data, nl::Int, f::Int; step_sign::Int, order::Symbol, ν=1.0 + 0im)
    nif = 4nl
    smap = Dict(if_var(f, i) => i - 1 for i in 1:nif)
    factors, νif = if_gaussian_factors(smap, data, nl, f)
    dense = embed(nif, factors[1])
    dense .*= ν * νif
    buf = similar(dense)
    for k in 2:length(factors)
        mul_monomials!(buf, dense, nif, factors[k])
        dense, buf = buf, dense
    end
    return integrate_out(dense, nif, collect(0:(nif - 1));
                         step_sign = step_sign, order = order)[1]
end

"""Step B：暴力 ∫Da[IF Gaussian × 接口测度] 对照 ν·Pf·exp(D)（系数级）。

在 (IF ∪ 接口轨迹腿) 的压缩空间 8N 上进行。返回 max|Δcoeff|。"""
function brute_force_elim_check(; order::Int, nl::Int, data, f::Int,
                                conv, elim_sign::Int)
    ν = data.normalization
    # 压缩空间：IF 变量 + 接口 x 腿（nl = 接口层数）
    cvars = Vector{Var}()
    for i in 1:(4nl)
        push!(cvars, if_var(f, i))
    end
    for t in 1:nl, br in (:+, :-), stat in (:bar, :eta)
        push!(cvars, partner_var(t, br, stat, f))
    end
    cmap = Dict(v => i - 1 for (i, v) in enumerate(cvars))
    L = length(cvars)
    factors = Monomials[]
    ifactors, _ = if_gaussian_factors(cmap, data, nl, f)
    append!(factors, ifactors)
    append!(factors, interface_measure_factors(cmap, nl, 1))
    dense = embed(L, factors[1])
    dense .*= ν
    buf = similar(dense)
    for k in 2:length(factors)
        mul_monomials!(buf, dense, L, factors[k])
        dense, buf = buf, dense
    end
    ifsites = [cmap[if_var(f, i)] for i in 1:(4nl)]
    elim = integrate_out(dense, L, ifsites; step_sign = conv.step_sign,
                         order = conv.order)
    # 解析侧：Pf·exp(D)
    pf = pf_mine(data, nl, f; step_sign = conv.step_sign, order = conv.order)
    Ginv = inv(data.G)
    dfactors = delta_eff_factors(cmap, Ginv, nl, f, elim_sign)
    ana = embed(L, dfactors[1])
    buf2 = similar(ana)
    for k in 2:length(dfactors)
        mul_monomials!(buf2, ana, L, dfactors[k])
        ana, buf2 = buf2, ana
    end
    ana .*= pf
    return maximum(abs.(elim .- ana)), pf
end

# ---------------------------------------------------------------------------
# 模型与主流程
# ---------------------------------------------------------------------------

const T_TOT = 0.4
const BETA = 1.3
const EPS_D = 0.7

"""约定旋钮：step_sign × order = 积分定向；closure_k = 迹闭合单项式符号。"""
struct Conv
    step_sign::Int
    order::Symbol
    closure_k::Int
end

single_models() =
    (up = FiniteBath([0.4;;], [0.23 + 0.11im], BETA),)
single_imp() = (eps = EPS_D, rho0 = ComplexF64[1 0; 0 0])

function passfail(ok::Bool)
    return ok ? "PASS" : "FAIL"
end

maxabsdiff(A, B) = maximum(abs.(A .- B))

"""Step A：约定校准（两扇区稠密副本，N=1 单自旋，order 1）。

返回 (两扇区收缩定向 conv_two, 全部通过组合数)。"""
function step_a_calibrate()
    up = single_models().up
    imp = single_imp()
    println("═"^72)
    println("Step A：约定校准（两扇区稠密副本 N=1 单自旋 order 1，T=$T_TOT）")
    println("═"^72)
    sched = trotter_schedule(1, 1, T_TOT)
    rho_ref = fock_chain_single(up, EPS_D, sched, imp.rho0)
    @printf("  Fock 参照：tr ρ = %.12f, ρ[2,2] = %.10f\n",
            real(tr(rho_ref)), real(rho_ref[2, 2]))
    quick = get(ENV, "FORMB_QUICK", "0") == "1"
    found = Conv[]
    combos = quick ? [(1, :inc, 1)] : [(ss, ord, ck) for ss in (1, -1),
                                                        ord in (:inc, :dec),
                                                        ck in (-1, 1)]
    for (ss, ord, ck) in combos
        conv = Conv(ss, ord, ck)
        res = dense_case(order = 1, N = 1, nflav = 1,
                         models = [(ifdata = influence_matrix_G(up, sched.cells),)],
                         imp = imp, mode = :twosector, conv = conv)
        errz = abs(res.Z - 1)
        errr = maxabsdiff(res.rho, rho_ref)
        @printf("  step_sign=%+d order=%-3s closure_k=%+d : Z-1 = %.3e, ρ-err = %.3e\n",
                ss, string(ord), ck, errz, errr)
        (errz < 1e-9 && errr < 1e-9) && push!(found, conv)
    end
    isempty(found) &&
        error("Step A：没有任何约定组合同时给出 Z=1 与 ρ=Fock —— 需进一步 bisect")
    @printf("  → %d 个组合通过；取第一个锁定。\n\n", length(found))
    return found[1]
end

"""Step B：消元代数暴力对照 + elim_sign 锁定。"""
function step_b_elimination(conv_a)
    up = single_models().up
    println("═"^72)
    println("Step B：消元代数（暴力 ∫Da 对照 ν·Pf·exp(D)）")
    println("═"^72)
    locked = nothing
    quick = get(ENV, "FORMB_QUICK", "0") == "1"
    grid = quick ? [(1, 1)] :
                  vcat([(o, n) for o in (1, 2) for n in (1, 2)], [(4, 1)])
    for (order, N) in grid
        sched = trotter_schedule(order, N, T_TOT / N)
        nl = length(sched.cells)          # 接口层数（order 4 = 3N）
        data = influence_matrix_G(up, sched.cells)
        for es in (-1, 1)
            err, pf = brute_force_elim_check(order = order, nl = nl, data = data,
                                             f = 1, conv = conv_a, elim_sign = es)
            @printf("  order %d N=%d  elim_sign=%+d : max|Δcoeff| = %.3e   Pf = %.6f%+.6fim\n",
                    order, N, es, err, real(pf), imag(pf))
            if err < 1e-9
                locked === nothing && (locked = es)
                locked == es || error("Step B：两个 elim_sign 都匹配（异常）")
            end
        end
    end
    locked === nothing &&
        error("Step B：暴力消元与解析 Δ_eff 不匹配 —— bisect D 的符号/配对")
    @printf("  → elim_sign = %+d 锁定（D = −½(G⁻¹[i,k]−G⁻¹[k,i])）。\n\n", locked)
    return locked
end

"""消元式校准：pf 定向 Ω_a × 收缩定向 Ω_x × closure_k 的网格。

消元式收缩的定向与两扇区收缩的定向相差一个排列符号（IF site 从积分序
中移除），故独立标定；pf 的定向必须与暴力消元一致（同一 Ω_a）。"""
function step_c_calibrate(elim_sign)
    up = single_models().up
    imp = single_imp()
    sched = trotter_schedule(1, 1, T_TOT)
    rho_ref = fock_chain_single(up, EPS_D, sched, imp.rho0)
    data = influence_matrix_G(up, sched.cells)
    found = nothing
    for ss_a in (1, -1), ord_a in (:inc, :dec),
        ss_x in (1, -1), ord_x in (:inc, :dec), ck in (-1, 1)
        pf = pf_mine(data, length(sched.cells), 1; step_sign = ss_a, order = ord_a)
        res = dense_case(order = 1, N = 1, nflav = 1, models = [(Ginv = inv(data.G),)],
                         imp = imp, mode = :elim,
                         conv = Conv(ss_x, ord_x, ck), elim_sign = elim_sign,
                         pf_scale = pf)
        if abs(res.Z - 1) < 1e-9 && maxabsdiff(res.rho, rho_ref) < 1e-9
            found = (pf_conv = Conv(ss_a, ord_a, ck), int_conv = Conv(ss_x, ord_x, ck))
            break
        end
    end
    found === nothing &&
        error("Step C：消元式收缩没有任何定向组合给出 Z=1 且 ρ=Fock —— bisect")
    @printf("  消元式锁定：pf 定向 = (step_sign %+d, order %s), ",
            found.pf_conv.step_sign, string(found.pf_conv.order))
    @printf("收缩定向 = (step_sign %+d, order %s), closure_k = %+d)\n\n",
            found.int_conv.step_sign, string(found.int_conv.order),
            found.int_conv.closure_k)
    return found
end

"""Step C：消元式收缩 vs Fock（N=1..4 × order 1/2/4，单自旋）。"""
function step_c_elim_contraction(elim_sign, pf_conv, int_conv)
    up = single_models().up
    imp = single_imp()
    println("═"^72)
    println("Step C：消元式单扇区收缩 vs Fock 离散链（单自旋）")
    println("═"^72)
    all_ok = true
    quick = get(ENV, "FORMB_QUICK", "0") == "1"
    # order 4 的 lattice 是 3N 层（dense L=4·3N+6）：dense 只到 N=1
    grid = quick ? [(1, 1), (1, 2)] :
                  vcat([(o, n) for o in (1, 2) for n in 1:4], [(4, 1)])
    for (order, N) in grid
        sched = trotter_schedule(order, N, T_TOT / N)
        rho_ref = fock_chain_single(up, EPS_D, sched, imp.rho0)
        data = influence_matrix_G(up, sched.cells)
        pf = pf_mine(data, length(sched.cells), 1; step_sign = pf_conv.step_sign,
                     order = pf_conv.order)
        res = dense_case(order = order, N = N, nflav = 1,
                         models = [(Ginv = inv(data.G),)], imp = imp,
                         mode = :elim, conv = int_conv, elim_sign = elim_sign,
                         pf_scale = pf)
        errz = abs(res.Z - 1)
        errr = maxabsdiff(res.rho, rho_ref)
        n_ex = real(rho_ref[2, 2]) / real(tr(rho_ref))
        n_pb = real(res.rho[2, 2]) / real(tr(res.rho))
        ok = errz < 1e-9 && errr < 1e-9
        all_ok &= ok
        @printf("  order %d N=%d : Z-1 = %.2e, ρ-err = %.2e, ⟨n⟩ = %.10f (Fock %.10f)  %s\n",
                order, N, errz, errr, n_pb, n_ex, passfail(ok))
    end
    println()
    return all_ok
end

"""Step D：两自旋消元式收缩 vs Fock vs IFGTEMPO GMPS（权威对齐）。"""
function step_d_two_spin(elim_sign, pf_conv, int_conv)
    up = FiniteBath([0.4 0.1; 0.1 0.3], [0.23 + 0.11im, 0.15 - 0.05im], BETA)
    dn = FiniteBath([0.38 0.08; 0.08 0.28], [0.19 - 0.07im, 0.12 + 0.04im], BETA)
    imp = TwoSpinImpurityModel(0.35, -0.22, 1.0)
    println("═"^72)
    println("Step D：两自旋消元式收缩 vs Fock vs IFGTEMPO GMPS（gtempo_density_matrix_direct）")
    println("═"^72)
    # 断言：我们的两 flavor 布局 == timelocallattice
    for N in (1, 2)
        smap, vars, L = build_sites(N, 2)
        lat = timelocallattice(N)
        for v in vars
            g = IFGTEMPO.index(lat, GrassmannVariable(v.typ, v.t, v.br, v.stat, v.f))
            g == smap[v] || error("lattice 布局与 timelocallattice 不一致：$v")
        end
    end
    println("  lattice 布局 == timelocallattice（N=1,2 全变量逐位断言通过）")
    all_ok = true
    quick = get(ENV, "FORMB_QUICK", "0") == "1"
    # 两自旋消元 L = 12 + 8·nl：order 4（nl=3N）N=1 已是 2^36，dense 不做
    grid = quick ? [(1, 1)] : [(o, n) for o in (1, 2) for n in (1, 2)]
    for (order, N) in grid
        dt = T_TOT / N
        sched = trotter_schedule(order, N, dt)
        rho_ref = fock_chain_two_spin(up, dn, imp, sched)
        path = build_impurity_path(up, dn, N, dt; order = order,
                                   impurity_model = imp, cutoff = 1e-12)
        rho_gmps = gtempo_density_matrix_direct(path.factors...)
        data_up = influence_matrix_G(up, sched.cells)
        data_dn = influence_matrix_G(dn, sched.cells)
        nl4 = length(sched.cells)
        pf = pf_mine(data_up, nl4, 1; step_sign = pf_conv.step_sign,
                     order = pf_conv.order) *
             pf_mine(data_dn, nl4, 2; step_sign = pf_conv.step_sign,
                     order = pf_conv.order)
        res = dense_case(order = order, N = N, nflav = 2,
                         models = [(Ginv = inv(data_up.G),), (Ginv = inv(data_dn.G),)],
                         imp = imp, mode = :elim, conv = int_conv,
                         elim_sign = elim_sign, pf_scale = pf)
        errz = abs(res.Z - 1)
        errf = maxabsdiff(res.rho, rho_ref)
        errg = maxabsdiff(res.rho, rho_gmps)
        errgf = maxabsdiff(rho_gmps, rho_ref)
        ok = errz < 1e-8 && errf < 1e-8 && errg < 1e-7
        all_ok &= ok
        @printf("  order %d N=%d : Z-1 = %.2e | 消元-Fock %.2e | 消元-GMPS %.2e | GMPS-Fock %.2e  %s\n",
                order, N, errz, errf, errg, errgf, passfail(ok))
    end
    println()
    return all_ok
end

function main()
    println("形态 B（消元式单扇区）稠密系数级验证")
    @printf("参数：T = %.2f, β = %.1f, ε_d = %.2f\n\n", T_TOT, BETA, EPS_D)
    conv_two = step_a_calibrate()
    elim_sign = step_b_elimination(conv_two)
    calib = step_c_calibrate(elim_sign)
    ok_c = step_c_elim_contraction(elim_sign, calib.pf_conv, calib.int_conv)
    ok_d = step_d_two_spin(elim_sign, calib.pf_conv, calib.int_conv)
    println("═"^72)
    println("总结：Step C（单自旋 N=1..4 × order 1/2/4）", passfail(ok_c),
            "；Step D（两自旋 GMPS 对齐）", passfail(ok_d))
    println("锁定约定：两扇区收缩 (step_sign $(conv_two.step_sign), order ",
            "$(conv_two.order), closure_k $(conv_two.closure_k))；消元式 pf 定向 ",
            "(step_sign $(calib.pf_conv.step_sign), order $(calib.pf_conv.order))；",
            "消元收缩 (step_sign $(calib.int_conv.step_sign), order ",
            "$(calib.int_conv.order), closure_k $(calib.int_conv.closure_k))；",
            "elim_sign = $elim_sign")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
