# =============================================================================
# imag_moment_solve.jl —— 秩-3 矩方程的精确求解（纯 row-exp 族 + 分辨率扩展族）
#
# 关键事实（已数值验证）：h 矩阵（order 2，全 εd）秩恰为 3 ⇒ Z 观测量 = 3 个
# 复矩 m_k(I) = z_k。纯 row-exp 族 (Δ∈C⁴, ν) 10 实参 vs 6 实约束 ⇒ 代数上应有
# 解；此前 LM 失败 = 坏盆地/冗余阻尼伪影。本脚本：SVD 提取秩-3 基，多起点
# 牛顿解矩方程。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_moment_solve.jl
# =============================================================================
using Printf
using LinearAlgebra
using Random
using Combinatorics: combinations
using GTEMPO
using Z2Tensors: NoTruncation

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(GT_ROOT, "benchmark", "exactif", "validation", "formb_dense.jl"))

const BETA = 1.0
const WB, AB = 1.0, 0.5
HB = ComplexF64[WB;;]
VB = ComplexF64[sqrt(AB) + 0.0im]
hcoupled() = begin
    m = length(VB)
    h = zeros(ComplexF64, m + 1, m + 1)
    for k in 1:m
        h[1, k+1] = conj(VB[k]); h[k+1, 1] = VB[k]; h[k+1, k+1] = HB[k, k]
    end
    h
end
annih_op(dim::Int, mode::Int) = begin
    op = zeros(ComplexF64, dim, dim)
    for n in 0:(dim - 1)
        (n >> mode) & 1 == 1 || continue
        sgn = isodd(count_ones(n & ((1 << mode) - 1))) ? -1 : 1
        op[(n ⊻ (1 << mode)) + 1, n + 1] = sgn
    end
    op
end
function fock_anchor(cells, gate_durs, eps_d)
    m = length(VB); nt = 1 + m; dim = 1 << nt
    ops = [annih_op(dim, k) for k in 0:(nt-1)]
    Hq = zeros(ComplexF64, dim, dim); h = hcoupled()
    for i in 1:nt, j in 1:nt
        iszero(h[i, j]) || (Hq .+= h[i, j] .* (ops[i]' * ops[j]))
    end
    Him = eps_d .* (ops[1]' * ops[1])
    Z = Matrix{ComplexF64}(I, dim, dim)
    for p in length(cells):-1:1
        W = exp(-cells[p] .* Hq) * exp(-gate_durs[p] .* Him)
        Z = W * Z
    end
    return tr(Z)
end
function ring_delta(cells::Vector{Float64})
    N = length(cells)
    h = hcoupled()
    mm = size(h, 1)
    Us = [exp(-cells[r] .* h) for r in 1:N]
    M = zeros(ComplexF64, mm * N, mm * N)
    for r in 1:N
        sgn = r == N ? +1.0 : -1.0
        for s2 in 1:mm, s1 in 1:mm
            M[(r-1)*mm + s2, (r-1)*mm + s1] += sgn * Us[r][s2, s1]
        end
        nxt = r == N ? 1 : r + 1
        for s in 1:mm
            M[(r-1)*mm + s, (nxt-1)*mm + s] += 1.0
        end
    end
    Iset = [(r-1)*mm + 1 for r in 1:N]
    Bset = [(r-1)*mm + s for r in 1:N for s in 2:mm]
    return -(M[Iset, Iset] - M[Iset, Bset] * (M[Bset, Bset] \ M[Bset, Iset]))
end

const NL = 2
const N = 2
const CELLS = fill(BETA / N, N)
const LAT = GTEMPO.GrassmannLattice(N=NL, δτ=BETA / NL, contour=:imag, bands=1)
const L = length(LAT)
const BARPOS(i) = LAT[GrassmannContourIndex(i + 1, :τ, 1, true)] - 1
const ETAPOS(j) = LAT[GrassmannContourIndex(j + 1, :τ, 1, false)] - 1
const ROWPOS = [(BARPOS(i), ETAPOS(j), BARPOS(i) > ETAPOS(j) ? -1.0 : 1.0) for i in 1:N, j in 1:N]

function mul_tables(a, b, L)
    monos = Tuple{Vector{Int}, ComplexF64}[]
    for idx in 1:length(b)
        b[idx] == 0 && continue
        push!(monos, ([p for p in 0:(L-1) if (idx - 1) >> p & 1 == 1], b[idx]))
    end
    return mul_monomials!(zeros(ComplexF64, 1 << L), a, L, Monomials(monos))
end
function if_table(Δ)
    acc = zeros(ComplexF64, 1 << L); acc[1] = 1.0
    for i in 1:N
        monos = Tuple{Vector{Int}, ComplexF64}[(Int[], 1.0 + 0.0im)]
        for j in 1:N
            iszero(Δ[i, j]) && continue
            pb, pe, sg = ROWPOS[i, j]
            push!(monos, (sort([pb, pe]), Δ[i, j] * sg))
        end
        acc = mul_monomials!(zeros(ComplexF64, 1 << L), acc, L, Monomials(monos))
    end
    return acc
end
function res_pair_table(cs::Vector{ComplexF64})
    acc = zeros(ComplexF64, 1 << L); acc[1] = 1.0
    for p in 1:N
        s1 = p + 1
        s2 = p == N ? 2 : p + 2
        pb = LAT[GrassmannContourIndex(s1, :τ, 1, true)] - 1
        pe = LAT[GrassmannContourIndex(s2, :τ, 1, false)] - 1
        sg = (pb > pe) ? -1.0 : 1.0
        acc = mul_monomials!(zeros(ComplexF64, 1 << L), acc, L,
                             Monomials([(Int[], 1.0 + 0.0im), (sort([pb, pe]), cs[p] * sg)]))
    end
    return acc
end

function main()
    # ---- 读 h 数据
    data = Dict{Tuple{Int,Float64}, Vector{ComplexF64}}()
    for (i, ln) in enumerate(eachline(joinpath(GT_ROOT, "benchmark", "exactif", "validation", "imag_h.csv")))
        i == 1 && continue
        f = split(strip(ln), ",")
        key = (parse(Int, f[1]), parse(Float64, f[2]))
        get!(data, key, zeros(ComplexF64, 32))[parse(Int, f[3])] = ComplexF64(parse(Float64, f[4]), parse(Float64, f[5]))
    end
    eds = sort([k[2] for k in keys(data) if k[1] == 2])
    # 偶集合列表（与 imag_h_grid.jl 的枚举一致）
    evensets = Vector{Vector{Int}}([Int[]])
    for k in (2, 4, 6), comb in combinations(0:(L-1), k)
        push!(evensets, collect(comb))
    end
    cfgidx = [1 + sum(1 << p for p in A; init=0) for A in evensets]
    Hm = hcat([data[(2, ed)] for ed in eds]...)          # 32 × nεd
    SVD = svd(Hm)
    rk = count(>(1e-8), SVD.S)
    println("h 秩 = $rk")
    U3 = SVD.U[:, 1:rk]                                   # 32×3 基
    # 锚的矩：z_k = (U3ᵀ-econ 投影后的坐标)——直接把锚向量投到基上
    zs = [fock_anchor(CELLS, CELLS, ed) for ed in eds]
    # （zs 是 nεd 维向量；投影见下方 Hmᵀ 行空间）
    # 注意：矩方程在“h-行空间”坐标下：m(I) := [Σ_b I_b h_ed(b)]_ed ∈ C^nεd，
    # 要求 m(I) = zs。投影：令 P = U3 (32×3 作用于 b-指标)——而 m(I) = Hmᵀ I ∈ C^nεd。
    # Hmᵀ I = Σ_k u_k·(U3列ᵀ... — 直接用 Hmᵀ 的 SVD 行空间：Hmᵀ = W Σ Vᵀ（nεd×32），
    # m(I) = Hmᵀ I；要求 = zs ∈ range(Hmᵀ)（验证 zs 残差！）
    Wt = svd(Hm')
    rkt = count(>(1e-8), Wt.S)
    println("Hmᵀ 秩 = $rkt")
    proj_zs = Wt.U * (Wt.U' * zs)
    println("锚 zs 在 Hmᵀ 行空间的投影残差 = ", norm(proj_zs - zs), "（应为 0）")

    # 矩方程：F(Δ, cs, ν) = Hmᵀ·I(Δ,cs)·ν − zs = 0（nεd 维，秩 rkt）
    function moment_resid(x; use_res::Bool)
        Δ = [ComplexF64(x[1], x[2]) ComplexF64(x[3], x[4]);
             ComplexF64(x[5], x[6]) ComplexF64(x[7], x[8])]
        itab = if_table(Δ)
        if use_res
            cs = [ComplexF64(x[9], x[10]), ComplexF64(x[11], x[12])]
            itab = mul_tables(itab, res_pair_table(cs), L)
        end
        m = [dot(data[(2, ed)], itab[cfgidx]) for ed in eds]
        return (use_res ? x[13] : x[9]) .* m - zs
    end
    function newton_solve(x0; use_res::Bool, iters=200)
        x = copy(x0)
        npar = use_res ? 13 : 9
        for it in 1:iters
            F = moment_resid(x; use_res=use_res)
            J = zeros(ComplexF64, length(F), npar)
            for p in 1:npar
                xp = copy(x); xp[p] += 1e-7
                xm = copy(x); xm[p] -= 1e-7
                J[:, p] = (moment_resid(xp; use_res=use_res) - moment_resid(xm; use_res=use_res)) ./ 2e-7
            end
            dx = real.(pinv(J' * J + 1e-10 * I, rtol=1e-10) * (-J' * F))
            x += dx
            norm(F) < 1e-13 && break
            norm(dx) < 1e-14 && break
        end
        return x, norm(moment_resid(x; use_res=use_res))
    end
    Δfull = ring_delta(CELLS)
    Random.seed!(42)
    for (tag, use_res) in [("纯 row-exp", false), ("row-exp×res", true)]
        best = (Inf, nothing)
        hits = 0
        for trial in 1:60
            x0 = use_res ?
                 vcat(vec(real.(Δfull)) .+ 0.5randn(4), zeros(4), 0.3 .* randn(2), [0.2 * (rand() - 0.5), 0.2 * (rand() - 0.5)], 1.0) :
                 vcat(vec(real.(Δfull)) .+ 0.5randn(4), zeros(4), [1.0])
            x, rf = newton_solve(x0; use_res=use_res)
            rf < best[1] && (best = (rf, copy(x)))
            rf < 1e-11 && (hits += 1; hits <= 3 && println("  [$tag] trial $trial 命中: Δ=", [round(real(v), digits=5) for v in x[1:8]]))
            hits >= 5 && break
        end
        println("[$(tag)] 命中数 = $(hits)，最优 |F| = ", best[1])
        if best[1] < 1e-9
            x = best[2]
            Δ = [ComplexF64(x[1], x[2]) ComplexF64(x[3], x[4]); ComplexF64(x[5], x[6]) ComplexF64(x[7], x[8])]
            println("  Δ = ", round.(Δ, digits=6))
            use_res && println("  c = ", [round(ComplexF64(x[9], x[10]), digits=6), round(ComplexF64(x[11], x[12]), digits=6)],
                               "  ν = ", round(x[13], digits=6))
            use_res || println("  ν = ", round(x[9], digits=6))
            println("  对照 D_full = ", round.(ring_delta(CELLS), digits=6))
        end
    end
end
main()
