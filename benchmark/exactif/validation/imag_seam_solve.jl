# =============================================================================
# imag_seam_solve.jl —— (N+1)×(N+1) 缝位影响矩阵判决
#
# 证据链：pinv-I* 的额外支撑全部锚在切片 1（缝位）腿上 ⇒ 所需 IF = 全
# (N+1)² 个 (bar_s, eta_s') 对（s, s' ∈ {缝,1..N}）。库 partialif 只铺 corr
# 节点行 ⇒ 缝行/列缺失 = 最终结构性缺块。
# 本脚本：3×3 Δ̃（18 实参）+ ν vs (Z + n) 矩行，多起点牛顿判决。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_seam_solve.jl
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
    Z = Matrix{ComplexF64}(I, dim, dim); Zn = Matrix{ComplexF64}(I, dim, dim)
    for p in length(cells):-1:1
        W = exp(-cells[p] .* Hq) * exp(-gate_durs[p] .* Him)
        Z = W * Z; Zn = W * Zn
    end
    Zn = ops[1]' * ops[1] * Zn
    return tr(Z), tr(Zn) / tr(Z)
end

const NL = 2
const N = 2
const CELLS = fill(BETA / N, N)
const LAT = GTEMPO.GrassmannLattice(N=NL, δτ=BETA / NL, contour=:imag, bands=1)
const L = length(LAT)
const EVENS = Vector{Vector{Int}}([Int[]])
for k in (2, 4, 6), comb in combinations(0:(L-1), k)
    push!(EVENS, collect(comb))
end
const CFGIDX = [1 + sum(1 << p for p in A; init=0) for A in EVENS]
# 3×3 切片位（切片 1=缝, 2, 3）：bar/eta 位置（0-based）
const SL = [1, 2, 3]                       # 切片编号
const BARP = [LAT[GrassmannContourIndex(s, :τ, 1, true)] - 1 for s in SL]
const ETAP = [LAT[GrassmannContourIndex(s, :τ, 1, false)] - 1 for s in SL]
const PR = [(BARP[a], ETAP[b], BARP[a] > ETAP[b] ? -1.0 : 1.0) for a in 1:3, b in 1:3]

function mul_tables(a, b, L)
    monos = Tuple{Vector{Int}, ComplexF64}[]
    for idx in 1:length(b)
        b[idx] == 0 && continue
        push!(monos, ([p for p in 0:(L-1) if (idx - 1) >> p & 1 == 1], b[idx]))
    end
    return mul_monomials!(zeros(ComplexF64, 1 << L), a, L, Monomials(monos))
end
"""全 3×3 IF 表：∏_a exp(Σ_b Δ̃[a,b]·θ_bar_a θ_eta_b)，行序 = 参数 order（缝先/后）。"""
function if3_table(x; seam_first::Bool)
    M = [ComplexF64(x[2*(3*(a-1)+b)-1], x[2*(3*(a-1)+b)]) for a in 1:3, b in 1:3]
    acc = zeros(ComplexF64, 1 << L); acc[1] = 1.0
    roworder = seam_first ? (1, 2, 3) : (2, 3, 1)
    for a in roworder
        monos = Tuple{Vector{Int}, ComplexF64}[(Int[], 1.0 + 0.0im)]
        for b in 1:3
            iszero(M[a, b]) && continue
            pb, pe, sg = PR[a, b]
            push!(monos, (sort([pb, pe]), M[a, b] * sg))
        end
        acc = mul_monomials!(zeros(ComplexF64, 1 << L), acc, L, Monomials(monos))
    end
    return acc
end
function load_csv(path; cols::Int)
    out = Dict{Tuple{Int,Float64}, Matrix{ComplexF64}}()
    for (i, ln) in enumerate(eachline(path))
        i == 1 && continue
        f = split(strip(ln), ",")
        key = (parse(Int, f[1]), parse(Float64, f[2]))
        ib = parse(Int, f[3])
        M = get!(out, key, zeros(ComplexF64, 32, cols))
        for c in 1:cols
            M[ib, c] = ComplexF64(parse(Float64, f[2 + 2c]), parse(Float64, f[3 + 2c]))
        end
    end
    return out
end

function main()
    hdata = load_csv(joinpath(GT_ROOT, "benchmark", "exactif", "validation", "imag_h.csv"); cols=1)
    gdata = load_csv(joinpath(GT_ROOT, "benchmark", "exactif", "validation", "imag_g.csv"); cols=3)
    eds = sort([k[2] for k in keys(gdata) if k[1] == 2])
    nz = [real(fock_anchor(CELLS, CELLS, ed)[2]) for ed in eds]

    function build_rows(variant::Int)
        rows = Vector{ComplexF64}[]; rhs = ComplexF64[]
        for (r, ed) in enumerate(eds)
            hz = hdata[(2, ed)][:, 1]
            gv = gdata[(2, ed)]
            ge, g1 = gv[:, 2], gv[:, 3]
            vrow = variant == 1 ? ge : variant == 2 ? (hz .- ge) :
                   variant == 3 ? g1 : (hz .- g1)
            push!(rows, hz); push!(rhs, fock_anchor(CELLS, CELLS, ed)[1])
            push!(rows, vrow .- nz[r] .* hz); push!(rhs, 0.0 + 0.0im)
        end
        return rows, rhs
    end
    function resid(x; variant::Int, seam_first::Bool)
        rows, rhs = build_rows(variant)
        itab = if3_table(x; seam_first=seam_first)
        ivec = itab[CFGIDX]
        return [x[19] * dot(rows[r], ivec) - rhs[r] for r in eachindex(rows)]
    end
    function newton(x0; variant::Int, seam_first::Bool, iters=400)
        x = copy(x0)
        for it in 1:iters
            F = resid(x; variant=variant, seam_first=seam_first)
            J = zeros(ComplexF64, length(F), 19)
            for p in 1:19
                xp = copy(x); xp[p] += 1e-7
                xm = copy(x); xm[p] -= 1e-7
                J[:, p] = (resid(xp; variant=variant, seam_first=seam_first) -
                           resid(xm; variant=variant, seam_first=seam_first)) ./ 2e-7
            end
            dx = real.(pinv(J' * J + 1e-11 * I, rtol=1e-10) * (-J' * F))
            x += dx
            norm(F) < 1e-13 && break
            norm(dx) < 1e-15 && break
        end
        return x, norm(resid(x; variant=variant, seam_first=seam_first))
    end
    Random.seed!(11)
    results = []
    for variant in (1, 2, 3, 4), seam_first in (true, false)
        best = (Inf, nothing); hits = 0
        for trial in 1:40
            x0 = vcat(0.6 .* randn(18), [0.5 + rand()])
            x, rf = newton(x0; variant=variant, seam_first=seam_first)
            rf < best[1] && (best = (rf, copy(x)))
            rf < 1e-11 && (hits += 1)
            hits >= 3 && break
        end
        @printf("variant=%d seam_first=%d : 命中 = %d，最优 |F| = %.3e\n",
                variant, seam_first, hits, best[1])
        push!(results, (best[1], variant, seam_first, best[2]))
    end
    sort!(results, by=t -> t[1])
    rf, variant, seam_first, x = results[1]
    println("\n全局最优：variant=$variant seam_first=$seam_first |F| = $rf")
    # 抛光：低阻尼多轮牛顿 + Z/n 残差拆分
    for round in 1:6
        x, rf = newton(x; variant=variant, seam_first=seam_first, iters=600)
        println("  抛光 round $round: |F| = ", rf)
        rf < 1e-12 && break
    end
    rowsA, rhsA = build_rows(variant)
    itab = if3_table(x; seam_first=seam_first)
    ivec = itab[CFGIDX]
    zres = [abs(x[19] * dot(rowsA[2r-1], ivec) - rhsA[2r-1]) for r in 1:length(eds)]
    nres = [abs(x[19] * dot(rowsA[2r], ivec) - rhsA[2r]) for r in 1:length(eds)]
    println("  Z-行残差: ", [round(v, sigdigits=2) for v in zres])
    println("  n-行残差: ", [round(v, sigdigits=2) for v in nres])
    if rf < 1e-9
        M = [ComplexF64(x[2*(3*(a-1)+b)-1], x[2*(3*(a-1)+b)]) for a in 1:3, b in 1:3]
        println("Δ̃ = ")
        display(round.(M, digits=5))
        println("ν = ", round(x[19], digits=6))
    end
end
main()
