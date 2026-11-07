# =============================================================================
# imag_final_solve.jl —— 组合矩方程求解（Z-行 + n-行）× (复/实 Δ) 模式判决
#
# 行：Z: m_z(ed; I) = z锚(ed)；n: m_g(ed; I) = n锚(ed)·m_z(ed; I)（线性化比值）。
# 模式：
#   a  Z-only  复Δ+ν   —— 期望 1 复维解族（Z 不唯一）
#   b  Z+n    复Δ+ν   —— 过定；命中 = 找到复合约定
#   c  Z+n    实Δ+ν   —— 强过定；命中 = Δ 实约定存在
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_final_solve.jl
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
const EVENS = Vector{Vector{Int}}([Int[]])
for k in (2, 4, 6), comb in combinations(0:(L-1), k)
    push!(EVENS, collect(comb))
end
const CFGIDX = [1 + sum(1 << p for p in A; init=0) for A in EVENS]

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
    # g-网格 εd（order 2）；z 用同一 εd 的 h 行
    eds = sort([k[2] for k in keys(gdata) if k[1] == 2])
    nz = [real(fock_anchor(CELLS, CELLS, ed)[2]) for ed in eds]
    Δfull = ring_delta(CELLS)
    println("g-网格 εd：", eds)

    """行向量组：variant ∈ {1:g_end, 2:1−g_end, 3:g_1, 4:1−g_1}。"""
    function build_system(mode::Symbol, variant::Int)
        rows = Vector{ComplexF64}[]; rhs = ComplexF64[]
        for (r, ed) in enumerate(eds)
            hz = hdata[(2, ed)][:, 1]
            gv = gdata[(2, ed)]            # 32×3（z, g_end, g_1）
            ge, g1 = gv[:, 2], gv[:, 3]
            vrow = variant == 1 ? ge : variant == 2 ? (hz .- ge) :
                   variant == 3 ? g1 : (hz .- g1)
            push!(rows, hz); push!(rhs, fock_anchor(CELLS, CELLS, ed)[1])
            if mode !== :Zonly
                push!(rows, vrow .- nz[r] .* hz); push!(rhs, 0.0 + 0.0im)
            end
        end
        return rows, rhs
    end

    function resid(x; mode::Symbol, realΔ::Bool, variant::Int=1)
        rows, rhs = build_system(mode, variant)
        if realΔ
            Δ = [ComplexF64(x[1], 0) ComplexF64(x[2], 0); ComplexF64(x[3], 0) ComplexF64(x[4], 0)]
            ν = x[5]
        else
            Δ = [ComplexF64(x[1], x[2]) ComplexF64(x[3], x[4]);
                 ComplexF64(x[5], x[6]) ComplexF64(x[7], x[8])]
            ν = x[9]
        end
        itab = if_table(Δ)
        ivec = itab[CFGIDX]
        return [ν * dot(rows[r], ivec) - rhs[r] for r in eachindex(rows)]
    end
    function newton(x0; mode::Symbol, realΔ::Bool, variant::Int=1, iters=300)
        x = copy(x0)
        npar = realΔ ? 5 : 9
        for it in 1:iters
            F = resid(x; mode=mode, realΔ=realΔ, variant=variant)
            J = zeros(ComplexF64, length(F), npar)
            for p in 1:npar
                xp = copy(x); xp[p] += 1e-7
                xm = copy(x); xm[p] -= 1e-7
                J[:, p] = (resid(xp; mode=mode, realΔ=realΔ, variant=variant) - resid(xm; mode=mode, realΔ=realΔ, variant=variant)) ./ 2e-7
            end
            dx = real.(pinv(J' * J + 1e-11 * I, rtol=1e-10) * (-J' * F))
            x += dx
            norm(F) < 1e-13 && break
            norm(dx) < 1e-15 && break
        end
        return x, norm(resid(x; mode=mode, realΔ=realΔ, variant=variant))
    end
    Random.seed!(7)
    for (tag, mode, realΔ, variant) in [("a Z-only 复Δ", :Zonly, false, 1),
                               ("b1 Z+n 复Δ var=g_end", :Zn, false, 1),
                               ("b2 Z+n 复Δ var=1−g_end", :Zn, false, 2),
                               ("b3 Z+n 复Δ var=g_1", :Zn, false, 3),
                               ("b4 Z+n 复Δ var=1−g_1", :Zn, false, 4),
                               ("c Z+n 实Δ var=g_end", :Zn, true, 1)]
        best = (Inf, nothing)
        hits = 0
        for trial in 1:40
            x0 = realΔ ?
                 vcat(vec(real.(Δfull)) .+ 0.5randn(4), [0.5 + rand()]) :
                 vcat(vec(real.(Δfull)) .+ 0.5randn(4), zeros(4), [0.5 + rand()])
            x, rf = newton(x0; mode=mode, realΔ=realΔ, variant=variant)
            rf < best[1] && (best = (rf, copy(x)))
            rf < 1e-11 && (hits += 1)
            hits >= 6 && break
        end
        println("\n[$(tag)] 命中 = $(hits)，最优 |F| = ", best[1])
        if best[1] < 1e-9
            x = best[2]
            if realΔ
                println("  Δ = ", [round(v, digits=6) for v in x[1:4]], "  ν = ", round(x[5], digits=6))
            else
                Δ = [ComplexF64(x[1], x[2]) ComplexF64(x[3], x[4]); ComplexF64(x[5], x[6]) ComplexF64(x[7], x[8])]
                println("  Δ = ", round.(Δ, digits=6), "  ν = ", round(x[9], digits=6))
                println("  |Im Δ| max = ", round(maximum(abs.(imag.(Δ))), digits=4))
            end
            println("  对照 D_full = ", round.(Δfull, digits=6))
        end
    end
end
main()
