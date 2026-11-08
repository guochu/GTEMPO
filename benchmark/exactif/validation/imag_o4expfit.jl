# =============================================================================
# imag_o4expfit.jl —— order-4 的 7×7 exp-form 影响矩阵拟合（干净行）
#
# 背景表示依赖发现后：seam_solve 的失败行作废。本脚本用 imag_n2o4ma_*.csv
# （nl 修复后的干净行）对 (nl+1)² = 7×7 exp-form 族（Δ̃ 98 实参 + ν）做牛顿拟合。
# 命中 ⇒ order-4 的 exp-form 闭合存在 → 构造 I-MPS → 端到端。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_o4expfit.jl
# =============================================================================
using Printf
using LinearAlgebra
using Random
using Combinatorics: combinations
using GTEMPO
using Z2Tensors: NoTruncation

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(GT_ROOT, "benchmark", "exactif", "validation", "formb_dense.jl"))
include(joinpath(GT_ROOT, "benchmark", "exactif", "src", "NativeExact.jl"))
using .NativeExact

const BETA = 1.0
const NL = 2
const WB, AB = 1.0, 0.5
HB = ComplexF64[WB;;]
VB = ComplexF64[sqrt(AB) + 0.0im]
hcoupled() = NativeExact._coupled_h(VB, HB)

annih_op_f(dim::Int, mode::Int) = begin
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
    ops = [annih_op_f(dim, k) for k in 0:(nt-1)]
    Hq = zeros(ComplexF64, dim, dim)
    h = hcoupled()
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

const SCHED = NativeExact.native_schedule(4, NL, BETA / NL)
const RLINKS = vcat(SCHED.links[1] + SCHED.links[end], SCHED.links[2:end-1])
const NL_INT = length(RLINKS)             # 6
const CELLS = SCHED.cells
const LAT = GTEMPO.GrassmannLattice(N=NL_INT, δτ=BETA / NL_INT, contour=:imag, bands=1)
const L = length(LAT)                     # 14
const K = LAT.k                           # 7
const EVENS = Vector{Vector{Int}}([Int[]])
for k in 2:2:L, comb in combinations(0:(L-1), k)
    push!(EVENS, collect(comb))
end
const NE = length(EVENS)                  # 8192
const CFGIDX = [1 + sum(1 << p for p in A; init=0) for A in EVENS]

# 7 切片位（切片 1=缝, 2..7）
const SL = collect(1:K)
const BARP = [LAT[GrassmannContourIndex(s, :τ, 1, true)] - 1 for s in SL]
const ETAP = [LAT[GrassmannContourIndex(s, :τ, 1, false)] - 1 for s in SL]
const PR = [(BARP[a], ETAP[b], BARP[a] > ETAP[b] ? -1.0 : 1.0) for a in 1:K, b in 1:K]

function load_csv(path; cols::Int)
    out = Dict{Float64, Matrix{ComplexF64}}()
    for (i, ln) in enumerate(eachline(path))
        i == 1 && continue
        f = split(strip(ln), ",")
        ed = parse(Float64, f[1])
        ib = parse(Int, f[2])
        M = get!(out, ed, zeros(ComplexF64, NE, cols))
        for c in 1:cols
            M[ib, c] = ComplexF64(parse(Float64, f[2 + 2c - 1]), parse(Float64, f[2 + 2c]))
        end
    end
    return out
end

"""7×7 exp-form IF 表：∏_a exp(Σ_b Δ̃[a,b]·θ_bar_a θ_eta_b)（行序 = a 升序）。"""
function if7_table(x::Vector{Float64})
    M = [ComplexF64(x[2*(K*(a-1)+b)-1], x[2*(K*(a-1)+b)]) for a in 1:K, b in 1:K]
    acc = zeros(ComplexF64, 1 << L); acc[1] = 1.0
    for a in 1:K
        monos = Tuple{Vector{Int}, ComplexF64}[(Int[], 1.0 + 0.0im)]
        for b in 1:K
            iszero(M[a, b]) && continue
            pb, pe, sg = PR[a, b]
            push!(monos, (sort([pb, pe]), M[a, b] * sg))
        end
        acc = mul_monomials!(zeros(ComplexF64, 1 << L), acc, L, Monomials(monos))
    end
    return acc
end

function main()
    outdir = joinpath(GT_ROOT, "benchmark", "exactif", "validation")
    hdata = load_csv(joinpath(outdir, "imag_n2o4ma_h1.csv"); cols=1)
    merge!(hdata, load_csv(joinpath(outdir, "imag_n2o4ma_h2.csv"); cols=1))
    gdata = load_csv(joinpath(outdir, "imag_n2o4ma_g1.csv"); cols=3)
    merge!(gdata, load_csv(joinpath(outdir, "imag_n2o4ma_g2.csv"); cols=3))
    heds = sort(collect(keys(hdata)))
    geds = sort(collect(keys(gdata)))
    @printf("h %d 点 / g %d 点\n", length(heds), length(geds))

    rows = Vector{Vector{ComplexF64}}(); rhs = ComplexF64[]; tags = String[]
    for ed in heds
        push!(rows, hdata[ed][:, 1]); push!(rhs, fock_anchor(CELLS, RLINKS, ed)[1])
        push!(tags, "Z@$ed")
    end
    for ed in geds
        hz = hdata[ed][:, 1]
        vend = gdata[ed][:, 2]
        nanc = real(fock_anchor(CELLS, RLINKS, ed)[2])
        push!(rows, hz .- vend .- nanc .* hz); push!(rhs, 0.0 + 0.0im)
        push!(tags, "n@$ed")
    end
    @printf("联合行 = %d（未知 98+1）\n", length(rows))

    function resid(x)
        itab = if7_table(x)
        ivec = itab[CFGIDX]
        return [x[2K*K+1] * dot(rows[r], ivec) - rhs[r] for r in eachindex(rows)]
    end
    function newton(x0; iters=60)
        x = copy(x0)
        for it in 1:iters
            F = resid(x)
            J = zeros(ComplexF64, length(F), 2K*K+1)
            for p in 1:(2K*K+1)
                xp = copy(x); xp[p] += 1e-6
                xm = copy(x); xm[p] -= 1e-6
                J[:, p] = (resid(xp) - resid(xm)) ./ 2e-6
            end
            dx = real.(pinv(J' * J + 1e-10 * I; rtol=1e-11) * (-J' * F))
            x += dx
            @printf("  iter %2d : |F| = %.3e\n", it, norm(F)); flush(stdout)
            norm(F) < 1e-12 && break
            norm(dx) < 1e-14 && break
        end
        return x, norm(resid(x))
    end
    Random.seed!(5)
    best = (Inf, nothing)
    for trial in 1:8
        x0 = vcat(0.3 .* randn(2K*K), [0.5 + rand()])
        x, rf = newton(x0)
        rf < best[1] && (best = (rf, copy(x)))
        rf < 1e-11 && break
    end
    @printf("\n最优 |F| = %.3e\n", best[1])
    if best[1] < 1e-9
        x = best[2]
        M = [ComplexF64(x[2*(K*(a-1)+b)-1], x[2*(K*(a-1)+b)]) for a in 1:K, b in 1:K]
        println("Δ̃(7×7) = ")
        display(round.(M, digits=6))
        @printf("ν = %.10f\n", x[2K*K+1])
    end
end
main()
